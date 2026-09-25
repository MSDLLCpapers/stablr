# predict_tos() propagates from the fitted value at time 0 and advances by each
# segment's DURATION, so absolute time is never read. A schedule not anchored at
# 0 was therefore evaluated as if it began at 0: t_start = 6, t_end = 12 returned
# y(6) = 88.69 instead of y(12) = 78.66 (11.3% high, reporting the product as
# fresher than it is), and t_start = 96, t_end = 102 returned that same y(6)
# = 88.69 rather than y(102) = 13.00. The returned row still carried the t_start
# and t_end it was handed, so the output asserted a timeline it had not used.
#
# The contiguity guard could not catch this: it compares starts[-1] against
# ends[-n], i.e. interior joins only, and is reached only when nrow > 1. The gap
# between 0 and the first t_start -- the one that rebases the entire schedule --
# was unchecked.
#
# t_start cannot be honoured instead of rejected. There is no state to propagate
# from at t = 6 except y0, and the schedule supplies no temperature for [0, 6),
# so the leading stretch must be given as its own segment. These tests pin the
# rejection AND the value the recommended two-segment form produces, so a fix
# that errors without giving a usable alternative would not pass.

origin_fit <- function() {
  g <- expand.grid(batch = c("A", "B"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- c(A = 100, B = 97)[g$batch] * exp(-k * g$time)
  suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                               intercept_cols = "batch"))
}

# Oracle: evaluate the fitted formula's RHS directly at absolute time t.
oracle <- function(fit, ...) {
  env <- list2env(as.list(stats::coef(fit)), parent = baseenv())
  eval(stats::formula(fit)[[3]],
       envir = data.frame(..., stringsAsFactors = FALSE), enclos = env)
}

# Segment rows only: these tests are about how the schedule's timeline is read,
# so the time-0 baseline row predict_tos() prepends is dropped here. The test at
# the bottom of the file checks it explicitly.
ptos <- function(fit, sch) {
  tos_segments(suppressWarnings(predict_tos(fit, sch, temp_unit = "K")))
}

test_that("a schedule not starting at 0 is rejected, not silently rebased", {
  fit <- origin_fit()

  err <- expect_error(
    ptos(fit, data.frame(t_start = 6, t_end = 12, temp = 298.15, batch = "A",
                         stringsAsFactors = FALSE))
  )
  # Previously returned 88.69204367 -- y(6), the duration -- with no warning.
  expect_match(conditionMessage(err), "begin at time 0")
  # The message must name the actual offending value and the remedy, since the
  # correct schedule is not obvious from "must begin at 0" alone.
  expect_match(conditionMessage(err), "6", fixed = TRUE)
  # expect_match(conditionMessage(err), "leading segment")
})

test_that("a far-future window is rejected rather than collapsed to its duration", {
  # t_start = 96 / t_end = 102 returned the same 88.69204367 as [6, 12]: the two
  # schedules are indistinguishable once only the duration is read, though the
  # true values differ by a factor of seven.
  fit <- origin_fit()
  expect_equal(oracle(fit, batch = "A", time = 102, tempk = 298.15),
               13.00287109, tolerance = 1e-6)

  err <- expect_error(
    ptos(fit, data.frame(t_start = 96, t_end = 102, temp = 298.15, batch = "A",
                         stringsAsFactors = FALSE))
  )
  expect_match(conditionMessage(err), "begin at time 0")
})

test_that("the documented two-segment form gives the value t_start could not", {
  # The remedy the error recommends must actually reproduce y(12). Without this
  # the fix would only be a refusal.
  fit <- origin_fit()
  sch <- data.frame(t_start = c(0, 6), t_end = c(6, 12),
                    temp = c(298.15, 298.15), batch = c("A", "A"),
                    stringsAsFactors = FALSE)
  p <- ptos(fit, sch)

  expect_equal(p$fit[1], oracle(fit, batch = "A", time = 6, tempk = 298.15),
               tolerance = 1e-8)
  expect_equal(p$fit[2], oracle(fit, batch = "A", time = 12, tempk = 298.15),
               tolerance = 1e-8)
  # The value the rebased single-segment schedule wrongly reported, and the one
  # it should have reported.
  expect_equal(p$fit[1], 88.69204367, tolerance = 1e-6)
  expect_equal(p$fit[2], 78.66278611, tolerance = 1e-6)

  # A genuine excursion: months 0-6 refrigerated, 6-12 at room temperature. The
  # second segment must continue from the first, not restart from y0.
  exc <- data.frame(t_start = c(0, 6), t_end = c(6, 12),
                    temp = c(278.15, 298.15), batch = c("A", "A"),
                    stringsAsFactors = FALSE)
  pe <- ptos(fit, exc)
  expect_lt(pe$fit[2], pe$fit[1])
  expect_gt(pe$fit[1], p$fit[1])          # colder first leg degrades less
  expect_lt(pe$fit[2], oracle(fit, batch = "A", time = 6, tempk = 298.15))
})

test_that("a leading gap is rejected even when every interior join is valid", {
  # The exact blind spot: starts[-1] == ends[-n] holds (9 == 9), so the interior
  # check passed and the schedule was evaluated as 0-3 then 3-6.
  fit <- origin_fit()
  err <- expect_error(
    ptos(fit, data.frame(t_start = c(6, 9), t_end = c(9, 12),
                         temp = c(298.15, 313.15), batch = c("A", "A"),
                         stringsAsFactors = FALSE))
  )
  expect_match(conditionMessage(err), "begin at time 0")
  # Must not be reported as an interior contiguity failure -- there is none.
  expect_no_match(conditionMessage(err), "contiguous")
})

test_that("a negative start time is rejected", {
  fit <- origin_fit()
  expect_error(
    ptos(fit, data.frame(t_start = -3, t_end = 12, temp = 298.15, batch = "A",
                         stringsAsFactors = FALSE)),
    "begin at time 0"
  )
})

test_that("the origin check runs after sorting, not on input row order", {
  # Rows may be supplied in any order, so a schedule that is anchored at 0 once
  # sorted must be accepted -- otherwise the check would reject valid input.
  fit <- origin_fit()
  sch <- data.frame(t_start = c(0, 6), t_end = c(6, 12),
                    temp = c(298.15, 313.15), batch = c("A", "A"),
                    stringsAsFactors = FALSE)
  expect_equal(ptos(fit, sch[c(2, 1), ])$fit, ptos(fit, sch)$fit)
})

test_that("anchored schedules are unaffected", {
  fit <- origin_fit()

  # Single segment.
  expect_equal(
    ptos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15, batch = "A",
                         stringsAsFactors = FALSE))$fit,
    oracle(fit, batch = "A", time = 12, tempk = 298.15), tolerance = 1e-8)

  # Zero-duration leading segment: t_start[1] is still 0.
  expect_no_error(
    ptos(fit, data.frame(t_start = c(0, 0, 6), t_end = c(0, 6, 12),
                         temp = c(298.15, 298.15, 313.15), batch = "A",
                         stringsAsFactors = FALSE)))

  # Fractional boundaries must not trip the tolerance.
  expect_no_error(
    ptos(fit, data.frame(t_start = c(0, 0.1, 0.3), t_end = c(0.1, 0.3, 0.7),
                         temp = c(298.15, 313.15, 298.15), batch = "A",
                         stringsAsFactors = FALSE)))

  # A t_start of floating-point dust is 0 for this purpose.
  expect_no_error(
    ptos(fit, data.frame(t_start = 1e-12, t_end = 12, temp = 298.15,
                         batch = "A", stringsAsFactors = FALSE)))
})

test_that("the returned time column describes the timeline actually used", {
  # The original complaint was that the output asserted a timeline it ignored.
  # Now that the timeline must start at 0, the labels are truthful: every row's
  # fit equals the fitted value at that row's absolute time -- the time-0
  # baseline row included, where that value is the fitted intercept.
  fit <- origin_fit()
  sch <- data.frame(t_start = c(0, 6, 12), t_end = c(6, 12, 18),
                    temp = 298.15, batch = "A", stringsAsFactors = FALSE)
  p <- suppressWarnings(predict_tos(fit, sch, temp_unit = "K"))

  # Baseline row, then one row per segment labelled by that segment's end.
  expect_equal(p$segment, c(0L, 1L, 2L, 3L))
  expect_equal(p$time, c(0, sch$t_end))
  for (i in seq_len(nrow(p))) {
    expect_equal(p$fit[i],
                 oracle(fit, batch = "A", time = p$time[i], tempk = 298.15),
                 tolerance = 1e-8, info = paste("row", i))
  }
})
