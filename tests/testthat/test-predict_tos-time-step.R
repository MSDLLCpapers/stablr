# predict_tos(time_step =) predicts within each segment, not only at its end. A
# schedule is a small number of long segments, so the boundary-only output could
# not be plotted as a trajectory: a 12-month storage leg contributed one point.
#
# An interior point at time t is the composite run over [0, t] -- every whole
# segment before t, then the segment t falls in truncated to end at t. The
# equivalence tests below are the ones that matter: a truncated trailing segment
# must give exactly what a schedule genuinely ENDING at t gives, since that is
# the same state reached by a different code path.

step_data <- function(noisy = FALSE) {
  g <- expand.grid(batch = c("A", "B"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- c(A = 100, B = 97)[g$batch] * exp(-k * g$time)
  if (noisy) {
    set.seed(1)
    g$y <- g$y + stats::rnorm(nrow(g), 0, 0.3)
  }
  g$yA <- g$y
  g$yB <- 100 - g$y
  g
}

step_fit <- function(noisy = FALSE) {
  suppressWarnings(fit_kinetic(step_data(noisy), kinetic_order = 1,
                               trend = "decreasing", intercept_cols = "batch"))
}

one_seg  <- data.frame(t_start = 0, t_end = 12, temp = 25, batch = "A",
                       stringsAsFactors = FALSE)
excursion <- data.frame(t_start = c(0, 12, 13), t_end = c(12, 13, 24),
                        temp = c(5, 25, 5), batch = "A",
                        stringsAsFactors = FALSE)

test_that("a step of 1 over a single 0-12 segment gives times 0 through 12", {
  p <- predict_tos(step_fit(), one_seg, time_step = 1)

  expect_equal(p$time, 0:12)
  # The baseline keeps segment 0; every interior point belongs to segment 1.
  expect_equal(p$segment, c(0L, rep(1L, 12)))
  expect_equal(unique(p$temp), 25)
  expect_equal(names(p),
               c("segment", "time", "temp", "variable", "fit", "lwr", "upr"))
})

test_that("time_step is the default-off addition it is meant to be", {
  # Omitting it must reproduce the boundary-only table exactly, including the
  # zero-duration-segment handling, so this is a pure extension.
  fit <- step_fit()
  p <- suppressWarnings(predict_tos(fit, excursion))

  expect_equal(p$time, c(0, 12, 13, 24))
  expect_equal(p$segment, 0:3)
  expect_equal(suppressWarnings(predict_tos(fit, excursion, time_step = NULL)), p)
})

test_that("an interior point equals a schedule that genuinely ends there", {
  # The load-bearing test. Truncating the trailing segment reaches the same state
  # as a schedule whose last segment really stops at t, so the two must agree to
  # the last bit -- fit and interval alike -- at every interior time, in every
  # segment, including across the two temperature changes.
  fit <- step_fit(noisy = TRUE)      # noise, so the interval has real width
  p <- suppressWarnings(predict_tos(fit, excursion, time_step = 1))

  ends_at <- function(t) {
    s <- excursion[excursion$t_start < t - 1e-12, , drop = FALSE]
    s$t_end[nrow(s)] <- t
    q <- suppressWarnings(predict_tos(fit, s))
    q[nrow(q), , drop = FALSE]
  }

  interior <- which(p$time > 0)
  expect_gt(length(interior), 20L)
  for (i in interior) {
    ref <- ends_at(p$time[i])
    expect_equal(p$fit[i], ref$fit, tolerance = 1e-10, info = paste("t =", p$time[i]))
    expect_equal(p$lwr[i], ref$lwr, tolerance = 1e-10, info = paste("t =", p$time[i]))
    expect_equal(p$upr[i], ref$upr, tolerance = 1e-10, info = paste("t =", p$time[i]))
    expect_equal(p$segment[i], ref$segment, info = paste("t =", p$time[i]))
    expect_equal(p$temp[i], ref$temp, info = paste("t =", p$time[i]))
  }
})

test_that("interior points match the fitted model on an isothermal schedule", {
  # An independent oracle is available when the temperature never changes: the
  # composite must agree with the fitted formula evaluated at absolute time t.
  fit <- step_fit()
  p <- predict_tos(fit, one_seg, time_step = 1)

  env <- list2env(as.list(stats::coef(fit)), parent = baseenv())
  oracle <- eval(stats::formula(fit)[[3]],
                 envir = data.frame(batch = "A", time = p$time, tempk = 298.15),
                 enclos = env)
  expect_equal(p$fit, as.numeric(oracle), tolerance = 1e-10)
})

test_that("segment boundaries survive a step that does not divide evenly", {
  # Dropping t_end would lose the excursion's start and end -- the times the
  # schedule was written to describe.
  p <- predict_tos(step_fit(), one_seg, time_step = 5)
  expect_equal(p$time, c(0, 5, 10, 12))

  # Every boundary of a multi-segment schedule is present, and a step wider than
  # a segment leaves that segment with its boundary alone.
  p2 <- suppressWarnings(predict_tos(step_fit(), excursion, time_step = 5))
  expect_true(all(c(0, 12, 13, 24) %in% p2$time))
  expect_equal(sum(p2$time == 13), 1L)
})

test_that("a boundary time is reported once, under the segment ending there", {
  # t = 12 ends segment 1 (at 5 C) and starts segment 2 (at 25 C). Attributing it
  # to both would duplicate the row and disagree about temp.
  p <- suppressWarnings(predict_tos(step_fit(), excursion, time_step = 1))

  expect_false(any(duplicated(p$time)))
  expect_equal(p$segment[p$time == 12], 1L)
  expect_equal(p$temp[p$time == 12], 5)
  expect_equal(p$segment[p$time == 13], 2L)
  expect_equal(p$temp[p$time == 13], 25)
  # Interior points of the third segment carry its own temperature.
  expect_equal(unique(p$temp[p$time > 13]), 5)
})

test_that("the returned table is ordered by time and monotone in fit", {
  p <- suppressWarnings(predict_tos(step_fit(), excursion, time_step = 1))
  expect_equal(p$time, sort(p$time))
  expect_true(all(diff(p$fit) < 0))
  expect_equal(rownames(p), as.character(seq_len(nrow(p))))
})

test_that("fractional steps and fractional boundaries do not drift", {
  # Floating-point seq() output must not produce a near-duplicate of t_end nor
  # spill a point into the neighbouring segment.
  fit <- step_fit()

  expect_equal(predict_tos(fit, data.frame(t_start = 0, t_end = 1, temp = 25,
                                           batch = "A"), time_step = 0.25)$time,
               c(0, 0.25, 0.5, 0.75, 1))

  sf <- data.frame(t_start = c(0, 0.1, 0.3), t_end = c(0.1, 0.3, 0.7),
                   temp = c(25, 40, 25), batch = "A", stringsAsFactors = FALSE)
  p <- predict_tos(fit, sf, time_step = 0.1)
  expect_equal(p$time, c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7))
  expect_equal(p$segment, c(0L, 1L, 2L, 2L, 3L, 3L, 3L, 3L))
})

test_that("a zero-duration segment keeps exactly one row", {
  # seq() over a zero-length segment yields only t_start, which belongs to the
  # previous row, so the segment would otherwise vanish from the output.
  sz <- data.frame(t_start = c(0, 0, 6), t_end = c(0, 6, 12),
                   temp = c(25, 25, 40), batch = "A", stringsAsFactors = FALSE)
  p <- suppressWarnings(predict_tos(step_fit(), sz, time_step = 3))

  expect_equal(p$segment, c(0L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(p$time, c(0, 0, 3, 6, 9, 12))
})

test_that("time_step must be a single positive number", {
  fit <- step_fit()
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "1", numeric(0))) {
    expect_error(predict_tos(fit, one_seg, time_step = bad),
                 "time_step must be a single positive number")
  }
})

test_that("an unfinishable grid is refused rather than attempted", {
  # Each point is a fresh composite build plus a Jacobian, so a mistyped step is
  # not a slow call but one that never returns.
  err <- expect_error(predict_tos(step_fit(), one_seg, time_step = 0.001))
  expect_match(conditionMessage(err), "12001 prediction times")
  expect_match(conditionMessage(err), "limit of 5000")
  # The message must name a step that works, not just refuse.
  expect_match(conditionMessage(err), "coarser time_step")
  expect_no_error(predict_tos(step_fit(), one_seg, time_step = 12 / 4999))
})

test_that("time_step works for the multi-response fitters", {
  fit <- suppressWarnings(fit_singledeg(step_data(), kinetic_order = 1,
                                        trend = "decreasing", Yinf = 0,
                                        intercept_cols = "batch"))
  pA <- predict_tos(fit, one_seg, var = "yA", time_step = 2)
  pB <- predict_tos(fit, one_seg, var = "yB", time_step = 2)

  expect_equal(pA$time, seq(0, 12, by = 2))
  expect_equal(pB$time, seq(0, 12, by = 2))
  expect_equal(unique(pA$variable), "yA")
  # A -> B conserves mass at every interior point, not just at the boundaries.
  expect_equal(pA$fit + pB$fit, rep(pA$fit[1] + pB$fit[1], nrow(pA)),
               tolerance = 1e-8)
})

test_that("time_step works for the fit_sb ODE path", {
  fit <- fit_sb(step_data(), rate_cols = "batch", intercept_cols = "batch",
                trend = "decreasing", Yinf = 0)
  p <- suppressWarnings(predict_tos(fit, excursion, time_step = 3))

  expect_equal(p$time, c(0, 3, 6, 9, 12, 13, 16, 19, 22, 24))
  expect_equal(p$segment, c(0L, 1L, 1L, 1L, 1L, 2L, 3L, 3L, 3L, 3L))
  expect_true(all(is.finite(p$fit)))
  expect_true(all(diff(p$fit) < 0))
  # The interior points must come off the sequential integration, so each must
  # match a schedule that genuinely ends there.
  s18 <- data.frame(t_start = c(0, 12, 13), t_end = c(12, 13, 19),
                    temp = c(5, 25, 5), batch = "A", stringsAsFactors = FALSE)
  ref <- suppressWarnings(predict_tos(fit, s18))
  expect_equal(p$fit[p$time == 19], ref$fit[nrow(ref)], tolerance = 1e-8)
})

test_that("attributes are unaffected by time_step", {
  # attr(, 'schedule') describes the schedule, not the output grid.
  fit <- step_fit()
  p <- suppressWarnings(predict_tos(fit, excursion, time_step = 1))
  s <- attr(p, "schedule")

  expect_equal(nrow(s), 3L)
  expect_equal(s$t_end, c(12, 13, 24))
  expect_equal(attr(p, "level"), 0.95)
  expect_equal(attr(p, "model_type"), "fit_kinetic")
})
