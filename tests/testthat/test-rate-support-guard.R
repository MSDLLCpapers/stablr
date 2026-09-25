# The rate-support guard -- unsupported_rate_cells() recording at fit time and
# check_rate_support_newdata() / check_observed_cells() enforcing at predict time
# -- had five ways to be got past or to stop a legitimate prediction. They share
# one root cause: the guard was written as a string comparison against a display
# label, and its "one timepoint" tier asked how the cell was sampled without
# asking whether k was actually unidentified.

# B1/B2 fully sampled; B3's sampling is the variable under test.
guard_data <- function(b3_times = 0, b3_temps = c(298.15, 313.15)) {
  g <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 313.15),
                   batch = c("B1", "B2"), stringsAsFactors = FALSE)
  g$y <- 100 * exp(-0.01 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15)) * g$time)
  b3 <- expand.grid(time = b3_times, tempk = b3_temps, batch = "B3",
                    stringsAsFactors = FALSE)
  b3$y <- 90 * exp(-0.01 * exp(-80000 / 8.314 * (1 / b3$tempk - 1 / 298.15)) * b3$time)
  rbind(g, b3)
}

guard_fit <- function(d = guard_data(), ...) {
  suppressWarnings(fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                               rate_cols = "batch", intercept_cols = "batch", ...))
}

# ── R03-01: the guard must key on the cell, not on its printed label ──────────

test_that("the guard survives a change of options(useFancyQuotes)", {
  # The recorded label is built with sQuote(), so a fit made with the default
  # fancy quotes recorded "(batch=U+2018C U+2019)" while predict() rebuilt
  # "(batch='C')". match() found nothing, the guard returned invisible(NULL), and
  # the flat forecast went out unchallenged. Nothing in the output says which
  # quoting style produced it, so the same script gave different safety on
  # different machines.
  fit <- guard_fit()
  nd  <- data.frame(batch = "B3", time = 24, tempk = 298.15)

  old <- options(useFancyQuotes = TRUE)
  on.exit(options(old), add = TRUE)
  for (q in list(TRUE, FALSE, "TeX")) {
    options(useFancyQuotes = q)
    expect_error(predict(fit, newdata = nd), "never estimated",
                 info = paste("useFancyQuotes =", format(q)))
  }
})

test_that("the recorded cell carries a quote-free matching key", {
  fit <- guard_fit()
  u <- attr(fit, "build")$settings$unsupported_rate_cells
  expect_true("key" %in% names(u))
  expect_equal(u$key, "B3")                       # raw value, no quoting
  expect_match(u$label, "B3")                     # label still there, for humans
})

test_that("a fit recorded before `key` existed still gets checked", {
  # Old saved fit objects hold only label/reason, so the label comparison has to
  # stay as a fallback rather than the guard silently opting out.
  legacy <- data.frame(label = cell_labels("B3", "batch"),
                       reason = "no elapsed time", stringsAsFactors = FALSE)
  expect_error(
    check_rate_support_newdata(data.frame(batch = "B3", time = 24), legacy,
                               "batch", "time"),
    "never estimated"
  )
})

# ── R03-04: the t == 0 exemption has to be two-sided ─────────────────────────

test_that("backward extrapolation off an unestimated rate is refused", {
  # The exemption was `tm > 0`, which reads as "only forward predictions use the
  # rate". exp(-k * t) is just as wrong at t = -6 as at t = +6, and a negative
  # time is a back-calculation to a manufacturing date -- a real use.
  fit <- guard_fit()
  expect_error(
    predict(fit, newdata = data.frame(batch = "B3", time = -6, tempk = 298.15)),
    "never estimated"
  )
})

test_that("t == 0 itself stays exempt", {
  # The rate is multiplied by t, so at t = 0 the answer is the fitted intercept
  # and no unidentified parameter enters.
  fit <- guard_fit()
  p <- predict(fit, newdata = data.frame(batch = "B3", time = 0, tempk = 298.15))
  expect_equal(unname(p), unname(coef(fit)["y0_batch_B3"]), tolerance = 1e-6)
})

# ── R03-06: "one timepoint" only matters when k is actually unidentified ──────

test_that("a pooled intercept identifies a single-timepoint rate", {
  # With no intercept_cols there is one y0 that B1/B2's timepoints pin down, so
  # B3's single observation is an equation in k alone. Flagging it made predict()
  # hard-refuse a rate the fitter estimates to the same relative precision as the
  # well-supported cells.
  d <- guard_data(6, b3_temps = 298.15)
  expect_equal(nrow(unsupported_rate_cells(d, "batch", "time",
                                           intercept_cols = character(),
                                           temp_var = "tempk")), 0L)
  # ... and the same design IS flagged once the intercept is exclusive to B3.
  expect_equal(nrow(unsupported_rate_cells(d, "batch", "time",
                                           intercept_cols = "batch",
                                           temp_var = "tempk")), 1L)
})

test_that("two temperatures at one timepoint identify the rate", {
  d <- guard_data(6, b3_temps = c(298.15, 313.15))
  expect_equal(nrow(unsupported_rate_cells(d, "batch", "time",
                                           intercept_cols = "batch",
                                           temp_var = "tempk")), 0L)
  # Without temp_var the two-temperature test cannot run, so the cell is treated
  # as single-temperature and flagged -- the conservative direction.
  expect_equal(nrow(unsupported_rate_cells(d, "batch", "time",
                                           intercept_cols = "batch")), 1L)
})

test_that("tier one is unconditional whatever the pooling or temperature spread", {
  # exp(-k * 0) == 1 for every k, so a cell with no elapsed time carries zero
  # information and no amount of pooling recovers it.
  d <- guard_data(0)                       # B3 at t = 0, two temperatures
  for (ic in list(character(), "batch")) {
    bad <- unsupported_rate_cells(d, "batch", "time", intercept_cols = ic,
                                  temp_var = "tempk")
    expect_equal(nrow(bad), 1L)
    expect_match(bad$reason, "no rows with elapsed time")
    expect_true(is.na(bad$obs_time))
  }
})

# ── W11-05: reproducing the fit's own rows is not extrapolation ───────────────

test_that("predict() on the fit's own training frame does not error", {
  # predict(fit) returned these rows without complaint while
  # predict(fit, newdata = <the same rows>) errored, so any workflow that passes
  # the model frame back explicitly -- diagnostics, residual plots, a refit
  # helper -- broke on a fit that was otherwise usable.
  d   <- guard_data(6, b3_temps = 298.15)   # B3 flagged: one time, one temp
  fit <- guard_fit(d)
  expect_equal(nrow(attr(fit, "build")$settings$unsupported_rate_cells), 1L)

  expect_silent(p_new <- predict(fit, newdata = d))
  expect_equal(unname(p_new), unname(predict(fit)), tolerance = 1e-8)

  # The exemption is only for the observed time. One month later is a forecast.
  expect_error(
    predict(fit, newdata = data.frame(batch = "B3", time = 7, tempk = 298.15)),
    "never estimated"
  )
})

test_that("the observed time is recorded so predict can recognise it", {
  d <- guard_data(6, b3_temps = 298.15)
  bad <- unsupported_rate_cells(d, "batch", "time", intercept_cols = "batch",
                                temp_var = "tempk")
  expect_equal(bad$obs_time, 6)
  expect_match(bad$reason, "only one distinct timepoint")
})

# ── R03-02: the single-global-rate branch has to be reachable ─────────────────

test_that("an unidentified GLOBAL rate is refused at predict time", {
  # With no rate_cols there are no grouping columns, so sanitize_map is empty and
  # apply_sanitize_map() returned before either guard ran. unsupported_rate_cells'
  # own "the single global rate" branch was therefore dead code through every
  # fitter: a fit whose every row sits at time 0 forecast two years of invented
  # degradation with nothing raised.
  d <- data.frame(time = 0, tempk = c(298.15, 298.15, 313.15, 313.15),
                  y = c(100, 99.8, 100, 99.5))
  fit <- suppressWarnings(fit_kinetic(d, kinetic_order = 1, trend = "decreasing"))
  expect_equal(length(attr(fit, "sanitize_map")), 0L)
  expect_match(attr(fit, "build")$settings$unsupported_rate_cells$label,
               "global rate")

  expect_error(predict(fit, newdata = data.frame(time = 24, tempk = 298.15)),
               "never estimated")
  # t == 0 still answers from the intercept.
  expect_silent(predict(fit, newdata = data.frame(time = 0, tempk = 298.15)))
})

# ── R03-07: an unobserved CELL at t == 0 is answered by the intercept ─────────

test_that("check_observed_cells exempts a t == 0 release prediction", {
  # (B,y) is never observed, so it has no rate parameter -- but the intercept is
  # keyed on batch alone, so y0_batch_B is fitted and is the whole answer at
  # t = 0. Refusing it blocked a release-time prediction while the sibling guard
  # allowed the identical situation, and the two messages disagreed about whether
  # t = 0 is safe.
  g <- expand.grid(batch = c("A", "B"), pkg = c("x", "y"), time = c(0, 3, 6, 12),
                   tempk = c(298.15, 313.15), stringsAsFactors = FALSE)
  g <- g[!(g$batch == "B" & g$pkg == "y"), ]
  g$y <- c(A = 100, B = 98)[g$batch] *
    exp(-0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15)) * g$time)
  fit <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                                      rate_cols = c("batch", "pkg"),
                                      intercept_cols = "batch"))

  p0 <- predict(fit, newdata = data.frame(batch = "B", pkg = "y", time = 0,
                                          tempk = 298.15))
  expect_equal(unname(p0), unname(coef(fit)["y0_batch_B"]), tolerance = 1e-6)

  # Any nonzero time still needs the rate the cell never had.
  for (tt in c(12, -12)) {
    expect_error(
      predict(fit, newdata = data.frame(batch = "B", pkg = "y", time = tt,
                                        tempk = 298.15)),
      "never observed"
    )
  }
})

test_that("check_observed_cells with no time_var exempts nothing", {
  # predict_tos() passes a schedule, which has no single time column and whose
  # segments always advance by a positive duration, so there is nothing to exempt.
  cells <- data.frame(batch = "A", pkg = "x", stringsAsFactors = FALSE)
  nd    <- data.frame(batch = "B", pkg = "y", time = 0, stringsAsFactors = FALSE)
  expect_error(check_observed_cells(nd, cells, c("batch", "pkg")),
               "never observed")
  expect_silent(check_observed_cells(nd, cells, c("batch", "pkg"),
                                     time_var = "time"))
})

# ── the fit-time warning and the predict-time refusal must agree ──────────────

test_that("a cell that warns at fit time is the cell refused at predict time", {
  # These were computed by separate calls with different arguments, so it was
  # possible to be warned about a cell that then predicted freely, or refused one
  # that never warned.
  d  <- guard_data(6, b3_temps = 298.15)
  ws <- character()
  fit <- withCallingHandlers(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                rate_cols = "batch", intercept_cols = "batch"),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  warned <- grep("insufficient data", ws, value = TRUE)[1]
  expect_match(warned, "B3")

  err <- expect_error(predict(fit, newdata = data.frame(batch = "B3", time = 24,
                                                        tempk = 298.15)))
  expect_match(conditionMessage(err), "B3")
  # Same reason string in both.
  expect_match(warned, "only one distinct timepoint")
  expect_match(conditionMessage(err), "only one distinct timepoint")
})

test_that("well-supported cells are unaffected by any of this", {
  fit <- guard_fit(guard_data(c(0, 6, 12)))
  expect_equal(nrow(attr(fit, "build")$settings$unsupported_rate_cells), 0L)
  for (b in c("B1", "B2", "B3")) {
    for (tt in c(-6, 0, 24)) {
      expect_silent(predict(fit, newdata = data.frame(batch = b, time = tt,
                                                      tempk = 313.15)))
    }
  }
})

# ── W09-02: predict.fit_sb must apply the same rate-support guard ─────────────

# fit_sb() takes the ODE path (make_group + interaction cell means), not
# apply_sanitize_map(), so it was the one fitter whose predict() skipped
# check_rate_support_newdata() entirely. Its k seed of 1/max(time) is a large
# positive rate, so the skipped guard produced a confident invented loss rather
# than the harmless flat line a zero seed would give.
sb_undersupported_data <- function() {
  g <- expand.grid(batch = c("B1", "B2"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40, 50), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  g$y <- 100 * exp(-0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15)) * g$time)
  g[!(g$batch == "B2" & g$time > 0), ]            # B2 observed only at time 0
}

test_that("fit_sb records the under-supported cell in its settings", {
  d <- sb_undersupported_data()
  fit <- suppressWarnings(fit_sb(d, intercept_cols = "batch", rate_cols = "batch",
                                 trend = "decreasing"))
  u <- attr(fit, "build")$settings$unsupported_rate_cells
  expect_equal(u$key, "B2")
  expect_match(u$reason, "no rows with elapsed time > 0")
})

test_that("predict.fit_sb errors on the under-supported extrapolation", {
  d  <- sb_undersupported_data()
  fit <- suppressWarnings(fit_sb(d, intercept_cols = "batch", rate_cols = "batch",
                                 trend = "decreasing"))
  nd <- data.frame(batch = "B2", time = c(0, 24), tempk = 298.15)

  # The forecast fit_sb used to return was its 1/max(time) seed, not an estimate.
  expect_equal(unname(fit$par["k_B2"]), 1 / max(d$time), tolerance = 1e-8)
  expect_error(predict(fit, newdata = nd), "never estimated from data")

  # Same refusal fit_kinetic gives on the identical frame -- the point of parity.
  fk <- suppressWarnings(fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                                     intercept_cols = "batch", rate_cols = "batch"))
  expect_error(predict(fk, newdata = nd), "never estimated from data")
})

test_that("predict.fit_sb still allows t == 0 and well-supported cells", {
  d  <- sb_undersupported_data()
  fit <- suppressWarnings(fit_sb(d, intercept_cols = "batch", rate_cols = "batch",
                                 trend = "decreasing"))

  # t == 0 release of the under-supported cell: the rate never enters.
  expect_silent(predict(fit, newdata = data.frame(batch = "B2", time = 0,
                                                  tempk = 298.15)))
  # A fully sampled cell predicts across time as before.
  expect_silent(predict(fit, newdata = data.frame(batch = "B1", time = c(0, 24),
                                                  tempk = 298.15)))
  # Predicting the training frame (newdata = NULL) is unaffected.
  expect_silent(predict(fit))
})
