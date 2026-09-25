# fit_sb integrates dalpha/dt from alpha = 1e-6, and deSolve::ode() places that
# initial condition at times[1]. predict_run() prepended t = 0 only when
# min(time) > 0 -- a strict test -- so a negative time left the initial condition
# at min(time), and y0 quietly stopped meaning Y(0).
#
# The prediction therefore stopped being a function of (group, temp, time):
# newdata time = c(0, 6, 12) gave 100.000, 56.9459, 32.4284 (matching the closed
# form to 6e-06), while c(-1, 0, 6, 12) gave 100.000, 91.0425, 51.8448, 29.5236 --
# the SAME three triples displaced by 8.96%, now off the closed form by 8.96%.
# At fit time the intercept absorbed the offset instead: data whose truth is
# y0 = 100 at t = 0 with one legitimate t = -3 reading returned y0_0 = 143.17,
# k_0 = 0.0525 (truth 0.02) and Ea = 37389 (truth 80000) with deviance 6711,
# reporting convergence while fitting neither the truth nor the data.
#
# Integrating backward is not available as a fix. Prepending 0 cannot work
# because sort(c(0, times)) returns the negative value to the front and the
# initial condition follows times[1], and deSolve cannot take a non-monotonic
# time vector. Nor is alpha < 0 representable: it is the extent of degradation on
# [0, 1] with alpha(0) ~ 0, growth_ode() clamps the state into that interval and
# residuals_fn() rejects alpha < -1e-6. So negative times are rejected, and
# re-referencing the axis -- which redefines y0 -- is left to the caller.

# n = 1, m = 0, Yinf = 0 gives dalpha/dt = k_T (1 - alpha), i.e. the closed form
# Y = y0 * exp(-k_T * t), so the ODE solution has an exact reference.
sb_kT <- function(tempk, k = 0.02, Ea = 80000) {
  k * exp(Ea / 8.314 * (1 / 298.15 - 1 / tempk))
}
sb_closed <- function(time, tempk, y0 = 100) y0 * exp(-sb_kT(tempk) * time)

sb_synth <- function(times = c(0, 3, 6, 12), temps_c = c(25, 40, 50)) {
  d <- expand.grid(batch = "B1", time = times, temp_c = temps_c,
                   stringsAsFactors = FALSE)
  d$tempk <- d$temp_c + 273.15
  d$y <- sb_closed(d$time, d$tempk)
  d
}

sb_fit <- function(times = c(0, 3, 6, 12)) {
  suppressWarnings(fit_sb(
    sb_synth(times = times), trend = "decreasing", Yinf = 0,
    start = c(Ea = 80000, n = 1, m = 0, k_0 = 0.02, y0_0 = 100),
    lower = c(Ea = 0, n = 1, m = 0, k_0 = 0, y0_0 = 0),
    upper = c(Ea = 1e6, n = 1, m = 0, k_0 = 100, y0_0 = 150)))
}

test_that("fit_sb rejects negative times in the fitting data", {
  # Previously fitted with y0_0 = 143.17 / k_0 = 0.0525 / Ea = 37389 and reported
  # convergence.
  err <- expect_error(sb_fit(times = c(-3, 0, 3, 6, 12)))
  expect_match(conditionMessage(err), "non-negative times")
  # Must name the column.
  expect_match(conditionMessage(err), "'time'")
})

test_that("predict.fit_sb rejects a negative time in newdata", {
  fit <- sb_fit()
  err <- expect_error(
    predict(fit, newdata = data.frame(time = c(-1, 0, 6, 12), tempk = 313.15))
  )
  expect_match(conditionMessage(err), "non-negative times")
  expect_match(conditionMessage(err), "newdata")
})

test_that("a prediction is a pure function of (group, temp, time)", {
  # The property the bug broke: adding a row must never move the other rows.
  fit <- sb_fit()

  base <- predict(fit, newdata = data.frame(time = c(0, 6, 12), tempk = 313.15))

  # An extra positive time interpolated in, and one extrapolated past the end.
  more <- predict(fit, newdata = data.frame(time = c(0, 3, 6, 12, 24),
                                            tempk = 313.15))
  expect_equal(as.numeric(more[c(1, 3, 4)]), as.numeric(base), tolerance = 1e-8)

  # An earlier POSITIVE time, the nearest legal analogue of the reported case:
  # t = 0 is no longer the minimum requested time but is still the origin.
  later <- predict(fit, newdata = data.frame(time = c(3, 6, 12), tempk = 313.15))
  expect_equal(as.numeric(later[c(2, 3)]), as.numeric(base[c(2, 3)]),
               tolerance = 1e-8)

  # Row order must not matter either.
  shuf <- predict(fit, newdata = data.frame(time = c(12, 0, 6), tempk = 313.15))
  expect_equal(as.numeric(shuf), as.numeric(base[c(3, 1, 2)]), tolerance = 1e-8)
})

test_that("predictions match the closed form, and the 8.96% shift is gone", {
  fit <- sb_fit()
  p  <- predict(fit, newdata = data.frame(time = c(0, 6, 12), tempk = 313.15))
  cf <- sb_closed(c(0, 6, 12), 313.15)

  expect_equal(as.numeric(p), cf, tolerance = 1e-5)
  # The values the report quotes as correct.
  expect_equal(as.numeric(p), c(100.0000, 56.9459, 32.4284), tolerance = 1e-4)
  # The displaced values must not be reachable: 91.0425 was y(6) reported at
  # t = 0's position after the origin moved to t = -1.
  expect_gt(abs(as.numeric(p)[2] - 91.0425), 30)
})

test_that("t = 0 is still prepended when every time is strictly positive", {
  # The prepend branch is the reason predictions are origin-anchored at all; it
  # must keep working when 0 is absent from the requested times.
  fit <- sb_fit()
  p <- predict(fit, newdata = data.frame(time = c(6, 12), tempk = 313.15))
  expect_equal(as.numeric(p), sb_closed(c(6, 12), 313.15), tolerance = 1e-5)
})

test_that("fit_sb still accepts data whose earliest time is 0 or positive", {
  # No false positives: the check must not reject ordinary designs.
  expect_s3_class(sb_fit(times = c(0, 3, 6, 12)), "fit_sb")
  expect_s3_class(sb_fit(times = c(3, 6, 12)), "fit_sb")
})

test_that("the closed-form fitters keep accepting negative times", {
  # They evaluate an algebraic expression in t, so y0 remains Y(0) and a negative
  # time neither misbehaves nor perturbs other rows. The fix must stay scoped to
  # the ODE path rather than becoming a package-wide restriction.
  g <- expand.grid(batch = c("A", "B"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- c(A = 100, B = 97)[g$batch] * exp(-k * g$time)
  fk <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                                     intercept_cols = "batch"))

  nd_pos <- data.frame(batch = "A", time = c(0, 6, 12), tempk = 313.15)
  nd_neg <- data.frame(batch = "A", time = c(-1, 0, 6, 12), tempk = 313.15)
  a <- suppressWarnings(predict(fk, newdata = nd_pos))
  b <- suppressWarnings(predict(fk, newdata = nd_neg))

  expect_equal(as.numeric(b)[2:4], as.numeric(a), tolerance = 1e-10)
  expect_gt(as.numeric(b)[1], as.numeric(b)[2])   # earlier time, less degraded
})

# ── the validator in isolation ────────────────────────────────────────────────

test_that("check_nonnegative_times reports the column and count", {
  df <- data.frame(time = c(0, 3, -2, 6, -5), y = 1)
  err <- expect_error(check_nonnegative_times(df, "time"))
  expect_match(conditionMessage(err), "2 negative values")
  expect_match(conditionMessage(err), "'time'")
})

test_that("check_nonnegative_times is silent on valid or absent input", {
  expect_silent(check_nonnegative_times(data.frame(time = c(0, 3, 6)), "time"))
  expect_silent(check_nonnegative_times(data.frame(time = c(3, 6)), "time"))
  # NA is the business of check_structural_nas(), not this check.
  expect_silent(check_nonnegative_times(data.frame(time = c(NA, 3)), "time"))
  # A missing or non-numeric column is reported by the checks that own it.
  expect_silent(check_nonnegative_times(data.frame(y = 1), "time"))
  expect_silent(check_nonnegative_times(data.frame(time = "a"), "time"))
  expect_silent(check_nonnegative_times(data.frame(time = 1), character()))
})

test_that("a non-default time_var is honoured", {
  # The message must name the caller's column, not a hardcoded "time".
  d <- sb_synth()
  names(d)[names(d) == "time"] <- "months"
  d$months[1] <- -1
  err <- expect_error(suppressWarnings(fit_sb(
    d, trend = "decreasing", Yinf = 0, time_var = "months",
    start = c(Ea = 80000, n = 1, m = 0, k_0 = 0.02, y0_0 = 100))))
  expect_match(conditionMessage(err), "'months'")
})
