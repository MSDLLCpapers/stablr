# The humidity term is exp(b * (rh - rh_ref)), so coding rh as a percent
# (40/60/75) while rh_ref stays a fraction (0.60) is an EXACT
# reparameterisation: b absorbs the factor of 100 and k_0 becomes the rate at
# 0.6% RH instead of 60% RH. sigma, residuals and AIC are bit-identical to the
# correctly scaled fit, so no goodness-of-fit diagnostic can detect it. The
# mismatch must be caught at input time.

make_rh_curve <- function(rh, divisor, Ea = 80000, R = 8.314,
                          k0 = 0.01, b = 2, y0 = 100) {
  d <- expand.grid(time = c(0, 3, 6, 12),
                   tempk = c(298.15, 313.15),
                   rh = rh)
  d$y <- y0 * exp(-k0 *
                    exp(-Ea / R * (1 / d$tempk - 1 / 298.15)) *
                    exp(b * (d$rh / divisor - 0.6)) * d$time)
  d
}

d_pct  <- make_rh_curve(c(40, 60, 75), 100)      # percent-coded
d_frac <- make_rh_curve(c(0.4, 0.6, 0.75), 1)    # fraction-coded

# Collect every warning a fit emits. A deliberately mis-scaled humidity fit does
# not converge (rh_ref = 60 against fraction data drives k_0 to ~5e4), so the
# convergence check fires alongside the scale warning under test; asserting on
# the collected set keeps the two from masking each other.
warnings_from <- function(expr) {
  ws <- character()
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

expect_scale_warning <- function(expr, pattern) {
  expect_gte(length(grep(pattern, warnings_from(expr))), 1)
}

test_that("percent-coded rh with a fractional rh_ref warns", {
  expect_warning(
    fit_kinetic(d_pct, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = expr_arrhenius_humidity()),
    "percent-coded humidity"
  )
})

test_that("the warning names the column, its range and the reference", {
  wmsg <- tryCatch(
    fit_kinetic(d_pct, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = expr_arrhenius_humidity()),
    warning = conditionMessage
  )
  expect_match(wmsg, "`rh`")
  expect_match(wmsg, "40")
  expect_match(wmsg, "75")
  expect_match(wmsg, "0\\.6")
  # It must tell the user how to fix it.
  expect_match(wmsg, "rh_ref = 60")
})

test_that("the silent misparameterisation is exactly what the warning claims", {
  # Documents the pathology the warning exists to catch: an exact
  # re-referencing, not a lack of fit.
  fit <- suppressWarnings(
    fit_kinetic(d_pct, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = expr_arrhenius_humidity())
  )
  cf <- coef(fit)
  # Analytic oracle: b -> b/100, k_0 -> k_0 * exp(-b * (0.6 - 0.006)).
  expect_equal(unname(cf[["b"]]), 2 / 100, tolerance = 1e-6)
  expect_equal(unname(cf[["k_0"]]), 0.01 * exp(-1.188), tolerance = 1e-6)
  # The fit looks flawless -- which is exactly why a range check is needed.
  expect_lt(sigma(fit), 1e-6)
})

test_that("fraction-coded rh with a percent rh_ref warns", {
  expect_scale_warning(
    fit_kinetic(d_frac, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = expr_arrhenius_humidity(rh_ref = 60)),
    "fraction-coded humidity"
  )
})

test_that("consistently scaled humidity does not warn", {
  # Fractions with a fractional reference.
  expect_no_warning(
    fit_frac <- fit_kinetic(d_frac, kinetic_order = 1, trend = "decreasing",
                            Yinf = 0,
                            temp_expr = expr_arrhenius_humidity())
  )
  expect_equal(unname(coef(fit_frac)[["k_0"]]), 0.01, tolerance = 1e-6)
  expect_equal(unname(coef(fit_frac)[["b"]]), 2, tolerance = 1e-6)

  # Percents with a percent reference: k_0 keeps its meaning as the rate at
  # the reference humidity, and b is simply per-percent.
  expect_no_warning(
    fit_pct <- fit_kinetic(d_pct, kinetic_order = 1, trend = "decreasing",
                           Yinf = 0,
                           temp_expr = expr_arrhenius_humidity(rh_ref = 60))
  )
  expect_equal(unname(coef(fit_pct)[["k_0"]]), 0.01, tolerance = 1e-6)
  expect_equal(unname(coef(fit_pct)[["b"]]), 2 / 100, tolerance = 1e-6)
})

test_that("models without a humidity term are unaffected", {
  expect_no_warning(
    fit_kinetic(d_frac, kinetic_order = 1, trend = "decreasing", Yinf = 0)
  )
})

test_that("check_rh_scale is a no-op on absent or unusable input", {
  expect_silent(check_rh_scale(d_frac, NULL, 0.6))
  expect_silent(check_rh_scale(d_frac, "rh", NULL))
  expect_silent(check_rh_scale(d_frac, "not_a_column", 0.6))
  expect_silent(check_rh_scale(data.frame(rh = NA_real_), "rh", 0.6))
  expect_silent(check_rh_scale(data.frame(rh = numeric(0)), "rh", 0.6))
  # Boundary: rh exactly at 1 is still a plausible fraction.
  expect_silent(check_rh_scale(data.frame(rh = c(0.4, 1)), "rh", 0.6))
})

test_that("fit_singledeg warns on a percent/fraction rh mismatch", {
  d <- d_pct
  d$yA <- d$y
  d$yB <- 100 - d$y

  expect_warning(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  yA_var = "yA", yB_var = "yB",
                  temp_expr = expr_arrhenius_humidity()),
    "percent-coded humidity"
  )
})

test_that("fit_parallel warns on a percent/fraction rh mismatch", {
  d <- d_pct
  d$yA <- d$y
  d$yB <- (100 - d$y) * 0.6
  d$yC <- (100 - d$y) * 0.4

  expect_warning(
    fit_parallel(d, kinetic_order = 1,
                 yA_var = "yA", yB_var = "yB", yC_var = "yC",
                 temp_expr = expr_arrhenius_humidity()),
    "percent-coded humidity"
  )
})
