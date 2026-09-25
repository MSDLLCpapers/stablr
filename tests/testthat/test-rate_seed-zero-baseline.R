# Regression tests for the rate_seed() zero-baseline defects.
#
# rate_seed() normalised the seed by y0 and clamped the denominator at 1e-8. When
# the CQA is ~0 at t = 0 the seed exploded to ~1e6-1e8, and gsl_nls returned that
# untouched seed as the "estimate" (isConv = FALSE, k == seed). The fix scales by
# the response range with a lower bound of 1 and bounds the seed to [0.01, 10].

R <- 8.314

test_that("rate_seed() stays bounded when the response is ~0 at release", {
  # Zero baseline: the old 1e-8 clamp sent this to ~1e8; it must now fall in the
  # documented [0.01, 10] band.
  g <- expand.grid(time = c(0, 6, 12, 24), tempk = c(298.15, 313.15),
                   stringsAsFactors = FALSE)
  g$y <- 100 * (1 - exp(-0.02 * g$time))
  g$y[g$time == 0] <- 0

  seed <- rate_seed(g, "y", "time")
  expect_gte(seed, 0.01)
  expect_lte(seed, 10)
})

test_that("rate_seed() falls back to the floor on a flat or degenerate series", {
  flat <- data.frame(time = c(0, 6, 12), y = c(0, 0, 0))
  expect_equal(rate_seed(flat, "y", "time"), 0.01)

  one_time <- data.frame(time = c(3, 3, 3), y = c(1, 2, 3))
  expect_equal(rate_seed(one_time, "y", "time"), 0.01)
})

test_that("A growth CQA with y(0) = 0 recovers truth", {
  # First-order approach to Yinf = 100 from 0; truth Ea = 80000, k_ref = 0.02.
  kref <- 0.02; Ea <- 80000; trefk <- 298.15
  kT <- function(Temp) kref * exp(-Ea / R * (1 / Temp - 1 / trefk))
  g <- expand.grid(time = c(0, 3, 6, 9, 12, 18, 24),
                   tempk = c(298.15, 308.15, 313.15),
                   stringsAsFactors = FALSE)
  g$y <- 100 * (1 - exp(-kT(g$tempk) * g$time))   # y(0) = 0 exactly

  fit <- suppressWarnings(
    fit_kinetic(g, kinetic_order = 1, trend = "increasing", Yinf = 100))
  co  <- coef(fit)

  expect_true(fit$convInfo$isConv)
  expect_equal(unname(co["Ea"]),  80000, tolerance = 1e-3)
  expect_equal(unname(co["k_0"]), 0.02,  tolerance = 1e-3)
  expect_lt(sum(residuals(fit)^2), 1e-6)
})

test_that("A impurity that is 0.0 at release recovers truth", {
  # Per-batch intercepts, zero release value; truth Ea = 90000, k_ref = 0.015.
  kref <- 0.015; Ea <- 90000; trefk <- 298.15
  kT <- function(Temp) kref * exp(-Ea / R * (1 / Temp - 1 / trefk))
  g <- expand.grid(time = c(0, 3, 6, 9, 12, 18, 24),
                   tempk = c(298.15, 308.15, 313.15),
                   batch = c("B1", "B2"),
                   stringsAsFactors = FALSE)
  g$y <- 100 * (1 - exp(-kT(g$tempk) * g$time))
  g$y[g$time == 0] <- 0                            # exact zero at release

  # The seed must not explode through the old 1e-8 clamp.
  expect_lte(rate_seed(g, "y", "time"), 10)

  fit <- suppressWarnings(
    fit_kinetic(g, kinetic_order = 1, trend = "increasing", Yinf = 100,
                intercept_cols = "batch", rate_cols = "batch", trefk = trefk))
  co  <- coef(fit)

  expect_true(fit$convInfo$isConv)
  expect_equal(unname(co["Ea"]), 90000, tolerance = 1e-2)
  expect_equal(unname(co[grep("^k", names(co))][1]), 0.015, tolerance = 1e-2)
  expect_lt(sum(residuals(fit)^2), 1e-6)
})
