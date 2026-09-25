# Regression tests for order-aware rate seeding in default_start().
#
#   X05-02 / X06-02  fit_parallel(kinetic_order = 2)
#   X05-03           fit_kinetic/fit_singledeg(kinetic_order = 2)
#
# rate_seed() returns a FRACTIONAL rate (units 1/time), which equals the rate
# constant k only at first order. For order n the closed form gives
# fractional_rate ~ k*y0^(n-1), so the order-blind seed was ~intercept^(n-1) too
# large and the fit collapsed to negative rate constants and absurd/negative Ea
# on noise-free data. default_start() now divides the seed by intercept^(n-1).

R <- 8.314

test_that("default_start() puts the rate seed on the scale of k for order 2", {
  # Fractional seed ~0.02 with an intercept of 100 must become ~2e-4 at order 2.
  df <- expand.grid(time = c(0, 6, 12, 24), tempk = c(298.15, 313.15),
                    stringsAsFactors = FALSE)
  df$y <- 1 / (1 / 100 + 2e-4 * df$time)   # order-2 decay from y0 = 100

  s1 <- default_start(c("Ea", "y0", "k_0"), df, R, "y", "decreasing", "time",
                      kinetic_order = 1)
  s2 <- default_start(c("Ea", "y0", "k_0"), df, R, "y", "decreasing", "time",
                      kinetic_order = 2)

  # The order-2 seed is the order-1 seed divided by the intercept (~100).
  expect_equal(unname(s2["k_0"]),
               unname(s1["k_0"]) / unname(s1["y0"]), tolerance = 1e-8)
  expect_lt(unname(s2["k_0"]), unname(s1["k_0"]))
})

test_that("X05-02/X06-02: fit_parallel(order 2) recovers truth on noise-free data", {
  trefk <- 278.15
  EaB <- 40000; EaC <- 34000; kB <- 1.2e-4; kC <- 6e-6; A0 <- 98
  gp <- expand.grid(time  = c(0, 3, 6, 9, 12, 18, 24),
                    tempk = c(278.15, 288.15, 298.15, 308.15),
                    stringsAsFactors = FALSE)
  kBapp <- kB * exp(-EaB / R * (1 / gp$tempk - 1 / trefk))
  kCapp <- kC * exp(-EaC / R * (1 / gp$tempk - 1 / trefk))
  ktot  <- kBapp + kCapp
  gp$yA <- 1 / (1 / A0 + ktot * gp$time)
  lost  <- A0 - gp$yA
  gp$yB <- lost * (kBapp / ktot)
  gp$yC <- lost * (kCapp / ktot)

  fit <- suppressWarnings(fit_parallel(gp, kinetic_order = 2, trefk = trefk))
  co  <- coef(fit)

  expect_true(fit$convInfo$isConv)
  expect_equal(unname(co["kB_0"]), kB,   tolerance = 1e-4)
  expect_equal(unname(co["EaB"]),  EaB,  tolerance = 1e-2)
  expect_lt(sum(residuals(fit)^2), 1e-6)
})

test_that("X05-03: fit_kinetic(order 2) recovers truth on noise-free data", {
  Ea <- 40000; kref <- 1.47e-4
  gk <- expand.grid(time  = c(0, 3, 6, 9, 12, 18, 24),
                    tempk = c(298.15, 308.15, 313.15),
                    batch = c("B1", "B2", "B3"),
                    stringsAsFactors = FALSE)
  kapp <- kref * exp(-Ea / R * (1 / gk$tempk - 1 / 298.15))
  gk$y <- 1 / (1 / 100 + kapp * gk$time)   # order-2 decay from y0 = 100

  fit <- suppressWarnings(
    fit_kinetic(gk, kinetic_order = 2, trend = "decreasing", Yinf = 0,
                intercept_cols = "batch"))
  co  <- coef(fit)

  expect_true(fit$convInfo$isConv)
  expect_gt(unname(co["Ea"]), 0)                        # no negative-Ea collapse
  expect_equal(unname(co["Ea"]),  Ea,   tolerance = 1e-2)
  expect_equal(unname(co["k_0"]), kref, tolerance = 1e-2)
  expect_lt(sum(residuals(fit)^2), 1e-6)
})

test_that("order-1 fits are unchanged by the order-aware seed", {
  # kinetic_order = 1 leaves the seed exactly as before (divide by intercept^0).
  g <- expand.grid(time = c(0, 3, 6, 9, 12, 18, 24),
                   tempk = c(298.15, 308.15, 313.15),
                   stringsAsFactors = FALSE)
  g$y <- 100 * exp(-0.02 * exp(-80000 / R * (1 / g$tempk - 1 / 298.15)) * g$time)

  fit <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing"))
  co  <- coef(fit)

  expect_true(fit$convInfo$isConv)
  expect_equal(unname(co["Ea"]),  80000, tolerance = 1e-2)
  expect_equal(unname(co["k_0"]), 0.02,  tolerance = 1e-3)
})
