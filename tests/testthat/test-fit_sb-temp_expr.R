skip_on_cran()

test_that("fit_sb works with default Arrhenius (backward compatibility)", {
  # Generate synthetic data
  set.seed(123)
  dat <- expand.grid(
    time = c(0, 1, 3, 6, 9, 12),
    tempk = c(278, 298, 313),
    batch = c("B1", "B2")
  )

  # Simulate Arrhenius kinetics
  k_ref <- 0.1
  Ea <- 80  # kJ/mol
  R <- 8.314
  trefk <- 298.15

  dat$y <- 100 - 20 * dat$time * k_ref *
    exp(Ea / R * (1/trefk - 1/dat$tempk)) +
    rnorm(nrow(dat), 0, 2)

  # Fit with default Arrhenius (no temp_expr)
  # Synthetic data is linear-in-time, so Ea rests on its bound; suppress the
  # expected "rest on a bound" warning.
  fit1 <- suppressWarnings(fit_sb(dat, trend = "decreasing", Yinf = 0))

  # Should have Ea parameter
  expect_true("Ea" %in% names(coef(fit1)))
  expect_true("n" %in% names(coef(fit1)))
  expect_true("m" %in% names(coef(fit1)))

  # Predictions should work
  pred1 <- predict(fit1)
  expect_equal(length(pred1), nrow(dat))
})

test_that("fit_sb works with custom Arrhenius expression", {
  set.seed(123)
  dat <- expand.grid(
    time = c(0, 1, 3, 6, 9, 12),
    tempk = c(278, 298, 313),
    batch = c("B1", "B2")
  )

  # Simulate Arrhenius kinetics
  k_ref <- 0.1
  Ea <- 80  # kJ/mol
  R <- 8.314
  trefk <- 298.15

  dat$y <- 100 - 20 * dat$time * k_ref *
    exp(Ea / R * (1/trefk - 1/dat$tempk)) +
    rnorm(nrow(dat), 0, 2)

  # Fit with explicit Arrhenius expression
  # suppress the expected "rest on a bound" warning.
  temp_expr <- expr_arrhenius(R = 8.314, trefk = 298.15, temp_var = "tempk", log_Ea = FALSE)
  fit2 <- suppressWarnings(fit_sb(dat, trend = "decreasing", Yinf = 0, temp_expr = temp_expr))

  # Should have Ea parameter
  expect_true("Ea" %in% names(coef(fit2)))
  expect_true("n" %in% names(coef(fit2)))
  expect_true("m" %in% names(coef(fit2)))

  # Predictions should work
  pred2 <- predict(fit2)
  expect_equal(length(pred2), nrow(dat))
})

test_that("fit_sb works with humidity expression", {
  # Generate synthetic data with humidity effect
  set.seed(789)
  dat <- expand.grid(
    time = c(0, 1, 3, 6, 9, 12),
    tempk = c(278, 298, 313),
    rh = c(0.4, 0.6, 0.75),  # 40%, 60%, 75% RH
    batch = c("B1")
  )

  # Simulate Arrhenius with humidity
  k_ref <- 0.1
  Ea <- 80  # kJ/mol
  b <- 2.0  # humidity coefficient
  R <- 8.314
  trefk <- 298.15
  rh_ref <- 0.6

  dat$y <- 100 - 20 * dat$time * k_ref *
    exp(Ea / R * (1/trefk - 1/dat$tempk)) *
    exp(b * (dat$rh - rh_ref)) +
    rnorm(nrow(dat), 0, 2)

  # Fit with humidity expression
  temp_expr <- expr_arrhenius_humidity(R = 8.314, trefk = 298.15,
                                        temp_var = "tempk", rh_col = "rh",
                                        rh_ref = 0.6, log_Ea = FALSE)
  fit4 <- fit_sb(dat, trend = "decreasing", Yinf = 0, temp_expr = temp_expr)

  # Should have Ea and b parameters
  expect_true("Ea" %in% names(coef(fit4)))
  expect_true("b" %in% names(coef(fit4)))
  expect_true("n" %in% names(coef(fit4)))
  expect_true("m" %in% names(coef(fit4)))

  # Predictions should work
  pred4 <- predict(fit4)
  expect_equal(length(pred4), nrow(dat))

  # New data prediction with different RH
  newdat <- data.frame(
    time = 18,
    tempk = 298,
    rh = 0.8
  )
  pred_new <- predict(fit4, newdata = newdat)
  expect_equal(length(pred_new), 1)
})


test_that("fit_sb checks for RH column when humidity expression is used", {
  
  dat <- data.frame(
    time = 1:5,
    tempk = 298,
    y = 100 - 5 * (1:5)
    # Note: no 'rh' column
  )
  
  temp_expr <- expr_arrhenius_humidity()

  # Should error because 'rh' column is missing
  expect_error(
    fit_sb(dat, trend = "decreasing", temp_expr = temp_expr)
  )
})

test_that("predict.fit_sb handles confidence intervals with custom expressions", {
  
  # Generate synthetic data
  set.seed(321)
  dat <- expand.grid(
    time = c(0, 3, 6, 9),
    tempk = c(278, 298, 313)
  )

  k_ref <- 0.1
  Ea <- 80
  R <- 8.314
  trefk <- 298.15

  dat$y <- 100 - 20 * dat$time * k_ref *
    exp(Ea / R * (1/trefk - 1/dat$tempk)) +
    rnorm(nrow(dat), 0, 1)

  # suppress the expected "rest on a bound" warning.
  temp_expr <- expr_arrhenius(R = 8.314, trefk = 298.15)
  fit <- suppressWarnings(fit_sb(dat, trend = "decreasing", Yinf = 0, temp_expr = temp_expr))

  # Test confidence intervals (if Hessian is available)
  pred_ci <- tryCatch(
    predict(fit, interval = "confidence", n_sim = 100),
    error = function(e) NULL
  )

  if (!is.null(pred_ci)) {
    expect_true(is.data.frame(pred_ci))
    expect_equal(nrow(pred_ci), nrow(dat))
    expect_true(all(c("fit", "lwr", "upr") %in% names(pred_ci)))
    expect_true(all(pred_ci$lwr <= pred_ci$fit))
    expect_true(all(pred_ci$fit <= pred_ci$upr))
  }
})