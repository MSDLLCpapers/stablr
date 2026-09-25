test_that("fit_kinetic accepts custom temp_expr parameter", {
  # Load test data
  data(evers_supp)
  test_data <- subset(evers_supp, cqa == "Purity" & !is.na(y))

  # Test 1: Default behavior (NULL temp_expr should use expr_arrhenius)
  fit_default <- fit_kinetic(
    test_data,
    kinetic_order = 1,
    trend = "decreasing",
    trefk = 298.15
  )

  # Should have Ea parameter
  expect_true("Ea" %in% names(coef(fit_default)))

  # Test 2: Custom Arrhenius with humidity (if humidity column exists)
  # Create synthetic humidity data for testing
  test_data_rh <- test_data
  set.seed(123)
  test_data_rh$rh <- runif(nrow(test_data_rh), 0.4, 0.8)

  custom_expr_rh <- expr_arrhenius_humidity(
    R = 8.314,
    trefk = 298.15,
    temp_var = "tempk",
    rh_col = "rh",
    rh_ref = 0.6,
    log_Ea = FALSE
  )

  fit_humidity <- fit_kinetic(
    test_data_rh,
    kinetic_order = 1,
    trend = "decreasing",
    temp_expr = custom_expr_rh
  )

  # Should have both Ea and b parameters
  expect_true("Ea" %in% names(coef(fit_humidity)))
  expect_true("b" %in% names(coef(fit_humidity)))

  # Test 3: Verify backward compatibility
  # Fit with old method (no temp_expr)
  fit_old <- fit_kinetic(
    test_data,
    kinetic_order = 1,
    trend = "decreasing",
    log_Ea = FALSE,
    trefk = 298.15
  )

  # Fit with new method but NULL temp_expr (should be identical)
  fit_new <- fit_kinetic(
    test_data,
    kinetic_order = 1,
    trend = "decreasing",
    log_Ea = FALSE,
    trefk = 298.15,
    temp_expr = NULL
  )

  # Coefficients should be the same
  expect_equal(coef(fit_old), coef(fit_new), tolerance = 1e-10)
})

