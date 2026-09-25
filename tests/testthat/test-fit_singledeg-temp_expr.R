test_that("fit_singledeg accepts custom temp_expr parameter", {

  # Load and prepare test data
  # data("nrcesds")
  test_data <- nrcesds |>
    tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  # Test 1: Default behavior (NULL temp_expr should use expr_arrhenius)
  fit_default <- fit_singledeg(
    df = test_data,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    trefk = 298.15
  )

  # Should have Ea parameter
  expect_true("Ea" %in% names(coef(fit_default)))

  # Test 2: Custom Arrhenius with humidity
  # Create synthetic humidity data for testing
  test_data_rh <- test_data
  set.seed(456)
  test_data_rh$rh <- runif(nrow(test_data_rh), 0.4, 0.8)

  custom_expr_rh <- expr_arrhenius_humidity(
    R = 8.314,
    trefk = 298.15,
    temp_var = "tempk",
    rh_col = "rh",
    rh_ref = 0.6,
    log_Ea = FALSE
  )

  fit_humidity <- fit_singledeg(
    df = test_data_rh,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    temp_expr = custom_expr_rh
  )

  # Should have both Ea and b parameters
  expect_true("Ea" %in% names(coef(fit_humidity)))
  expect_true("b" %in% names(coef(fit_humidity)))

  # Test 3: Verify backward compatibility
  # Fit with old method (no temp_expr)
  fit_old <- fit_singledeg(
    df = test_data,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    log_Ea = FALSE,
    trefk = 298.15
  )

  # Fit with new method but NULL temp_expr (should be identical)
  fit_new <- fit_singledeg(
    df = test_data,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    log_Ea = FALSE,
    trefk = 298.15,
    temp_expr = NULL
  )

  # Coefficients should be the same
  expect_equal(coef(fit_old), coef(fit_new), tolerance = 1e-10)
})

test_that("fit_singledeg with log-scale humidity parameters works correctly", {

  # data("nrcesds")
  test_data <- nrcesds |>
    tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  # Test with log_Ea = TRUE in humidity expression
  custom_expr_log_ea <- expr_arrhenius_humidity(
    R = 8.314,
    trefk = 298.15,
    temp_var = "tempk",
    rh_col = "rh",
    rh_ref = 0.6,
    log_Ea = TRUE
  )

  # Add synthetic humidity data
  test_data_rh <- test_data
  test_data_rh$rh <- sample(c(0.60,0.2),size=nrow(test_data_rh),replace=TRUE)

  fit_log_ea <- fit_singledeg(
    df = test_data_rh,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    temp_expr = custom_expr_log_ea
  )

  # Should have Ea_log and b parameters
  expect_true("Ea_log" %in% names(coef(fit_log_ea)))
  expect_true("b" %in% names(coef(fit_log_ea)))
  expect_false("Ea" %in% names(coef(fit_log_ea)))
})