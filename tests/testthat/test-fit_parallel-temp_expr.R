test_that("fit_parallel maintains backward compatibility with default Arrhenius", {

  # Use built-in data
  # data("upsec")
  upsec_wide <- upsec |>
    tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  # Fit with default (implicit Arrhenius)
  fit_default <- fit_parallel(
    df = upsec_wide,
    kinetic_order = 1,
    yA_var = "Total Main",
    yB_var = "HMWS",
    yC_var = "LMWS"
  )

  # Fit with explicit temp_expr = NULL (should be identical)
  fit_explicit_null <- fit_parallel(
    df = upsec_wide,
    kinetic_order = 1,
    temp_expr = NULL,
    yA_var = "Total Main",
    yB_var = "HMWS",
    yC_var = "LMWS"
  )

  # Check that parameter names are identical
  expect_equal(names(coef(fit_default)), names(coef(fit_explicit_null)))

  # Check that coefficients are very close (allowing for numerical differences)
  expect_equal(coef(fit_default), coef(fit_explicit_null), tolerance = 1e-8)

  # Check that default includes EaB and EaC parameters
  param_names <- names(coef(fit_default))
  expect_true("EaB" %in% param_names)
  expect_true("EaC" %in% param_names)
})

test_that("fit_parallel works with expr_arrhenius function", {

  # data("upsec")
  upsec_wide <- upsec |>
    tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  # Fit with expr_arrhenius as function
  fit_arrhenius <- fit_parallel(
    df = upsec_wide,
    kinetic_order = 1,
    temp_expr = expr_arrhenius(),
    yA_var = "Total Main",
    yB_var = "HMWS",
    yC_var = "LMWS"
  )

  # Check parameter names
  param_names <- names(coef(fit_arrhenius))
  expect_true("EaB" %in% param_names)
  expect_true("EaC" %in% param_names)

  # Should NOT have log versions when log_Ea = FALSE
  expect_false("EaB_log" %in% param_names)
  expect_false("EaC_log" %in% param_names)
})

test_that("fit_parallel works with expr_arrhenius_humidity function", {

  # Coherent parallel + humidity data generated from the package's own model
  # form, so the fit genuinely converges and recovers the known parameters --
  # the previous version used independent random yA/yB/yC, which has no Arrhenius
  # structure and drove the fit to a spurious negative Ea.
  set.seed(123)
  R <- 8.314; Tref <- 298.15; rh_ref <- 0.60
  EaB <- 80000; EaC <- 70000
  krefB <- 0.010; krefC <- 0.006
  bB <- 2.0; bC <- 1.2
  A0 <- 100
  df_humidity <- expand.grid(time = seq(0, 24, length.out = 9),
                             tempk = c(278.15, 298.15, 313.15),
                             rh = c(0.40, 0.60, 0.80))
  kB <- krefB * exp(EaB / R * (1 / Tref - 1 / df_humidity$tempk)) *
    exp(bB * (df_humidity$rh - rh_ref))
  kC <- krefC * exp(EaC / R * (1 / Tref - 1 / df_humidity$tempk)) *
    exp(bC * (df_humidity$rh - rh_ref))
  kA <- kB + kC
  df_humidity$yA <- A0 * exp(-kA * df_humidity$time) +
    rnorm(nrow(df_humidity), 0, 0.3)
  df_humidity$yB <- (kB / kA) * A0 * (1 - exp(-kA * df_humidity$time)) +
    rnorm(nrow(df_humidity), 0, 0.15)
  df_humidity$yC <- (kC / kA) * A0 * (1 - exp(-kA * df_humidity$time)) +
    rnorm(nrow(df_humidity), 0, 0.15)

  # Fit with expr_arrhenius_humidity as function
  fit_humidity <- fit_parallel(
    df = df_humidity,
    kinetic_order = 1,
    temp_expr = expr_arrhenius_humidity(),
    yA_var = "yA",
    yB_var = "yB",
    yC_var = "yC"
  )

  # The humidity path must produce Ea AND b parameters for both branches.
  param_names <- names(coef(fit_humidity))
  expect_true("EaB" %in% param_names)
  expect_true("EaC" %in% param_names)
  expect_true("bB" %in% param_names)
  expect_true("bC" %in% param_names)

  # ...and recover the generating values, confirming the path actually fits.
  expect_true(fit_humidity$convInfo$isConv)
  cf <- coef(fit_humidity)
  expect_equal(unname(cf["EaB"]), EaB, tolerance = 0.02)
  expect_equal(unname(cf["EaC"]), EaC, tolerance = 0.02)
  expect_equal(unname(cf["bB"]),  bB,  tolerance = 0.05)
  expect_equal(unname(cf["bC"]),  bC,  tolerance = 0.05)
})



