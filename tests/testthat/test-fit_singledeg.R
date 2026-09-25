test_that("fit_singledeg function works with second-order CE-SDS CQAs.", {
  # Major and Minor
  #data("nrcesds")
  df <- nrcesds |> dplyr::filter(!is.na(y))

  df_wide <- df |>
    tidyr::pivot_wider(
      id_cols = c("temp", "tempk", "time", "batch"),
      names_from = "cqa",
      values_from = "y"
    ) |>
    dplyr::filter(!is.na(`NR Major Peak`) & !is.na(`NR Minor Peaks`))

  # Stablr - singledeg - Major and Minor
  stablr_cesdsnr_singledeg <- stablr::fit_singledeg(
    df = df_wide,
    kinetic_order = 2,
    trend = "decreasing",
    log_Ea = FALSE,
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    trefk = 278.15,
    Yinf = 100
  )
  # No stablr
  start <- c(kref = NA, Ea = NA, A0 = NA, B0 = NA)
  singledeg_fxn <- function(params, data) {
    kref <- params[1]
    Ea <- params[2]
    A0 <- params[3]
    B0 <- params[4]
    tempk <- data$tempk
    time <- data$time
    k <- kref * exp(Ea / 8.314 * (( 1 / 278.15)-( 1 / tempk)))
    `NR Major Peak` <- A0 / (1 + k * time * A0)
    `NR Minor Peaks` <- B0 + (A0 - `NR Major Peak`)
    vec <- c(`NR Major Peak`,`NR Minor Peaks`)
    return(c(vec))
  }
  vec_fin <- c(df_wide$`NR Major Peak`, df_wide$`NR Minor Peaks`)
  mod_cesdsnr_singledeg <- gslnls::gsl_nls(
    fn = singledeg_fxn,
    y = vec_fin,
    start = start,
    data = df_wide
  )
  
  # Extract coefs
  stablr_coef_nr <- coef(stablr_cesdsnr_singledeg)
  other_coef_nr <- coef(mod_cesdsnr_singledeg)
  
  # testthat:: fxns
  testthat::expect_equal(
    stablr_coef_nr[["k_0"]],
    other_coef_nr[["kref"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_nr[["Ea"]],
    other_coef_nr[["Ea"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_nr[["A0"]],
    other_coef_nr[["A0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_nr[["B0"]],
    other_coef_nr[["B0"]],
    tolerance = 1e-5
  )
})

test_that("predict.fit_singledeg works with tolerance intervals", {
  
  # Prepare data
  # data("nrcesds")
  df <- nrcesds |> dplyr::filter(!is.na(y))
  
  df_wide <- df |>
    tidyr::pivot_wider(
      id_cols = c("temp", "tempk", "time", "batch"),
      names_from = "cqa",
      values_from = "y"
    ) |>
    dplyr::filter(!is.na(`NR Major Peak`) & !is.na(`NR Minor Peaks`))
  
  # Fit model
  fit <- stablr::fit_singledeg(
    df = df_wide,
    kinetic_order = 2,
    trend = "decreasing",
    log_Ea = FALSE,
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    trefk = 278.15,
    Yinf = 100
  )
  
  # Create prediction grid
  pred_grid <- data.frame(
    time = seq(0, 12, by = 3),
    tempk = 278.15,
    batch = df_wide$batch[1]
  )
  
  # Test confidence intervals
  ci <- predict(fit, newdata = pred_grid, interval = "confidence", level = 0.95)
  expect_s3_class(ci, "data.frame")
  expect_true(all(c("fit_yA", "lwr_yA", "upr_yA", "fit_yB", "lwr_yB", "upr_yB") %in% names(ci)))
  
  # Test prediction intervals
  pi <- predict(fit, newdata = pred_grid, interval = "prediction", level = 0.95)
  expect_s3_class(pi, "data.frame")
  expect_true(all(c("fit_yA", "lwr_yA", "upr_yA", "fit_yB", "lwr_yB", "upr_yB") %in% names(pi)))
  
  # Test two-sided tolerance intervals
  ti_two <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.90, level = 0.95, side = "two")
  )
  expect_s3_class(ti_two, "data.frame")
  expect_true(all(c("fit_yA", "lwr_yA", "upr_yA", "fit_yB", "lwr_yB", "upr_yB") %in% names(ti_two)))
  
  # Test one-sided lower tolerance intervals
  ti_lower <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.90, level = 0.95, side = "lower")
  )
  expect_s3_class(ti_lower, "data.frame")
  expect_true(all(ti_lower$upr_yA == Inf))
  expect_true(all(ti_lower$upr_yB == Inf))
  
  # Test one-sided upper tolerance intervals
  ti_upper <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.90, level = 0.95, side = "upper")
  )
  expect_s3_class(ti_upper, "data.frame")
  expect_true(all(ti_upper$lwr_yA == -Inf))
  expect_true(all(ti_upper$lwr_yB == -Inf))
  
  # Verify tolerance intervals with high coverage are wider than prediction intervals
  # Use 95% coverage to ensure intervals are wider
  ti_high <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.95, level = 0.95, side = "two")
  )
  expect_true(all((ti_high$upr_yA - ti_high$lwr_yA) > (ci$upr_yA - ci$lwr_yA)))
  
  # Test error handling
  expect_error(
    predict(fit, newdata = pred_grid, interval = "tolerance", coverage = 1.5),
    "'coverage' must be between 0 and 1"
  )
  expect_error(
    predict(fit, newdata = pred_grid, interval = "tolerance", level = 1.5),
    "'level' must be between 0 and 1"
  )
})
