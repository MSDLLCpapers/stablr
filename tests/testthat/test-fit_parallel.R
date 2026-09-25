test_that("fit_parallel function works with suitable data (e.g., UPSEC).", {
    
  ### fit_parallel with SEC Monomer, HMW, and LMW data ###
  data("upsec")
  df <- upsec |> dplyr::filter(cqa == "Total Main" | cqa == "HMWS" | cqa == "LMWS")
  df_wide <- df |>
    tidyr::pivot_wider(
      id_cols = c("batch", "temp", "time", "tempk"),
      names_from = "cqa",
      values_from = "y"
    )
  
  # Stablr - sec parallel
  stablr_sec_parallel <- fit_parallel(
    df = df_wide,
    kinetic_order = 1,
    log_Ea = FALSE,
    yA_var = "Total Main",
    yB_var = "HMWS",
    yC_var = "LMWS",
    trefk = 278.15
  )
  # No stablr - sec parallel
  start <- c(kBref = NA, kCref = NA, EaB = NA, EaC = NA, A0 = NA, B0 = NA, C0 = NA)
  parallel_fxn <- function(params, data) {
    kBref <- params[1]
    kCref <- params[2]
    EaB <- params[3]
    EaC <- params[4]
    A0 <- params[5]
    B0 <- params[6]
    C0 <- params[7]
    tempk <- data$tempk
    time <- data$time
    kB <- kBref * exp(EaB / 8.314 * ((1 / 278.15) - (1 / tempk)))
    kC <- kCref * exp(EaC / 8.314 * ((1 / 278.15) - (1 / tempk)))
    kApp <- kB + kC
    `Total Main` <- A0 * exp(-kApp * time)
    HMWS <- B0 + (kB / kApp) * (A0 - `Total Main`)
    LMWS <- C0 + (kC / kApp) * (A0 - `Total Main`)
    vec <- c(`Total Main`, HMWS, LMWS)
    vec
  }
  vec_fin <- c(df_wide$`Total Main`, df_wide$HMWS, df_wide$LMWS)
  no_stablr_sec_parallel <- gslnls::gsl_nls(
    fn = parallel_fxn,
    y = vec_fin,
    start = start,
    data = df_wide
  )
  # Extract coefs
  stablr_coef_sec <- coef(stablr_sec_parallel)
  other_coef_sec <- coef(no_stablr_sec_parallel)
  
  # testthat:: fxns
  testthat::expect_equal(
    stablr_coef_sec[["kB_0"]],
    other_coef_sec[["kBref"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_sec[["kC_0"]],
    other_coef_sec[["kCref"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_sec[["EaB"]],
    other_coef_sec[["EaB"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_sec[["EaC"]],
    other_coef_sec[["EaC"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_sec[["A0"]],
    other_coef_sec[["A0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_sec[["B0"]],
    other_coef_sec[["B0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    stablr_coef_sec[["C0"]],
    other_coef_sec[["C0"]],
    tolerance = 1e-5
  )
})

test_that("predict.fit_parallel works with tolerance intervals", {
    
  # Prepare data
  data("upsec")
  df <- upsec |> dplyr::filter(cqa == "Total Main" | cqa == "HMWS" | cqa == "LMWS")
  df_wide <- df |>
    tidyr::pivot_wider(
      id_cols = c("batch", "temp", "time", "tempk"),
      names_from = "cqa",
      values_from = "y"
    )
  
  # Fit model
  fit <- stablr::fit_parallel(
    df = df_wide,
    kinetic_order = 1,
    log_Ea = FALSE,
    yA_var = "Total Main",
    yB_var = "HMWS",
    yC_var = "LMWS",
    trefk = 278.15
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
  expect_true(all(c("fit_yA", "lwr_yA", "upr_yA", "fit_yB", "lwr_yB", "upr_yB", 
                    "fit_yC", "lwr_yC", "upr_yC") %in% names(ci)))
  
  # Test prediction intervals
  pi <- predict(fit, newdata = pred_grid, interval = "prediction", level = 0.95)
  expect_s3_class(pi, "data.frame")
  expect_true(all(c("fit_yA", "lwr_yA", "upr_yA", "fit_yB", "lwr_yB", "upr_yB",
                    "fit_yC", "lwr_yC", "upr_yC") %in% names(pi)))
  
  # Test two-sided tolerance intervals
  ti_two <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.90, level = 0.95, side = "two")
  )
  expect_s3_class(ti_two, "data.frame")
  expect_true(all(c("fit_yA", "lwr_yA", "upr_yA", "fit_yB", "lwr_yB", "upr_yB",
                    "fit_yC", "lwr_yC", "upr_yC") %in% names(ti_two)))
  
  # Test one-sided lower tolerance intervals
  ti_lower <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.90, level = 0.95, side = "lower")
  )
  expect_s3_class(ti_lower, "data.frame")
  expect_true(all(ti_lower$upr_yA == Inf))
  expect_true(all(ti_lower$upr_yB == Inf))
  expect_true(all(ti_lower$upr_yC == Inf))
  
  # Test one-sided upper tolerance intervals
  ti_upper <- suppressWarnings(
    predict(fit, newdata = pred_grid, interval = "tolerance", 
            coverage = 0.90, level = 0.95, side = "upper")
  )
  expect_s3_class(ti_upper, "data.frame")
  expect_true(all(ti_upper$lwr_yA == -Inf))
  expect_true(all(ti_upper$lwr_yB == -Inf))
  expect_true(all(ti_upper$lwr_yC == -Inf))
  
  # Verify tolerance intervals with high coverage are wider than confidence intervals
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

test_that("intercept prefixes do not mangle level values containing 'y0'", {
  set.seed(1)

  make_df <- function(levels) {
    d <- expand.grid(batch = levels, tempk = c(298.15, 313.15),
                     time = c(0, 3, 6, 9, 12), stringsAsFactors = FALSE)
    d$A <- 100 - 0.4 * d$time + rnorm(nrow(d), 0, 0.2)
    d$B <- 0.3 * d$time + rnorm(nrow(d), 0, 0.2)
    d$C <- 0.1 * d$time + rnorm(nrow(d), 0, 0.2)
    d
  }

  # Renaming y0 -> A0/B0/C0 on the built expression used to rewrite the quoted
  # level literals too, leaving the formula referencing parameters that do not
  # exist.
  # (The data has no temperature dependence, so Ea is unidentified and can drift
  # slightly negative; that is out of scope here, so suppress the sign warning.)
  d <- make_df(c("Ly0t1", "Ly0t2"))
  expect_no_error(
    suppressWarnings(fit_singledeg(d, yA_var = "A", yB_var = "B", trend = "decreasing",
                  Yinf = 0, intercept_cols = "batch"))
  )
  expect_no_error(
    suppressWarnings(fit_parallel(d, yA_var = "A", yB_var = "B", yC_var = "C",
                 Yinf = 0, intercept_cols = "batch"))
  )

  # Worse: levels that collide with the new prefix produced duplicate indicator
  # terms, silently pooling groups instead of erroring. Each level must keep its
  # own intercept, recovering the true values.
  d2 <- expand.grid(batch = c("y0", "A0", "B0"), tempk = c(298.15, 313.15),
                    time = c(0, 3, 6, 9, 12), stringsAsFactors = FALSE)
  d2$A <- ifelse(d2$batch == "y0", 100, ifelse(d2$batch == "A0", 70, 50)) -
    0.4 * d2$time + rnorm(nrow(d2), 0, 0.2)
  d2$B <- 0.3 * d2$time + rnorm(nrow(d2), 0, 0.2)

  fit <- suppressWarnings(fit_singledeg(d2, yA_var = "A", yB_var = "B", trend = "decreasing",
                       Yinf = 0, intercept_cols = "batch"))
  ints <- coef(fit)[c("A0_batch_y0", "A0_batch_A0", "A0_batch_B0")]
  expect_equal(unname(ints), c(100, 70, 50), tolerance = 1)

  # The grouping column name itself may contain "y0".
  d3 <- make_df(c("P", "Q"))
  names(d3)[names(d3) == "batch"] <- "y0lot"
  expect_no_error(
    suppressWarnings(fit_parallel(d3, yA_var = "A", yB_var = "B", yC_var = "C",
                 Yinf = 0, intercept_cols = "y0lot"))
  )
})

# ---------------------------------------------------------------------------
# The parallel model no longer exposes `trend` (increasing produced negative
# concentrations). The argument is gone, and `Yinf` is a documented no-op that
# warns on a non-zero value.
# ---------------------------------------------------------------------------
test_that("fit_parallel and build_parallel reject a `trend` argument", {
  data("upsec")
  df_wide <- upsec |> tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  expect_error(
    fit_parallel(df_wide, kinetic_order = 1, trend = "increasing",
                 yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS"),
    "unused argument"
  )
  expect_error(
    fit_parallel(df_wide, kinetic_order = 1, trend = "decreasing",
                 yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS"),
    "unused argument"
  )
  expect_error(
    build_parallel(df_wide, kinetic_order = 1, trend = "increasing",
                   yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS",
                   temp_fml = expr_arrhenius()$formula),
    "unused argument"
  )
})

test_that("build_parallel warns and ignores a non-zero Yinf", {
  data("upsec")
  df_wide <- upsec |> tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  expect_warning(
    build_parallel(df_wide, kinetic_order = 1, Yinf = 50,
                   yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS",
                   temp_fml = expr_arrhenius()$formula),
    "Yinf.*not used"
  )

  # The warning does not change the model: Yinf = 50 fits identically to Yinf = 0.
  b0  <- build_parallel(df_wide, kinetic_order = 1, Yinf = 0,
                        yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS",
                        temp_fml = expr_arrhenius()$formula)
  b50 <- suppressWarnings(
    build_parallel(df_wide, kinetic_order = 1, Yinf = 50,
                   yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS",
                   temp_fml = expr_arrhenius()$formula)
  )
  pars <- stats::setNames(rep(0.01, length(b0$param_names)), b0$param_names)
  expect_equal(b0$fn(pars, df_wide), b50$fn(pars, df_wide))
})
