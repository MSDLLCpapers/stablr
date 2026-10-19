test_that("fit_kinetic works for orders and trends.", {
  # All calls below are namespaced (dplyr::, stablr::, gslnls::, testthat::), so
  # tidyverse need not be attached. Loading it here only emitted "package built
  # under R 4.4.3" startup warnings.

  ### Zero-order decreasing - Purity ###
  # data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  # stablr function - zero order - purity
  stablr_zero_purity <- suppressWarnings(stablr::fit_kinetic(
    df,
    kinetic_order = 0,
    trend = "decreasing",
    log_Ea = TRUE,
    trefk = 278.15
  ))

  # no stablr - zero order - purity
  estimate_kkf <- function(
    df, 
    trefk = 278.15, 
    R = 8.314, 
    start = NULL
  ){
    # check data contains required variables
    stopifnot(all(c("tempk", "y", "time") %in% colnames(df)))
    #  check that tempk > 0
    stopifnot(all(df$tempk > 0))
    # starting parameters
    if (is.null(start)) start <- c(kref = NA, Ea = NA, A0 = NA)
    ### Fit model using NLS regression ###
    mod <- gslnls::gsl_nls(
      fn = as.formula(
        "y ~ A0 - kref * exp(exp(Ea) / R * (1 / trefk - 1 / tempk)) * time"
      ),
      data = df,
      start = start
    )
    mod
  }

  nostablr_zero_purity <- suppressWarnings(estimate_kkf(df))

  # Coefs - zero order - purity
  coefs_stablr_zero_purity <- coef(stablr_zero_purity)
  coefs_nostablr_zero_purity <- coef(nostablr_zero_purity)

  # testthat - zero order - purity
  testthat::expect_equal(
    coefs_nostablr_zero_purity[["kref"]],
    coefs_stablr_zero_purity[["k_0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_zero_purity[["Ea"]],
    coefs_stablr_zero_purity[["Ea_log"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_zero_purity[["A0"]],
    coefs_stablr_zero_purity[["y0"]],
    tolerance = 1e-5
  )

  ### First-order decreasing - purity ###
  # data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  # stablr function - first order - purity
  stablr_first_purity <- suppressWarnings(stablr::fit_kinetic(
    df,
    kinetic_order = 1,
    trend = "decreasing",
    log_Ea = FALSE,
    trefk = 278.15
  ))

  # no stablr - first order - purity
  estimate_firstorder <- function(
    df,
    trefk = 278.15,
    R = 8.314,
    start = NULL
  ) {
    # check data contains required variables
    stopifnot(all(c("tempk", "y", "time") %in% colnames(df)))

    #  check that tempk > 0
    stopifnot(all(df$tempk > 0))

    # starting parameters
    if (is.null(start)) start <- c(kref = NA, Ea = NA, A0 = NA)

    firstorder_fxn <- function(params, data) {
      kref <- params[1]
      ea <- params[2]
      a0 <- params[3]

      tempk <- data$tempk
      time <- data$time

      k <- kref * exp(ea / R * ((1 / trefk) - (1 / tempk)))
      y_a <- a0 * exp(-k * time)
      y_a
    }

    mod <- gslnls::gsl_nls(
      fn = firstorder_fxn,
      y = df$y,
      start = start,
      data = df
    )
    mod
  }

  nostablr_first_purity <- suppressWarnings(estimate_firstorder(df))

  # Coefs - first order - purity
  coefs_stablr_first_purity <- coef(stablr_first_purity)
  coefs_nostablr_first_purity <- coef(nostablr_first_purity)

  # testthat - first order - monomer
  testthat::expect_equal(
    coefs_nostablr_first_purity[["kref"]],
    coefs_stablr_first_purity[["k_0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_first_purity[["Ea"]],
    coefs_stablr_first_purity[["Ea"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_first_purity[["A0"]],
    coefs_stablr_first_purity[["y0"]],
    tolerance = 1e-5
  )

  ### Zero-order Increasing - HMW ###
  #data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "HMW" & !is.na(y))

  # stablr function - zero order - HMW
  stablr_zero_hmw <- suppressWarnings(stablr::fit_kinetic(
    df,
    kinetic_order = 0,
    trend = "increasing",
    log_Ea = TRUE,
    trefk = 278.15
  ))

  # no stablr - zero order - hmw
  estimate_kkf <- function(
    df,
    trefk = 278.15,
    R = 8.314,
    start = NULL
  ) {
    # check data contains required variables
    stopifnot(all(c("tempk", "y", "time") %in% colnames(df)))
    #  check that tempk > 0
    stopifnot(all(df$tempk > 0))
    # starting parameters
    if (is.null(start)) start <- c(kref = NA, Ea = NA, A0 = NA)
    ### Fit model using NLS regression ###
    mod <- gslnls::gsl_nls(
      fn = as.formula(
        "y ~ A0 + kref * exp(exp(Ea) / R * (1 / trefk - 1 / tempk)) * time"
      ),
      data = df,
      start = start
    )
    mod
  }

  nostablr_zero_hmw <- suppressWarnings(estimate_kkf(df))

  # Coefs - zero order - hmw
  coefs_stablr_zero_hmw <- coef(stablr_zero_hmw)
  coefs_nostablr_zero_hmw <- coef(nostablr_zero_hmw)

  # testthat - zero order - monomer
  testthat::expect_equal(
    coefs_nostablr_zero_hmw[["kref"]],
    coefs_stablr_zero_hmw[["k_0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_zero_hmw[["Ea"]],
    coefs_stablr_zero_hmw[["Ea_log"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_zero_hmw[["A0"]],
    coefs_stablr_zero_hmw[["y0"]],
    tolerance = 1e-5
  )

  ### First-order decreasing - hmw ###
  # data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "HMW" & !is.na(y))

  # stablr function - first order - hmw
  # HMW trends upward, but we deliberately fit a decreasing model to compare
  # numerics against the hand-rolled oracle below, so silence the expected
  # trend-mismatch warning.
  stablr_first_hmw <- suppressWarnings(stablr::fit_kinetic(
    df,
    kinetic_order = 1,
    trend = "decreasing",
    log_Ea = FALSE,
    trefk = 278.15
  ))

  # no stablr - first order - monomer
  estimate_firstorder <- function(
    df,
    trefk = 278.15,
    R = 8.314,
    start = NULL
  ) {
    # check data contains required variables
    stopifnot(all(c("tempk", "y", "time") %in% colnames(df)))

    #  check that tempk > 0
    stopifnot(all(df$tempk > 0))

    # starting parameters
    if (is.null(start)) start <- c(kref = NA, Ea = NA, A0 = NA)

    firstorder_fxn <- function(params, data) {
      kref <- params[1]
      ea <- params[2]
      a0 <- params[3]

      tempk <- data$tempk
      time <- data$time

      k <- kref * exp(ea / R * ((1 / trefk) - (1 / tempk)))
      y_a <- a0 * exp(-k * time)
      y_a
    }

    mod <- gslnls::gsl_nls(
      fn = firstorder_fxn,
      y = df$y,
      start = start,
      data = df
    )
    mod
  }

  nostablr_first_hmw <- suppressWarnings(estimate_firstorder(df))

  # Coefs - first order - hmw
  coefs_stablr_first_hmw <- coef(stablr_first_hmw)
  coefs_nostablr_first_hmw <- coef(nostablr_first_hmw)

  # testthat - first order - hmw
  testthat::expect_equal(
    coefs_nostablr_first_hmw[["kref"]],
    coefs_stablr_first_hmw[["k_0"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_first_hmw[["Ea"]],
    coefs_stablr_first_hmw[["Ea"]],
    tolerance = 1e-5
  )
  testthat::expect_equal(
    coefs_nostablr_first_hmw[["A0"]],
    coefs_stablr_first_hmw[["y0"]],
    tolerance = 1e-5
  )
})
