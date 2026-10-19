skip_on_cran()

# Single constant-temperature segment used throughout.
.seg_time <- 24      # duration (same time units as the fitted data)
.seg_TT   <- 298.15  # constant temperature, in Kelvin (temp_unit = "K")

# ---------------------------------------------------------------------------
# fit_kinetic  (evers_supp)
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_kinetic (order 1, decreasing)", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  fit <- fit_kinetic(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    intercept_cols = "batch"
  )

  batch <- "Combi_F1-PM1"
  ref <- as.numeric(
    predict(fit, newdata = data.frame(time = .seg_time, tempk = .seg_TT, batch = batch))
  )

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT, batch = batch)
  tos <- predict_tos(fit, schedule, temp_unit = "K")

  expect_equal(tos_segments(tos)$fit, ref, tolerance = 1e-4)
})

test_that("predict_tos matches predict.fit_kinetic (order 2, increasing)", {
  # EXPECTED TO FAIL until the kinetic composite honors kinetic_order / trend.
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "HMW" & !is.na(y))

  fit <- fit_kinetic(
    df,
    kinetic_order  = 2,
    trend          = "increasing",
    Yinf           = 100,
    intercept_cols = "batch"
  )

  batch <- "Combi_F1-PM1"
  ref <- as.numeric(
    predict(fit, newdata = data.frame(time = .seg_time, tempk = .seg_TT, batch = batch))
  )

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT, batch = batch)
  tos <- predict_tos(fit, schedule, temp_unit = "K")

  expect_equal(tos_segments(tos)$fit, ref, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# fit_singledeg  (nrcesds)
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_singledeg (yA and yB)", {
  # EXPECTED TO FAIL until the singledeg composite multiplies the segment
  # rate by its duration dt.
  data("nrcesds")
  df <- nrcesds |> tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  fit <- fit_singledeg(
    df,
    kinetic_order  = 2,
    trend          = "decreasing",
    yA_var         = "NR Major Peak",
    yB_var         = "NR Minor Peaks",
    intercept_cols = "batch",
    trefk          = 313.15
  )

  batch <- 2
  ref <- predict(fit, newdata = data.frame(time = .seg_time, tempk = .seg_TT, batch = batch))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT, batch = batch)
  tos_A <- predict_tos(fit, schedule, temp_unit = "K", var = "yA")
  tos_B <- predict_tos(fit, schedule, temp_unit = "K", var = "yB")

  expect_equal(tos_segments(tos_A)$fit, ref$fit_yA, tolerance = 1e-4)
  expect_equal(tos_segments(tos_B)$fit, ref$fit_yB, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# fit_parallel  (upsec)
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_parallel (yA, yB, yC)", {
  data("upsec")
  df <- upsec |> tidyr::pivot_wider(names_from = "cqa", values_from = "y")

  fit <- fit_parallel(
    df,
    kinetic_order  = 1,
    yA_var         = "Total Main",
    yB_var         = "HMWS",
    yC_var         = "LMWS",
    intercept_cols = "batch"
  )

  batch <- 2
  ref <- predict(fit, newdata = data.frame(time = .seg_time, tempk = .seg_TT, batch = batch))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT, batch = batch)
  tos_A <- predict_tos(fit, schedule, temp_unit = "K", var = "yA")
  tos_B <- predict_tos(fit, schedule, temp_unit = "K", var = "yB")
  tos_C <- predict_tos(fit, schedule, temp_unit = "K", var = "yC")

  expect_equal(tos_segments(tos_A)$fit, ref$fit_yA, tolerance = 1e-4)
  expect_equal(tos_segments(tos_B)$fit, ref$fit_yB, tolerance = 1e-4)
  expect_equal(tos_segments(tos_C)$fit, ref$fit_yC, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# fit_sb  (evers_supp)
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_sb", {
  # Slightly looser tolerance: the SB path integrates the ODE numerically in
  # both the reference prediction and the composite propagation.
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  fit <- fit_sb(
    df,
    rate_cols      = "batch",
    intercept_cols = "batch",
    trend          = "decreasing",
    Yinf           = 0
  )

  batch <- "Combi_F1-PM1"
  ref <- as.numeric(
    predict(fit, newdata = data.frame(time = .seg_time, tempk = .seg_TT, batch = batch))
  )

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT, batch = batch)
  tos <- predict_tos(fit, schedule, temp_unit = "K")

  expect_equal(tos_segments(tos)$fit, ref, tolerance = 1e-3)
})