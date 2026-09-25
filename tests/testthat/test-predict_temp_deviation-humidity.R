# =============================================================================
# predict_tos() with a humidity-modified temperature expression
# -----------------------------------------------------------------------------
# When a model is fitted with expr_arrhenius_humidity(), the temperature
# expression references a data column (`rh`) that is neither a parameter nor a
# grouping variable. The composite propagation must therefore read that column
# from the schedule, per segment.
#
# The tightest available check is the constant-condition equivalence used in
# test-predict_temp_deviation-equivalence.R: a single-segment schedule holding
# one (temp, rh) pair for duration t is the same experiment as predicting at
# (time = t, tempk = T, rh = RH), so the point predictions must agree.
# =============================================================================


.seg_time <- 24
.seg_TT   <- 298.15
.seg_RH   <- 0.7

.hum_expr <- function(...) expr_arrhenius_humidity(rh_col = "rh", rh_ref = 0.6, ...)

# ---------------------------------------------------------------------------
# fit_kinetic
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_kinetic with a humidity expression", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))
  set.seed(123)
  df$rh <- runif(nrow(df), 0.4, 0.8)

  fit <- fit_kinetic(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    intercept_cols = "batch",
    temp_expr      = .hum_expr()
  )
  expect_true(all(c("Ea", "b") %in% names(coef(fit))))

  batch <- "Combi_F1-PM1"
  ref <- as.numeric(predict(fit, newdata = data.frame(
    time = .seg_time, tempk = .seg_TT, batch = batch, rh = .seg_RH
  )))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT,
                         batch = batch, rh = .seg_RH)
  tos <- tos_segments(predict_tos(fit, schedule, temp_unit = "K"))

  expect_equal(tos$fit, ref, tolerance = 1e-6)
})

test_that("predict_tos propagates a humidity schedule segment by segment", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))
  set.seed(123)
  df$rh <- runif(nrow(df), 0.4, 0.8)

  fit <- fit_kinetic(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    intercept_cols = "batch",
    temp_expr      = .hum_expr()
  )
  batch <- "Combi_F1-PM1"

  # Two segments at the same temperature but different RH. The durations are
  # deliberately unequal: with equal durations the first-order exponent
  # k * sum(arr(rh_i) * dt) is symmetric under swapping the RH values, so the
  # swap check below would pass even if RH were ignored per segment.
  schedule <- data.frame(
    t_start = c(0, 12),
    t_end   = c(12, 48),
    temp    = c(.seg_TT, .seg_TT),
    batch   = batch,
    rh      = c(0.45, 0.80)
  )
  tos <- tos_segments(predict_tos(fit, schedule, temp_unit = "K"))
  expect_equal(nrow(tos), 2L)
  expect_true(all(is.finite(tos$fit)))

  # Segment 1 alone must equal the corresponding two-segment cumulative row.
  seg1 <- tos_segments(predict_tos(fit, schedule[1, ], temp_unit = "K"))
  expect_equal(tos$fit[1], seg1$fit, tolerance = 1e-8)

  # RH must actually influence the prediction: pairing the high RH with the
  # long segment degrades more than pairing it with the short one.
  schedule_swapped <- schedule
  schedule_swapped$rh <- rev(schedule$rh)
  tos_swapped <- tos_segments(predict_tos(fit, schedule_swapped, temp_unit = "K"))
  expect_gt(abs(tos$fit[2] - tos_swapped$fit[2]), 1e-6)
})

test_that("predict_tos errors informatively when the schedule omits rh", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))
  set.seed(123)
  df$rh <- runif(nrow(df), 0.4, 0.8)

  fit <- fit_kinetic(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    intercept_cols = "batch",
    temp_expr      = .hum_expr()
  )

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT,
                         batch = "Combi_F1-PM1")
  expect_error(predict_tos(fit, schedule, temp_unit = "K"), "rh")

  # Present but NA is equally unusable.
  schedule$rh <- NA_real_
  expect_error(predict_tos(fit, schedule, temp_unit = "K"), "NA")
})

test_that("predict_tos still works for models fitted without humidity", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  fit <- fit_kinetic(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    intercept_cols = "batch"
  )

  batch <- "Combi_F1-PM1"
  ref <- as.numeric(predict(fit, newdata = data.frame(
    time = .seg_time, tempk = .seg_TT, batch = batch
  )))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT,
                         batch = batch)
  tos <- tos_segments(predict_tos(fit, schedule, temp_unit = "K"))

  expect_equal(tos$fit, ref, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# fit_singledeg
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_singledeg with a humidity expression", {
  data("nrcesds")
  df <- nrcesds |> tidyr::pivot_wider(names_from = "cqa", values_from = "y")
  set.seed(9)
  df$rh <- runif(nrow(df), 0.4, 0.8)

  fit <- fit_singledeg(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    yA_var         = "NR Major Peak",
    yB_var         = "NR Minor Peaks",
    intercept_cols = "batch",
    trefk          = 313.15,
    temp_expr      = .hum_expr()
  )
  expect_true(all(c("Ea", "b") %in% names(coef(fit))))

  batch <- 2
  ref <- predict(fit, newdata = data.frame(
    time = .seg_time, tempk = .seg_TT, batch = batch, rh = .seg_RH
  ))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT,
                         batch = batch, rh = .seg_RH)
  tos_A <- tos_segments(predict_tos(fit, schedule, temp_unit = "K", var = "yA"))
  tos_B <- tos_segments(predict_tos(fit, schedule, temp_unit = "K", var = "yB"))

  expect_equal(tos_A$fit, ref$fit_yA, tolerance = 1e-6)
  expect_equal(tos_B$fit, ref$fit_yB, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# fit_parallel  (parameters are branch-suffixed: Ea -> EaB/EaC, b -> bB/bC,
#                but the `rh` data column is shared by both branches)
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_parallel with a humidity expression", {
  data("upsec")
  df <- upsec |> tidyr::pivot_wider(names_from = "cqa", values_from = "y")
  set.seed(11)
  df$rh <- runif(nrow(df), 0.4, 0.8)

  fit <- fit_parallel(
    df,
    kinetic_order  = 1,
    yA_var         = "Total Main",
    yB_var         = "HMWS",
    yC_var         = "LMWS",
    intercept_cols = "batch",
    temp_expr      = .hum_expr()
  )
  expect_true(all(c("EaB", "bB", "EaC", "bC") %in% names(coef(fit))))

  batch <- 2
  ref <- predict(fit, newdata = data.frame(
    time = .seg_time, tempk = .seg_TT, batch = batch, rh = .seg_RH
  ))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT,
                         batch = batch, rh = .seg_RH)

  for (v in c("yA", "yB", "yC")) {
    tos <- tos_segments(predict_tos(fit, schedule, temp_unit = "K", var = v))
    expect_equal(tos$fit, ref[[paste0("fit_", v)]], tolerance = 1e-6)
  }
})

# ---------------------------------------------------------------------------
# fit_sb
# ---------------------------------------------------------------------------
test_that("predict_tos matches predict.fit_sb with a humidity expression", {
  set.seed(789)
  dat <- expand.grid(
    time  = c(0, 1, 3, 6, 9, 12),
    tempk = c(278, 298, 313),
    rh    = c(0.4, 0.6, 0.75)
  )
  dat$y <- 100 - 20 * dat$time * 0.1 *
    exp(80 / 8.314 * (1 / 298.15 - 1 / dat$tempk)) *
    exp(2 * (dat$rh - 0.6)) +
    rnorm(nrow(dat), 0, 2)

  fit <- fit_sb(dat, trend = "decreasing", Yinf = 0, temp_expr = .hum_expr())
  expect_true(all(c("Ea", "b") %in% names(coef(fit))))

  ref <- as.numeric(predict(fit, newdata = data.frame(
    time = .seg_time, tempk = 298, rh = .seg_RH
  )))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = 298,
                         rh = .seg_RH)
  tos <- tos_segments(predict_tos(fit, schedule, temp_unit = "K"))

  # The SB path integrates the ODE numerically on both sides.
  expect_equal(tos$fit, ref, tolerance = 1e-3)
})

# ---------------------------------------------------------------------------
# Custom column name / reference value must be honoured, not hard-coded.
# ---------------------------------------------------------------------------
test_that("predict_tos honours a custom humidity column name", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))
  set.seed(321)
  df$humidity <- runif(nrow(df), 0.4, 0.8)

  fit <- fit_kinetic(
    df,
    kinetic_order  = 1,
    trend          = "decreasing",
    intercept_cols = "batch",
    temp_expr      = expr_arrhenius_humidity(rh_col = "humidity", rh_ref = 0.75)
  )

  batch <- "Combi_F1-PM1"
  ref <- as.numeric(predict(fit, newdata = data.frame(
    time = .seg_time, tempk = .seg_TT, batch = batch, humidity = .seg_RH
  )))

  schedule <- data.frame(t_start = 0, t_end = .seg_time, temp = .seg_TT,
                         batch = batch, humidity = .seg_RH)
  tos <- tos_segments(predict_tos(fit, schedule, temp_unit = "K"))

  expect_equal(tos$fit, ref, tolerance = 1e-6)

  # And it must complain about `humidity`, not the default `rh`.
  expect_error(
    predict_tos(fit, schedule[, setdiff(names(schedule), "humidity")],
                temp_unit = "K"),
    "humidity"
  )
})
