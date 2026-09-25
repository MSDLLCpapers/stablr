test_that("predict_tos function works with fit_kinetic models.", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  # First-order
  fit1 <- fit_kinetic(df,
    kinetic_order = 1,
    trend = "decreasing",
    intercept_cols = "batch"
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(278.15, 278.15, 298.15, 278.15),
    batch = "Combi_F1-PM1"
  )

  expect_no_error(
    predict_tos(
      fit = fit1,
      schedule = schedule,
      temp_unit = "K"
    )
  )

  schedule_nobatch <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(278.15, 278.15, 298.15, 278.15)
  )

  expect_error(
    predict_tos(
      fit = fit1,
      schedule = schedule_nobatch,
      temp_unit = "K"
    )
  )

  df <- evers_supp |> dplyr::filter(cqa == "HMW" & !is.na(y))

  # Second-order
  fit2 <- fit_kinetic(df,
    kinetic_order = 2,
    trend = "increasing",
    Yinf = 100,
    intercept_cols = "batch"
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(5, 5, 40, 5),
    batch = "Combi_F1-PM2"
  )

  expect_no_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C"
    )
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(5, 5, 40, 5)
  )

  expect_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C"
    )
  )
})

test_that("predict_tos function works with fit_singledeg models.", {
  data("nrcesds")
  df <- nrcesds |> tidyr::pivot_wider(names_from = "cqa",values_from = "y")
  df$container <- sample(c("container 1","container 2"),nrow(df),replace=TRUE) # mock rate_cols column
  
  # First-order
  fit1 <- fit_singledeg(df,
    kinetic_order = 1,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    intercept_cols = "batch",
    rate_cols = "container",
    trefk = 313.15
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(278.15, 278.15, 313.15, 278.15),
    batch = 2,
    container = "container 2"
  )

  expect_no_error(
    predict_tos(
      fit = fit1,
      schedule = schedule,
      temp_unit = "K",
      var = "yB"
    )
  )

  schedule_nobatch_nocontainer <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(278.15, 278.15, 298.15, 278.15)
  )

  # No intercept_cols or rate_cols specified in schedule.
  expect_error(
    predict_tos(
      fit = fit1,
      schedule = schedule_nobatch_nocontainer,
      temp_unit = "K"
    )
  )

  # Second-order
  fit2 <- fit_singledeg(df,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    intercept_cols = "batch",
    rate_cols = "container",
    trefk = 313.15
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(5, 5, 40, 5),
    batch = 3,
    container = "container 1"
  )

  expect_no_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C",
      var = "yA"
    )
  )

  # No var specified
  expect_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C"
    )
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(5, 5, 40, 5)
  )

  # No intercept_cols or rate_cols specified in schedule.
  expect_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C",
      var = "yA"
    )
  )
})

############# HI
test_that("predict_tos function works with fit_parallel models.", {
  data("nrcesds")
  df <- nrcesds |> tidyr::pivot_wider(names_from = "cqa",values_from = "y")
  df$container <- sample(c("container 1","container 2"),nrow(df),replace=TRUE) # mock rate_cols column
  
  # First-order
  fit1 <- fit_singledeg(df,
    kinetic_order = 1,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    intercept_cols = "batch",
    rate_cols = "container",
    trefk = 313.15
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(278.15, 278.15, 313.15, 278.15),
    batch = 2,
    container = "container 2"
  )

  expect_no_error(
    predict_tos(
      fit = fit1,
      schedule = schedule,
      temp_unit = "K",
      var = "yB"
    )
  )

  schedule_nobatch_nocontainer <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(278.15, 278.15, 298.15, 278.15)
  )

  # No intercept_cols or rate_cols specified in schedule.
  expect_error(
    predict_tos(
      fit = fit1,
      schedule = schedule_nobatch_nocontainer,
      temp_unit = "K"
    )
  )

  # Second-order
  fit2 <- fit_singledeg(df,
    kinetic_order = 2,
    trend = "decreasing",
    yA_var = "NR Major Peak",
    yB_var = "NR Minor Peaks",
    intercept_cols = "batch",
    rate_cols = "container",
    trefk = 313.15
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(5, 5, 40, 5),
    batch = 3,
    container = "container 1"
  )

  expect_no_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C",
      var = "yA"
    )
  )

  # No var specified
  expect_error(
    predict_tos(
      fit = fit2,
      schedule = schedule,
      temp_unit = "C"
    )
  )

  schedule <- data.frame(
    t_start = c(0, 0, 24, 30),
    t_end = c(0, 24, 30, 36),
    temp = c(5, 5, 40, 5)
  )

  # No intercept_cols or rate_cols specified in schedule.
  expect_error(
    predict_tos(
      fit = fit2,
      schedule = schedule_nobatch,
      temp_unit = "C",
      var = "yA"
    )
  )
})


test_that("predict_tos requires a contiguous schedule.", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  fit1 <- fit_kinetic(df,
    kinetic_order = 1,
    trend = "decreasing",
    intercept_cols = "batch"
  )

  base <- data.frame(
    t_start = c(0, 24),
    t_end = c(24, 30),
    temp = c(278.15, 298.15),
    batch = "Combi_F1-PM1"
  )

  expect_no_error(predict_tos(fit1, base, temp_unit = "K"))

  # Gap: segment 1 ends at 10 but segment 2 starts at 20, so the propagation
  # would silently omit that stretch of the thermal history.
  gapped <- base
  gapped$t_end[1] <- 10
  gapped$t_start[2] <- 20
  expect_error(
    predict_tos(fit1, gapped, temp_unit = "K"),
    "contiguous"
  )

  # Overlap: the 12-24 window would be counted twice.
  overlapping <- base
  overlapping$t_start[2] <- 12
  expect_error(
    predict_tos(fit1, overlapping, temp_unit = "K"),
    "contiguous"
  )

  # Zero-duration leading segments remain valid.
  with_zero <- data.frame(
    t_start = c(0, 0, 24),
    t_end = c(0, 24, 30),
    temp = c(278.15, 278.15, 298.15),
    batch = "Combi_F1-PM1"
  )
  expect_no_error(predict_tos(fit1, with_zero, temp_unit = "K"))

  # Rows need not be supplied in order; sorting happens internally.
  expect_equal(
    predict_tos(fit1, base[c(2, 1), ], temp_unit = "K")$fit,
    predict_tos(fit1, base, temp_unit = "K")$fit
  )

  # Fractional segment boundaries must not trip the tolerance.
  fractional <- data.frame(
    t_start = c(0, 0.1, 0.3),
    t_end = c(0.1, 0.3, 0.7),
    temp = c(278.15, 298.15, 278.15),
    batch = "Combi_F1-PM1"
  )
  expect_no_error(predict_tos(fit1, fractional, temp_unit = "K"))
})


test_that("predict_tos filters by var for single-response models.", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  fit1 <- fit_kinetic(df,
    kinetic_order = 1,
    trend = "decreasing",
    intercept_cols = "batch"
  )

  schedule <- data.frame(
    t_start = c(0, 24),
    t_end = c(24, 30),
    temp = c(278.15, 298.15),
    batch = "Combi_F1-PM1"
  )

  unfiltered <- predict_tos(fit1, schedule, temp_unit = "K")

  # fit_kinetic labels its single response "y"; requesting it explicitly must
  # return those rows rather than silently filtering everything out.
  filtered <- predict_tos(fit1, schedule, var = "y", temp_unit = "K")
  expect_equal(nrow(filtered), nrow(unfiltered))
  expect_equal(unique(filtered$variable), "y")
  expect_equal(filtered$fit, unfiltered$fit)

  expect_error(predict_tos(fit1, schedule, var = "yZ", temp_unit = "K"))
})


test_that("predict_tos rejects schedule levels absent from the fit.", {
  data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  fit1 <- fit_kinetic(df,
    kinetic_order = 1,
    trend = "decreasing",
    intercept_cols = "batch"
  )

  schedule <- data.frame(
    t_start = c(0, 24),
    t_end = c(24, 30),
    temp = c(278.15, 298.15),
    batch = "NotARealBatch"
  )

  # An unmapped level must not become NA and propagate into the prediction.
  expect_error(
    predict_tos(fit1, schedule, temp_unit = "K"),
    "unseen levels"
  )
})