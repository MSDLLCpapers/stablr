test_that("predict_tos rejects an unseen grouping level on the fit_sb path", {
  # data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))
  df$batch <- paste("Batch", df$batch)
  df$pkg   <- paste("Container-", df$pkg, sep = "")

  fit <- fit_sb(df = df, trend = "decreasing", Yinf = 0, intercept_cols = "batch", rate_cols = "pkg")

  # fit_sb stores no sanitize_map, so predict_tos used to sanitize the schedule
  # and skip the unseen-level guard: an unfitted level matched no indicator and
  # every column of the returned row came back NA with no condition signalled.
  sched_bad <- data.frame(t_start = 0, t_end = 12, temp = 40,
                          batch = "NOT_A_BATCH", pkg = "Container-PM2")
  testthat::expect_error(
    predict_tos(fit, sched_bad, temp_unit = "C"),
    "unseen levels in column 'batch'"
  )

  # The same schedule errors on the fit_kinetic path; the two must agree.
  fit_k <- fit_kinetic(df = df, kinetic_order = 1, trend = "decreasing",
                       Yinf = 0, intercept_cols = "batch", rate_cols = "pkg")
  testthat::expect_error(
    predict_tos(fit_k, sched_bad, temp_unit = "C"),
    "unseen levels in column 'batch'"
  )

  # Levels that sanitize() rewrites must still be accepted and predict finitely.
  sched_ok <- data.frame(t_start = 0, t_end = 12, temp = 40,
                         batch = "Batch Combi_F1-PM1", pkg = "Container-PM2")
  res <- predict_tos(fit, sched_ok, temp_unit = "C")
  testthat::expect_true(all(is.finite(res$fit)))
})
