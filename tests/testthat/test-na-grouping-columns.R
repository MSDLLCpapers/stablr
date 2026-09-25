# Only the OUTCOME was checked for NA; nothing checked the columns that define
# the model structure. Each fitter then failed differently and none said "NA":
#
#   fit_kinetic  -- nls() dropped the NA rows via na.action, fitting 8 of 12 rows
#     with df.residual = 2, while STILL emitting a phantom `k_batch_NA`
#     parameter for the dropped level. sigma came back 4.2e-10 off a subset.
#   fit_sb       -- did NOT drop them: an NA group label matches no y0/k, so
#     those rows returned residuals equal to the whole observed value (55-89).
#     Deviance was 2.1e4 instead of ~1.6e-27 and Ea was dragged from 105956 to
#     116923, with no error, no warning and no NA in the output to notice.
#   fit_singledeg / fit_parallel -- aborted with "missing values returned by
#     'fn' at starting values", blaming the starting values rather than the data.
#
# NA grouping values are highly plausible from an upstream join.

na_data <- function(na_time = 6) {
  d <- expand.grid(time = c(0, 3, 6), tempk = c(298.15, 313.15),
                   batch = c("B1", "B2"), stringsAsFactors = FALSE)
  d$y  <- 100 * exp(-0.01 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15)) * d$time)
  d$yA <- d$y
  d$yB <- 100 - d$y
  d$yB2 <- (100 - d$y) * 0.6
  d$yC <- (100 - d$y) * 0.4
  if (!is.null(na_time)) d$batch[d$time == na_time] <- NA_character_
  d
}

test_that("fit_kinetic rejects NA in a grouping column", {
  # Previously: no error, 8 of 12 rows fitted, plus a phantom k_batch_NA.
  d <- na_data()
  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch", rate_cols = "batch"),
    "NA values in df column"
  )
})

test_that("fit_sb rejects NA in a grouping column", {
  # Previously: converged with deviance 2.1e4 built from phantom residuals.
  d <- na_data()
  expect_error(
    fit_sb(d, intercept_cols = "batch", rate_cols = "batch", trend = "decreasing", Yinf = 0),
    "NA values in df column"
  )
})

test_that("fit_singledeg rejects NA in a grouping column", {
  d <- na_data()
  expect_error(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  yA_var = "yA", yB_var = "yB",
                  intercept_cols = "batch", rate_cols = "batch"),
    "NA values in df column"
  )
})

test_that("fit_parallel rejects NA in a grouping column", {
  d <- na_data()
  expect_error(
    fit_parallel(d, kinetic_order = 1, yA_var = "yA", yB_var = "yB2",
                 yC_var = "yC", intercept_cols = "batch", rate_cols = "batch"),
    "NA values in df column"
  )
})

test_that("the error names the column, the row count and a remedy", {
  d <- na_data()
  msg <- tryCatch(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch", rate_cols = "batch"),
    error = conditionMessage
  )
  expect_match(msg, "'batch' \\(4 rows\\)")
  expect_match(msg, "unknown")   # the suggested explicit label
  # The old failure mode named nothing at all; the message must explain why the
  # rows cannot simply be ignored.
  expect_match(msg, "no fitted parameter")
})

test_that("row count is pluralised correctly", {
  d <- na_data(na_time = NULL)
  d$batch[1] <- NA_character_
  msg <- tryCatch(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch"),
    error = conditionMessage
  )
  expect_match(msg, "\\(1 row\\)")
  expect_no_match(msg, "\\(1 rows\\)")
})

test_that("NA is rejected in a rate_cols-only column", {
  # The column need not be an intercept column to corrupt the fit.
  d <- na_data()
  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                rate_cols = "batch"),
    "NA values in df column"
  )
})

test_that("NA in time, temperature and humidity is rejected", {
  d <- na_data(na_time = NULL)

  d_t <- d; d_t$time[2] <- NA_real_
  expect_error(fit_kinetic(d_t, kinetic_order = 1, trend = "decreasing",
                           intercept_cols = "batch"), "'time'")

  d_k <- d; d_k$tempk[2] <- NA_real_
  expect_error(fit_kinetic(d_k, kinetic_order = 1, trend = "decreasing",
                           intercept_cols = "batch"), "'tempk'")

  d_h <- d
  d_h$rh <- rep(c(0.4, 0.75), length.out = nrow(d_h))
  d_h$rh[2] <- NA_real_
  expect_error(
    fit_kinetic(d_h, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch",
                temp_expr = expr_arrhenius_humidity()),
    "'rh'"
  )
})

test_that("several affected columns are all reported", {
  d <- na_data()
  d$tempk[1] <- NA_real_
  msg <- tryCatch(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch", rate_cols = "batch"),
    error = conditionMessage
  )
  expect_match(msg, "batch")
  expect_match(msg, "tempk")
})

test_that("a column is named once even when used for both roles", {
  # 'batch' is passed as both intercept_cols and rate_cols; the old fit_sb
  # message listed it twice ("found NA in: batch, batch").
  d <- na_data()
  msg <- tryCatch(
    fit_sb(d, intercept_cols = "batch", rate_cols = "batch", trend = "decreasing", Yinf = 0),
    error = conditionMessage
  )
  expect_equal(length(gregexpr("'batch'", msg)[[1]]), 1L)
})

test_that("an NA in an unused column does not block the fit", {
  # Only columns that define the model matter; the check must not be a blanket
  # complete.cases() over the whole data frame.
  d <- na_data(na_time = NULL)
  d$site <- NA_character_
  expect_silent(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch", rate_cols = "batch")
  )
})

test_that("clean data is unaffected and the outcome check still fires", {
  d <- na_data(na_time = NULL)
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     intercept_cols = "batch", rate_cols = "batch")
  # All 12 rows used -- no silent drop.
  expect_equal(nobs(fit), 12L)

  d$y[1] <- NA_real_
  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch", rate_cols = "batch"),
    "NAs in the outcome"
  )
})

test_that("check_structural_nas ignores absent and unnamed columns", {
  d <- na_data(na_time = NULL)
  expect_silent(check_structural_nas(d, "not_a_column"))
  expect_silent(check_structural_nas(d, c("", NA_character_)))
  expect_silent(check_structural_nas(d, character()))
})

test_that("predict names an NA grouping value instead of calling it unseen", {
  # Previously: "unseen levels in column 'batch': NA", which sends the user
  # looking for a level they never fitted rather than at a failed join.
  d <- na_data(na_time = NULL)
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     intercept_cols = "batch", rate_cols = "batch")
  msg <- tryCatch(
    predict(fit, newdata = data.frame(batch = NA_character_, time = 3,
                                      tempk = 313.15)),
    error = conditionMessage
  )
  expect_match(msg, "NA values in newdata column")
  expect_match(msg, "join")
  expect_no_match(msg, "unseen")
})

test_that("predict on clean newdata still works", {
  d <- na_data(na_time = NULL)
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     intercept_cols = "batch", rate_cols = "batch")
  p <- predict(fit, newdata = data.frame(batch = "B1", time = 3,
                                         tempk = 313.15))
  expect_true(is.finite(p))
  expect_lt(p, 100)
})
