# Temperature and time are arithmetic inputs, but only fit_sb checked their type.
#
# predict_tos() applied `schedule$temp + 273.15` with no check, so a factor temp
# warned "'+' not meaningful for factors", produced NA, and then let the function
# run to completion and return a normal-looking one-row data frame with
# fit/se_mean/se_pred/lwr/upr all NA -- easy to lose in a loop over products.
# Under temp_unit = "K" the conversion does no arithmetic at all, so the factor
# survived into the Arrhenius term and warned about '/' once per numDeriv
# evaluation (18 times for a single segment), naming an operator even further
# from the cause. A CHARACTER temp raised "non-numeric argument to binary
# operator" instead, so two spellings of one mistake behaved differently.
#
# The same class of error reached the closed-form fitters as "missing values
# returned by 'fn' at starting values", which reads as a starting-value problem.
# (That is an error, so the earlier fit-time bug -- silent coercion to integer
# level codes, converging against temperatures of 1-3 Kelvin -- was already gone;
# only the diagnosis was missing.)
#
# Factor temps arise naturally: read.csv() with stringsAsFactors, a
# group_by()/summarise() table, or a factor-coded condition column. The check
# refuses rather than coercing, because as.numeric() on a factor yields level
# codes and as.character()-then-numeric would quietly rescue genuinely
# mislabelled data too.

num_fit <- function() {
  g <- expand.grid(batch = c("A", "B"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- c(A = 100, B = 97)[g$batch] * exp(-k * g$time)
  suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                               intercept_cols = "batch"))
}

num_data <- function() {
  g <- expand.grid(batch = c("A", "B"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- c(A = 100, B = 97)[g$batch] * exp(-k * g$time)
  g
}

sched <- function(temp, unit = "C") {
  data.frame(t_start = 0, t_end = 12, temp = temp, batch = "A",
             stringsAsFactors = FALSE)
}

test_that("predict_tos rejects a factor temp instead of returning an NA row", {
  fit <- num_fit()

  err <- expect_error(
    suppressWarnings(predict_tos(fit, sched(factor(25)), temp_unit = "C"))
  )
  expect_match(conditionMessage(err), "must be numeric")
  # Must name the column and its actual type.
  expect_match(conditionMessage(err), "'temp' (factor)", fixed = TRUE)
  expect_match(conditionMessage(err), "schedule")
})

test_that("predict_tos rejects a factor temp under temp_unit = 'K' too", {
  # The K branch does no arithmetic at conversion, so the factor used to reach
  # the Arrhenius term and warn about '/' repeatedly instead of '+' once.
  fit <- num_fit()
  err <- expect_error(
    suppressWarnings(predict_tos(fit, sched(factor(298.15)), temp_unit = "K"))
  )
  expect_match(conditionMessage(err), "'temp' (factor)", fixed = TRUE)
})

test_that("a factor temp produces no warning-and-NA path at all", {
  # The original signal was a base R warning plus a returned all-NA row. Neither
  # may happen now: the call must fail outright, with no bare '+' warning.
  fit <- num_fit()
  expect_error(predict_tos(fit, sched(factor(25)), temp_unit = "C"))

  w <- tryCatch(predict_tos(fit, sched(factor(25)), temp_unit = "C"),
                error = function(e) "errored",
                warning = function(w) conditionMessage(w))
  expect_identical(w, "errored")
})

test_that("factor and character temps are rejected the same way", {
  # The asymmetry was the tell: character errored, factor returned NA.
  fit <- num_fit()
  f_err <- expect_error(suppressWarnings(
    predict_tos(fit, sched(factor(25)), temp_unit = "C")))
  c_err <- expect_error(suppressWarnings(
    predict_tos(fit, sched("25"), temp_unit = "C")))

  expect_match(conditionMessage(f_err), "must be numeric")
  expect_match(conditionMessage(c_err), "must be numeric")
  expect_match(conditionMessage(c_err), "'temp' (character)", fixed = TRUE)
})

test_that("factor t_start and t_end are rejected as well", {
  # Both are differenced to get each segment's duration, so the same type error
  # applies; it previously surfaced as an NA-comparison failure in the t_end <
  # t_start guard.
  fit <- num_fit()
  s <- sched(25)
  s$t_end <- factor(12)
  err <- expect_error(suppressWarnings(predict_tos(fit, s, temp_unit = "C")))
  expect_match(conditionMessage(err), "'t_end' (factor)", fixed = TRUE)

  s2 <- sched(25)
  s2$t_start <- factor(0)
  err2 <- expect_error(suppressWarnings(predict_tos(fit, s2, temp_unit = "C")))
  expect_match(conditionMessage(err2), "'t_start' (factor)", fixed = TRUE)
})

test_that("a numeric schedule still predicts, in both units", {
  # No false positives, and the value the report expected.
  fit <- num_fit()
  p_c <- tos_segments(predict_tos(fit, sched(25), temp_unit = "C"))
  p_k <- tos_segments(predict_tos(fit, sched(298.15), temp_unit = "K"))

  expect_equal(p_c$fit, 78.662786, tolerance = 1e-6)
  expect_equal(p_c$fit, p_k$fit, tolerance = 1e-10)
  expect_false(any(is.na(c(p_c$fit, p_c$lwr, p_c$upr))))
  # Integer columns are numeric and must not be rejected.
  expect_no_error(predict_tos(fit, sched(25L), temp_unit = "C"))
})

test_that("all four fitters name the column for a factor temperature", {
  # fit_kinetic previously said "missing values returned by 'fn' at starting
  # values"; fit_sb already had its own message. All four must agree now.
  d <- num_data()
  d$tempk <- factor(d$tempk)

  k_err <- expect_error(suppressWarnings(fit_kinetic(
    d, kinetic_order = 1, trend = "decreasing", intercept_cols = "batch")))
  expect_match(conditionMessage(k_err), "'tempk' (factor)", fixed = TRUE)

  s_err <- expect_error(suppressWarnings(fit_sb(
    d, trend = "decreasing", Yinf = 0, start = c(Ea = 80000, n = 1, m = 0, k_0 = 0.02, y0_0 = 100))))
  expect_match(conditionMessage(s_err), "'tempk' (factor)", fixed = TRUE)

  # A factor TIME column is the same class of error.
  d2 <- num_data()
  d2$time <- factor(d2$time)
  t_err <- expect_error(suppressWarnings(fit_kinetic(
    d2, kinetic_order = 1, trend = "decreasing", intercept_cols = "batch")))
  expect_match(conditionMessage(t_err), "'time' (factor)", fixed = TRUE)
})

test_that("grouping columns may still be factors", {
  # Only the arithmetic columns are constrained; a factor batch is ordinary and
  # must keep working, at fit time and in a schedule.
  d <- num_data()
  d$batch <- factor(d$batch)
  fit <- suppressWarnings(fit_kinetic(d, kinetic_order = 1,
                                      trend = "decreasing",
                                      intercept_cols = "batch"))
  s <- data.frame(t_start = 0, t_end = 12, temp = 25, batch = factor("A"))
  expect_no_error(predict_tos(fit, s, temp_unit = "C"))
})

# ── the validator in isolation ────────────────────────────────────────────────

test_that("check_numeric_cols names every offending column and its class", {
  df <- data.frame(a = 1, b = factor(2), c = "3", d = 4,
                   stringsAsFactors = FALSE)
  err <- expect_error(check_numeric_cols(df, c("a", "b", "c", "d")))
  expect_match(conditionMessage(err), "'b' (factor)", fixed = TRUE)
  expect_match(conditionMessage(err), "'c' (character)", fixed = TRUE)
  # Valid columns must not be listed.
  expect_no_match(conditionMessage(err), "'a'", fixed = TRUE)
  expect_no_match(conditionMessage(err), "'d'", fixed = TRUE)
})

test_that("check_numeric_cols only mentions level codes when a factor is present", {
  df <- data.frame(x = "3", stringsAsFactors = FALSE)
  err <- expect_error(check_numeric_cols(df, "x"))
  expect_no_match(conditionMessage(err), "level codes")
})

test_that("check_numeric_cols is silent on numeric, absent or empty input", {
  expect_silent(check_numeric_cols(data.frame(x = 1, y = 2L), c("x", "y")))
  expect_silent(check_numeric_cols(data.frame(x = 1), "nonexistent"))
  expect_silent(check_numeric_cols(data.frame(x = 1), character()))
  expect_silent(check_numeric_cols(data.frame(x = 1), NULL))
  expect_silent(check_numeric_cols(data.frame(x = 1), ""))
  # NA in a numeric column is check_structural_nas()' business, not this one.
  expect_silent(check_numeric_cols(data.frame(x = c(1, NA)), "x"))
})
