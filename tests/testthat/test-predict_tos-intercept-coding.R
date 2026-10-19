two_factor_data <- function(orient_levels = c("inv", "upr"),
                            effects = c(inv = 0, upr = -5)) {
  g <- expand.grid(batch = c("A", "B"), orient = orient_levels,
                   time = c(0, 3, 6, 12), temp_c = c(25, 40),
                   stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- (c(A = 100, B = 97)[g$batch] + effects[g$orient]) * exp(-k * g$time)
  g
}

two_factor_fit <- function(...) {
  suppressWarnings(fit_kinetic(two_factor_data(...), kinetic_order = 1,
                               trend = "decreasing",
                               intercept_cols = c("batch", "orient")))
}

# Evaluate the fitted formula's RHS directly in base R at the target
# point. Independent of every expression the composite builder assembles.
oracle <- function(fit, ...) {
  env <- list2env(as.list(stats::coef(fit)), parent = baseenv())
  eval(stats::formula(fit)[[3]], envir = data.frame(..., stringsAsFactors = FALSE),
       enclos = env)
}

test_that("predict_tos matches the fitted model at a NON-reference level", {
  fit <- two_factor_fit()
  # The fit recovers truth, so the defect cannot be blamed on estimation.
  expect_equal(unname(coef(fit)["y0_orient_upr"]), -5, tolerance = 1e-6)

  p <- predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15,
                                   batch = "A", orient = "upr"),
                   temp_unit = "K")
  p <- tos_segments(p)   # drop the time-0 baseline row

  expect_equal(p$fit, oracle(fit, batch = "A", orient = "upr", time = 12,
                             tempk = 298.15),
               tolerance = 1e-8)
  # Previously 78.66278611: high by exactly the dropped y0_orient_upr = -5.
  expect_equal(p$fit, 74.72964680, tolerance = 1e-6)
})

test_that("predict_tos matches the fitted model at the reference level", {
  # This level was already correct and must stay correct.
  fit <- two_factor_fit()
  p <- predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15,
                                   batch = "A", orient = "inv"),
                   temp_unit = "K")
  p <- tos_segments(p)   # drop the time-0 baseline row
  expect_equal(p$fit, oracle(fit, batch = "A", orient = "inv", time = 12,
                             tempk = 298.15),
               tolerance = 1e-8)
})

test_that("the interval covers the fitted value on data with noise", {
  g <- two_factor_data()
  set.seed(42)
  g$y <- g$y + stats::rnorm(nrow(g), 0, 0.3)
  fit <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                                      intercept_cols = c("batch", "orient")))

  p <- predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15,
                                   batch = "A", orient = "upr"),
                   temp_unit = "K")
  p <- tos_segments(p)   # drop the time-0 baseline row
  target <- oracle(fit, batch = "A", orient = "upr", time = 12, tempk = 298.15)

  expect_equal(p$fit, target, tolerance = 1e-8)
  expect_gt(p$upr - p$lwr, 0)     # a real SE, so the interval has width
  expect_lt(p$lwr, target)
  expect_gt(p$upr, target)
})

test_that("every batch x orient combination matches the fitted model", {
  fit <- two_factor_fit()
  for (b in c("A", "B")) {
    for (o in c("inv", "upr")) {
      p <- predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 313.15,
                                       batch = b, orient = o),
                       temp_unit = "K")
      p <- tos_segments(p)   # drop the time-0 baseline row
      expect_equal(p$fit, oracle(fit, batch = b, orient = o, time = 12,
                                 tempk = 313.15),
                   tolerance = 1e-8,
                   info = paste("batch", b, "orient", o))
    }
  }
})

test_that("a multi-segment schedule matches the fitted model in every segment", {
  # Every cumulative sub-schedule must evaluate the full fitted intercept, so
  # check a non-reference level of the second intercept column at both segments.
  fit <- two_factor_fit(orient_levels = c("inv", "upr", "xyz"),
                        effects = c(inv = 0, upr = -5, xyz = -8))
  expect_equal(unname(coef(fit)["y0_orient_upr"]), -5, tolerance = 1e-6)

  sch <- data.frame(t_start = c(0, 6), t_end = c(6, 12),
                    temp = c(298.15, 298.15), batch = c("A", "A"),
                    orient = c("upr", "upr"), stringsAsFactors = FALSE)
  p <- predict_tos(fit, sch, temp_unit = "K")
  p <- tos_segments(p)   # drop the time-0 baseline row

  expect_equal(p$fit[1], oracle(fit, batch = "A", orient = "upr", time = 6,
                                tempk = 298.15), tolerance = 1e-8)
  expect_equal(p$fit[2], oracle(fit, batch = "A", orient = "upr", time = 12,
                                tempk = 298.15), tolerance = 1e-8)
})

test_that("a grouping column that changes between segments is rejected (V01-16)", {
  fit <- two_factor_fit(orient_levels = c("inv", "upr", "xyz"),
                        effects = c(inv = 0, upr = -5, xyz = -8))

  sch <- data.frame(t_start = c(0, 6), t_end = c(6, 12),
                    temp = c(298.15, 298.15), batch = c("A", "A"),
                    orient = c("upr", "xyz"), stringsAsFactors = FALSE)
  expect_error(predict_tos(fit, sch, temp_unit = "K"),
               "orient.*vary across segments")

  # The message must name every offending column, not just the first.
  sch2 <- transform(sch, batch = c("A", "B"))
  expect_error(predict_tos(fit, sch2, temp_unit = "K"), "'batch', 'orient'")
})

test_that("a single intercept column is unaffected", {
  g <- two_factor_data()
  fit <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                                      intercept_cols = "batch"))
  p <- predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15,
                                   batch = "A"), temp_unit = "K")
  p <- tos_segments(p)   # drop the time-0 baseline row
  expect_equal(p$fit, oracle(fit, batch = "A", time = 12, tempk = 298.15),
               tolerance = 1e-8)
})

test_that("rate_cols are also evaluated from the fitted expression", {
  g <- two_factor_data()
  fit <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                                      intercept_cols = c("batch", "orient"),
                                      rate_cols = "batch"))
  p <- predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 313.15,
                                   batch = "B", orient = "upr"),
                   temp_unit = "K")
  p <- tos_segments(p)   # drop the time-0 baseline row
  expect_equal(p$fit, oracle(fit, batch = "B", orient = "upr", time = 12,
                             tempk = 313.15),
               tolerance = 1e-8)
})

# ── expr_intercept(): a single-level later factor is legitimate input ─────────

test_that("expr_intercept tolerates a later factor observed at one level", {
  one_term <- "((batch == 'A') * y0_batch_A)"

  one_row <- data.frame(batch = "A", orient = "upr", stringsAsFactors = FALSE)
  ip <- expr_intercept(one_row, c("batch", "orient"))
  expect_equal(as.character(ip$expr), one_term)
  expect_equal(as.character(ip$params), "y0_batch_A")

  two_rows <- data.frame(batch = c("A", "A"), orient = c("upr", "upr"),
                         stringsAsFactors = FALSE)
  expect_equal(as.character(expr_intercept(two_rows, c("batch", "orient"))$expr),
               one_term)
})

# ── the backstop assertion ───────────────────────────────────────────────────

test_that("check_composite_params errors when a fitted parameter is unreachable", {
  err <- expect_error(
    check_composite_params(c("((batch == 'A') * y0_batch_A)", "k_0"),
                           c("Ea", "y0_batch_A", "y0_orient_upr", "k_0"))
  )
  expect_match(conditionMessage(err), "y0_orient_upr")
  expect_match(conditionMessage(err), "Ea")
  # Must say why it cannot be caught by inspecting the output.
  expect_match(conditionMessage(err), "silently wrong")
})

test_that("check_composite_params is silent when every parameter is covered", {
  expect_silent(check_composite_params(
    c("((batch == 'A') * y0_batch_A) + ((orient == 'upr') * y0_orient_upr)",
      "k_0", "exp(Ea / 8.314 * (1/298.15 - 1/tempk))"),
    c("Ea", "y0_batch_A", "y0_orient_upr", "k_0")
  ))
  expect_silent(check_composite_params(character(), c("Ea")))
  expect_silent(check_composite_params("k_0", character()))
})

test_that("a fit without the captured expressions is not rebuilt", {
  # predict_tos() only ever reuses the fitted intercept/rate expressions; there is
  # no schedule-rebuild fallback to silently re-derive reference coding from.
  fit <- two_factor_fit()
  rm("int_part", envir = environment(stats::formula(fit)))

  expect_error(
    predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15,
                                batch = "A", orient = "upr"),
                temp_unit = "K"),
    "int_part"
  )
})
