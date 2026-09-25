# Temperature-expression parameter detection, and parameter/data name clashes.
#
# Two failure modes are covered:
#
#   1. Parameters were detected with fixed greps (grepl("\\bb\\s*\\*", ...)), so
#      an algebraically identical expression with the operands swapped, or any
#      renamed parameter, was dropped from param_names and never got a start
#      value.
#   2. The residual-function fitters evaluate with
#      eval_tidy(expr, data = data, env = penv); the data mask outranks `env`,
#      so a data column sharing a parameter's name silently shadows it.

# ── fixtures ───────────────────────────────────────────────────────────────
sim_humidity <- function(kref = 0.02, Ea = 80000, b = 1.5) {
  g <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15),
                   rh = c(0.4, 0.6, 0.75))
  k <- kref * exp(-Ea / 8.314 * (1 / g$tempk - 1 / 298.15)) *
    exp(b * (g$rh - 0.6))
  g$y  <- 100 * exp(-k * g$time)
  g$yA <- g$y
  g$yB <- 2 + (100 - g$yA)
  g$yC <- 0
  g
}

as_temp_expr <- function(f, rh_col = "rh", rh_ref = 0.6) {
  list(formula = f, R = 8.314, trefk = 298.15, temp_var = "tempk",
       log_Ea = FALSE, rh_col = rh_col, rh_ref = rh_ref)
}

ARR_SWAPPED <- "exp(Ea / 8.314 * (1/298.15 - 1/tempk)) * exp((rh - 0.6) * b)"
ARR_RENAMED <- "exp(Eact / 8.314 * (1/298.15 - 1/tempk)) * exp(b * (rh - 0.6))"

# ── 1. parameter detection ─────────────────────────────────────────────────

test_that("temp_expr_params is insensitive to operand order", {
  canonical <- as.character(expr_arrhenius_humidity()$formula)

  expect_setequal(temp_expr_params(canonical,   "tempk", "rh"), c("Ea", "b"))
  expect_setequal(temp_expr_params(ARR_SWAPPED, "tempk", "rh"), c("Ea", "b"))
})

test_that("temp_expr_params finds renamed and log-scale parameters", {
  expect_setequal(temp_expr_params(ARR_RENAMED, "tempk", "rh"), c("Eact", "b"))
  expect_equal(temp_expr_params("exp(exp(Ea_log) / 8.314 * (1/298.15 - 1/tempk))",
                                "tempk"), "Ea_log")
  expect_equal(temp_expr_params("exp(Ea / 8.314 * (1/298.15 - 1/tempk))",
                                "tempk"), "Ea")
})

test_that("temp_expr_params excludes covariates, temp_var, literals and funs", {
  p <- temp_expr_params(
    "exp(Ea / 8.314 * (1/298.15 - 1/TK)) * sqrt(abs(log(hum))) * pi",
    temp_var = "TK", covar_cols = "hum")
  expect_equal(p, "Ea")

  expect_equal(temp_expr_params(NULL, "tempk"), character())
  expect_equal(temp_expr_params(mean, "tempk"), character())
})

test_that("a column named like a parameter does not reclassify it as data", {
  # Deciding "is this a parameter?" by "is there a column of that name?" would
  # be circular and would drop the parameter instead of reporting the clash.
  expect_setequal(temp_expr_params(ARR_SWAPPED, "tempk", "rh"), c("Ea", "b"))
})

test_that("build_kinetic collects swapped and renamed temp_expr parameters", {
  g <- sim_humidity()

  b1 <- build_kinetic(g, kinetic_order = 1, trend = "decreasing",
                      temp_fml = ARR_SWAPPED, rh_col = "rh", rh_ref = 0.6)
  expect_true(all(c("Ea", "b") %in% b1$param_names))

  b2 <- build_kinetic(g, kinetic_order = 1, trend = "decreasing",
                      temp_fml = ARR_RENAMED, rh_col = "rh", rh_ref = 0.6)
  expect_true(all(c("Eact", "b") %in% b2$param_names))
})

test_that("build_singledeg agrees with build_parallel on the same expression", {
  d <- sim_humidity()

  bs <- build_singledeg(d, temp_fml = ARR_SWAPPED, rh_col = "rh", rh_ref = 0.6)
  bp <- build_parallel(d,  temp_fml = ARR_SWAPPED, rh_col = "rh", rh_ref = 0.6)

  expect_true(all(c("Ea", "b") %in% bs$param_names))
  # parallel suffixes per branch
  expect_true(all(c("EaB", "bB", "EaC", "bC") %in% bp$param_names))
})

test_that("fit_kinetic recovers truth from a swapped-operand expression", {
  g <- sim_humidity()
  m <- fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                   temp_expr = as_temp_expr(ARR_SWAPPED),
                   start = c(Ea = 80000, b = 0, y0 = 100, k_0 = 0.01))
  cc <- coef(m)
  expect_equal(unname(cc[["Ea"]]),  80000, tolerance = 1e-4)
  expect_equal(unname(cc[["b"]]),   1.5,   tolerance = 1e-4)
  expect_equal(unname(cc[["y0"]]),  100,   tolerance = 1e-4)
  expect_equal(unname(cc[["k_0"]]), 0.02,  tolerance = 1e-4)
})

test_that("fit_kinetic recovers truth from a renamed activation energy", {
  g <- sim_humidity()
  m <- fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                   temp_expr = as_temp_expr(ARR_RENAMED),
                   start = c(Eact = 80000, b = 0, y0 = 100, k_0 = 0.01))
  expect_equal(unname(coef(m)[["Eact"]]), 80000, tolerance = 1e-4)
  expect_equal(unname(coef(m)[["b"]]),    1.5,   tolerance = 1e-4)
})

# ── 2. parameter / data-column clashes ─────────────────────────────────────

test_that("fit_singledeg errors when a column shadows a parameter", {
  g  <- sim_humidity()
  te <- expr_arrhenius_humidity()

  # baseline: the fit is exact without the stray column
  m <- fit_singledeg(g, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                     yA_var = "yA", yB_var = "yB", temp_expr = te)
  expect_equal(unname(coef(m)[["b"]]),  1.5,   tolerance = 1e-4)
  expect_equal(unname(coef(m)[["Ea"]]), 80000, tolerance = 1e-4)

  for (col in c("b", "Ea", "k_0", "A0", "B0")) {
    bad <- g
    bad[[col]] <- 99
    expect_error(
      fit_singledeg(bad, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                    yA_var = "yA", yB_var = "yB", temp_expr = te),
      "collide with model parameter name",
      info = paste("column", col)
    )
  }
})

test_that("fit_kinetic and fit_parallel also reject shadowing columns", {
  g  <- sim_humidity()
  te <- expr_arrhenius_humidity()

  bad <- g; bad$Ea <- 99
  expect_error(fit_kinetic(bad, kinetic_order = 1, trend = "decreasing",
                           temp_expr = te),
               "collide with model parameter name")

  bad2 <- g; bad2$A0 <- 99
  expect_error(build_parallel(bad2, temp_fml = ARR_SWAPPED,
                              rh_col = "rh", rh_ref = 0.6),
               "collide with model parameter name")
})

test_that("the clash error names the offending column", {
  g <- sim_humidity(); g$b <- 99
  msg <- tryCatch(
    fit_singledeg(g, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                  yA_var = "yA", yB_var = "yB",
                  temp_expr = expr_arrhenius_humidity()),
    error = conditionMessage)
  expect_match(msg, "b")
  expect_match(msg, "mask the parameter")
})

test_that("intercept and rate parameters are covered by the clash check", {
  g <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15),
                   batch = c("A", "B"), stringsAsFactors = FALSE)
  k <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
  g$y <- 100 * exp(-k * g$time)

  # y0_batch_A is a generated intercept parameter name
  bad <- g
  bad$y0_batch_A <- 1
  expect_error(fit_kinetic(bad, kinetic_order = 1, trend = "decreasing",
                           intercept_cols = "batch"),
               "collide with model parameter name")
})

test_that("clean data with a plain rh covariate still fits", {
  g <- sim_humidity()
  expect_silent({
    m <- fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                     temp_expr = expr_arrhenius_humidity())
  })
  expect_equal(unname(coef(m)[["b"]]), 1.5, tolerance = 1e-4)
})
