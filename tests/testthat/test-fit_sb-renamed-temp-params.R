test_that("fit_sb reads parameters from the expression, like the closed-form fitters (V10-03)", {
  sb_kT <- function(tempk, k = 0.02, Ea = 80000)
    k * exp(Ea / 8.314 * (1 / 298.15 - 1 / tempk))

  set.seed(1)
  d <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15),
                   rh = c(0.4, 0.7))
  d$y <- 100 * exp(-sb_kT(d$tempk) * d$time * exp(-1.2 * (d$rh - 0.6))) +
         rnorm(nrow(d), 0, 0.2)

  # A renamed humidity coefficient (b -> bh) must fit, not error on a missing
  # object, and must not inject a phantom "b".
  te_h <- expr_arrhenius_humidity(rh_col = "rh", rh_ref = 0.6)
  te_h$formula <- gsub("\\bb\\b", "bh", as.character(te_h$formula), perl = TRUE)
  f_h <- suppressWarnings(fit_sb(d, temp_expr = te_h, trend = "decreasing", Yinf = 0))
  expect_true("bh" %in% names(f_h$par))
  expect_false("b" %in% names(f_h$par))

  # A renamed activation energy (Ea -> Eact) must fit and must NOT create a
  # phantom "Ea" beside it.
  d2 <- subset(d, rh == 0.4)
  te_e <- expr_arrhenius()
  te_e$formula <- gsub("\\bEa\\b", "Eact", as.character(te_e$formula), perl = TRUE)
  f_e <- suppressWarnings(fit_sb(d2, temp_expr = te_e, trend = "decreasing", Yinf = 0))
  expect_true("Eact" %in% names(f_e$par))
  expect_false("Ea" %in% names(f_e$par))

  # The renamed energy is still seeded and bounded as an activation energy, so it
  # recovers the truth rather than resting at a covariate's tight bounds; and it
  # agrees with fit_kinetic on the same expression.
  expect_gt(unname(f_e$par["Eact"]), 1e4)
  f_k <- suppressWarnings(
    fit_kinetic(d2, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = te_e))
  expect_equal(unname(f_e$par["Eact"]), unname(coef(f_k)["Eact"]),
               tolerance = 0.05)
})
