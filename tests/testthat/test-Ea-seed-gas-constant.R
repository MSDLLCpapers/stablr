# The humidity fitters need an explicit Ea seed, and that seed used to be a
# hard-coded 80000 -- a J/mol value. Fitting the same physics on the cal/mol
# scale (R = 1.987, Ea ~ 19120) started 4.2x too high and diverged to a
# negative activation energy and a negative rate constant, warning only about a
# "singular gradient matrix". The seed must now scale with the user's R.

# Same degradation curve expressed on either energy scale.
Ea_si  <- 80000
R_si   <- 8.314
R_cal  <- 1.987
Ea_cal <- Ea_si / R_si * R_cal   # 19119.56 cal/mol

make_curve <- function(Ea, R, k0 = 0.01, b = 2, y0 = 100, rh_ref = 0.6) {
  d <- expand.grid(time = c(0, 3, 6, 12),
                   tempk = c(298.15, 313.15),
                   rh = c(0.4, 0.6, 0.75))
  d$y <- y0 * exp(-k0 *
                    exp(-Ea / R * (1 / d$tempk - 1 / 298.15)) *
                    exp(b * (d$rh - rh_ref)) * d$time)
  d
}

test_that("default_Ea_seed scales the seed with the gas constant", {
  expect_equal(default_Ea_seed(R_si), Ea_si)
  expect_equal(default_Ea_seed(R_cal), Ea_cal)
  # Seeded temperature sensitivity Ea/R is unit-free.
  expect_equal(default_Ea_seed(R_cal) / R_cal, default_Ea_seed(R_si) / R_si)

  expect_error(default_Ea_seed(0), "positive")
  expect_error(default_Ea_seed(-8.314), "positive")
  expect_error(default_Ea_seed(NA_real_), "positive")
  expect_error(default_Ea_seed(c(8.314, 1.987)), "single")
})

test_that("fit_kinetic recovers cal/mol humidity parameters without a hand seed", {
  d <- make_curve(Ea_cal, R_cal)

  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                     R = R_cal, temp_expr = expr_arrhenius_humidity(R = R_cal))
  cf <- coef(fit)

  expect_equal(unname(cf[["Ea"]]), Ea_cal, tolerance = 1e-4)
  expect_equal(unname(cf[["b"]]), 2, tolerance = 1e-4)
  expect_equal(unname(cf[["k_0"]]), 0.01, tolerance = 1e-4)
  # The old failure signature: a physically impossible negative rate constant.
  expect_gt(cf[["k_0"]], 0)
  expect_gt(cf[["Ea"]], 0)
})

test_that("fit_kinetic SI humidity behaviour is unchanged", {
  d <- make_curve(Ea_si, R_si)

  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                     temp_expr = expr_arrhenius_humidity())
  cf <- coef(fit)

  expect_equal(unname(cf[["Ea"]]), Ea_si, tolerance = 1e-4)
  expect_equal(unname(cf[["b"]]), 2, tolerance = 1e-4)
  expect_equal(unname(cf[["k_0"]]), 0.01, tolerance = 1e-4)
})

test_that("SI and cal/mol fits agree on the same physics", {
  d_si  <- make_curve(Ea_si, R_si)
  d_cal <- make_curve(Ea_cal, R_cal)

  fit_si <- fit_kinetic(d_si, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                        temp_expr = expr_arrhenius_humidity())
  fit_cal <- fit_kinetic(d_cal, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                         R = R_cal, temp_expr = expr_arrhenius_humidity(R = R_cal))

  # Ea/R is the scale-invariant quantity.
  expect_equal(coef(fit_si)[["Ea"]] / R_si,
               coef(fit_cal)[["Ea"]] / R_cal, tolerance = 1e-4)
  expect_equal(coef(fit_si)[["b"]], coef(fit_cal)[["b"]], tolerance = 1e-4)
  expect_equal(coef(fit_si)[["k_0"]], coef(fit_cal)[["k_0"]], tolerance = 1e-4)
})

test_that("the log_Ea humidity seed also scales with the gas constant", {
  d <- make_curve(Ea_cal, R_cal)

  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                     R = R_cal,
                     temp_expr = expr_arrhenius_humidity(R = R_cal, log_Ea = TRUE))
  cf <- coef(fit)

  expect_equal(exp(unname(cf[["Ea_log"]])), Ea_cal, tolerance = 1e-3)
  expect_gt(cf[["k_0"]], 0)
})

test_that("fit_singledeg seeds Ea on the user's energy scale", {
  d <- make_curve(Ea_cal, R_cal)
  d$yA <- d$y
  d$yB <- 100 - d$y

  fit <- fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                       yA_var = "yA", yB_var = "yB",
                       temp_expr = expr_arrhenius_humidity(R = R_cal))
  cf <- coef(fit)

  expect_equal(unname(cf[["Ea"]]), Ea_cal, tolerance = 1e-2)
  expect_gt(cf[["Ea"]], 0)
  expect_true(all(cf[grep("^k_", names(cf))] > 0))
})

test_that("fit_parallel seeds EaB/EaC on the user's energy scale", {
  d <- make_curve(Ea_cal, R_cal)
  d$yA <- d$y
  d$yB <- (100 - d$y) * 0.6
  d$yC <- (100 - d$y) * 0.4

  fit <- fit_parallel(d, kinetic_order = 1,
                      yA_var = "yA", yB_var = "yB", yC_var = "yC",
                      temp_expr = expr_arrhenius_humidity(R = R_cal))
  cf <- coef(fit)

  for (nm in intersect(c("EaB", "EaC"), names(cf))) {
    expect_equal(unname(cf[[nm]]), Ea_cal, tolerance = 0.05 * Ea_cal)
    expect_gt(cf[[nm]], 0)
  }
})
