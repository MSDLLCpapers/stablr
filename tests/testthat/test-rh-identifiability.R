# The humidity term is exp(b * (rh - rh_ref)). When rh never varies that factor
# is a single constant multiplying the rate, so b is perfectly aliased with the
# rate parameters and is not identifiable. Optimisers do not fail on this --
# they leave b at its seed. fit_sb returned b = 0.1 (its hard-coded start) with
# ZERO warnings; fit_kinetic returned b = 0 with NaN standard errors. Either
# could be published as a measured moisture sensitivity.

make_curve <- function(rh, Ea = 80000, R = 8.314, k0 = 0.01, b = 2, y0 = 100) {
  d <- expand.grid(time = c(0, 3, 6, 12),
                   tempk = c(298.15, 313.15),
                   rh = rh)
  d$y <- y0 * exp(-k0 *
                    exp(-Ea / R * (1 / d$tempk - 1 / 298.15)) *
                    exp(b * (d$rh - 0.6)) * d$time)
  d
}

d_const <- make_curve(0.6)                 # single humidity chamber
d_vary  <- make_curve(c(0.4, 0.6, 0.75))   # identifiable

# Collect every warning a fit emits. A constant-humidity fit is genuinely
# singular, so the optimiser adds its own gradient warnings; asserting on the
# collected set keeps those from masking (or being mistaken for) the
# identifiability warning under test.
warnings_from <- function(expr) {
  ws <- character()
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

expect_id_warning <- function(expr) {
  ws <- warnings_from(expr)
  expect_gte(length(grep("not identifiable", ws)), 1)
}

expect_no_id_warning <- function(expr) {
  ws <- warnings_from(expr)
  expect_length(grep("not identifiable", ws), 0)
}

test_that("fit_sb warns that b is unidentifiable under constant humidity", {
  expect_id_warning(
    fit_sb(d_const, trend = "decreasing", Yinf = 0, temp_expr = expr_arrhenius_humidity())
  )
})

test_that("fit_sb no longer reports its hard-coded seed silently", {
  # The regression: b came back as exactly 0.1, the seed in fit_sb, with no
  # warning of any kind. The estimate is still unidentifiable -- what changed
  # is that the user is now told.
  ws <- warnings_from(
    fit_sb(d_const, trend = "decreasing", Yinf = 0, temp_expr = expr_arrhenius_humidity())
  )
  wmsg <- grep("not identifiable", ws, value = TRUE)[1]
  expect_match(wmsg, "`b`")
  expect_match(wmsg, "starting value, not an estimate")
  expect_match(wmsg, "0\\.6")

  fit <- suppressWarnings(
    fit_sb(d_const, trend = "decreasing", Yinf = 0, temp_expr = expr_arrhenius_humidity())
  )
  # Documents that the returned b is still the seed, hence the warning.
  expect_equal(unname(coef(fit)[["b"]]), 0.1, tolerance = 1e-8)
})

test_that("fit_kinetic warns under constant humidity", {
  expect_id_warning(
    fit_kinetic(d_const, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = expr_arrhenius_humidity())
  )
})

test_that("fit_singledeg warns under constant humidity", {
  d <- d_const
  d$yA <- d$y
  d$yB <- 100 - d$y

  expect_id_warning(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  yA_var = "yA", yB_var = "yB",
                  temp_expr = expr_arrhenius_humidity())
  )
})

test_that("fit_parallel warns under constant humidity", {
  d <- d_const
  d$yA <- d$y
  d$yB <- (100 - d$y) * 0.6
  d$yC <- (100 - d$y) * 0.4

  expect_id_warning(
    fit_parallel(d, kinetic_order = 1,
                 yA_var = "yA", yB_var = "yB", yC_var = "yC",
                 temp_expr = expr_arrhenius_humidity())
  )
})

test_that("varying humidity does not trigger the identifiability warning", {
  expect_no_id_warning(fit_kinetic(d_vary, kinetic_order = 1,
                                   trend = "decreasing", Yinf = 0,
                                   temp_expr = expr_arrhenius_humidity()))
  expect_no_id_warning(fit_sb(d_vary, trend = "decreasing", Yinf = 0,
                              temp_expr = expr_arrhenius_humidity()))

  # And b is genuinely recovered when it is identifiable.
  fit <- suppressWarnings(
    fit_kinetic(d_vary, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                temp_expr = expr_arrhenius_humidity())
  )
  expect_equal(unname(coef(fit)[["b"]]), 2, tolerance = 1e-4)
})

test_that("models with no humidity term are unaffected", {
  expect_no_id_warning(
    fit_kinetic(d_vary, kinetic_order = 1, trend = "decreasing", Yinf = 0)
  )
})

test_that("check_rh_identifiable flags near-constant humidity too", {
  # Floating-point equality would miss this, but it is just as unidentifiable.
  expect_warning(
    check_rh_identifiable(data.frame(rh = c(0.6, 0.6 + 1e-12)), "rh"),
    "not identifiable"
  )
})

test_that("check_rh_identifiable is a no-op on absent or unusable input", {
  expect_silent(check_rh_identifiable(d_const, NULL))
  expect_silent(check_rh_identifiable(d_const, "not_a_column"))
  expect_silent(check_rh_identifiable(data.frame(rh = 0.6), "rh"))
  expect_silent(check_rh_identifiable(data.frame(rh = numeric(0)), "rh"))
  expect_silent(check_rh_identifiable(d_vary, "rh"))
})

test_that("check_rh_identifiable stays silent when b is not a model parameter", {
  # A custom temperature expression may reference rh without a `b` coefficient;
  # there is then no unidentifiable parameter to warn about.
  expect_silent(
    check_rh_identifiable(d_const, "rh", "exp(Ea / 8.314 * (1/298.15 - 1/tempk))")
  )
  expect_warning(
    check_rh_identifiable(d_const, "rh",
                          "exp(Ea / 8.314 * (1/298.15 - 1/tempk)) * exp(b * (rh - 0.6))"),
    "not identifiable"
  )
})
