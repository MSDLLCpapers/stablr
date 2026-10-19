skip_on_cran()

# fit_sb evaluates its Arrhenius factor inside an ODE, in an environment it
# populates by hand.  It used to bind the temperature under the literal "tempk"
# while expr_arrhenius() had written the caller's `temp_var` into the formula, so
# any name but the lowercase default died with "object '<temp_var>' not found".
# The fix is to bind under `temp_var`; these tests pin the behaviour by checking
# that the *name* of the column cannot change the numbers.

# One data frame carrying the same Kelvin temperatures under several names, so a
# fit under each is comparable coefficient-for-coefficient.
sb_data <- function() {
  g  <- expand.grid(time = c(0, 3, 6, 9, 12, 18), T = c(278.15, 298.15, 313.15),
                    batch = c("b1", "b2"), stringsAsFactors = FALSE)
  k  <- 0.05 * exp(70000 / 8.314 * (1 / 298.15 - 1 / g$T))
  y0 <- ifelse(g$batch == "b1", 99, 98)
  df <- data.frame(
    time  = g$time,
    batch = g$batch,
    y     = y0 - 20 * (1 - exp(-k * g$time)),
    stringsAsFactors = FALSE
  )
  # Same values, four names.  "tempk" is the default; the others are not.
  for (nm in c("tempk", "tempK", "myT", "Temperature_K")) df[[nm]] <- g$T
  df
}

fit_with <- function(df, temp_var, ...) {
  fit_sb(df, outcome_var = "y", intercept_cols = "batch",
         temp_var = temp_var, time_var = "time", trend = "decreasing", ...)
}

test_that("a non-default temp_var fits at all", {
  df <- sb_data()

  # Each of these raised "object '<temp_var>' not found" from inside deSolve.
  for (nm in c("tempK", "myT", "Temperature_K")) {
    testthat::expect_no_error(fit_with(df, nm))
  }
})

test_that("the temperature column's name does not change the estimates", {
  df  <- sb_data()
  ref <- coef(fit_with(df, "tempk"))

  for (nm in c("tempK", "myT", "Temperature_K")) {
    testthat::expect_equal(coef(fit_with(df, nm)), ref)
  }
})

test_that("predict.fit_sb agrees across temp_var names, for every interval", {
  df <- sb_data()
  nd <- data.frame(time = c(0, 6, 24), batch = "b1", stringsAsFactors = FALSE)

  for (int in c("none", "confidence", "prediction", "tolerance")) {
    ref <- NULL
    for (nm in c("tempk", "tempK", "myT")) {
      nd_nm <- nd
      nd_nm[[nm]] <- 298.15
      got <- predict(fit_with(df, nm), newdata = nd_nm, interval = int)
      if (is.null(ref)) ref <- got else testthat::expect_equal(got, ref)
    }
  }
})

test_that("predict_tos agrees across temp_var names", {
  df    <- sb_data()
  sched <- data.frame(t_start = c(0, 3), t_end = c(3, 12), temp = c(25, 40),
                      batch = "b1", stringsAsFactors = FALSE)

  ref <- NULL
  for (nm in c("tempk", "tempK", "myT")) {
    got <- predict_tos(fit_with(df, nm), sched)
    if (is.null(ref)) ref <- got else testthat::expect_equal(got, ref)
  }
  # Guard against the comparison passing on a degenerate answer: the profile
  # must actually degrade across the two segments.
  testthat::expect_true(ref$fit[2] < ref$fit[1])
})

test_that("a custom temp_expr's temp_var is honoured too", {
  # Here the name reaches build_sb() through temp_expr rather than the argument,
  # which is the path resolve_temp_expr() overwrites.
  df <- sb_data()
  te <- expr_arrhenius(temp_var = "myT")

  fit <- fit_with(df, "tempk", temp_expr = te)
  testthat::expect_equal(fit_with(df, "myT") |> coef(), coef(fit))
})

test_that("a temp_var named after an ODE parameter is rejected, not silently fitted", {
  # growth_ode() takes k/n/m and the temperature out of one named vector, so a
  # column called "n" overwrote the Weibull exponent.  The old behaviour was the
  # dangerous kind: it converged, to Ea in the tens instead of tens of thousands.
  df <- sb_data()

  # "Ea" is caught one step earlier than k/n/m: expr_arrhenius() itself refuses
  # to name a covariate after a parameter of its own expression (V10-02), so the
  # rejection never reaches build_alpha_ode()'s reserved-name check. Both
  # messages are accepted here -- what matters is that neither name fits.
  for (nm in c("k", "n", "m", "Ea")) {
    d <- df
    d[[nm]] <- d$tempk
    testthat::expect_error(
      fit_with(d, nm),
      "share a name with a model parameter|names a parameter in this expression"
    )
  }
})

test_that("log_Ea's parameter name is reserved as well", {
  df <- sb_data()
  df$Ealog <- df$tempk

  testthat::expect_error(
    fit_with(df, "Ealog", temp_expr = expr_arrhenius(temp_var = "Ealog",
                                                     log_Ea = TRUE)),
    "share a name with a model parameter"
  )
})

test_that("every spelling the log_Ea rename recognises is reserved", {
  for (nm in c("Ealog", "log_Ea", "logEa")) {
    df <- sb_data()
    df[[nm]] <- df$tempk
    testthat::expect_error(
      fit_with(df, nm, temp_expr = expr_arrhenius(temp_var = nm, log_Ea = TRUE)),
      "share a name with a model parameter"
    )
  }
})

test_that("the log_Ea rename leaves a column that merely contains a spelling alone", {
  df <- sb_data()
  df$Ea_logK <- df$tempk

  ref <- coef(fit_with(df, "tempk",
                       temp_expr = expr_arrhenius(temp_var = "tempk", log_Ea = TRUE)))
  got <- coef(fit_with(df, "Ea_logK",
                       temp_expr = expr_arrhenius(temp_var = "Ea_logK", log_Ea = TRUE)))
  testthat::expect_equal(got, ref)
  testthat::expect_true("Ealog" %in% names(got))
})

test_that("the humidity column is checked against the same reserved names", {
  # rh_col rides in the same parameter vector, so it has the same exposure.
  df <- sb_data()
  df$rh <- 60
  df$m  <- 60

  testthat::expect_error(
    fit_sb(df, outcome_var = "y", intercept_cols = "batch", time_var = "time",
           trend = "decreasing",
           temp_expr = expr_arrhenius_humidity(temp_var = "tempk",
                                               rh_col = "m", rh_ref = 60)),
    "share a name with a model parameter"
  )
})

test_that("humidity fits still work under a non-default temp_var", {
  # The RH path binds a second value into the same environment; make sure the
  # temperature rename did not disturb it. Three temperatures keep Ea
  # identifiable so the fit lands in the interior rather than on the Ea = 0
  # bound (a degenerate two-temperature design admits an Ea -> 0 optimum where
  # the humidity term absorbs the temperature effect).
  g  <- expand.grid(time = c(0, 3, 6, 9, 12), T = c(298.15, 308.15, 313.15),
                    rh = c(40, 75))
  k  <- 0.05 * exp(70000 / 8.314 * (1 / 298.15 - 1 / g$T) + 0.02 * (g$rh - 60))
  df <- data.frame(time = g$time, myT = g$T, humid = g$rh,
                   y = 99 - 20 * (1 - exp(-k * g$time)))
  df$tempk <- df$myT

  te_def <- expr_arrhenius_humidity(temp_var = "tempk", rh_col = "humid",
                                    rh_ref = 60)
  te_alt <- expr_arrhenius_humidity(temp_var = "myT", rh_col = "humid",
                                    rh_ref = 60)

  fit_def <- fit_sb(df, outcome_var = "y", time_var = "time", trend = "decreasing",
                    temp_expr = te_def)
  fit_alt <- fit_sb(df, outcome_var = "y", time_var = "time", trend = "decreasing",
                    temp_expr = te_alt)

  testthat::expect_equal(coef(fit_alt), coef(fit_def))
  testthat::expect_true("b" %in% names(coef(fit_alt)))
})
