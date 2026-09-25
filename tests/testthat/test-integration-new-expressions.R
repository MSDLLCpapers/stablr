# Integration tests for new expression functions
# These tests verify that expr_arrhenius_humidity()
# works correctly with the existing fitting infrastructure

test_that("expr_arrhenius_humidity adds humidity dependence", {
  # Standard Arrhenius
  arr_expr <- expr_arrhenius()

  # Arrhenius with humidity
  arr_hum_expr <- expr_arrhenius_humidity()

  # Humidity expression should have both temperature and humidity terms
  expect_true(grepl("exp\\(Ea", arr_hum_expr[["formula"]]))
  expect_true(grepl("exp\\(b", arr_hum_expr[["formula"]]))
  expect_true(grepl("\\brh\\b", arr_hum_expr[["formula"]]))

  # Standard Arrhenius should not have humidity terms
  expect_false(grepl("\\bb\\b", arr_expr[["formula"]]))
  expect_false(grepl("\\brh\\b", arr_expr[["formula"]]))
})

test_that("New expressions are compatible with build functions", {

  # Get standard expression components
  # data("evers_supp")
  df <- evers_supp |> dplyr::filter(cqa == "Purity")

  int_part <- expr_intercept(df, "pkg")
  rate_part <- expr_rate(df, "pkg")

  # Standard Arrhenius
  arr_expr <- expr_arrhenius(log_Ea = TRUE)

  # Arrhenius with humidity (would need rh column in data)
  hum_expr <- expr_arrhenius_humidity(log_Ea = TRUE)

  # Both temperature expressions should be valid
  expect_type(arr_expr[["formula"]], "character")
  expect_type(hum_expr[["formula"]], "character")

  # They should all be parseable R code
  expect_no_error(parse(text = arr_expr[["formula"]]))
  expect_no_error(parse(text = hum_expr[["formula"]]))

  # Build example formula (as build_kinetic would do)
  formula_arr <- paste("y ~", int_part$expr, "-", rate_part$expr, "*", arr_expr[["formula"]], "* time")

  # Should parse as valid formula
  expect_s3_class(as.formula(formula_arr), "formula")
})

test_that("Parameter names for new expressions are distinct", {
  
  # Each expression uses different parameter names
  arr_expr <- expr_arrhenius()
  hum_expr <- expr_arrhenius_humidity()

  # Standard Arrhenius uses Ea
  expect_true(grepl("\\bEa\\b", arr_expr[["formula"]]))
  expect_false(grepl("\\bb\\b", arr_expr[["formula"]]))

  # Humidity model uses Ea and b
  expect_true(grepl("\\bEa\\b", hum_expr[["formula"]]))
  expect_true(grepl("\\bb\\b", hum_expr[["formula"]]))
})

test_that("New expressions handle reference values correctly", {
  
  # Humidity model with custom reference humidity
  hum_expr <- expr_arrhenius_humidity(rh_ref = 0.75)
  expect_true(grepl("0.75", hum_expr[["formula"]], fixed = TRUE))

  # Custom column names
  hum_expr_custom <- expr_arrhenius_humidity(
    temp_var = "T_kelvin",
    rh_col = "RH_frac"
  )
  expect_true(grepl("T_kelvin", hum_expr_custom[["formula"]]))
  expect_true(grepl("RH_frac", hum_expr_custom[["formula"]]))
})
