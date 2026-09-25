test_that("expr_arrhenius_humidity returns correctly formed expression", {
  # Default parameters
  expr <- expr_arrhenius_humidity()
  expected <- "exp(Ea / 8.314 * (1/298.15 - 1/tempk)) * exp(b * (rh - 0.6))"
  expect_equal(expr[["formula"]], expected)

  # Check expression is parseable
expect_no_error(parse(text = expr[["formula"]]))
})

test_that("expr_arrhenius_humidity handles log_Ea option", {
  # With log_Ea = TRUE
  expr <- expr_arrhenius_humidity(log_Ea = TRUE)
  expected <- "exp(exp(Ea_log) / 8.314 * (1/298.15 - 1/tempk)) * exp(b * (rh - 0.6))"
  expect_equal(expr[["formula"]], expected)

  # Check expression is parseable
  expect_no_error(parse(text = expr[["formula"]]))
})

test_that("expr_arrhenius_humidity respects custom parameters", {
  # Custom R constant
  expr <- expr_arrhenius_humidity(R = 10)
  expect_true(grepl("/ 10 *\\*", expr[["formula"]]))

  # Custom reference temperature
  expr <- expr_arrhenius_humidity(trefk = 300)
  expect_true(grepl("1/300", expr[["formula"]]))

  # Custom temperature variable
  expr <- expr_arrhenius_humidity(temp_var = "temperature")
  expect_true(grepl("1/temperature", expr[["formula"]]))

  # Custom humidity column
  expr <- expr_arrhenius_humidity(rh_col = "humidity")
  expected <- "exp(Ea / 8.314 * (1/298.15 - 1/tempk)) * exp(b * (humidity - 0.6))"
  expect_equal(expr[["formula"]], expected)

  # Custom reference humidity
  expr <- expr_arrhenius_humidity(rh_ref = 0.75)
  expected <- "exp(Ea / 8.314 * (1/298.15 - 1/tempk)) * exp(b * (rh - 0.75))"
  expect_equal(expr[["formula"]], expected)
})

test_that("expr_arrhenius_humidity combines all options correctly", {
  # All custom parameters at once
  expr <- expr_arrhenius_humidity(
    R = 10,
    trefk = 300,
    temp_var = "T",
    rh_col = "RH",
    rh_ref = 0.80,
    log_Ea = TRUE
  )

  expected <- "exp(exp(Ea_log) / 10 * (1/300 - 1/T)) * exp(b * (RH - 0.8))"
  expect_equal(expr[["formula"]], expected)

  # Check expression is parseable
  expect_no_error(parse(text = expr[["formula"]]))
})

test_that("expr_arrhenius_humidity expression structure is correct", {
  expr <- expr_arrhenius_humidity()

  # Should contain Arrhenius term
  expect_true(grepl("exp\\(Ea / .* \\* \\(1/.* - 1/.*\\)\\)", expr[["formula"]]))

  # Should contain humidity term
  expect_true(grepl("exp\\(b \\* \\(.*\\)\\)", expr[["formula"]]))

  # Should have multiplication between terms
  expect_true(grepl("\\) \\* exp\\(", expr[["formula"]]))

  # Should reference humidity sensitivity parameter 'b'
  expect_true(grepl("\\bb\\b", expr[["formula"]]))
})

test_that("expr_arrhenius_humidity handles edge cases", {
  # Zero reference humidity
  expr <- expr_arrhenius_humidity(rh_ref = 0)
  expect_true(grepl("\\(rh - 0\\)", expr[["formula"]]))

  # Reference humidity of 1 (100%)
  expr <- expr_arrhenius_humidity(rh_ref = 1)
  expect_true(grepl("\\(rh - 1\\)", expr[["formula"]]))

  # Very high reference temperature
  expr <- expr_arrhenius_humidity(trefk = 1000)
  expect_true(grepl("1/1000", expr[["formula"]]))

  # All expressions should be parseable
  expect_no_error(parse(text = expr_arrhenius_humidity(rh_ref = 0)[["formula"]]))
  expect_no_error(parse(text = expr_arrhenius_humidity(rh_ref = 1)[["formula"]]))
  expect_no_error(parse(text = expr_arrhenius_humidity(trefk = 1000)[["formula"]]))
})