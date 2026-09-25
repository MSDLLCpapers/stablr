test_that("tolerance_bounds_km validates sigma and d_squared", {
  # A negative sigma silently returned an inverted interval (lwr > upr); a
  # negative d_squared fed sqrt() an NaN that propagated into the bounds. Both
  # must error rather than emit a wrong-but-plausible tolerance bound.
  expect_error(
    stablr:::tolerance_bounds_km(0, 0.1, sigma = -1, df = 9,
                                 coverage = 0.95, level = 0.95),
    "sigma")
  expect_error(
    stablr:::tolerance_bounds_km(0, d_squared = -1, sigma = 1, df = 9,
                                 coverage = 0.95, level = 0.95),
    "d_squared")
  expect_error(
    stablr:::tolerance_bounds_km(0, 0.1, sigma = NA_real_, df = 9,
                                 coverage = 0.95, level = 0.95),
    "sigma")
  expect_error(
    stablr:::tolerance_bounds_km(0, d_squared = Inf, sigma = 1, df = 9,
                                 coverage = 0.95, level = 0.95),
    "d_squared")

  # A valid call is unaffected and returns an ordered interval.
  ok <- stablr:::tolerance_bounds_km(0, 0.1, sigma = 1, df = 9,
                                     coverage = 0.95, level = 0.95)
  expect_true(all(ok$lwr <= ok$upr))
})
