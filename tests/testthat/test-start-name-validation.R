test_that("fitters reject unknown start names and report missing ones", {
  df <- evers_supp |> dplyr::filter(cqa == "Purity" & !is.na(y))

  # A valid start for a first-order fit_kinetic on this data.
  good <- c(Ea = 80000, k_0 = 0.01, y0 = 100)

  # Extra (unrecognized) name is rejected rather than silently fitted.
  expect_error(
    fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                start = c(good, bogus = 7)),
    "Unknown starting value name"
  )

  # Missing name is still reported.
  expect_error(
    fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                start = c(Ea = 80000, k_0 = 0.01)),
    "Missing starting values"
  )

  # fit_sb rejects an unknown start name via its name-aligned bookkeeping.
  expect_error(
    suppressWarnings(
      fit_sb(df, intercept_cols = "batch", trend = "decreasing", Yinf = 0,
             start = c(Ea = 80000, n = 1, m = 0.1, k_0 = 0.01, bogus = 7))
    ),
    "does not match model parameters|Unknown"
  )
})
