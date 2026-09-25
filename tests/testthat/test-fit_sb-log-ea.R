test_that("log_Ea is a reparameterisation: Ea agrees with the linear path (X18-03)", {
  sb_kT <- function(tempk, k = 0.02, Ea = 80000)
    k * exp(Ea / 8.314 * (1 / 298.15 - 1 / tempk))

  # Increasing-trend fixtures across a range of true Ea, including the high
  # region (>1.2e5) where the old log_Ea seed of log(100) ran to the guard rail.
  for (Ea_true in c(80000, 1.2e5, 1.5e5)) {
    set.seed(1)
    d <- expand.grid(time = c(0, 3, 6, 12, 24),
                     tempk = c(298.15, 313.15, 323.15))
    d$y <- 10 + 90 * (1 - exp(-sb_kT(d$tempk, k = 0.02, Ea = Ea_true) * d$time))

    f_lin <- suppressWarnings(
      fit_sb(d, trend = "increasing", Yinf = 100, log_Ea = FALSE))
    f_log <- suppressWarnings(
      fit_sb(d, trend = "increasing", Yinf = 100, log_Ea = TRUE))

    ea_lin <- unname(coef(f_lin)[["Ea"]])
    ea_log <- exp(unname(coef(f_log)[["Ealog"]]))

    # Both parameterisations must recover the same activation energy, and it must
    # not rest on the log upper rail exp(log(1e6)) = 1e6.
    expect_equal(ea_log, ea_lin, tolerance = 1e-3)
    expect_lt(ea_log, 1e6 * (1 - 1e-6))
  }
})
