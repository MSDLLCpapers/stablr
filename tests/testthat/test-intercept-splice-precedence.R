# expr_intercept() joins its terms with " + " and does not enclose the sum, so a
# bare splice into a subtraction ("Yinf - (Yinf - t1 + t2)") let R's precedence
# flip the sign of every term after the first. The optimiser absorbed the flip,
# leaving SSR, fitted() and predict() correct while coef() reported negated
# intercepts and predict_tos() -- which sums the terms as written -- disagreed.

make_increasing <- function(seed = 11) {
  set.seed(seed)
  y0 <- c(B1 = 0.83, B2 = 0.58, B3 = 0.56)
  d <- expand.grid(batch = names(y0), tempk = c(298.15, 308.15, 313.15, 323.15),
                   time = c(0, 1, 3, 6, 9, 12), stringsAsFactors = FALSE)
  k <- 3e-4 * exp(-92000 / 8.314 * (1 / d$tempk - 1 / 298.15))
  d$y <- 100 - (100 - y0[d$batch]) * exp(-k * d$time) + stats::rnorm(nrow(d), 0, 0.01)
  d
}

test_that("trend='increasing' reports intercepts with the sign of an nls oracle", {
  d <- make_increasing()
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "increasing", Yinf = 100,
                     intercept_cols = "batch")
  cf <- coef(fit)[c("y0_batch_B1", "y0_batch_B2", "y0_batch_B3")]

  # All three levels rise from a positive baseline, so no coefficient may be
  # negative -- the flipped fit reported -0.578 and -0.559 for B2 and B3.
  testthat::expect_true(all(cf > 0))

  # Explicit-indicator oracle: same model, parenthesisation forced by construction.
  d$i1 <- d$batch == "B1"; d$i2 <- d$batch == "B2"; d$i3 <- d$batch == "B3"
  orc <- stats::nls(
    y ~ 100 - (100 - (i1 * p1 + i2 * p2 + i3 * p3)) *
      exp(-k0 * exp(Ea / 8.314 * (1 / 298.15 - 1 / tempk)) * time),
    data = d,
    start = list(p1 = .8, p2 = .6, p3 = .6, k0 = 3e-4, Ea = 9e4)
  )
  testthat::expect_equal(unname(cf), unname(coef(orc)[c("p1", "p2", "p3")]),
                         tolerance = 1e-6)
})

test_that("reported intercepts equal the model's own t=0 prediction on every branch", {
  set.seed(7)
  y0 <- c(B1 = 95, B2 = 92, B3 = 90)
  gens <- list(
    "0-increasing"  = function(a, k, t) a + k * t,
    "0-decreasing"  = function(a, k, t) a - k * t,
    "1-increasing"  = function(a, k, t) 100 - (100 - a) * exp(-k * t),
    "1-decreasing"  = function(a, k, t) a * exp(-k * t),
    "2-increasing"  = function(a, k, t) 100 - ((100 - a)^-1 - k * t * -1)^-1,
    "2-decreasing"  = function(a, k, t) (a^-1 - k * t * -1)^-1
  )
  for (nm in names(gens)) {
    parts <- strsplit(nm, "-", fixed = TRUE)[[1]]
    ko <- as.integer(parts[1]); tr <- parts[2]

    d <- expand.grid(batch = names(y0), tempk = c(298.15, 308.15, 313.15),
                     time = c(0, 3, 6, 9, 12), stringsAsFactors = FALSE)
    k <- 2e-3 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15))
    d$y <- gens[[nm]](y0[d$batch], k, d$time) + stats::rnorm(nrow(d), 0, 0.02)

    fit <- fit_kinetic(d, kinetic_order = ko, trend = tr, Yinf = 100,
                       intercept_cols = "batch")
    cf <- coef(fit)[paste0("y0_batch_", names(y0))]
    at0 <- predict(fit, newdata = data.frame(batch = names(y0), time = 0,
                                             tempk = 298.15))
    testthat::expect_equal(unname(cf), unname(at0), tolerance = 1e-6,
                           info = paste("branch", nm))
  }
})

test_that("predict_tos agrees with predict on an increasing multi-level fit", {
  d <- make_increasing()
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "increasing", Yinf = 100,
                     intercept_cols = "batch")
  for (bb in c("B1", "B2", "B3")) {
    p_direct <- predict(fit, newdata = data.frame(batch = bb, time = 12,
                                                 tempk = 298.15))
    p_tos <- tos_segments(
      predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 298.15,
                                  batch = bb), temp_unit = "K"))$fit
    testthat::expect_equal(unname(p_tos), unname(p_direct), tolerance = 1e-8,
                           info = paste("batch", bb))
  }
})
