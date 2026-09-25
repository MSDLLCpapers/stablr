# predict_tos() used to accept several malformed or implausible schedules and
# return a normal-looking row rather than failing: NA in a numeric column gave a
# full row of NA fit/se/lwr/upr with no condition, an out-of-range `level` gave a
# zero-width, infinite or inverted interval, a non-physical temperature gave
# fit = 0 with se = 0, and an omitted grouping column surfaced only as
# "object 'batch' not found" from deep inside the composite.

fit_for_guards <- function(seed = 31) {
  set.seed(seed)
  d <- expand.grid(batch = c("B1", "B2"), tempk = c(298.15, 313.15, 323.15),
                   time = c(0, 3, 6, 9, 12), stringsAsFactors = FALSE)
  d$y <- 98 * exp(-4e-3 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15)) * d$time) +
    stats::rnorm(nrow(d), 0, 0.05)
  fit_kinetic(d, kinetic_order = 1, trend = "decreasing", intercept_cols = "batch")
}

sched <- function(...) data.frame(...)

test_that("predict_tos rejects malformed schedules instead of returning NA rows", {
  fit <- fit_for_guards()

  testthat::expect_error(
    predict_tos(fit, sched(t_start = numeric(0), t_end = numeric(0),
                           temp = numeric(0), batch = character(0))),
    "no rows"
  )
  testthat::expect_error(
    predict_tos(fit, sched(t_start = 0, t_end = 6, temp = NA_real_, batch = "B1")),
    "NA in column\\(s\\): temp"
  )
  testthat::expect_error(
    predict_tos(fit, sched(t_start = 0, t_end = NA_real_, temp = 25, batch = "B1")),
    "NA in column\\(s\\): t_end"
  )
  # A grouping column the fit used must be supplied for every segment.
  testthat::expect_error(
    predict_tos(fit, sched(t_start = 0, t_end = 6, temp = 25)),
    "missing grouping column\\(s\\).*batch"
  )
})

test_that("predict_tos requires level strictly inside (0, 1)", {
  fit <- fit_for_guards()
  s <- sched(t_start = 0, t_end = 6, temp = 25, batch = "B1")
  # 0 gave a zero-width interval, 1 gave +/-Inf, -0.5 inverted it (lwr > upr),
  # 1.5 gave NaN behind a bare "NaNs produced" warning.
  for (bad in list(0, 1, -0.5, 1.5, NA_real_, c(0.9, 0.95), "0.95")) {
    testthat::expect_error(predict_tos(fit, s, level = bad),
                           "level must be a single number")
  }
  testthat::expect_no_error(predict_tos(fit, s, level = 0.90))
})

test_that("predict_tos rejects non-physical and flags implausible temperatures", {
  fit <- fit_for_guards()

  # tempk <= 0 is an error: the Arrhenius factor stays computable there, so this
  # previously returned fit = 0 with se = 0.
  testthat::expect_error(
    predict_tos(fit, sched(t_start = 0, t_end = 6, temp = -400, batch = "B1")),
    "non-physical absolute temperature"
  )

  # Celsius passed while temp_unit = "K" reads 25 as 25 K, where the rate is ~0,
  # so the excursion reports essentially no degradation. Warn, don't stop.
  testthat::expect_warning(
    predict_tos(fit, sched(t_start = 0, t_end = 24, temp = 25, batch = "B1"),
                temp_unit = "K"),
    "outside the usual stability range"
  )

  # Legitimate cold storage must still be usable, so the band cannot be an error:
  # -80 C freezer and liquid nitrogen are real stability conditions.
  testthat::expect_warning(
    r <- predict_tos(fit, sched(t_start = 0, t_end = 24, temp = -80, batch = "B1")),
    "outside the usual stability range"
  )
  testthat::expect_true(all(is.finite(r$fit)))

  # Ordinary conditions warn about nothing.
  testthat::expect_no_warning(
    predict_tos(fit, sched(t_start = c(0, 6), t_end = c(6, 9), temp = c(5, 25),
                           batch = "B1"))
  )
})

test_that("predict_tos fails on a non-finite prediction rather than returning NaN", {
  # The fit_sb composite integrates its ODE numerically; a stiff segment returned
  # alpha = NaN, which reached the output as fit = NaN with a NaN interval and no
  # condition signalled. No input check can anticipate this, so it is caught on
  # the way out.
  set.seed(4)
  d <- expand.grid(batch = c("B1", "B2"), tempk = c(298.15, 313.15, 323.15),
                   time = c(0, 3, 6, 9, 12, 24, 48), stringsAsFactors = FALSE)
  ks <- 0.5 * exp(-60000 / 8.314 * (1 / d$tempk - 1 / 298.15))
  d$y <- 100 * exp(-ks * d$time) + stats::rnorm(nrow(d), 0, 0.05)
  fsb <- suppressWarnings(fit_sb(
    d, rate_cols = "batch", intercept_cols = "batch", trend = "decreasing", Yinf = -1,
    start = c(Ea = 60000, n = 1, m = 0, k_B1 = .5, k_B2 = .5,
              y0_B1 = 100, y0_B2 = 100)
  ))

  res <- suppressWarnings(try(
    predict_tos(fsb, sched(t_start = 0, t_end = 200, temp = 60, batch = "B1")),
    silent = TRUE
  ))
  # Either it evaluates finitely or it errors -- what it must never do is hand
  # back a NaN fit.
  if (inherits(res, "try-error")) {
    testthat::expect_match(conditionMessage(attr(res, "condition")),
                           "not finite")
  } else {
    testthat::expect_true(all(is.finite(res$fit)))
  }
})
