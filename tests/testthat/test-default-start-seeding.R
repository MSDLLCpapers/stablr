# Every parameter now gets an explicit starting value.
#
# Previously the fitters only seeded parameters when a humidity coefficient was
# present; otherwise `start` was an all-NA vector, which gsl_nls() treats as a
# request to self-start -- and self-start sets every NA parameter to 1. With an
# intercept of 1 against a response of 80-98, and k = 1 against times up to 24
# weeks, the initial curve is exp(-24) ~ 4e-11: eleven orders of magnitude below
# the data. Levenberg-Marquardt climbs out of that with three parameters, but not
# once `rate_cols` splits the rate into one parameter per batch. It then failed
# two different ways:
#
#   fit_kinetic    -> hard error, "Missing value or an infinity produced when
#                     evaluating the model" (k goes negative, exp(-k*t) overflows,
#                     numericDeriv aborts) before a single iteration completes.
#   fit_singledeg  -> silently returned isConv = TRUE with RSS = NaN and rate
#                     constants of 2.2e6, 5.5e6 and -1.5e5. Worse than the error:
#                     nothing in the object marks it as garbage.
#
# Both are start-point failures, not model failures: the same models fit to a
# LOWER RSS than their submodels once seeded.

pur <- function() subset(evers_supp, cqa == "Purity" & !is.na(y))

upsec_wide <- function() {
  uw <- tidyr::pivot_wider(upsec, names_from = "cqa", values_from = "y")
  keep <- c("Total Main", "HMWS", "LMWS", "time", "tempk", "batch")
  uw[stats::complete.cases(uw[, keep]), ]
}

test_that("fit_kinetic with rate_cols fits instead of erroring", {
  # Was: Error, "Missing value or an infinity produced when evaluating the model".
  m <- fit_kinetic(pur(), kinetic_order = 1, trend = "decreasing",
                   rate_cols = "batch")
  expect_s3_class(m, "fit_kinetic")
  expect_true(m$convInfo$isConv)
  expect_true(is.finite(sum(residuals(m)^2)))
})

test_that("splitting the rate improves the fit rather than breaking it", {
  # The submodel is nested in the split model, so the split model's RSS must be
  # no worse. That it errored at all was the tell that the start point, not the
  # model, was the problem.
  d <- pur()
  m0 <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing")
  m1 <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                    rate_cols = "batch")
  m2 <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                    rate_cols = "batch", intercept_cols = "batch")

  expect_lte(sum(residuals(m1)^2), sum(residuals(m0)^2) + 1e-8)
  expect_lte(sum(residuals(m2)^2), sum(residuals(m1)^2) + 1e-8)
})

test_that("fit_singledeg with rate_cols no longer returns NaN as a converged fit", {
  # Was: isConv = TRUE, RSS = NaN, k_batch_* in the millions including negatives.
  m <- fit_singledeg(upsec_wide(), kinetic_order = 1, trend = "decreasing",
                     yA_var = "Total Main", yB_var = "HMWS",
                     rate_cols = "batch")
  expect_true(m$convInfo$isConv)
  expect_true(is.finite(sum(residuals(m)^2)))

  k <- coef(m)[grep("^k", names(coef(m)))]
  expect_true(all(is.finite(k)))
  expect_true(all(k > 0))          # a degradation rate cannot be negative
  expect_true(all(abs(k) < 1))     # and 2.2e6 per week is not a rate
})

test_that("fit_parallel with rate_cols converges to a finite RSS", {
  m <- fit_parallel(upsec_wide(), kinetic_order = 1,
                    yA_var = "Total Main", yB_var = "HMWS", yC_var = "LMWS",
                    rate_cols = "batch")
  expect_true(m$convInfo$isConv)
  expect_true(is.finite(sum(residuals(m)^2)))
})

test_that("a user-supplied start still wins", {
  # The seeding must be a default, not an override: someone who supplies starts
  # has a reason.
  d <- pur()
  bad <- c(Ea = 1e5, y0 = 90, k_0 = 0.02)
  m <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", start = bad)
  expect_true(m$convInfo$isConv)

  # An impossible start is reported as such rather than being quietly replaced.
  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                start = c(Ea = 8e4, wrong_name = 90, k_0 = 0.01))
  )
})

test_that("default_start seeds Ea on the scale implied by R", {
  # Same requirement as the Ea-seed bug: the seed is an energy, so it must track
  # the gas constant the user supplied.
  p <- c("Ea", "y0", "k_0")
  expect_equal(default_start(p, R = 8.314)[["Ea"]], default_Ea_seed(8.314))
  expect_equal(default_start(p, R = 1.987)[["Ea"]], default_Ea_seed(1.987))
  expect_gt(default_start(p, R = 8.314)[["Ea"]],
            default_start(p, R = 1.987)[["Ea"]])
})

test_that("default_start puts log_Ea on the log scale", {
  s <- default_start(c("Ea_log", "y0", "k_0"), R = 8.314)
  expect_equal(s[["Ea_log"]], log(default_Ea_seed(8.314)))
  # "Ea_log" must not also be caught by the plain-Ea branch.
  expect_false(s[["Ea_log"]] == default_Ea_seed(8.314))
})

test_that("default_start recognises the branch-suffixed parallel names", {
  p <- c("kB_0", "kC_0", "EaB", "EaC", "bB", "bC", "A0", "B0", "C0")
  s <- default_start(p, R = 8.314)
  expect_true(all(is.finite(s[c("kB_0", "kC_0", "EaB", "EaC", "bB", "bC")])))
  expect_equal(unname(s[c("EaB", "EaC")]), rep(default_Ea_seed(8.314), 2))
  expect_equal(unname(s[c("bB", "bC")]), c(0, 0))
  # No data frame, so the rate seed falls back to the historical 0.01 default.
  expect_equal(unname(s[c("kB_0", "kC_0")]), c(0.01, 0.01))

  s2 <- default_start(c("Ea_logB", "Ea_logC", "kB_0"), R = 8.314)
  expect_equal(unname(s2[c("Ea_logB", "Ea_logC")]),
               rep(log(default_Ea_seed(8.314)), 2))
})

test_that("default_start seeds the humidity coefficient at zero", {
  # 0 means "no moisture effect until the data says otherwise".
  s <- default_start(c("Ea", "b", "y0", "k_0"), R = 8.314)
  expect_equal(s[["b"]], 0)
})

test_that("default_start seeds the rate from the observed decay", {
  # A fixed k = 0.01 is blind to the time scale. exp(-0.01 * t) is flat
  # over an ordinary window, so the Ea gradient at the seed is ~0 and the fit
  # walks off to a spurious minimum. Seed from the fractional change per unit
  # time instead: |y_end - y0| / (y0 * span).
  d <- data.frame(y = c(100, 60, 30), tempk = 298.15, time = c(0, 0.5, 1))
  expect_equal(default_start(c("Ea", "y0", "k_0"), d, 8.314, "y",
                             "decreasing")[["k_0"]], abs(30 - 100) / (100 * 1))
  # An explicit k_seed still overrides the data-driven default.
  expect_equal(default_start(c("Ea", "y0", "k_0"), d, 8.314, "y",
                             "decreasing", k_seed = 0.05)[["k_0"]], 0.05)
  # No usable response or time column falls back to the historical 0.01.
  expect_equal(default_start(c("Ea", "y0", "k_0"))[["k_0"]], 0.01)
})

test_that("rate_seed reads fast decay and floors slow decay at 0.01", {
  # Fast decay (response nearly vanishes) needs a large seed; a near-flat series
  # keeps the historical 0.01 floor, since its rate genuinely is tiny and it
  # fitted fine there.
  fast <- data.frame(y = c(100, 1), time = c(0, 1))
  slow <- data.frame(y = c(100, 99.5), time = c(0, 6))
  expect_equal(rate_seed(fast, "y", "time"), 0.99)
  expect_equal(rate_seed(slow, "y", "time"), 0.01)
  # Fallbacks to 0.01 whenever no estimate is possible.
  expect_equal(rate_seed(data.frame(time = c(0, 4)), "nope", "time"), 0.01)
  expect_equal(rate_seed(NULL, "y", "time"), 0.01)
})

test_that("fit_kinetic recovers Ea when the true rate is far above 0.01", {
  # With a fixed 0.01 seed this fit converged to Ea ~ 1e6 at
  # more than double the RSS of the true minimum. The time-based seed reaches
  # the good minimum by default.
  R <- 8.314; Tref <- 298.15; Ea <- 90000; k_ref <- 3.0
  set.seed(1)
  d <- expand.grid(time = c(0, .25, .5, .75, 1),
                   tempk = c(25, 40, 50, 60) + 273.15)
  kT <- k_ref * exp(-(Ea / R) * (1 / d$tempk - 1 / Tref))
  d$y <- 100 * exp(-kT * d$time) + rnorm(nrow(d), 0, 0.3)

  m <- fit_kinetic(d, outcome_var = "y", temp_var = "tempk", time_var = "time",
                   trend = "decreasing", kinetic_order = 1)
  expect_equal(coef(m)[["Ea"]], Ea, tolerance = 0.05)
})

test_that("intercepts seed from the response at the earliest timepoint", {
  # The intercept IS the response at t = 0, so read it off the data.
  d <- data.frame(y = c(98, 90, 80), tempk = 298.15, time = 0:2)
  expect_equal(default_start(c("Ea", "y0", "k_0"), d, 8.314, "y",
                             "decreasing")[["y0"]], 98)
  # ...and the same value regardless of `trend`, because time settles it.
  expect_equal(default_start(c("Ea", "y0", "k_0"), d, 8.314, "y",
                             "increasing")[["y0"]], 98)
})

test_that("intercept_seed ignores a trend that disagrees with the data", {
  # This is the bug the time-based seed fixes. fit_parallel() defaults to
  # trend = "increasing", so a DECAYING yA was seeded at its final value: with
  # yA running 100 -> 46.8, the old rule seeded A0 = 46.8. From there EaB ran to
  # -1.0e7 and the fit returned isConv = TRUE with an NaN residual sum of
  # squares -- non-convergence dressed up as success.
  d <- data.frame(y = c(100, 70, 46.8), tempk = 298.15, time = c(0, 6, 12))
  expect_equal(intercept_seed(d, "y", "time", "increasing"), 100)
  expect_equal(intercept_seed(d, "y", "time", "decreasing"), 100)
})

test_that("intercept_seed averages replicates at the earliest timepoint", {
  # With the intercept split by a grouping column, every group starts from the
  # pooled baseline rather than from one arbitrary row.
  d <- data.frame(y = c(98, 100, 80, 60), time = c(0, 0, 5, 5))
  expect_equal(intercept_seed(d, "y", "time", "decreasing"), 99)
})

test_that("intercept_seed falls back to the trend's extreme without a time column", {
  d <- data.frame(y = c(80, 90, 98))
  expect_equal(intercept_seed(d, "y", "time", "decreasing"), 98)
  expect_equal(intercept_seed(d, "y", "time", "increasing"), 80)
  # Same fallback when the time column exists but is unusable.
  d2 <- data.frame(y = c(80, 90, 98), time = c(NA, NA, NA))
  expect_equal(intercept_seed(d2, "y", "time", "decreasing"), 98)
})

test_that("intercept_seed returns NA when it cannot seed", {
  expect_true(is.na(intercept_seed(NULL, "y")))
  expect_true(is.na(intercept_seed(data.frame(z = 1), "y")))
  expect_true(is.na(intercept_seed(data.frame(y = 1), NULL)))
  expect_true(is.na(intercept_seed(data.frame(y = c(NA_real_, NA_real_)), "y")))
})

test_that("intercept_seed skips non-finite responses at t = 0", {
  d <- data.frame(y = c(NA, 100, 70), time = c(0, 0, 6))
  expect_equal(intercept_seed(d, "y", "time", "decreasing"), 100)
})

test_that("fit_parallel seeds a decaying yA to a finite fit", {
  # Regression on the seeding fix, end to end: the intercept A0 is read from the
  # response at the earliest time (not inferred from a trend), so a decaying yA
  # is seeded at its starting value and the fit stays finite.
  d <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15))
  d$yA <- 100 * exp(-0.01 * exp(-80000 / 8.314 *
                                  (1 / d$tempk - 1 / 298.15)) * d$time)
  d$yB <- (100 - d$yA) * 0.6
  d$yC <- (100 - d$yA) * 0.4

  m <- fit_parallel(d, kinetic_order = 1,
                    yA_var = "yA", yB_var = "yB", yC_var = "yC")
  expect_true(is.finite(sum(residuals(m)^2)))
  expect_false(is.na(coef(m)[["A0"]]))
})

test_that("every intercept parameter is seeded when the design splits them", {
  d <- pur()
  m <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                   intercept_cols = "batch", rate_cols = "batch")
  s <- default_start(attr(m, "build")$param_names, d, 8.314, "y", "decreasing")
  expect_false(anyNA(s))
  # Nothing may be left for gsl_nls self-start to fill in with 1.
  expect_true(all(is.finite(s)))
})

test_that("default_start falls back gracefully without data", {
  # Intercepts stay NA (self-start) rather than erroring when the response is
  # unavailable -- there is no better guess to make.
  s <- default_start(c("Ea", "y0", "k_0"))
  expect_true(is.na(s[["y0"]]))
  expect_true(is.finite(s[["Ea"]]))
  expect_true(is.finite(s[["k_0"]]))

  s2 <- default_start(c("Ea", "y0"), data.frame(z = 1), 8.314, "y")
  expect_true(is.na(s2[["y0"]]))

  s3 <- default_start(c("Ea", "y0"), data.frame(y = c(NA_real_, NA_real_)),
                      8.314, "y")
  expect_true(is.na(s3[["y0"]]))

  expect_length(default_start(character()), 0L)
})

test_that("default_start returns names aligned with the parameters it was given", {
  p <- c("Ea", "y0_batch_1", "y0_batch_2", "k_batch_1", "k_batch_2")
  s <- default_start(p, pur(), 8.314, "y", "decreasing")
  expect_named(s, p)
})
