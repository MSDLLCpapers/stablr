decay_arm <- function(temp_c, times, k) {
  data.frame(temp_c = temp_c, tempk = temp_c + 273.15, time = times,
             y = 100 * exp(-k * times))
}

rise_arm <- function(temp_c, times, k) {
  data.frame(temp_c = temp_c, tempk = temp_c + 273.15, time = times,
             y = 100 - 40 * exp(-k * times))
}

# Cold arm followed 36 months (barely moves); hot arm followed ~1 week (steep).
# Every arm decreases, but pooled cov(time, y) is positive.
imbalanced_decay <- function() {
  rbind(
    decay_arm(5,  c(0, 6, 12, 18, 24, 36), k = 0.0004),
    decay_arm(25, c(0, 3, 6, 12),          k = 0.004),
    decay_arm(40, c(0, 1, 3, 6),           k = 0.03),
    decay_arm(50, c(0, 0.1, 0.23),         k = 0.7)
  )
}

test_that("the setup actually has the sign-flipping imbalance", {
  d <- imbalanced_decay()
  arm_slopes <- vapply(split(d, d$temp_c),
                       function(a) unname(stats::coef(stats::lm(y ~ time, a))[2]),
                       numeric(1))
  expect_true(all(arm_slopes < 0))          # unambiguously decreasing
  expect_gt(stats::cov(d$time, d$y), 0)     # yet pooled cov flips positive
})

test_that("fit_sb fits imbalanced data: no check reads the pooled slope", {
  d <- imbalanced_decay()
  # fit_kinetic (takes trend from the caller) is the behavioural oracle: it fits.
  expect_no_error(fit_kinetic(d, kinetic_order = 1, trend = "decreasing"))
  # Y0 (~100) and Yinf (0) bracket every observation, so the geometry check passes
  # regardless of the pooled-covariance sign flip.
  f <- fit_sb(d, trend = "decreasing", Yinf = 0)
  expect_s3_class(f, "fit_sb")
  expect_equal(length(fitted(f)), nrow(d))
})

test_that("an unspecified Yinf takes the asymptote the trend implies", {
  dec <- fit_sb(imbalanced_decay(), trend = "decreasing")
  expect_equal(attr(dec, "build")$settings$Yinf, 0)
  expect_equal(attr(dec, "build")$settings$trend, "decreasing")

  inc <- fit_sb(rbind(rise_arm(25, c(0, 3, 6, 9, 12), 0.05),
                      rise_arm(40, c(0, 3, 6, 9, 12), 0.15)),
                trend = "increasing")
  expect_equal(attr(inc, "build")$settings$Yinf, 100)
  expect_equal(attr(inc, "build")$settings$trend, "increasing")
})

test_that("a Yinf inside the data warns but still fits", {
  d <- rbind(decay_arm(25, c(0, 3, 6, 9, 12), 0.05),
             decay_arm(40, c(0, 3, 6, 9, 12), 0.15))

  w <- capture_warnings(f <- fit_sb(d, trend = "decreasing", Yinf = 50))
  expect_match(w, "does not lie below every observation", all = FALSE)
  expect_s3_class(f, "fit_sb")
  expect_no_warning(fit_sb(d, trend = "decreasing", Yinf = 0))

  r <- rbind(rise_arm(25, c(0, 3, 6, 9, 12), 0.05),
             rise_arm(40, c(0, 3, 6, 9, 12), 0.15))
  w <- capture_warnings(f <- fit_sb(r, trend = "increasing", Yinf = 80))
  expect_match(w, "does not lie above every observation", all = FALSE)
  expect_s3_class(f, "fit_sb")
  expect_no_warning(fit_sb(r, trend = "increasing", Yinf = max(r$y) + 10))
})

test_that("supplying the default Yinf is the same as omitting it", {
  d <- rbind(decay_arm(25, c(0, 3, 6, 9, 12), 0.05),
             decay_arm(40, c(0, 3, 6, 9, 12), 0.15))
  d$y[which.min(d$y)] <- -0.05

  # The sub-zero reading sits just past Yinf = 0, which warns but still fits.
  expect_no_error(suppressWarnings(fit_sb(d, trend = "decreasing", Yinf = 0)))
  expect_no_error(suppressWarnings(fit_sb(d, trend = "decreasing")))
})

test_that("a time-zero replicate above the group mean stays reachable", {
  d <- rbind(decay_arm(25, c(0, 3, 6, 9, 12), 0.05),
             decay_arm(40, c(0, 3, 6, 9, 12), 0.15))
  d <- rbind(d, transform(d[d$time == 0, ], y = 97))
  d$y[d$temp_c == 25 & d$time == 3] <- 99.5

  expect_no_error(fit_sb(d, trend = "decreasing", Yinf = 0))
})

test_that("a defaulted Yinf on a backwards trend warns about the geometry", {
  d <- rbind(decay_arm(25, c(0, 3, 6, 9, 12), 0.05),
             decay_arm(40, c(0, 3, 6, 9, 12), 0.15))
  w <- capture_warnings(fit_sb(d, trend = "increasing"))
  expect_match(w, "outside the range spanned by Y0 and Yinf", all = FALSE)
})

test_that("a Yinf that cuts through the data on the far side of Y0 warns", {
  # Balanced decay from 100 down to ~1; Y0 ~ 100. Yinf = 200 puts the reachable
  # interval at [100, 200], but the data runs down to ~1.
  d <- rbind(decay_arm(25, c(0, 3, 6, 9, 12), 0.05),
             decay_arm(40, c(0, 3, 6, 9, 12), 0.15))
  w <- capture_warnings(fit_sb(d, trend = "increasing", Yinf = 200))
  expect_match(w, "outside the range spanned by Y0 and Yinf", all = FALSE)
  expect_no_warning(fit_sb(d, trend = "decreasing", Yinf = 0))  # default Yinf brackets the data
})

test_that("increasing data with Yinf below Y0 warns", {
  over_arm <- function(temp_c, times, k) {
    data.frame(temp_c = temp_c, tempk = temp_c + 273.15, time = times,
               y = 100 * exp(k * times))               # rises above Y0 ~ 100
  }
  d <- rbind(over_arm(25, c(0, 3, 6, 9, 12), 0.02),
             over_arm(40, c(0, 3, 6, 9, 12), 0.05))
  # Yinf = 0 gives reachable interval [0, 100], but the data climbs past 100.
  w <- capture_warnings(fit_sb(d, trend = "decreasing", Yinf = 0))
  expect_match(w, "outside the range spanned by Y0 and Yinf", all = FALSE)
  # A Yinf above the highest observation is warning-free.
  expect_no_warning(fit_sb(d, trend = "increasing", Yinf = max(d$y) + 10))
})

test_that("reachability is judged per intercept group, not on a pooled Y0", {
  # Two batches with very different starting levels. A single pooled Y0 (~75)
  # would wrongly reject B1's readings near 100; per-group Y0 must not.
  mk <- function(batch, y0) rbind(
    within(decay_arm(25, c(0, 3, 6, 9, 12), 0.05), { y <- y * y0 / 100; batch <- batch }),
    within(decay_arm(40, c(0, 3, 6, 9, 12), 0.15), { y <- y * y0 / 100; batch <- batch })
  )
  d <- rbind(mk("B1", 100), mk("B2", 50))   # B1 starts ~100, B2 starts ~50

  # Yinf = 0 is below both groups' Y0, so both are fully bracketed -> no warning.
  expect_no_warning(fit_sb(d, intercept_cols = "batch", trend = "decreasing",
                           Yinf = 0))
})
