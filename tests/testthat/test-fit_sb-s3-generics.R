# W09-04 / W11-02: the model-summary generics must exist for fit_sb and return
# the same quantities they do for the three gsl_nls-based classes.
#
# fit_sb inherits from minpack.lm's nls.lm, which ships none of fitted / logLik /
# nobs / AIC / BIC, and stats::sigma.default turns the missing nobs() into a
# silent NaN rather than an error. The gsl_nls fitters inherit working methods,
# so a class-agnostic bake-off loop -- sum((y - fitted(fit))^2), sigma(fit),
# AIC(fit) -- silently recorded RSS 0, RMSE NaN and AIC NA for fit_sb alone,
# ranking it first on RSS and dropping it from AIC.

sb_g <- function(seed = 3) {
  g <- expand.grid(time = c(0, 3, 6, 12), temp_c = c(25, 40, 50),
                   stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  set.seed(seed)
  g$y <- 100 * exp(-k * g$time) + stats::rnorm(nrow(g), 0, 0.3)
  g
}

test_that("fitted() has one value per observation and reproduces residuals()", {
  g <- sb_g()
  f <- fit_sb(g, trend = "decreasing")

  expect_length(fitted(f), nrow(g))
  # residuals() == observed - fitted(), the identity all four classes hold.
  expect_equal(residuals(f), g$y - fitted(f))
  # The silent-zero idiom now equals deviance() instead of collapsing to 0.
  expect_equal(sum((g$y - fitted(f))^2), deviance(f))
  expect_gt(sum((g$y - fitted(f))^2), 0)
})

test_that("nobs / logLik / AIC / BIC / sigma match the Gaussian closed form", {
  g <- sb_g()
  f <- fit_sb(g, trend = "decreasing")

  n   <- nrow(g)
  p   <- f$n_free
  RSS <- deviance(f)
  ll  <- -n / 2 * (log(2 * pi) + log(RSS / n) + 1)

  expect_equal(nobs(f), n)
  expect_equal(as.numeric(logLik(f)), ll)
  expect_equal(attr(logLik(f), "df"), p + 1L)      # +1 for sigma^2
  expect_equal(attr(logLik(f), "nobs"), n)
  expect_equal(AIC(f), -2 * ll + 2 * (p + 1))
  expect_equal(BIC(f), -2 * ll + log(n) * (p + 1))
  expect_equal(sigma(f), sqrt(RSS / (n - p)))
  expect_false(is.na(sigma(f)))                    # was a silent NaN
})

test_that("the generics agree with a gsl_nls-based fitter on the same data", {
  g  <- sb_g()
  fs <- fit_sb(g, trend = "decreasing")
  fk <- fit_kinetic(g, kinetic_order = 1, trend = "decreasing")

  # Same closed form drives both, so the identities that hold for fit_kinetic
  # must now hold for fit_sb too.
  for (f in list(fs, fk)) {
    expect_equal(length(fitted(f)), nobs(f))
    expect_equal(sum((g$y - fitted(f))^2), deviance(f), tolerance = 1e-8)
    expect_equal(sigma(f), sqrt(deviance(f) / df.residual(f)), tolerance = 1e-8)
    expect_equal(AIC(f), -2 * as.numeric(logLik(f)) + 2 * attr(logLik(f), "df"))
  }
})

test_that("a pinned parameter is not counted as estimated by the new generics", {
  # n and m pinned via lower == upper: par has 7 entries but only 5 are free.
  # df.residual / sigma / logLik-df must use n_free, as vcov() and summary() do,
  # not the full par length that sigma.default would have used.
  d <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 308.15),
                   batch = c("B1", "B2"), stringsAsFactors = FALSE)
  set.seed(1)
  d$y <- pmax(100 * exp(-0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / d$tempk)) *
                          d$time) + stats::rnorm(nrow(d), 0, 0.3), 1)
  pinned <- suppressWarnings(fit_sb(
    d, rate_cols = "batch", intercept_cols = "batch", trend = "decreasing", Yinf = 0,
    start = c(Ea = 100, n = 1, m = 0.1, k_B1 = 1/12, k_B2 = 1/12,
              y0_B1 = 100, y0_B2 = 100),
    lower = c(Ea = 0,   n = 1, m = 0, k_B1 = 0,   k_B2 = 0,   y0_B1 = 50,  y0_B2 = 50),
    upper = c(Ea = 1e6, n = 1, m = 0, k_B1 = 100, k_B2 = 100, y0_B1 = 150, y0_B2 = 150)
  ))
  expect_equal(names(which(pinned$fixed)), c("n", "m"))

  n <- nobs(pinned)
  expect_equal(df.residual(pinned), n - pinned$n_free)
  # sigma keyed on n_free (matches df.residual), not on the full par length.
  expect_equal(sigma(pinned), sqrt(deviance(pinned) / (n - pinned$n_free)))
  expect_equal(attr(logLik(pinned), "df"), pinned$n_free + 1L)
})

test_that("summary() enriches the non-positive-definite Hessian error", {
  # summary.nls.lm() Cholesky-factorises the Hessian; a non-positive-definite
  # one aborts with an opaque LAPACK message. summary.fit_sb() must catch that
  # and add the modelling context and remedies, while leaving a healthy summary
  # untouched.
  g <- sb_g()
  f <- fit_sb(g, trend = "decreasing")

  # A healthy fit summarises without error.
  expect_error(invisible(utils::capture.output(summary(f))), NA)

  # Force the leading minor negative so the Cholesky in summary.nls.lm() fails.
  f_bad <- f
  h <- f_bad$hessian
  h[1, 1] <- -abs(h[1, 1])
  f_bad$hessian <- h

  expect_error(summary(f_bad), "leading minor of order 1 is not positive")
  expect_error(summary(f_bad), "Hessian is not positive definite")
  expect_error(summary(f_bad), "widen")
})
