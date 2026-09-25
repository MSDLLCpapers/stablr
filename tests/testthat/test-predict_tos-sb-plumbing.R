# predict_tos() on the fit_sb path diverged from the fit in two ways that no
# input check could catch, because both live in how the composite is assembled
# rather than in what the caller supplied.
#
# V02-10: coef.nls.lm() returns the full parameter vector (fit_sb reconstructs it
# so that pinned parameters print), but vcov.fit_sb() covers only the free ones.
# The delta-method product J V J' was therefore n x p against n_free x n_free and
# aborted with a bare "non-conformable arguments" naming neither the fit nor the
# parameter -- so any fit_sb with a parameter pinned via lower == upper could not
# be used with predict_tos() at all.
#
# V02-16: the composite's own deSolve::ode() call omitted the rtol/atol that
# build_alpha_ode()'s predict_run() integrates with, leaving deSolve's defaults.
# The default atol = 1e-6 is exactly the initial state (alpha = 1e-6), so the
# absolute error tolerance equalled the state itself at the start of every
# segment, and the answer depended on how the schedule was sliced.

sb_data <- function(seed = 1) {
  d <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 308.15),
                   batch = c("B1", "B2"), stringsAsFactors = FALSE)
  set.seed(seed)
  d$y <- pmax(100 * exp(-0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / d$tempk)) *
                          d$time) + stats::rnorm(nrow(d), 0, 0.3), 1)
  d
}

sb_fit <- function(...) {
  fit_sb(sb_data(), rate_cols = "batch", intercept_cols = "batch", trend = "decreasing", Yinf = 0, ...)
}

test_that("predict_tos works on a fit_sb with a pinned parameter (V02-10)", {
  # n and m pinned via lower == upper, the documented way to fix a parameter.
  pinned <- sb_fit(
    start = c(Ea = 100, n = 1, m = 0.1, k_B1 = 1/12, k_B2 = 1/12,
              y0_B1 = 100, y0_B2 = 100),
    lower = c(Ea = 0,   n = 1, m = 0, k_B1 = 0,   k_B2 = 0,   y0_B1 = 50,  y0_B2 = 50),
    upper = c(Ea = 1e6, n = 1, m = 0, k_B1 = 100, k_B2 = 100, y0_B1 = 150, y0_B2 = 150)
  )
  expect_equal(names(which(pinned$fixed)), c("n", "m"))
  # The mismatch the fix has to absorb.
  expect_lt(ncol(stats::vcov(pinned)), length(stats::coef(pinned)))

  sch <- data.frame(t_start = c(0, 6), t_end = c(6, 12), temp = c(5, 25),
                    batch = "B1")
  p <- tos_segments(predict_tos(pinned, sch))

  expect_equal(nrow(p), 2L)
  expect_true(all(is.finite(c(p$fit, p$lwr, p$upr))))
  # A pinned parameter is a constant, so it must contribute no variance: the
  # interval stays finite and of strictly positive width, since the free
  # parameters still contribute theirs.
  expect_true(all(p$upr - p$lwr > 0))
  expect_true(all(p$lwr < p$fit & p$fit < p$upr))

  # summary.fit_sb() cat()s its fixed-parameter block; predict_tos() must not
  # leak that to the console mid-prediction.
  expect_silent(predict_tos(pinned, sch))
})

test_that("a fit with no pinned parameter is unaffected by the alignment (V02-10)", {
  free <- sb_fit()
  expect_false(any(free$fixed))
  sch <- data.frame(t_start = 0, t_end = 12, temp = 25, batch = "B1")
  p <- tos_segments(predict_tos(free, sch))
  expect_equal(ncol(stats::vcov(free)), length(stats::coef(free)))
  expect_true(all(is.finite(c(p$fit, p$lwr, p$upr))))
  expect_gt(p$upr - p$lwr, 0)
})

test_that("a constant hold gives the same answer however it is sliced (V02-16)", {
  fit <- sb_fit()
  # Same physics, different segmentation: 25 C from 0 to 12 in 1, 2, 4 and 12
  # equal segments. At deSolve's defaults these spanned 78.73818 to 78.73881.
  finals <- vapply(c(1, 2, 4, 12), function(nseg) {
    br <- seq(0, 12, length.out = nseg + 1)
    sch <- data.frame(t_start = br[-length(br)], t_end = br[-1], temp = 25,
                      batch = "B1")
    p <- predict_tos(fit, sch)
    p$fit[nrow(p)]
  }, numeric(1))

  expect_equal(diff(range(finals)), 0, tolerance = 1e-5)

  # And the unsplit hold must agree with predict.fit_sb(), which solves the same
  # trajectory in one go at rtol = 1e-8 / atol = 1e-10.
  reference <- predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                batch = "B1"))
  expect_equal(finals, rep(as.numeric(reference), 4), tolerance = 1e-6)
})
