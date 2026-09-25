test_that("fit_sb delta CI is not inflated by the FD step in a refrigerated-reference regime (V11-02)", {
  # Refrigerated reference (5 C) with an accelerated 40-50 C arm: the Arrhenius
  # factor on the hot arm is ~190x the reference, so an *absolute* finite-
  # difference step on the tiny rate k_0 is a huge relative perturbation that
  # attenuates the k_0 Jacobian column and destroys the negative Ea-k
  # cross-covariance cancellation, inflating the delta CI several-fold. A
  # *relative* step (theta_j * 1e-4) must reproduce the fine-step Jacobian.
  set.seed(1)
  trefk <- 278.15
  d <- expand.grid(time = c(0, 6, 12, 24, 36),
                   tempk = c(278.15, 313.15, 323.15))
  arr <- exp(110000 / 8.314 * (1 / trefk - 1 / d$tempk))
  d$y <- 100 * exp(-5e-5 * arr * d$time) + rnorm(nrow(d), 0, 0.3)

  f  <- suppressWarnings(fit_sb(d, temp_expr = expr_arrhenius(trefk = trefk),
                                trend = "decreasing", Yinf = 0))
  nd <- data.frame(time = 36, tempk = 313.15)

  # Reconstruct the delta-method CI width with a much finer relative step, using
  # the fit's own predict() -- this is the exact-Jacobian oracle.
  oracle_width <- function(step) {
    coef_est  <- f$par
    fit_vals  <- as.numeric(predict(f, newdata = nd))
    resid     <- f$fvec
    n_free    <- length(coef_est)
    sigma2    <- sum(resid^2) / (length(resid) - n_free)
    vcov_free <- sigma2 * stablr:::invert_hessian(f$hessian)
    G <- matrix(NA_real_, length(fit_vals), n_free)
    for (j in seq_len(n_free)) {
      pp  <- coef_est
      th  <- abs(pp[j])
      eps <- if (th > 0) th * step else step
      pp[j] <- pp[j] + eps
      ot <- f; ot$par <- pp
      G[, j] <- (as.numeric(predict(ot, newdata = nd)) - fit_vals) / eps
    }
    v <- rowSums((G %*% vcov_free) * G)
    2 * stats::qt(0.975, length(resid) - n_free) * sqrt(v)
  }

  ci <- suppressWarnings(
    predict(f, newdata = nd, interval = "confidence", method = "delta"))
  reported <- diff(as.numeric(ci[, c("lwr", "upr")]))

  # The old absolute-step defect made this ratio ~5.46; the relative step keeps
  # it at ~1. Allow 3% for the coarser production step vs the fine oracle.
  # (The refrigerated reference barely degrades, so Ea is near-unidentified and
  # the Hessian is rank-deficient -- the prediction CI is still accurate, but
  # invert_hessian() warns; suppressed here since this test is about the step.)
  expect_equal(reported / suppressWarnings(oracle_width(1e-7)), 1,
               tolerance = 0.03)
})
