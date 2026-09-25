test_that("invert_hessian verifies against the original Hessian, not the preconditioned one (V11-06/W02b-03)", {
  # A well-conditioned SPD matrix is inverted exactly and silently.
  set.seed(1)
  A <- crossprod(matrix(rnorm(25), 5, 5)) + diag(5)
  Ai <- invert_hessian(A)
  expect_lt(max(abs(A %*% Ai - diag(5))), 1e-7)

  # Diagonal preconditioning D^{-1} H D^{-1} can make a rank-deficient Hessian's
  # correlation form pass a Cholesky-inverse check while the rescaled inverse is
  # off by orders of magnitude on H itself. Here H = D C D with C well-
  # conditioned enough that chol2inv(chol(S)) reconstructs S to ~3e-11 (the old
  # check on S would accept), but the ~1e18 diagonal spread amplifies the
  # rescaling error to ~1e-2 on H. The acceptance test must be on H, so this must
  # warn and fall back to the pseudo-inverse rather than silently return it.
  C <- matrix(0.3, 4, 4); diag(C) <- 1
  C[1, 2] <- C[2, 1] <- 0.999999          # near-collinear pair
  D <- sqrt(c(1e-6, 1e-1, 1e3, 1e12))     # ~18 orders of magnitude (Ea..k_0)
  H <- outer(D, D) * C
  H <- (H + t(H)) / 2

  # Confirm the construction: the old S-check would have passed silently.
  dsc <- sqrt(abs(diag(H)))
  S   <- H / outer(dsc, dsc)
  expect_lt(max(abs(S %*% chol2inv(chol(S)) - diag(4))), 1e-7)

  expect_warning(Hi <- invert_hessian(H), "singular|rank def")
  expect_true(is.matrix(Hi))
})

test_that("predict.fit_sb warns rather than silently inflating a rank-deficient CI (V11-06)", {
  # Refrigerated reference that barely degrades: only k*arr is identified on the
  # accelerated arm, so Ea is unidentified and the Hessian is rank-deficient.
  set.seed(7)
  d <- expand.grid(time = c(0, 3, 6, 9, 12, 18, 24, 36),
                   tempk = c(313.15, 313.65))          # 0.5 C apart
  kT <- 0.05 * exp(110000 / 8.314 * (1 / 313.15 - 1 / d$tempk))
  d$y <- 100 * exp(-kT * d$time) + rnorm(nrow(d), 0, 0.3)

  f <- suppressWarnings(fit_sb(d, trend = "decreasing", Yinf = 0))
  nd <- data.frame(time = c(12, 24), tempk = 313.15)

  expect_warning(
    predict(f, newdata = nd, interval = "confidence", method = "delta"),
    "singular|rank def")
})
