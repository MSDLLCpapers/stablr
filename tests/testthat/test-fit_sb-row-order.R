test_that("fit_sb is invariant to input row order (H11-09)", {
  sb_kT <- function(tempk, k = 0.01, Ea = 80000)
    k * exp(Ea / 8.314 * (1 / 298.15 - 1 / tempk))

  set.seed(7)
  d <- expand.grid(batch = c("B1", "B2"), time = c(0, 1, 3, 6, 9),
                   temp_c = c(25, 40, 50), stringsAsFactors = FALSE)
  d$tempk <- d$temp_c + 273.15
  y0v <- c(B1 = 98, B2 = 60)
  d$y <- y0v[d$batch] * exp(-sb_kT(d$tempk) * d$time) + rnorm(nrow(d), 0, 0.05)

  f_ref <- fit_sb(d, intercept_cols = "batch", trend = "decreasing", Yinf = 0)

  set.seed(99)
  f_shuf <- fit_sb(d[sample(nrow(d)), ], intercept_cols = "batch", trend = "decreasing", Yinf = 0)

  # A permuted extract of the same table must give the identical fit: same
  # estimates, same deviance, and the same reported standard errors.
  expect_identical(unname(f_shuf$par), unname(f_ref$par))
  expect_identical(f_shuf$deviance, f_ref$deviance)
  expect_identical(unname(sqrt(diag(vcov(f_shuf)))),
                   unname(sqrt(diag(vcov(f_ref)))))
})
