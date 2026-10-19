skip_on_cran()

# An identical call means a different MODEL in fit_sb than in the closed-form
# fitters.
#
#   fit_kinetic/fit_singledeg/fit_parallel -- expr_intercept() reference coding:
#     |A| + |B| - 1 additive intercepts (3 for a 2x2).
#   fit_sb -- make_group() = interaction(sep = "||"), full cell means:
#     |A| x |B| intercepts (4 for a 2x2).
#
# On 32 rows that is 8 parameters / 24 residual df versus 11 / 21. Neither coding
# is wrong, so both are kept; the difference is documented rather than warned
# about. These tests pin the divergence so it cannot change silently.

coding_data <- function() {
  d <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15),
                   batch = c("B1", "B2"), orient = c("inv", "upr"),
                   stringsAsFactors = FALSE)
  d$y <- 98 * exp(-0.01 * exp(-80000 / 8.314 *
                                (1 / d$tempk - 1 / 298.15)) * d$time)
  d
}

# Non-additive cell intercepts: cell means can fit these, reference coding cannot.
interaction_data <- function() {
  d <- coding_data()
  y0 <- c("B1.inv" = 98, "B2.inv" = 95, "B1.upr" = 90, "B2.upr" = 99)
  d$y <- y0[paste(d$batch, d$orient, sep = ".")] *
    exp(-0.01 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15)) * d$time)
  d
}

test_that("the two paths really do fit different parameter counts", {
  # Documents the divergence rather than asserting it away: this is the current
  # behaviour of both fitters.
  d <- coding_data()
  # Synthetic data has no real batch/orient effect, so the extra parameters are
  # unidentifiable and gsl_nls warns of a singular gradient; irrelevant to the
  # parameter counts asserted here.
  k <- suppressWarnings(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = c("batch", "orient"),
                rate_cols = c("batch", "orient")))
  s <- fit_sb(d, intercept_cols = c("batch", "orient"),
              rate_cols = c("batch", "orient"), trend = "decreasing", Yinf = 0)

  expect_length(grep("^y0", names(coef(k))), 3L)   # additive
  expect_length(grep("^y0", names(s$par)), 4L)     # cell means
  expect_length(coef(k), 8L)
  expect_length(s$par, 11L)
  expect_equal(df.residual(k), 24L)
})

test_that("one intercept column is coded identically by both paths", {
  # A single column is cell means in both paths, so the counts agree.
  d <- coding_data()
  k <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                   intercept_cols = "batch")
  s <- fit_sb(d, intercept_cols = "batch", trend = "decreasing", Yinf = 0)
  expect_equal(length(grep("^y0", names(coef(k)))),
               length(grep("^y0", names(s$par))))
})

test_that("crossing the columns by hand makes the two paths agree", {
  # The remedy the documentation recommends must actually work.
  d <- coding_data()
  d$cell <- interaction(d$batch, d$orient, drop = TRUE)

  k <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                   intercept_cols = "cell", rate_cols = "cell")
  s <- fit_sb(d, intercept_cols = "cell", rate_cols = "cell",
              trend = "decreasing", Yinf = 0)
  expect_length(grep("^y0", names(coef(k))), 4L)
  expect_length(grep("^y0", names(s$par)), 4L)
})

test_that("the coding difference changes the fit, not just the names", {
  # With a real batch x orient interaction, cell means recover all four
  # intercepts and reference coding structurally cannot -- which is why this is
  # a substantive divergence and not a naming convention.
  d <- interaction_data()
  k <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                   intercept_cols = c("batch", "orient"),
                   rate_cols = c("batch", "orient"))
  s <- fit_sb(d, intercept_cols = c("batch", "orient"),
              rate_cols = c("batch", "orient"), trend = "decreasing", Yinf = 0)
  expect_gt(sum(residuals(k)^2), 1)          # cannot reach the cell means
  expect_lt(sum(s$fvec^2), 1e-6)             # fits them essentially exactly

  # And the additive coefficients are distorted enough to reverse the ordering
  # of the two batches within orient = "inv" (true: B1 98 > B2 95).
  cf <- coef(k)
  expect_gt(cf[["y0_batch_B2"]], cf[["y0_batch_B1"]])
})

test_that("the other closed-form fitters reference-code multi-column intercepts", {
  # fit_singledeg and fit_parallel share expr_intercept(), so a 2x2 gives the
  # same 3 additive intercepts as fit_kinetic, not the 4 cell means of fit_sb.
  d <- coding_data()
  d$yA <- d$y
  d$yB <- 100 - d$y
  d$yC <- (100 - d$y) * 0.4

  # Unidentifiable extra parameters may draw a singular-gradient warning from
  # gsl_nls; not what this test is about.
  sd <- suppressWarnings(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  yA_var = "yA", yB_var = "yB",
                  intercept_cols = c("batch", "orient")))
  pl <- suppressWarnings(
    fit_parallel(d, kinetic_order = 1,
                 yA_var = "yA", yB_var = "yB", yC_var = "yC",
                 intercept_cols = c("batch", "orient")))
  expect_length(grep("^A0", names(coef(sd))), 3L)
  expect_length(grep("^A0", names(coef(pl))), 3L)
})
