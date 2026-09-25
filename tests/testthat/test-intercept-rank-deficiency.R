# expr_intercept() uses reference coding for 2+ intercept_cols: the first column
# keeps every level, later columns drop their first. That removes the linear
# dependency for a CROSSED design but not for a NESTED one -- if batches A,B only
# appear at site S1 and C,D only at S2, site is a function of batch and its
# indicator is already a sum of batch indicators. The design then emits 5
# parameters for 4 distinct cells (rank 4). The fit reproduced the data, and the
# per-cell sums were right, but the individual coefficients were an arbitrary
# ridge point that changed with the order of intercept_cols -- e.g.
# y0_batch_C = 121.91 next to y0_site_S2 = -29.91 for a cell whose true
# intercept is 92 -- flagged only by a generic "singular gradient matrix".

nested_data <- function() {
  b2s <- c(A = "S1", B = "S1", C = "S2", D = "S2")
  y0  <- c(A = 100, B = 96, C = 92, D = 88)
  d <- expand.grid(time = c(0, 3, 6, 9, 12),
                   tempk = c(298.15, 313.15),
                   batch = c("A", "B", "C", "D"),
                   stringsAsFactors = FALSE)
  d$site <- unname(b2s[d$batch])
  d$y <- unname(y0[d$batch]) *
    exp(-0.01 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15)) * d$time)
  d
}

# Crossed and genuinely additive: y0 = 100 - 4*(batch B) - 2*(orient upr).
crossed_data <- function() {
  d <- expand.grid(time = c(0, 3, 6, 9, 12),
                   tempk = c(298.15, 313.15),
                   batch = c("A", "B"),
                   orient = c("inv", "upr"),
                   stringsAsFactors = FALSE)
  y0 <- 100 + ifelse(d$batch == "B", -4, 0) + ifelse(d$orient == "upr", -2, 0)
  d$y <- y0 * exp(-0.01 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15)) * d$time)
  d
}

test_that("expr_intercept rejects a nested intercept design", {
  d <- nested_data()
  expect_error(expr_intercept(d, c("batch", "site")), "rank deficient")
})

test_that("the rejection does not depend on column order", {
  # The old behaviour produced a DIFFERENT parameter set and different numbers
  # for each order; both are equally unidentifiable, so both must be rejected.
  d <- nested_data()
  expect_error(expr_intercept(d, c("batch", "site")), "rank deficient")
  expect_error(expr_intercept(d, c("site", "batch")), "rank deficient")
})

test_that("the error reports the parameter count and rank", {
  d <- nested_data()
  msg <- tryCatch(expr_intercept(d, c("batch", "site")),
                  error = conditionMessage)
  expect_match(msg, "5 parameters")   # 4 batches + 1 non-reference site
  expect_match(msg, "rank 4")         # only 4 distinct cells
})

test_that("the fitters surface the nested design error", {
  d <- nested_data()
  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = c("batch", "site")),
    "rank deficient"
  )
})

test_that("a crossed design still uses reference coding and fits cleanly", {
  d <- crossed_data()

  expect_silent(ei <- expr_intercept(d, c("batch", "orient")))
  # batch keeps both levels; orient drops its reference ("inv").
  expect_equal(as.character(ei$params),
               c("y0_batch_A", "y0_batch_B", "y0_orient_upr"))

  ws <- character()
  fit <- withCallingHandlers(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = c("batch", "orient")),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("singular gradient", ws), 0)

  cf <- coef(fit)
  expect_equal(unname(cf[["y0_batch_A"]]), 100, tolerance = 1e-4)
  expect_equal(unname(cf[["y0_batch_B"]]), 96, tolerance = 1e-4)
  expect_equal(unname(cf[["y0_orient_upr"]]), -2, tolerance = 1e-4)
})

test_that("the finest-grained column recovers every nested cell intercept", {
  # The remedy the error message recommends must actually work.
  d <- nested_data()
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     intercept_cols = "batch")
  cf <- coef(fit)

  expect_equal(unname(cf[["y0_batch_A"]]), 100, tolerance = 1e-4)
  expect_equal(unname(cf[["y0_batch_B"]]), 96, tolerance = 1e-4)
  expect_equal(unname(cf[["y0_batch_C"]]), 92, tolerance = 1e-4)
  expect_equal(unname(cf[["y0_batch_D"]]), 88, tolerance = 1e-4)
})

test_that("single-column and empty intercept designs are unaffected", {
  d <- nested_data()
  expect_silent(expr_intercept(d, "batch"))
  expect_silent(expr_intercept(d, "site"))
  expect_silent(expr_intercept(d, character()))
  expect_silent(expr_intercept(d, NULL))
})

test_that("a duplicated intercept column is rejected as rank deficient", {
  # The most degenerate nesting: a column nested in itself.
  d <- nested_data()
  d$batch2 <- d$batch
  expect_error(expr_intercept(d, c("batch", "batch2")), "rank deficient")
})

test_that("partially crossed designs of full rank are accepted", {
  # Every batch/orient combination present except one cell: still full rank, so
  # the check must not reject it.
  d <- crossed_data()
  d <- d[!(d$batch == "B" & d$orient == "upr"), ]
  expect_silent(expr_intercept(d, c("batch", "orient")))
})

test_that("the rank check applies to prefixed multi-response builders", {
  # build_singledeg()/build_parallel() pass prefix = "A0"/"B0"/"C0"; the same
  # collinearity applies there.
  d <- nested_data()
  expect_error(expr_intercept(d, c("batch", "site"), prefix = "A0"),
               "rank deficient")
  expect_silent(expr_intercept(crossed_data(), c("batch", "orient"),
                               prefix = "B0"))
})
