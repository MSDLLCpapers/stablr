# expr_rate()/make_group() enumerate whatever rate cells are OBSERVED, with no
# check that a cell has enough distinct times to identify a rate. A batch present
# only at time == 0 contributes residuals mathematically independent of its k --
# exp(-k * 0) == 1 for every k -- so the parameter is completely unconstrained.
# The optimiser left it at 0, the worst possible default: k = 0 means "never
# degrades", and predict() would extrapolate a flat 24-month forecast from a
# parameter never estimated from data. The singularity also turned the WHOLE
# covariance matrix NaN, so one badly-supported batch destroyed the standard
# errors of every other parameter. The only signal was a generic "singular
# gradient matrix" warning naming nothing.

# B1/B2 fully sampled; B3's sampling is the variable under test.
support_data <- function(b3_times = 0, include_b3 = TRUE,
                         b3_temps = c(298.15, 313.15)) {
  g <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 313.15),
                   batch = c("B1", "B2"), stringsAsFactors = FALSE)
  g$y <- 100 * exp(-0.01 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15)) * g$time)
  if (!include_b3) return(g)

  b3 <- expand.grid(time = b3_times, tempk = b3_temps, batch = "B3",
                    stringsAsFactors = FALSE)
  b3$y <- 90 * exp(-0.01 * exp(-80000 / 8.314 * (1 / b3$tempk - 1 / 298.15)) * b3$time)
  rbind(g, b3)
}

fit_support <- function(d) {
  suppressWarnings(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                rate_cols = "batch", intercept_cols = "batch")
  )
}

# Collect every warning a fit emits. An under-supported design is genuinely
# singular, so the optimiser adds its own gradient warnings; asserting on the
# collected set keeps those from masking the support warning under test.
warnings_from <- function(expr) {
  ws <- character()
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

support_warnings <- function(d) {
  warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                rate_cols = "batch", intercept_cols = "batch")
  )
}

expect_support_warning <- function(expr) {
  expect_gte(length(grep("insufficient data", warnings_from(expr))), 1)
}

test_that("a single-timepoint batch is named in a fit-time warning", {
  expect_support_warning(
    fit_kinetic(support_data(0), kinetic_order = 1, trend = "decreasing",
                rate_cols = "batch", intercept_cols = "batch")
  )

  msg <- grep("insufficient data", support_warnings(support_data(0)), value = TRUE)[1]
  # The old failure gave only "singular gradient matrix", naming nothing.
  expect_match(msg, "B3")
  expect_match(msg, "no rows with elapsed time > 0")
})

test_that("predict refuses to extrapolate an unidentified rate", {
  m <- fit_support(support_data(0))
  expect_error(
    predict(m, newdata = data.frame(batch = "B3", time = 24, tempk = 313.15)),
    "never estimated from data"
  )
  # The regression: this used to return 90 -- a flat, zero-degradation forecast.
  expect_error(
    predict(m, newdata = data.frame(batch = "B3", time = 24, tempk = 313.15)),
    "B3"
  )
})

test_that("predicting the unsupported cell at time 0 is still allowed", {
  # At t = 0 the rate does not enter the prediction, so there is nothing unsound.
  m <- fit_support(support_data(0))
  expect_equal(
    unname(predict(m, newdata = data.frame(batch = "B3", time = 0, tempk = 313.15))),
    90, tolerance = 1e-3
  )
})

test_that("well-supported batches still predict normally in the same model", {
  m <- fit_support(support_data(0))
  expect_silent(
    p <- predict(m, newdata = data.frame(batch = "B1", time = 24, tempk = 313.15))
  )
  expect_lt(p, 100)   # genuinely degrading, not a flat line
})

test_that("a batch sampled at several times triggers no warning and predicts", {
  d <- support_data(c(0, 6, 12))
  expect_length(grep("insufficient data", support_warnings(d)), 0)

  m <- fit_support(d)
  p <- predict(m, newdata = data.frame(batch = "B3", time = 24, tempk = 313.15))
  expect_lt(p, 90)
})

test_that("a two-timepoint batch is treated as supported", {
  # Two distinct times is the minimum that carries rate information.
  d <- support_data(c(0, 12))
  expect_length(grep("insufficient data", support_warnings(d)), 0)
})

test_that("a batch at one nonzero time and one temperature is flagged", {
  # One timepoint, one temperature, and an intercept exclusive to B3: y0 and k
  # trade off with nothing to break the tie. The fit puts k at 0.0187 against a
  # truth of 0.01 and every standard error in the model goes NaN.
  d <- support_data(6, b3_temps = 298.15)
  msg <- grep("insufficient data", support_warnings(d), value = TRUE)[1]
  expect_match(msg, "only one distinct timepoint")
  expect_match(msg, "B3")
})

test_that("one timepoint at TWO temperatures is supported", {
  # Ea is shared, so the two temperatures differ only through a known factor:
  # two observations, two B3-specific unknowns (y0 and k), both identified. This
  # was flagged and then hard-refused at predict() time, on a rate the fit
  # recovers exactly.
  d <- support_data(6)   # b3_temps defaults to both
  expect_length(grep("insufficient data", support_warnings(d)), 0)

  m <- fit_support(d)
  expect_equal(unname(coef(m)["k_batch_B3"]), 0.01, tolerance = 1e-4)
  expect_false(any(is.nan(summary(m)$coefficients[, 2])))
  expect_silent(predict(m, newdata = data.frame(batch = "B3", time = 24,
                                                tempk = 313.15)))
})

test_that("one timepoint with a POOLED intercept is supported", {
  # No intercept_cols, so every batch shares one y0 that B1/B2's five timepoints
  # pin down. B3's single observation then determines its k on its own.
  d <- support_data(6, b3_temps = 298.15)
  expect_equal(nrow(unsupported_rate_cells(d, "batch", "time",
                                           intercept_cols = character(),
                                           temp_var = "tempk")), 0L)

  ws <- warnings_from(
    m <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     rate_cols = "batch")
  )
  expect_length(grep("insufficient data", ws), 0)
  expect_silent(predict(m, newdata = data.frame(batch = "B3", time = 24,
                                                tempk = 313.15)))
})

test_that("fully sampled data is unaffected", {
  d <- support_data(include_b3 = FALSE)
  expect_length(grep("insufficient data", support_warnings(d)), 0)

  m <- fit_support(d)
  expect_silent(predict(m, newdata = data.frame(batch = "B1", time = 24,
                                                tempk = 313.15)))
  # Standard errors are finite -- the NaN contagion is absent.
  expect_false(any(is.nan(summary(m)$coefficients[, 2])))
})

test_that("unsupported_rate_cells identifies the offending cells", {
  d <- support_data(0)
  bad <- unsupported_rate_cells(d, "batch", "time")
  expect_equal(nrow(bad), 1L)
  expect_match(bad$label, "B3")
  expect_match(bad$reason, "no rows with elapsed time")

  expect_equal(nrow(unsupported_rate_cells(support_data(c(0, 6, 12)),
                                           "batch", "time")), 0L)
})

test_that("unsupported_rate_cells handles a single global rate", {
  # No rate_cols: the whole data set is one cell, and an all-time-0 data set
  # cannot identify even the global rate.
  d <- support_data(include_b3 = FALSE)
  expect_equal(nrow(unsupported_rate_cells(d, character(), "time")), 0L)

  d0 <- d[d$time == 0, ]
  bad <- unsupported_rate_cells(d0, character(), "time")
  expect_equal(nrow(bad), 1L)
  expect_match(bad$label, "global rate")
})

test_that("unsupported_rate_cells is a no-op on unusable input", {
  d <- support_data(0)
  expect_equal(nrow(unsupported_rate_cells(d, "batch", "not_a_column")), 0L)
  expect_equal(nrow(unsupported_rate_cells(d, "not_a_column", "time")), 0L)
})

test_that("fit_singledeg warns about an under-supported rate cell", {
  d <- support_data(0)
  d$yA <- d$y
  d$yB <- 100 - d$y

  expect_support_warning(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  yA_var = "yA", yB_var = "yB", rate_cols = "batch")
  )
})

test_that("fit_parallel warns about an under-supported rate cell", {
  d <- support_data(0)
  d$yA <- d$y
  d$yB <- (100 - d$y) * 0.6
  d$yC <- (100 - d$y) * 0.4

  expect_support_warning(
    fit_parallel(d, kinetic_order = 1, yA_var = "yA", yB_var = "yB",
                 yC_var = "yC", rate_cols = "batch")
  )
})
