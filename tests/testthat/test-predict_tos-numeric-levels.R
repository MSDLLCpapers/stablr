# All-digit grouping levels must be compared as strings in the composite
# propagation, exactly as in the fitted formula.
#
# format_level() emits an all-digit level unquoted, and build_kinetic() repairs
# that by running quote_numeric_levels() over the ASSEMBLED RHS -- not over the
# intercept and rate parts, which stay unquoted in the formula's closure. So
# .build_composite_fn_kinetic() cannot simply reuse those parts: it has to
# re-apply the quoting itself.
#
# Without it the composite emits (batch == 0000803274). R coerces the level to
# the numeric 803274, which never equals the character "0000803274", so the
# indicator is FALSE on every row and BOTH the intercept and the rate drop out.
# The prediction is 0 -- for stability work, 0% of the attribute remaining, the
# worst possible wrong answer -- and since every gradient is zero too, the
# delta-method interval collapses to zero width and the function reports
# certainty. Levels without leading zeros round-trip through the coercion
# unharmed, so the failure is invisible unless a level actually has them.

digit_level_fit <- function(lv, rate_cols = "batch",
                            intercept_cols = "batch") {
  g <- expand.grid(batch = lv, time = c(0, 3, 6, 12), temp_c = c(25, 40),
                   stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y <- stats::setNames(c(100, 99), lv)[g$batch] * exp(-k * g$time)
  suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                               intercept_cols = intercept_cols,
                               rate_cols = rate_cols))
}

# Oracle: evaluate the fitted formula's RHS directly. The fitted RHS is quoted,
# so this is independent of whatever the composite builder assembles.
oracle <- function(fit, ...) {
  env <- list2env(as.list(stats::coef(fit)), parent = baseenv())
  eval(stats::formula(fit)[[3]],
       envir = data.frame(..., stringsAsFactors = FALSE), enclos = env)
}

# Segment rows only: these tests are about the level quoting in the propagation,
# not about the time-0 baseline row predict_tos() prepends.
tos1 <- function(fit, level, temp = 298.15, t_end = 12) {
  tos_segments(suppressWarnings(predict_tos(
    fit, data.frame(t_start = 0, t_end = t_end, temp = temp, batch = level,
                    stringsAsFactors = FALSE),
    temp_unit = "K")))
}

test_that("a leading-zero all-digit level predicts the fitted value", {
  lv  <- c("0000803274", "0000803279")
  fit <- digit_level_fit(lv)
  p   <- tos1(fit, lv[1])

  expect_equal(p$fit, oracle(fit, batch = lv[1], time = 12, tempk = 298.15),
               tolerance = 1e-8)
  expect_equal(p$fit, 78.66278611, tolerance = 1e-6)
  # The old failure mode, stated directly: not 0% remaining.
  expect_gt(p$fit, 1)
})

test_that("the stored expression parts are quoted at construction", {
  # This test used to assert the opposite -- that the parts were unquoted, and
  # that quote_numeric_levels() repaired them afterwards -- and said that if the
  # builders ever quoted the parts themselves the extra step could go. They now
  # do: sanitized_level() quotes every comparison target as it is written,
  # because sanitized_col() coerces the column to character and a bare numeric
  # target cannot match a string. quote_numeric_levels() is still applied to the
  # assembled right-hand sides but is a no-op on an already-quoted target.
  #
  # Matched on the comparison alone: the column reference is wrapped in a
  # sanitising gsub() (see sanitized_col()), so "batch ==" is no longer adjacent
  # text. What matters is that the leading zeros survive, in both places.
  lv  <- c("0000803274", "0000803279")
  fit <- digit_level_fit(lv)
  stored <- get("int_part", envir = environment(stats::formula(fit)))$expr

  expect_match(as.character(stored), "== '0000803274'", fixed = TRUE)
  expect_false(grepl("== 0000803274", as.character(stored), fixed = TRUE))
  # ...and the fitted formula compares against the quoted level too.
  expect_match(paste(deparse(stats::formula(fit)[[3]]), collapse = ""),
               "== \"0000803274\"", fixed = TRUE)
})

test_that("both the intercept and the rate survive the coercion", {
  # The unquoted comparison zeroed the intercept AND the rate. Fit each role on
  # its own so a fix that only repaired one of them is caught.
  lv <- c("0000803274", "0000803279")

  int_only <- digit_level_fit(lv, rate_cols = character())
  p_int <- tos1(int_only, lv[1])
  expect_equal(p_int$fit, oracle(int_only, batch = lv[1], time = 12,
                                 tempk = 298.15), tolerance = 1e-8)
  expect_gt(p_int$fit, 1)

  rate_only <- digit_level_fit(lv, intercept_cols = character())
  p_rate <- tos1(rate_only, lv[1])
  expect_equal(p_rate$fit, oracle(rate_only, batch = lv[1], time = 12,
                                  tempk = 298.15), tolerance = 1e-8)
  expect_gt(p_rate$fit, 1)
})

test_that("levels with and without leading zeros agree", {
  # Same digits minus the leading zeros round-trip through numeric coercion
  # unharmed, which is what made the bug invisible. Both must now be right.
  lead <- digit_level_fit(c("0000803274", "0000803279"))
  bare <- digit_level_fit(c("803274", "803279"))

  expect_equal(tos1(lead, "0000803274")$fit, tos1(bare, "803274")$fit,
               tolerance = 1e-8)
})

test_that("short all-digit levels are handled", {
  lv  <- c("01", "02")
  fit <- digit_level_fit(lv)
  p   <- tos1(fit, "01")
  expect_equal(p$fit, oracle(fit, batch = "01", time = 12, tempk = 298.15),
               tolerance = 1e-8)
  expect_gt(p$fit, 1)
})

test_that("predict_tos agrees with predict() on a constant-temperature schedule", {
  # Independent cross-check: over a single isothermal segment the composite must
  # reproduce the ordinary prediction path, which never rebuilt the expression.
  lv  <- c("0000803274", "0000803279")
  fit <- digit_level_fit(lv)

  nd <- data.frame(batch = lv[1], time = 12, tempk = 298.15,
                   stringsAsFactors = FALSE)
  direct <- suppressWarnings(stats::predict(fit, newdata = nd))

  expect_equal(tos1(fit, lv[1])$fit, as.numeric(direct)[1], tolerance = 1e-8)
})

test_that("a multi-segment schedule keeps the level quoted in every segment", {
  lv  <- c("0000803274", "0000803279")
  fit <- digit_level_fit(lv)

  sch <- data.frame(t_start = c(0, 6), t_end = c(6, 12),
                    temp = c(298.15, 313.15), batch = c(lv[1], lv[1]),
                    stringsAsFactors = FALSE)
  p <- tos_segments(suppressWarnings(predict_tos(fit, sch, temp_unit = "K")))

  expect_equal(p$fit[1], oracle(fit, batch = lv[1], time = 6, tempk = 298.15),
               tolerance = 1e-8)
  expect_true(all(p$fit > 1))
  # Warmer second segment must degrade further, not reset to zero.
  expect_lt(p$fit[2], p$fit[1])
})
