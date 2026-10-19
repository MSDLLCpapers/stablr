# predict_tos() returns a time-0 baseline row ahead of the segment rows, so the
# returned table describes the whole trajectory from the start of the schedule
# rather than beginning at the end of the first segment. The baseline is the
# composite evaluated over zero elapsed time -- which reduces to the model's
# fitted intercept -- carried through the same delta-method interval as every
# other row, so it is directly plottable alongside them.

base_data <- function() {
  g <- expand.grid(batch = c("A", "B"), time = c(0, 3, 6, 12),
                   temp_c = c(25, 40), stringsAsFactors = FALSE)
  g$tempk <- g$temp_c + 273.15
  k <- 0.02 * exp(-80000 / 8.314 * (1 / g$tempk - 1 / 298.15))
  g$y  <- c(A = 100, B = 97)[g$batch] * exp(-k * g$time)
  g$yA <- g$y
  g$yB <- 100 - g$y
  g
}

base_fit <- function() {
  suppressWarnings(fit_kinetic(base_data(), kinetic_order = 1,
                               trend = "decreasing", intercept_cols = "batch"))
}

three_seg <- function(batch = "A") {
  data.frame(t_start = c(0, 12, 13), t_end = c(12, 13, 24),
             temp = c(5, 25, 5), batch = batch, stringsAsFactors = FALSE)
}

test_that("the output carries exactly the documented columns", {
  p <- predict_tos(base_fit(), three_seg())
  expect_equal(names(p),
               c("segment", "time", "temp", "variable", "fit", "lwr", "upr"))
})

test_that("a k-segment schedule returns k + 1 rows, the first at time 0", {
  p <- predict_tos(base_fit(), three_seg())

  expect_equal(nrow(p), 4L)
  expect_equal(p$segment, 0:3)
  # time is each segment's end; segment 0 sits at the origin.
  expect_equal(p$time, c(0, 12, 13, 24))
  # The baseline row takes the first segment's temperature, that being the
  # condition the unit is in when the schedule starts.
  expect_equal(p$temp, c(5, 5, 25, 5))
})

test_that("the baseline row is the model's fitted intercept", {
  fit <- base_fit()
  p <- predict_tos(fit, three_seg())

  # Two independent references: the ordinary prediction path at time 0, and the
  # fitted coefficient itself.
  # At time 0 the prediction is the intercept whatever the temperature, so this
  # reference is taken inside the fitted range to avoid an extrapolation warning.
  at0 <- as.numeric(predict(fit, newdata = data.frame(batch = "A", time = 0,
                                                     tempk = 298.15)))
  expect_equal(p$fit[1], at0, tolerance = 1e-8)
  expect_equal(p$fit[1], unname(coef(fit)["y0_batch_A"]), tolerance = 1e-8)
  expect_equal(p$fit[1], 100, tolerance = 1e-6)   # the value the data was built on
})

test_that("the baseline row follows the unit, not the schedule", {
  # The intercept is a property of the grouping level, so a different batch moves
  # it and a different schedule does not.
  fit <- base_fit()

  expect_equal(predict_tos(fit, three_seg("B"))$fit[1],
               unname(coef(fit)["y0_batch_B"]), tolerance = 1e-8)
  expect_equal(predict_tos(fit, three_seg("B"))$fit[1], 97, tolerance = 1e-6)

  # Same unit, unrelated schedule (and a different leading temperature): the
  # baseline fit is unchanged, since no time has elapsed at that row.
  other <- data.frame(t_start = 0, t_end = 6, temp = 40, batch = "A",
                      stringsAsFactors = FALSE)
  expect_equal(predict_tos(fit, other)$fit[1],
               predict_tos(fit, three_seg())$fit[1], tolerance = 1e-8)
})

test_that("the baseline row carries a real prediction interval", {
  # It must be usable in the same plot as the segment rows, so it needs a finite
  # interval of positive width that brackets the fit -- not a placeholder.
  g <- base_data()
  set.seed(42)
  g$y <- g$y + stats::rnorm(nrow(g), 0, 0.3)
  fit <- suppressWarnings(fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                                      intercept_cols = "batch"))
  b <- tos_baseline(predict_tos(fit, three_seg()))

  expect_equal(nrow(b), 1L)
  expect_true(all(is.finite(c(b$fit, b$lwr, b$upr))))
  expect_lt(b$lwr, b$fit)
  expect_gt(b$upr, b$fit)
})

test_that("the trajectory is monotone through the baseline row", {
  # The point of the row: the series starts at the intercept and degrades from
  # there, with no discontinuity between the baseline and the first segment.
  p <- predict_tos(base_fit(), three_seg())
  expect_true(all(diff(p$fit) < 0))
})

test_that("each variable of a multi-response fit gets its own baseline row", {
  fit <- suppressWarnings(fit_singledeg(base_data(), kinetic_order = 1,
                                        trend = "decreasing", Yinf = 0,
                                        intercept_cols = "batch"))
  ref <- predict(fit, newdata = data.frame(batch = "A", time = 0,
                                           tempk = 298.15))

  pA <- predict_tos(fit, three_seg(), var = "yA")
  pB <- predict_tos(fit, three_seg(), var = "yB")

  expect_equal(nrow(pA), 4L)
  expect_equal(nrow(pB), 4L)
  expect_equal(unique(pA$variable), "yA")
  expect_equal(unique(pB$variable), "yB")
  expect_equal(pA$fit[1], as.numeric(ref$fit_yA), tolerance = 1e-8)
  expect_equal(pB$fit[1], as.numeric(ref$fit_yB), tolerance = 1e-8)
})

test_that("the fit_sb path gets a baseline row too", {
  # Noise keeps m off its lower bound and log_Ea puts Ea on the scale of the
  # other parameters. On exact data with linear Ea the Hessian is so badly
  # conditioned that solve() in vcov() fails on some platforms (Windows).
  d <- base_data()
  set.seed(1)
  d$y <- d$y + stats::rnorm(nrow(d), 0, 0.3)
  fit <- fit_sb(d, rate_cols = "batch", intercept_cols = "batch",
                trend = "decreasing", Yinf = 0, log_Ea = TRUE)
  p <- predict_tos(fit, three_seg())

  expect_equal(nrow(p), 4L)
  expect_equal(p$segment, 0:3)
  # fit_sb starts its state at alpha = 1e-6 rather than exactly 0, in the fit and
  # in the composite alike, so the baseline sits that far off the intercept.
  expect_equal(p$fit[1], unname(coef(fit)["y0_A"]), tolerance = 1e-4)
  expect_true(all(diff(p$fit) < 0))
})

test_that("a single-segment schedule still returns two rows", {
  p <- predict_tos(base_fit(),
                   data.frame(t_start = 0, t_end = 12, temp = 25, batch = "A",
                              stringsAsFactors = FALSE))
  expect_equal(nrow(p), 2L)
  expect_equal(p$segment, 0:1)
  expect_equal(p$time, c(0, 12))
})
