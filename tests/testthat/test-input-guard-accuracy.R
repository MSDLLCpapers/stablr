# Guards that used to refuse legitimate data, accept invalid data, or name the
# wrong cause. One test per defect, each asserting the behaviour the guard was
# missing rather than the wording it happens to use.

pur <- function() subset(evers_supp, cqa == "Purity" & !is.na(y))

warnings_from <- function(expr) {
  ws <- character()
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

test_that("a duplicated `start` name is refused, not silently doubled", {
  # Previously: coef() came back with two entries of the same name and DIFFERENT
  # values (A0 = 100 alongside A0 = 86.6 from fit_singledeg), the parameter count
  # was inflated so sigma and the residual df were wrong, and the only signal was
  # a generic "singular gradient matrix".
  expect_error(
    fit_kinetic(pur(), kinetic_order = 1, trend = "decreasing",
                start = c(Ea = 80000, k_0 = 0.01, y0 = 100, y0 = 50)),
    "Duplicated name"
  )

  d <- pur(); d$yA <- d$y; d$yB <- 100 - d$y
  expect_error(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  start = c(Ea = 80000, k_0 = 0.01, A0 = 100, B0 = 1, A0 = 5)),
    "Duplicated name"
  )
})

test_that("a character `start` is refused instead of being discarded", {
  # gsl_nls drops a non-numeric start and self-starts, so quoting the numbers
  # meant the supplied values were ignored and the fit reported nothing.
  expect_error(
    fit_kinetic(pur(), kinetic_order = 1, trend = "decreasing",
                start = c(Ea = "80000", k_0 = "0.01", y0 = "100")),
    "must be numeric"
  )
  # An all-NA start is a self-start request, which default_start() itself emits.
  expect_no_error(
    suppressWarnings(fit_kinetic(pur(), kinetic_order = 1, trend = "decreasing",
                                 start = c(Ea = NA, k_0 = NA, y0 = NA)))
  )
})

test_that("levels differing only by surrounding whitespace are one level", {
  # Previously refused outright as a punctuation collision, so a batch column
  # carrying 'Lot-1' and 'Lot-1 ' from a spreadsheet export could not be fitted
  # at all. Internal punctuation must still collide -- see
  # test-predict-sanitize-values.R.
  d <- pur()
  b1 <- unique(d$batch)[1]
  rows <- which(d$batch == b1)
  d$batch[rows[seq(1, length(rows), by = 2)]] <- paste0(b1, " ")

  fit <- expect_no_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch")
  )
  # One parameter for the one real level, not two and not none.
  expect_length(grep("^y0_batch_", names(coef(fit))), length(unique(d$batch)) - 1L)

  # And predicting from either spelling works.
  nd <- data.frame(time = 12, tempk = 298.15, batch = c(b1, paste0(b1, " ")))
  expect_no_error(predict(fit, nd))
})

test_that("the trend guard uses the within-temperature slope", {
  # Pooling across temperatures cancelled the trend out: on this data the pooled
  # correlation is -0.07, under the guard's own 0.1 floor, while every
  # within-temperature correlation is below -0.94. The reversed trend therefore
  # went unwarned and died inside the optimiser instead.
  d <- pur()
  expect_lt(abs(stats::cor(d$time, d$y)), 0.1)

  ws <- warnings_from(try(
    fit_kinetic(d, kinetic_order = 1, trend = "increasing"), silent = TRUE))
  expect_gte(length(grep("Declared trend", ws)), 1)

  # The correct declaration stays quiet.
  expect_length(grep("Declared trend", warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing"))), 0)
})

test_that("the trend guard requires every temperature to contradict `trend`", {
  # Previously the verdict was the median of the within-temperature
  # correlations, so a minority of arms could carry it: four rising arms
  # outvoted one genuinely falling arm and the guard fired on data that does not
  # unambiguously disagree. Only unanimity is enough evidence.
  arm <- function(Tk, slope, sigma = 0.2) {
    tt <- rep(c(0, 3, 6, 9, 12, 18, 24), 3)
    data.frame(tempk = Tk, time = tt,
               y = 98 + slope * tt + stats::rnorm(length(tt), 0, sigma))
  }
  temps <- c(278.15, 298.15, 303.15, 310.15, 313.15)
  mk <- function(slopes) do.call(rbind, Map(arm, temps, slopes))

  set.seed(42)
  rising  <- mk(c(0.3, 0.4, 0.5, 0.6, 0.7))
  falling <- mk(-c(0.3, 0.4, 0.5, 0.6, 0.7))
  mixed   <- mk(c(0.3, 0.4, 0.5, 0.6, -0.5))   # one arm disagrees
  flat    <- mk(c(0.3, 0.4, 0.5, 0.6, 0))      # one arm carries no direction

  fires <- function(d, trend) {
    length(grep("Declared trend",
                warnings_from(check_trend_direction(d, "y", "time", trend, "tempk")))) > 0
  }

  # Unanimous contradiction, in both directions.
  expect_true(fires(rising, "decreasing"))
  expect_true(fires(falling, "increasing"))

  # Unanimous agreement.
  expect_false(fires(falling, "decreasing"))
  expect_false(fires(rising, "increasing"))

  # Not unanimous: one dissenting arm, and one near-flat arm, each suppress it.
  expect_false(fires(mixed, "decreasing"))
  expect_false(fires(flat, "decreasing"))
})

test_that("`weights` and `subset` reach gsl_nls", {
  # Every spelling used to abort with "..1 used in an incorrect context, no ...
  # to look in", because gsl_nls re-evaluates them in the caller's frame.
  d <- pur(); d$w <- ifelse(d$tempk > 300, 2, 1)

  mw <- expect_no_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing", weights = w))
  m  <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing")
  expect_false(isTRUE(all.equal(coef(mw)[["Ea"]], coef(m)[["Ea"]])))

  ms <- expect_no_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing", subset = time < 12))
  expect_lt(length(residuals(ms)), nrow(d))

  # print() must not deparse the whole model frame.
  expect_lt(length(grep("Combi", capture.output(print(mw)))), 1L)
})

test_that("a non-syntactic grouping column name is named, not left to the parser", {
  # Previously: "<text>:1:41: unexpected symbol", pointing into the generated
  # formula and naming neither the column nor the reason.
  d <- pur(); d[["test site"]] <- rep(c("A", "B"), length.out = nrow(d))
  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "test site"),
    "not valid R names"
  )
})

test_that("the single-temperature warning is silent when nothing estimates Ea", {
  # A temp_expr with no free parameter has nothing unidentifiable, so warning
  # about "an arbitrary, start-driven Ea" named a cause that did not exist: the
  # fit's coef() is y0 and k_0 only.
  iso <- subset(pur(), tempk == sort(unique(pur()$tempk))[3])
  flat <- expr_arrhenius()
  flat$formula <- "1"

  ws <- warnings_from(fit_kinetic(iso, kinetic_order = 1, trend = "decreasing",
                                 temp_expr = flat))
  expect_length(grep("not identifiable", ws), 0)

  # With an Ea to estimate it still warns, and names it.
  ws2 <- warnings_from(fit_kinetic(iso, kinetic_order = 1, trend = "decreasing"))
  msg <- grep("not identifiable", ws2, value = TRUE)[1]
  expect_false(is.na(msg))
  expect_match(msg, "Ea")
})
