# Every humidity guard ran at fit time only. check_rh_scale(), check_numeric_cols()
# and check_rh_identifiable() were all called from the four fitters and from
# nowhere else, so a model fitted correctly could still be handed a percent-coded
# rh or a factor tempk at predict time and answer without complaint. The two
# halves of one contract disagreed about what counts as acceptable input.

rh_te <- expr_arrhenius_humidity()          # rh_ref = 0.60, fraction scale

rh_data <- function(rh_levels = c(0.40, 0.60, 0.75), sd = 0, seed = 1) {
  g <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 313.15),
                   rh = rh_levels, stringsAsFactors = FALSE)
  k <- 0.01 * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk)) *
    exp(2 * (g$rh - 0.60))
  set.seed(seed)
  g$y <- 100 * exp(-k * g$time) + stats::rnorm(nrow(g), 0, sd)
  g
}

rh_fit <- function(d = rh_data(), ...) {
  suppressWarnings(fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                               temp_expr = rh_te, ...))
}

# Truth, for checking the guards did not disturb correct answers.
rh_truth <- function(time, rh, tempk = 298.15) {
  100 * exp(-0.01 * exp(80000 / 8.314 * (1 / 298.15 - 1 / tempk)) *
              exp(2 * (rh - 0.60)) * time)
}

# ── W07-02: a percent-coded rh at predict time ───────────────────────────────

test_that("percent-coded rh in newdata is refused, not answered with 0", {
  # `b` is already fixed by the fit, so exp(2 * (60 - 0.60)) is e^119 ~ 1e51. The
  # rate is scaled by that, exp(-k * t) underflows and predict() returned exactly
  # 0 -- a forecast of total degradation, which is plausible enough for a
  # stability study that nothing looks wrong.
  fit <- rh_fit()
  expect_equal(unname(coef(fit)["b"]), 2, tolerance = 1e-3)

  err <- expect_error(predict(fit, newdata = data.frame(time = 12,
                                                        tempk = 298.15, rh = 60)))
  expect_match(conditionMessage(err), "percent-coded")
  # The reference is the humidity the model was FITTED on, not rh_ref -- see the
  # percent-fitted test below for why rh_ref is the wrong yardstick here.
  expect_match(conditionMessage(err), "fitted on humidity ranging 0.4 to 0.75")
  # The message must say what the old behaviour was, or the reader assumes the
  # prediction merely errored rather than that it used to be silently wrong.
  expect_match(conditionMessage(err), "exactly 0")
})

test_that("percent-coded rh in a predict_tos schedule is refused", {
  # rh reaches predict_tos() as a plain temperature-expression covariate, so it
  # met no guard at all on this path.
  fit <- rh_fit()
  expect_error(
    predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 25, rh = 60)),
    "percent-coded"
  )
})

test_that("a fraction-coded rh against a percent-scale fit is also refused", {
  # The mirror image: rh_ref = 60 means the model wants percent, so 0.60 is the
  # mistake. Both directions produce a wrong number, so both must be caught.
  d   <- transform(rh_data(), rh = rh * 100)
  fit <- suppressWarnings(fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                                      temp_expr = expr_arrhenius_humidity(rh_ref = 60)))
  expect_error(
    predict(fit, newdata = data.frame(time = 12, tempk = 298.15, rh = 0.60)),
    "fraction-coded"
  )
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                 rh = 60)))
})

test_that("a model fitted on percent rh accepts percent rh at predict time", {
  # The predict-side check must judge newdata against the scale the model was
  # FITTED on, not against rh_ref. Pairing percent rh with rh_ref = 0.60 is an
  # exact reparameterisation -- internally consistent, which is why the fit only
  # warns -- so percent rh is the CORRECT input for that model's predictions, and
  # comparing it against rh_ref would refuse the one coding that works.
  d   <- transform(rh_data(), rh = rh * 100)
  fit <- rh_fit(d)                                    # rh_ref = 0.60, percent data
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                  rh = 60)))
  # Predictions still track the reparameterised truth: b absorbed the 1/100.
  p <- predict(fit, newdata = data.frame(time = 12, tempk = 298.15, rh = 75))
  expect_equal(unname(p), rh_truth(12, 0.75), tolerance = 1e-3)
  # And the fraction coding is now the mistake, on this model.
  expect_error(predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                 rh = 0.60)),
               "fraction-coded")
})

test_that("the fit-time check still WARNS rather than erroring", {
  # A mis-scaled fit is an exact reparameterisation: internally consistent, just
  # differently parameterised, and refusing it would break existing analyses. It
  # is the prediction that is wrong, so only the predict side is fatal.
  d <- transform(rh_data(), rh = rh * 100)
  expect_warning(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing", temp_expr = rh_te),
    "percent-coded"
  )
})

# ── W07-07: predictor column types at predict time ───────────────────────────

test_that("a factor rh or tempk in newdata is refused with a named message", {
  # Previously returned NA behind nothing but base R's "'-' not meaningful for
  # factors" -- which names neither the column nor predict() -- and the NA then
  # travelled on as a normal-looking prediction.
  fit <- rh_fit()
  for (col in c("rh", "tempk", "time")) {
    nd <- data.frame(time = 12, tempk = 298.15, rh = 0.60)
    nd[[col]] <- factor(nd[[col]])
    err <- expect_error(predict(fit, newdata = nd), info = col)
    expect_match(conditionMessage(err), col, info = col)
    expect_match(conditionMessage(err), "must be numeric", info = col)
  }
})

test_that("a character rh is refused too", {
  fit <- rh_fit()
  expect_error(
    predict(fit, newdata = data.frame(time = 12, tempk = 298.15, rh = "0.60")),
    "must be numeric"
  )
})

test_that("a factor rh in a predict_tos schedule is refused", {
  fit <- rh_fit()
  expect_error(
    predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 25,
                                rh = factor(0.60))),
    "must be numeric"
  )
})

# ── W07-03: the identifiability tolerance is a design threshold ───────────────

test_that("a humidity spread too small to estimate b is warned about", {
  # The tolerance was 1e-8, which only caught humidity constant to machine
  # precision. A spread of 1e-4 -- a rounding artefact in a nominally
  # single-humidity study -- passed, and b came back -12.7 against a truth of 2
  # with the RH-shifted forecast 15% wrong and nothing warned.
  for (spread in c(1e-6, 1e-4, 1e-3)) {
    d <- rh_data(c(0.60, 0.60 + spread), sd = 0.3)
    expect_warning(
      fit_kinetic(d, kinetic_order = 1, trend = "decreasing", temp_expr = rh_te),
      "not identifiable",
      info = format(spread)
    )
  }
})

test_that("the warning distinguishes a seed from a noise-fitted b", {
  # With no spread the optimiser cannot move b and returns the seed verbatim;
  # with a little spread it chases noise and returns something large. Calling the
  # second one "the starting value" would send the user hunting for a seed
  # nowhere near what they got.
  #
  # A perfectly constant rh aliases b with the rate, so nls also emits an
  # (expected) "singular gradient matrix" warning; muffle only that one so the
  # identifiability message is the thing under test.
  muffle_singular <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("singular gradient", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  }
  const <- muffle_singular(expect_warning(
    fit_kinetic(rh_data(0.60), kinetic_order = 1, trend = "decreasing",
                temp_expr = rh_te),
    "is constant"
  ))
  expect_match(conditionMessage(const), "starting value, not an estimate")

  narrow <- muffle_singular(expect_warning(
    fit_kinetic(rh_data(c(0.60, 0.6001), sd = 0.3), kinetic_order = 1,
                trend = "decreasing", temp_expr = rh_te),
    "spans only"
  ))
  expect_match(conditionMessage(narrow), "fitted to noise")
  expect_match(conditionMessage(narrow), "below the")
})

test_that("a genuinely informative humidity spread is not warned about", {
  # One percentage point of RH is the threshold; real designs span 35 of them.
  for (rh in list(c(0.60, 0.61), c(0.40, 0.60, 0.75))) {
    d  <- rh_data(rh, sd = 0.3)
    ws <- character()
    withCallingHandlers(
      fit_kinetic(d, kinetic_order = 1, trend = "decreasing", temp_expr = rh_te),
      warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }
    )
    expect_length(grep("not identifiable", ws), 0)
  }
})

test_that("the threshold scales with the humidity coding", {
  # 0.01 on the fraction scale and 0.75 on a percent scale topping out at 75 are
  # the same physical spread, so neither coding may be privileged.
  d_frac <- rh_data(c(0.60, 0.605), sd = 0.3)
  d_pct  <- transform(d_frac, rh = rh * 100)
  te_pct <- expr_arrhenius_humidity(rh_ref = 60)

  expect_warning(fit_kinetic(d_frac, kinetic_order = 1, trend = "decreasing",
                             temp_expr = rh_te), "not identifiable")
  expect_warning(fit_kinetic(d_pct, kinetic_order = 1, trend = "decreasing",
                             temp_expr = te_pct), "not identifiable")
})

# ── the guards must not disturb correct humidity models ──────────────────────

test_that("a correctly coded humidity model predicts exactly as before", {
  fit <- rh_fit()
  for (rh in c(0.40, 0.60, 0.75)) {
    p <- predict(fit, newdata = data.frame(time = 12, tempk = 298.15, rh = rh))
    expect_equal(unname(p), rh_truth(12, rh), tolerance = 1e-6)
  }
  # Extrapolating to an unobserved humidity is still allowed -- it is a modelling
  # judgement, not a coding error, and the scale is right.
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                  rh = 0.90)))
  # predict_tos() on the same fit, across a temperature excursion.
  p <- tos_segments(predict_tos(fit, data.frame(t_start = c(0, 6), t_end = c(6, 12),
                                                temp = c(25, 40), rh = 0.60)))
  expect_equal(nrow(p), 2L)
  expect_true(all(is.finite(c(p$fit, p$lwr, p$upr))))
  expect_lt(p$fit[2], p$fit[1])
})

test_that("a model with no humidity term is untouched", {
  # settings$rh_col is NULL there, so both new checks must no-op rather than
  # reach into a column that does not exist.
  fit <- suppressWarnings(fit_kinetic(rh_data(), kinetic_order = 1,
                                      trend = "decreasing"))
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 298.15)))
  expect_silent(predict_tos(fit, data.frame(t_start = 0, t_end = 12, temp = 25)))
  # An rh column present but unused by the model is not the guard's business.
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                 rh = 60)))
})

test_that("fit_sb's humidity path gets the same predict-side guards", {
  d   <- rh_data(sd = 0.3)
  fit <- suppressWarnings(fit_sb(d, temp_expr = rh_te, trend = "decreasing", Yinf = 0))
  expect_error(
    predict(fit, newdata = data.frame(time = 12, tempk = 298.15, rh = 60)),
    "percent-coded"
  )
  expect_error(
    predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                      rh = factor(0.60))),
    "must be numeric"
  )
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                                 rh = 0.60)))
})
