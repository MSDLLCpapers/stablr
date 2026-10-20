# Two silent failures behind the "multi-column rate_cols explosion" report.
#
# 1. gsl_nls() records its exit status in convInfo (isConv, stopMessage,
#    finIter), but nothing looked at it. A fit that hit the iteration limit came
#    back with the same class as a converged one, so coef(), summary(), predict()
#    and every interval treated the last-visited parameter vector as an estimate,
#    with no warning and nothing in the printed output to notice. That is what
#    makes an RSS *worse than a nested submodel's* reportable: for a correctly
#    converged superset model it is a mathematical impossibility, so it is proof
#    of non-convergence -- but only if someone checks.
#
# 2. check_rate_support() covers cells with too few TIMEPOINTS, but a cell can
#    have many timepoints and still be unidentifiable if they all sit at one
#    TEMPERATURE: Ea is shared across cells while k_ref is per-cell, so any Ea
#    can be absorbed by rescaling that cell's k_ref. Crossing a second column
#    into rate_cols multiplies cells and shatters ones that had spanned the whole
#    temperature range into fragments that do not.

hmw <- function() subset(evers_supp, cqa == "HMW" & !is.na(y))

# fac2 aligned with temperature: every (batch, fac2) cell sits at one tempk.
one_temp_data <- function() {
  d <- hmw()
  d$fac2 <- ifelse(d$tempk >= 313.15, "L1", "L2")
  d
}

# fac2 interleaved within batch: cells still span the full temperature range.
spanning_data <- function() {
  d <- hmw()
  d <- d[order(d$batch, d$tempk, d$time), ]
  d$fac2 <- ifelse(ave(seq_len(nrow(d)), d$batch,
                       FUN = function(i) seq_along(i) %% 2) == 0, "L1", "L2")
  d
}

warnings_from <- function(expr) {
  ws <- character()
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

test_that("a fit that hits the iteration limit warns", {
  # Previously: returned silently with isConv = FALSE and no signal at all.
  ws <- warnings_from(
    fit_kinetic(hmw(), kinetic_order = 1, trend = "increasing",
                rate_cols = "batch", intercept_cols = "batch",
                control = list(maxiter = 2))
  )
  expect_gte(length(grep("did not converge", ws)), 1)
})

test_that("the non-convergence warning reports the reason", {
  msg <- grep("did not converge", warnings_from(
    fit_kinetic(hmw(), kinetic_order = 1, trend = "increasing",
                rate_cols = "batch", intercept_cols = "batch",
                control = list(maxiter = 2))
  ), value = TRUE)[1]

  expect_match(msg, "max number of iterations")
  expect_match(msg, "after 2 iterations")
})

test_that("the fit object is still returned so the user can inspect it", {
  # Warn, do not error: the user may want convInfo to diagnose the failure.
  m <- suppressWarnings(
    fit_kinetic(hmw(), kinetic_order = 1, trend = "increasing",
                rate_cols = "batch", intercept_cols = "batch",
                control = list(maxiter = 2))
  )
  expect_s3_class(m, "fit_kinetic")
  expect_false(m$convInfo$isConv)
})

test_that("a converged fit emits no convergence warning", {
  ws <- warnings_from(
    fit_kinetic(hmw(), kinetic_order = 1, trend = "increasing",
                rate_cols = "batch", intercept_cols = "batch")
  )
  expect_length(grep("did not converge", ws), 0)
})

test_that("a stall AT a stationary point does not warn", {
  # gsl_nls reports isConv = FALSE ("iteration is not making progress") when the
  # step search stalls because it is already at the solution. Here rh varies but
  # the model has no humidity term, so RSS 411.9 is the global optimum -- every
  # start below lands on it. isConv alone would cry wolf on ordinary lack of fit,
  # so the check confirms with the scaled gradient J'r.
  d <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15),
                   rh = c(0.4, 0.6, 0.75))
  d$y <- 100 * exp(-0.01 * exp(-80000 / 8.314 * (1 / d$tempk - 1 / 298.15)) *
                     exp(2 * (d$rh - 0.6)) * d$time)

  # An all-NA start is a gsl_nls self-start request; it stalls after one
  # iteration. This used to be what fit_kinetic passed by default, and it is
  # still the shortest way to provoke the stall.
  selfstart <- c(Ea = NA_real_, y0 = NA_real_, k_0 = NA_real_)
  m <- suppressWarnings(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                start = selfstart)
  )
  # Whether the self-start stalls (isConv = FALSE) is platform-dependent: some
  # toolchains (e.g. clang/flang on Debian) converge outright. Either way, no
  # non-convergence warning may be emitted at the optimum.
  expect_type(m$convInfo$isConv, "logical")
  expect_length(grep("did not converge", warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                start = selfstart)
  )), 0)                                    # ...but it is at the solution

  # A supplied start reaches the same RSS and reports success, confirming that
  # the stalled point really was the optimum.
  m2 <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0,
                    start = c(Ea = 80000, y0 = 100, k_0 = 0.01))
  expect_true(m2$convInfo$isConv)
  expect_equal(sum(residuals(m)^2), sum(residuals(m2)^2), tolerance = 1e-6)

  # And the seeded default now converges outright on the same data.
  m3 <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", Yinf = 0)
  expect_true(m3$convInfo$isConv)
  expect_equal(sum(residuals(m3)^2), sum(residuals(m2)^2), tolerance = 1e-6)
})

test_that("check_convergence is a no-op on a converged or convInfo-less object", {
  m <- fit_kinetic(hmw(), kinetic_order = 1, trend = "increasing",
                   rate_cols = "batch", intercept_cols = "batch")
  expect_silent(check_convergence(m, 13L))
  expect_silent(check_convergence(list(), 5L))
  expect_silent(check_convergence(list(convInfo = list()), 5L))
})

test_that("single-temperature rate cells are named in a warning", {
  ws <- warnings_from(
    fit_kinetic(one_temp_data(), kinetic_order = 1, trend = "increasing",
                rate_cols = c("batch", "fac2"), intercept_cols = "batch")
  )
  msg <- grep("only one", ws, value = TRUE)[1]
  expect_false(is.na(msg))
  expect_match(msg, "6 of 12 rate cell")
  expect_match(msg, "tempk")
  expect_match(msg, "fac2")
  # Must explain the Ea/k_ref confounding, not just count cells.
  expect_match(msg, "activation energy is shared")
  expect_match(msg, "residual sum")
})

test_that("cells that span temperatures do not warn", {
  # Multi-column rate_cols is legitimate when the cells retain temperature
  # spread; the check must not fire on every crossed design.
  ws <- warnings_from(
    fit_kinetic(spanning_data(), kinetic_order = 1, trend = "increasing",
                rate_cols = c("batch", "fac2"), intercept_cols = "batch")
  )
  expect_length(grep("only one", ws), 0)
})

test_that("crossing a factor does not make the fit worse than its submodel", {
  # The report's core symptom: a strict superset model with HIGHER RSS. Both
  # models must converge and the superset must not be worse.
  d <- spanning_data()
  m1 <- fit_kinetic(d, kinetic_order = 1, trend = "increasing",
                    rate_cols = "batch", intercept_cols = "batch")
  m2 <- fit_kinetic(d, kinetic_order = 1, trend = "increasing",
                    rate_cols = c("batch", "fac2"), intercept_cols = "batch")

  expect_true(m1$convInfo$isConv)
  expect_true(m2$convInfo$isConv)
  expect_lte(sum(residuals(m2)^2), sum(residuals(m1)^2) + 1e-8)
})

test_that("a single rate column aligned with temperature is caught too", {
  # One column carries the whole problem when it IS the temperature arm: every
  # cell of `fac2` sits at one tempk, so Ea is exactly unidentified. This used to
  # be exempt on the grounds that one column cannot shatter the design.
  d <- one_temp_data()
  msg <- grep("only one", warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "increasing",
                rate_cols = "fac2", intercept_cols = "batch")
  ), value = TRUE)[1]
  expect_false(is.na(msg))
  expect_match(msg, "aligned with")
  expect_match(msg, "Drop it from")

  # A single column whose levels each span the temperature range is fine.
  expect_length(grep("only one", warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "increasing",
                rate_cols = "batch", intercept_cols = "batch")
  )), 0)

  # A global rate has nothing to confound.
  expect_length(grep("only one", warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "increasing",
                intercept_cols = "batch")
  )), 0)
})

test_that("check_rate_temp_support handles degenerate input", {
  d <- one_temp_data()
  expect_silent(check_rate_temp_support(d, c("batch", "fac2"), "not_a_col"))
  expect_silent(check_rate_temp_support(d, c("batch", "nope"), "tempk"))
  expect_silent(check_rate_temp_support(d, character(), "tempk"))
  expect_silent(check_rate_temp_support(d, "batch", "tempk"))

  # An isothermal study has no Ea information anywhere: a different problem, and
  # not one this check should claim to have found.
  iso <- d[d$tempk == d$tempk[1], ]
  expect_silent(check_rate_temp_support(iso, c("batch", "fac2"), "tempk"))
})

test_that("the returned labels identify the affected cells", {
  labs <- suppressWarnings(
    check_rate_temp_support(one_temp_data(), c("batch", "fac2"), "tempk")
  )
  expect_length(labs, 6L)
  expect_true(all(grepl("batch=", labs)))
  expect_true(all(grepl("fac2=", labs)))
})

test_that("the other fitters also check convergence and temperature support", {
  d <- one_temp_data()
  d$yA <- d$y
  d$yB <- 100 - d$y

  ws <- warnings_from(
    fit_singledeg(d, kinetic_order = 1, trend = "increasing",
                  yA_var = "yA", yB_var = "yB",
                  rate_cols = c("batch", "fac2"), intercept_cols = "batch")
  )
  expect_gte(length(grep("only one", ws)), 1)

  ws2 <- warnings_from(
    fit_singledeg(d, kinetic_order = 1, trend = "increasing",
                  yA_var = "yA", yB_var = "yB",
                  rate_cols = "batch", intercept_cols = "batch",
                  control = list(maxiter = 2))
  )
  expect_gte(length(grep("did not converge", ws2)), 1)
})

# Negative activation energy
#
# Genuinely decreasing data, but the HOTTER arm degrades LESS -> the only
# least-squares solution is a negative Ea. The trend guard stays silent because
# the trend really is "decreasing".
inverted_arrhenius_data <- function() {
  set.seed(3)
  tempsK <- c(5, 25, 40) + 273.15
  times  <- c(0, 3, 6, 9, 12, 18, 24, 36)
  drop   <- c(`278.15` = 0.05, `298.15` = 0.03, `313.15` = 0.012)
  do.call(rbind, lapply(c("B1", "B2"), function(b) {
    g <- expand.grid(time = times, tempk = tempsK)
    r <- drop[as.character(g$tempk)]
    g$y  <- 100 * exp(-r * g$time) + rnorm(nrow(g), 0, 0.15)
    g$yA <- g$y
    g$yB <- 100 - g$y
    g$batch <- b
    g
  }))
}

test_that("a fit that converges to a negative Ea warns", {
  d <- inverted_arrhenius_data()
  msg <- grep("activation energy is negative", warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch")
  ), value = TRUE)[1]
  expect_false(is.na(msg))
  # The message must name the offending parameter and its value.
  expect_match(msg, "activation energy is negative")
  expect_match(msg, "Ea = -")
})

test_that("a normal positive-Ea fit emits no sign warning", {
  set.seed(1)
  d <- expand.grid(time = c(0, .25, .5, .75, 1),
                   tempk = c(25, 40, 50, 60) + 273.15)
  kT <- 3.0 * exp(-(90000 / 8.314) * (1 / d$tempk - 1 / 298.15))
  d$y <- 100 * exp(-kT * d$time) + rnorm(nrow(d), 0, 0.3)
  ws <- warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing")
  )
  expect_length(grep("activation energy is negative", ws), 0)
})

test_that("log_Ea excludes the sign check (Ea is structurally positive)", {
  # Ea_log is on the log scale, so a negative Ea_log just means Ea < 1 -- it is
  # not a negative activation energy and must not warn.
  d <- inverted_arrhenius_data()
  ws <- warnings_from(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch", log_Ea = TRUE)
  )
  expect_length(grep("activation energy is negative", ws), 0)
})

test_that("fit_singledeg also warns on a negative Ea", {
  d <- inverted_arrhenius_data()
  ws <- warnings_from(
    fit_singledeg(d, kinetic_order = 1, trend = "decreasing",
                  yA_var = "yA", yB_var = "yB", intercept_cols = "batch")
  )
  expect_gte(length(grep("activation energy is negative", ws)), 1)
})

test_that("check_activation_energy_sign matches only plain, negative Ea", {
  fired <- function(x) {
    ws <- character()
    withCallingHandlers(check_activation_energy_sign(x),
      warning = function(w) { ws <<- c(ws, conditionMessage(w));
                              invokeRestart("muffleWarning") })
    length(ws) > 0
  }
  expect_true(fired(c(Ea = -100, k_0 = 0.1)))
  expect_true(fired(c(EaB = -5, EaC = 100)))     # branch-suffixed parallel name
  expect_false(fired(c(Ea = 90000, k_0 = 0.1)))
  expect_false(fired(c(Ea_log = -12, k_0 = 0.1)))  # kinetic log spelling
  expect_false(fired(c(Ealog = -3, k = 0.1)))      # fit_sb log spelling
  expect_false(fired(c(Ea = NA_real_)))
  expect_false(fired(stats::setNames(-1, NULL)))
})

