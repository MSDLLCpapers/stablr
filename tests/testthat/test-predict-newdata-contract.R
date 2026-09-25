# The predict() newdata contract.
#
# Three defects shared one root cause: the guarantees a caller depends on lived
# only in the stablr predict.* methods, so anything that reached the underlying
# object by another route silently lost them.
#
#   * W11-01 -- the indicator expressions compared a raw grouping column against
#     a sanitize()d level value. predict.fit_kinetic() sanitized newdata first,
#     so the fitters were shielded; the inherited stats::predict.nls(), reached
#     whenever the fit_kinetic class is lost, was not. Every indicator was FALSE
#     for every level and the sum came back as exactly 0 -- total degradation,
#     the worst available wrong answer for a stability forecast, with no warning.
#     Fixed in the expression itself (see sanitized_col()) rather than in the
#     method, so it holds on whichever path evaluates it.
#   * H05-04 -- a newdata column sharing a name with a fitted parameter outranks
#     the estimate during evaluation, silently replacing the fit.
#   * W04-05 -- extrapolating outside the fitted temperature range is exponential
#     in 1/T, so it moves the answer by a large factor while staying numerically
#     well behaved, and said nothing.

# ── helpers ────────────────────────────────────────────────────────────────

# Levels carry a "-", which sanitize() strips: this is what made the raw column
# and the baked-in level value disagree. Three temperatures so temp_range has an
# interior, and 298.15 is included so trefk is not an endpoint.
TEMPS <- c(298.15, 308.15, 313.15)

sim_kin <- function(k = 0.02, Ea = 80000, y0 = c(`A-1` = 100, B = 90)) {
  do.call(rbind, lapply(names(y0), function(b) {
    g  <- expand.grid(time = c(0, 3, 6, 12), tempk = TEMPS)
    kk <- k * exp(Ea / 8.314 * (1 / 298.15 - 1 / g$tempk))
    data.frame(batch = b, time = g$time, tempk = g$tempk,
               y = y0[[b]] * exp(-kk * g$time))
  }))
}

sim_single <- function() {
  d <- sim_kin(y0 = c(`A-1` = 100, B = 100))
  d$yA <- d$y
  d$yB <- 100 - d$y
  d$y  <- NULL
  d
}

sim_par <- function() {
  do.call(rbind, lapply(c("A-1", "B"), function(b) {
    g  <- expand.grid(time = c(0, 1, 3, 6, 9, 12), tempk = TEMPS)
    kB <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    kC <- 0.01 * exp(70000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    kA <- kB + kC
    yA <- 100 * exp(-kA * g$time)
    data.frame(batch = b, time = g$time, tempk = g$tempk, yA = yA,
               yB = (kB / kA) * (100 - yA), yC = (kC / kA) * (100 - yA))
  }))
}

# ── W11-01: the expression must be correct off the stablr method ───────────

test_that("a fit that loses its stablr class still predicts correctly", {
  # The trigger is the loss of the class, not serialization as such: an RDS
  # round-trip auto-loads the stablr namespace (the `build` attribute holds
  # closures whose environments live there), so dispatch survives it. Any layer
  # that rebuilds or strips the class -- a caching layer, a manual
  # stats::predict.nls() call -- does not.
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = "batch", intercept_cols = "batch")
  nd  <- data.frame(batch = "A-1", time = 12, tempk = 298.15)

  expected <- 100 * exp(-0.02 * 12)
  expect_equal(unname(predict(fit, newdata = nd)), expected, tolerance = 1e-4)

  bare <- fit
  class(bare) <- setdiff(class(fit), "fit_kinetic")
  got <- unname(predict(bare, newdata = nd))

  # The specific wrong answer this replaces. Asserted explicitly because 0 is
  # not an implausible-looking stability result -- it reads as total loss.
  expect_false(isTRUE(all.equal(got, 0)))
  expect_equal(got, expected, tolerance = 1e-4)
})

test_that("an RDS round-trip predicts correctly", {
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = "batch")
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(fit, f)

  nd <- data.frame(batch = "A-1", time = 12, tempk = 308.15)
  expect_equal(unname(predict(readRDS(f), newdata = nd)),
               unname(predict(fit, newdata = nd)), tolerance = 1e-8)
})

test_that("built expressions evaluate on raw data, not only sanitized data", {
  # expr_intercept() and expr_rate() are exported, so a caller composing a model
  # by hand never passes through the fitters' sanitization step at all.
  d <- data.frame(batch = c("A-1", "B"), stringsAsFactors = FALSE)

  ie <- expr_intercept(d, "batch")
  re <- expr_rate(d, "batch")

  ev <- function(x) {
    p <- stats::setNames(rep(7, length(x$params)), x$params)
    rlang::eval_tidy(rlang::parse_expr(x$expr), data = d,
                     env = list2env(as.list(p), parent = baseenv()))
  }
  # Every row must be selected by exactly one indicator. Summing to 0 is the
  # failure; each row picking up 7 is the contract.
  expect_equal(ev(ie), c(7, 7))
  expect_equal(ev(re), c(7, 7))

  # ... and equally on data that is already sanitized, since sanitize() is
  # idempotent and both sides of the comparison go through it.
  d2 <- data.frame(batch = c("A1", "B"), stringsAsFactors = FALSE)
  ev2 <- function(x) {
    p <- stats::setNames(rep(7, length(x$params)), x$params)
    rlang::eval_tidy(rlang::parse_expr(x$expr), data = d2,
                     env = list2env(as.list(p), parent = baseenv()))
  }
  expect_equal(ev2(ie), c(7, 7))
  expect_equal(ev2(re), c(7, 7))
})

test_that("a numeric grouping column with decimals still matches its rows", {
  # Sanitizing the column means comparing character, so the target must be a
  # quoted string. It briefly was not: sanitize() strips the "." out of 0.4,
  # leaving the all-digit "04", which format_level() emitted bare as the number
  # 4 -- and "04" == 4 is FALSE, so every indicator was FALSE and the expression
  # evaluated to 0. Same failure mode as W11-01, entering through the comparison
  # target rather than the column. See sanitized_level().
  d <- data.frame(rh = c(0.4, 0.75))
  for (x in list(expr_intercept(d, "rh"), expr_rate(d, "rh"))) {
    p <- stats::setNames(rep(7, length(x$params)), x$params)
    expect_equal(
      rlang::eval_tidy(rlang::parse_expr(x$expr), data = d,
                       env = list2env(as.list(p), parent = baseenv())),
      c(7, 7))
  }

  # Leading zeros must survive too -- an unquoted 0000803274 would read as
  # 803274. This is what quote_numeric_levels() used to repair after the fact.
  z <- data.frame(batch = c("0000803274", "0000803279"),
                  stringsAsFactors = FALSE)
  zr <- expr_rate(z, "batch")
  expect_match(zr$expr, "'0000803274'", fixed = TRUE)
  p <- stats::setNames(rep(7, length(zr$params)), zr$params)
  expect_equal(rlang::eval_tidy(rlang::parse_expr(zr$expr), data = z,
                                env = list2env(as.list(p),
                                               parent = baseenv())),
               c(7, 7))
})

test_that("a numeric grouping column fits and predicts correctly end to end", {
  g <- do.call(rbind, lapply(c(0.4, 0.75), function(rh) {
    gg <- expand.grid(time = c(0, 3, 6, 12), tempk = TEMPS)
    k  <- (if (rh > 0.5) 0.03 else 0.01) *
      exp(80000 / 8.314 * (1 / 298.15 - 1 / gg$tempk))
    data.frame(rh = rh, time = gg$time, tempk = gg$tempk,
               y = 100 * exp(-k * gg$time))
  }))
  m <- fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                   rate_cols = "rh")
  expect_lt(max(abs(stats::residuals(m))), 1e-8)

  nd <- data.frame(rh = 0.75, time = 12, tempk = 298.15)
  expect_equal(unname(predict(m, newdata = nd)), 100 * exp(-0.03 * 12),
               tolerance = 1e-6)
  # and off the stablr method, where newdata is not pre-sanitized
  bare <- m
  class(bare) <- setdiff(class(m), "fit_kinetic")
  expect_equal(unname(predict(bare, newdata = nd)),
               100 * exp(-0.03 * 12), tolerance = 1e-6)
})

# ── H05-04: a newdata column may not shadow a fitted parameter ────────────

test_that("a newdata column named after a fitted parameter is refused", {
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = "batch", intercept_cols = "batch")
  pars <- attr(fit, "build")$param_names
  base <- data.frame(batch = "A-1", time = 12, tempk = 313.15)

  # Every parameter, not a subset: the data mask outranks the parameter
  # environment for all of them.
  expect_true(all(c("Ea", "k_batch_A1", "y0_batch_A1") %in% pars))
  for (p in pars) {
    nd <- base
    nd[[p]] <- 999
    expect_error(predict(fit, newdata = nd), p, fixed = TRUE, info = p)
  }
})

test_that("the shadowing refusal is uniform across all four fitters", {
  base <- data.frame(batch = "A-1", time = 12, tempk = 313.15)

  fits <- list(
    kinetic   = fit_kinetic(sim_kin(), kinetic_order = 1,
                            trend = "decreasing", rate_cols = "batch"),
    singledeg = fit_singledeg(sim_single(), kinetic_order = 1,
                              trend = "decreasing", rate_cols = "batch"),
    parallel  = fit_parallel(sim_par(), kinetic_order = 1,
                             rate_cols = "batch"),
    # fit_sb assembles its own named parameter vector rather than evaluating an
    # expression against a data mask, so it is structurally immune. The check is
    # applied there anyway: the contract should read the same for all four, not
    # "the same, except in fit_sb".
    sb        = fit_sb(sim_kin(), trend = "decreasing", rate_cols = "batch")
  )

  for (nm in names(fits)) {
    # The shadowing name is read off each fit rather than assumed. Hardcoding
    # `Ea` gave a false negative for fit_parallel, which renames its temperature
    # parameters per branch (EaB / EaC) and so has no parameter called `Ea` to
    # shadow -- the check was right not to fire, and the test was wrong.
    pars <- attr(fits[[nm]], "build")$param_names
    expect_true(length(pars) > 0L, info = nm)
    for (p in pars) {
      nd <- base
      nd[[p]] <- 999
      expect_error(predict(fits[[nm]], newdata = nd),
                   "share a name with a fitted parameter",
                   info = paste(nm, p))
      expect_error(predict(fits[[nm]], newdata = nd), p, fixed = TRUE,
                   info = paste(nm, p))
    }
  }
})

test_that("an unrelated extra newdata column is still allowed", {
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = "batch")
  nd  <- data.frame(batch = "A-1", time = 12, tempk = 298.15,
                    analyst = "kf", lot_note = 1)
  expect_silent(p <- predict(fit, newdata = nd))
  expect_length(p, 1L)
})

# ── W04-05: temperature extrapolation must not be silent ──────────────────

test_that("all four fitters record the fitted temperature range", {
  fits <- list(
    fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing"),
    fit_singledeg(sim_single(), kinetic_order = 1, trend = "decreasing"),
    fit_parallel(sim_par(), kinetic_order = 1),
    fit_sb(sim_kin(), trend = "decreasing")
  )
  for (f in fits) {
    expect_equal(attr(f, "build")$settings$temp_range, range(TEMPS))
  }
})

test_that("extrapolating outside the fitted temperature range warns", {
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = "batch")
  nd <- function(tk) data.frame(batch = "A-1", time = 12, tempk = tk)

  # Both directions. The upward direction is the dangerous one -- 20 K above the
  # fitted maximum returns a number indistinguishable from total degradation --
  # and was the direction the original report did not mention.
  expect_warning(predict(fit, newdata = nd(278.15)), "extrapolat")
  expect_warning(predict(fit, newdata = nd(333.15)), "extrapolat")
  expect_warning(predict(fit, newdata = nd(333.15)), "298.15 to 313.15")

  # Inside the range, including at the endpoints, it must stay quiet: a warning
  # on ordinary interpolation would be trained away within a day.
  for (tk in c(298.15, 305, 313.15)) {
    expect_silent(predict(fit, newdata = nd(tk)))
  }
  # A newdata mixing interior and exterior temperatures warns once, naming both.
  expect_warning(predict(fit, newdata = nd(c(278.15, 305, 333.15))),
                 "as low as 278.15.*as high as 333.15")
})

test_that("the extrapolation warning does not suppress the prediction", {
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing")
  nd  <- data.frame(time = 12, tempk = 333.15)
  p   <- suppressWarnings(predict(fit, newdata = nd))
  expect_length(p, 1L)
  expect_true(is.finite(p))
})

test_that("a fit recorded before temp_range existed does not warn", {
  # settings is carried on the fitted object, so an object built by an older
  # version has no temp_range. That must skip the check, not guess a range or
  # error -- otherwise the fix breaks every stored fit.
  fit <- fit_kinetic(sim_kin(), kinetic_order = 1, trend = "decreasing")
  old <- fit
  b <- attr(old, "build")
  b$settings$temp_range <- NULL
  attr(old, "build") <- b

  expect_silent(predict(old, newdata = data.frame(time = 12, tempk = 333.15)))
})

test_that("fit_sb warns on temperature extrapolation too", {
  fit <- fit_sb(sim_kin(), trend = "decreasing")
  expect_warning(predict(fit, newdata = data.frame(time = 12, tempk = 333.15)),
                 "extrapolat")
  expect_silent(predict(fit, newdata = data.frame(time = 12, tempk = 305)))
})
