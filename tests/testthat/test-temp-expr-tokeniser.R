# Temperature-expression tokenising and parameter naming
#
# `temp_expr_params()` decides which symbols in a temperature expression are
# free parameters and which are data columns. It used to do that by splitting
# the string on non-alphanumerics and subtracting a whitelist of known function
# names, which misclassified in both directions:
#
#   * function names, numeric literals in scientific notation and logical
#     constants were read as parameters (V10-05);
#   * any data column the expression referenced other than temp_var/rh_col was
#     read as a parameter too, so a legitimate custom expression either could
#     not be written at all (V09-12, V10-10) or was silently mis-fitted
#     (R03-03).
#
# The fix is to classify with R's own parser -- all.vars() on the parsed
# expression -- and to have the expression *declare* its data columns via
# `covar_cols`, which resolve_temp_expr() threads through to the builders.

sim_moist <- function(k = 0.02, Ea = 80000, bm = 1.5, ref = 0.6) {
  d <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 308.15, 313.15),
                   moist = c(0.4, 0.6, 0.75))
  d$y <- 100 * exp(-k * exp(Ea / 8.314 * (1 / 298.15 - 1 / d$tempk)) *
                     exp(bm * (d$moist - ref)) * d$time)
  d
}

sim_parallel_moist <- function() {
  d <- expand.grid(time = c(0, 1, 3, 6, 9, 12), tempk = c(298.15, 313.15),
                   moist = c(0.4, 0.75))
  kB <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / d$tempk)) *
    exp(1.5 * (d$moist - 0.6))
  kC <- 0.01 * exp(70000 / 8.314 * (1 / 298.15 - 1 / d$tempk)) *
    exp(1.0 * (d$moist - 0.6))
  kA <- kB + kC
  d$yA <- 100 * exp(-kA * d$time)
  d$yB <- (kB / kA) * (100 - d$yA)
  d$yC <- (kC / kA) * (100 - d$yA)
  d
}

moist_expr <- function(...) {
  utils::modifyList(
    list(formula  = paste("exp(Ea/8.314 * (1/298.15 - 1/tempk)) *",
                          "exp(bm * (moist - 0.6))"),
         R = 8.314, trefk = 298.15, temp_var = "tempk", log_Ea = FALSE,
         covar_cols = "moist"),
    list(...))
}

# ── V10-05: phantom parameters from the whitelist tokeniser ────────────────

test_that("function names, scientific-notation literals and constants are not parameters", {
  # Each of these used to leak a phantom parameter into the fitted object: a
  # free coefficient with no gradient, inflating sigma and deflating residual df.
  expect_equal(temp_expr_params("exp(Ea/8.314 * (1/298.15 - 1/tempk))", "tempk"),
               "Ea")
  # "1e-9" split into "1e" and "9" on the "-"; "1e" survived as a parameter.
  expect_equal(temp_expr_params("exp(Ea/8.314 * (1e-9 - 1/tempk))", "tempk"),
               "Ea")
  expect_equal(temp_expr_params("exp(Ea * 2.5e-3 / tempk)", "tempk"), "Ea")
  # Functions outside the whitelist. The whitelist cannot be completed: the set
  # of valid R functions is open-ended.
  expect_equal(temp_expr_params("tanh(Ea / tempk)", "tempk"), "Ea")
  expect_equal(temp_expr_params("pmax(Ea / tempk, 0)", "tempk"), "Ea")
  expect_equal(temp_expr_params("ifelse(TRUE, Ea / tempk, 0)", "tempk"), "Ea")
})

test_that("a syntactically invalid expression is reported as such", {
  # Previously an unparseable expression tokenised happily and the failure
  # surfaced much later, from inside the optimiser.
  expect_error(temp_expr_params("exp(Ea / tempk", "tempk"),
               "not valid R")
})

# ── V10-07: a dotted temp_var must not invent a phantom `b` ────────────────

test_that("a temp_var containing a word-boundary b is not read as a parameter", {
  # The mirror of V10-05: the old fit_sb tokeniser was grepl("\\bb\\b", expr),
  # and "." is a non-word character, so a temperature column named `temp.b`
  # matched the humidity-coefficient pattern. build_alpha_ode then added a free
  # parameter `b` to a plain-Arrhenius model with no humidity term and no
  # humidity data, reported it at its hard-coded seed 0.1, and made the Hessian
  # exactly singular so summary()/vcov()/confint() all aborted. all.vars() reads
  # `temp.b` as one symbol, so the parser classifies it as the temperature
  # column, not as `Ea` plus a stray `b`.
  expr <- expr_arrhenius(temp_var = "temp.b")$formula
  expect_equal(temp_expr_params(expr, "temp.b"), "Ea")
})

test_that("fit_sb on a dotted temp_var fits plain Arrhenius with no phantom b", {
  set.seed(3L)
  d <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15),
                   stringsAsFactors = FALSE)
  d$y <- 100 * exp(-0.05 * exp(80000 / 8.314 * (1 / 298.15 - 1 / d$tempk)) * d$time) +
    stats::rnorm(nrow(d), 0, 0.3)
  d$temp.b <- d$tempk

  fit <- suppressWarnings(fit_sb(d, trend = "decreasing", temp_expr = expr_arrhenius(temp_var = "temp.b")))

  # Exactly the five parameters of a humidity-free model -- no invented `b`.
  expect_setequal(names(stats::coef(fit)), c("Ea", "n", "m", "k_0", "y0_0"))
  expect_false("b" %in% names(stats::coef(fit)))
  expect_null(attr(fit, "build")$settings$rh_col)

  # The phantom used to consume a residual df and make inference die.
  expect_equal(df.residual(fit), nrow(d) - 5L)
  V <- vcov(fit)
  expect_equal(dim(V), c(5L, 5L))
  expect_true(all(is.finite(V)))
})

# ── V09-12 / V10-10: a declared covariate is data, not a parameter ─────────

test_that("fit_kinetic accepts a custom temp_expr with a declared covariate", {
  d   <- sim_moist()
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     temp_expr = moist_expr())

  cf <- stats::coef(fit)
  expect_false("moist" %in% names(cf))          # a column, not a parameter
  expect_equal(unname(cf["Ea"]), 80000, tolerance = 1e-4)
  expect_equal(unname(cf["bm"]),   1.5, tolerance = 1e-4)
  expect_equal(unname(cf["k_0"]), 0.02, tolerance = 1e-4)
})

test_that("fit_sb accepts a custom temp_expr with a declared covariate", {
  # The ODE path bound only rh_col into the evaluation environment and seeded
  # only a coefficient literally named `b`, so this failed inside deSolve with
  # "object 'moist' not found".
  fit <- fit_sb(sim_moist(), trend = "decreasing", temp_expr = moist_expr())

  cf <- stats::coef(fit)
  expect_true(all(c("Ea", "bm") %in% names(cf)))
  expect_equal(unname(cf["Ea"]), 80000, tolerance = 1e-3)
  expect_equal(unname(cf["bm"]),   1.5, tolerance = 1e-3)
})

test_that("fit_parallel does not suffix a declared covariate onto its branches", {
  # build_parallel() renames every temperature parameter to <p>B / <p>C. A
  # column read as a parameter therefore became two free branch parameters.
  b <- build_parallel(sim_parallel_moist(), temp_fml = moist_expr()$formula,
                      covar_cols = "moist")
  expect_false(any(c("moist", "moistB", "moistC") %in% b$param_names))
  expect_true(all(c("EaB", "bmB", "EaC", "bmC") %in% b$param_names))
})

# ── R03-03: an UNdeclared covariate must be refused, not absorbed ──────────

test_that("fit_parallel refuses an undeclared covariate instead of fitting it", {
  # The silent failure this replaces: `moist` was read as a parameter, renamed
  # to moistB/moistC, fitted freely, and the data column it named stopped being
  # referenced at all -- yielding a converged fit with a plausible RSS.
  d <- sim_parallel_moist()

  expect_error(
    fit_parallel(d, temp_expr = moist_expr(covar_cols = NULL)),
    "moist"
  )
})

# ── V09-11: names must be checked on the union of all blocks ───────────────

test_that("a temp_expr carrying its own rate parameter is rejected", {
  # expr_rate() and expr_intercept() each check only their own names, so a `k_0`
  # supplied by the temperature expression produced two entries called `k_0` in
  # one start vector.
  d <- sim_moist()
  te <- moist_expr(formula = paste("k_0 * exp(Ea/8.314 * (1/298.15 - 1/tempk))",
                                   "* exp(bm * (moist - 0.6))"))

  expect_error(
    fit_kinetic(d, kinetic_order = 1, trend = "decreasing", temp_expr = te),
    "k_0"
  )
})

# ── V10-02: a covariate may not be named after a parameter ─────────────────

test_that("a temp_var or rh_col named after a parameter is refused at construction", {
  # rlang::eval_tidy()'s data mask outranks the parameter environment, so the
  # column would shadow the parameter and the model would quietly lose it.
  expect_error(expr_arrhenius(temp_var = "Ea"), "Ea")
  expect_error(expr_arrhenius(temp_var = "Ea_log", log_Ea = TRUE), "Ea_log")
  expect_error(expr_arrhenius_humidity(rh_col = "b"), "\\bb\\b")
  expect_error(expr_arrhenius_humidity(temp_var = "Ea"), "Ea")

  # ordinary names are unaffected
  expect_silent(expr_arrhenius(temp_var = "myT"))
  expect_silent(expr_arrhenius_humidity(temp_var = "myT", rh_col = "hum"))
})

# ── V10-04: an omitted list element must not NULL the argument ─────────────

test_that("a temp_expr that omits an element falls back to the argument", {
  d <- sim_moist()

  # expr_arrhenius() returns five elements and expr_arrhenius_humidity() seven,
  # so a hand-built list is easy to under-populate. Overwriting from it
  # unconditionally replaced the caller's value with NULL: omitting trefk or
  # log_Ea fitted silently with NULL in settings, omitting temp_var failed deep
  # inside the optimiser.
  for (drop in c("R", "trefk", "log_Ea", "temp_var")) {
    te <- moist_expr()
    te[[drop]] <- NULL
    fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                       temp_expr = te)
    s <- attr(fit, "build")$settings
    expect_false(is.null(s$trefk), info = drop)
    expect_false(is.null(s$R), info = drop)
    expect_equal(s$temp_var, "tempk", info = drop)
  }
})

test_that("a temp_expr with no formula, or not a list, is rejected clearly", {
  expect_error(resolve_temp_expr(temp_expr = "exp(Ea/tempk)"),
               "must be a list")
  expect_error(resolve_temp_expr(temp_expr = list(R = 8.314)),
               "no `formula` element")
  expect_error(resolve_temp_expr(temp_expr = moist_expr(covar_cols = 1:2)),
               "character vector")
})

test_that("resolve_temp_expr reports covar_cols as the union with rh_col", {
  te <- resolve_temp_expr(expr_arrhenius_humidity(rh_col = "hum"))
  expect_equal(te$covar_cols, "hum")
  expect_equal(resolve_temp_expr(moist_expr())$covar_cols, "moist")
  expect_null(resolve_temp_expr()$covar_cols)
})

# ── predict-time contract for a declared covariate ────────────────────────

test_that("predict requires the declared covariate rather than falling back", {
  # Without covar_cols in required_vars, newdata missing `moist` silently
  # returned predictions for the *fit* data -- 45 values from a 1-row newdata.
  fit <- fit_kinetic(sim_moist(), kinetic_order = 1, trend = "decreasing",
                     temp_expr = moist_expr())

  expect_error(predict(fit, newdata = data.frame(time = 12, tempk = 298.15)),
               "moist")
  p <- predict(fit, newdata = data.frame(time = 12, tempk = 298.15,
                                         moist = 0.75))
  expect_length(p, 1L)
  expect_equal(unname(p),
               100 * exp(-0.02 * exp(1.5 * (0.75 - 0.6)) * 12),
               tolerance = 1e-3)
})

test_that("a non-numeric declared covariate is refused", {
  d <- sim_moist()
  d$moist <- factor(d$moist)
  expect_error(fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                           temp_expr = moist_expr()),
               "moist")
})

# ── R03-08: dead helper removed ───────────────────────────────────────────

test_that("temp_expr_data_cols is gone", {
  # Dead code whose roxygen named build_parallel() as its caller, which made it
  # look like the live classifier while build_parallel() used the tokeniser.
  expect_false(exists("temp_expr_data_cols",
                      where = asNamespace("stablr"), inherits = FALSE))
})
