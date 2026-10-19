# Parameter-name construction tests
#
# Parameter names are built by joining a prefix, grouping column name(s) and
# level value(s) with "_".  Two independent failure modes are covered here:
#
#   1. The recipe is not injective, so distinct (column, level) pairs can
#      collapse onto one name and silently share a coefficient. The collapse is
#      resolved by suffixing the duplicates with ".1", ".2", ... rather than by
#      refusing the design: "_" is legal in both column names and sanitized
#      level values, so a perfectly ordinary crossed design can trip the recipe
#      and used to be rejected outright (V09-16). Note that sanitize() strips
#      ".", which is what makes it a safe disambiguator.
#   2. Rewriting a *built* expression string to re-prefix parameters cannot tell
#      a parameter stem from a column name or a quoted level value.

# ── helpers ────────────────────────────────────────────────────────────────
sim_kinetic <- function(cells, y0, kref = 0.02, Ea = 80000,
                        times = c(0, 1, 2, 3, 6, 9, 12),
                        temps = c(298.15, 313.15)) {
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    g <- expand.grid(time = times, tempk = temps)
    k <- kref * exp(Ea / 8.314 * (1 / 298.15 - 1 / g$tempk))
    cbind(cells[i, , drop = FALSE], row.names = NULL,
          data.frame(time = g$time, tempk = g$tempk,
                     y = y0[i] * exp(-k * g$time)))
  }))
}

sim_parallel <- function(col, lvls, kref,
                         times = c(0, 1, 3, 6, 9, 12),
                         temps = c(298.15, 313.15)) {
  do.call(rbind, lapply(seq_along(lvls), function(i) {
    g  <- expand.grid(time = times, tempk = temps)
    kB <- kref[i] * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    kC <- kref[i] * 0.5 * exp(70000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    kA <- kB + kC
    yA <- 100 * exp(-kA * g$time)
    out <- data.frame(lv = lvls[i], time = g$time, tempk = g$tempk, yA = yA,
                      yB = (kB / kA) * (100 - yA), yC = (kC / kA) * (100 - yA))
    names(out)[1] <- col
    out
  }))
}

# ── 1. name collisions ─────────────────────────────────────────────────────

test_that("expr_intercept disambiguates colliding (column, level) pairs", {
  df <- expand.grid(f = c("a_b", "z"), f_a = c("b", "a"),
                    stringsAsFactors = FALSE)

  # column 'f' level 'a_b' and column 'f_a' level 'b' both want y0_f_a_b
  ie <- expr_intercept(df, c("f", "f_a"))

  # The point of the fix: one parameter per source, never two sources sharing a
  # coefficient, and never a refusal.
  expect_equal(length(ie$params), length(unique(ie$params)))
  expect_true(all(c("y0_f_a_b.1", "y0_f_a_b.2") %in% ie$params))

  # Every parameter must actually appear in the expression that will be
  # evaluated -- a renamed parameter left out of the expression would be seeded
  # and reported but have no effect on the fit.
  for (p in ie$params) expect_true(grepl(p, ie$expr, fixed = TRUE))
})

test_that("a collision-prone design fits through the fit_* entry points", {
  cells <- expand.grid(f = c("a_b", "z"), f_a = c("b", "a"),
                       stringsAsFactors = FALSE)
  d <- sim_kinetic(cells, y0 = c(89, 92, 97, 100))

  # This test is about parameter NAMING with multi-column intercept_cols.
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     intercept_cols = c("f", "f_a"))
  nms <- names(stats::coef(fit))
  expect_equal(length(nms), length(unique(nms)))
  expect_true(all(c("y0_f_a_b.1", "y0_f_a_b.2") %in% nms))
})

test_that("expr_rate disambiguates cells that collide via the '_' collapse", {
  df <- expand.grid(f = c("a_b", "a"), f_a = c("b", "b_b"),
                    stringsAsFactors = FALSE)

  # (f='a_b', f_a='b') and (f='a', f_a='b_b') both want k_f_f_a_a_b_b
  re <- expr_rate(df, c("f", "f_a"))

  expect_equal(length(re$params), length(unique(re$params)))
  expect_true(all(c("k_f_f_a_a_b_b.1", "k_f_f_a_a_b_b.2") %in% re$params))

  # Each suffixed parameter must be guarded by its own indicator, so the two
  # cells that collided are still told apart at evaluation time. Matched with
  # `[^+]*` standing in for the column reference, which is not what this test is
  # about. `[^+]*` rather than `.*` keeps each match inside a single additive
  # term.
  expect_match(re$expr,
               "== 'a_b'\\) \\* \\([^+]*== 'b'\\)\\) \\* k_f_f_a_a_b_b\\.1")
  expect_match(re$expr,
               "== 'a'\\) \\* \\([^+]*== 'b_b'\\)\\) \\* k_f_f_a_a_b_b\\.2")
})

test_that("no false positives: distinct names are accepted and fit correctly", {
  cells <- expand.grid(batch = c("A", "B"), orient = c("inv", "upr"),
                       stringsAsFactors = FALSE)
  y0 <- c(89, 92, 97, 100)
  d  <- sim_kinetic(cells, y0 = y0)

  expect_silent(expr_intercept(d, c("batch", "orient")))
  expect_silent(expr_rate(d, c("batch", "orient")))

  # Two intercept_cols: reference coding, which is what this test is about.
  m <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                   intercept_cols = c("batch", "orient"))
  expect_length(unique(names(coef(m))), length(names(coef(m))))

  # reference coding must reproduce every per-cell intercept
  nd <- cbind(cells, time = 0, tempk = 298.15)
  expect_equal(unname(predict(m, newdata = nd)), y0, tolerance = 1e-4)
})

test_that("single-factor cell-means coding is unaffected", {
  df <- data.frame(batch = c("A", "B"), stringsAsFactors = FALSE)
  ie <- expr_intercept(df, "batch")
  expect_equal(ie$params, c("y0_batch_A", "y0_batch_B"))
})

# ── 2. branch prefixes must not touch columns or level values ──────────────

test_that("expr_rate honours prefix without mangling columns or levels", {
  df <- data.frame(k_lot = c("x", "y"), stringsAsFactors = FALSE)
  rp <- expr_rate(df, "k_lot", prefix = "kB")

  expect_equal(rp$params, c("kB_k_lot_x", "kB_k_lot_y"))
  # The column reference in the indicator must survive verbatim -- the prefix
  # substitution must not rewrite `k_lot` into `kB_lot`. Asserted on the symbols
  # of the parsed expression rather than on its text: a substring test cannot
  # tell the column apart from the parameter `kB_k_lot_x`, which contains
  # "k_lot" too.
  syms <- all.vars(parse(text = rp$expr))
  expect_true("k_lot" %in% syms)
  expect_false("kB_lot" %in% syms)

  # every parameter named in params must actually appear in expr
  for (p in rp$params) expect_match(rp$expr, p, fixed = TRUE)
})

test_that("expr_rate prefix leaves quoted level values verbatim", {
  df <- data.frame(grp = c("k_1", "lotk_2"), stringsAsFactors = FALSE)
  rp <- expr_rate(df, "grp", prefix = "kB")

  expect_match(rp$expr, "'k_1'", fixed = TRUE)
  expect_match(rp$expr, "'lotk_2'", fixed = TRUE)
  expect_equal(rp$params, c("kB_grp_k_1", "kB_grp_lotk_2"))
})

test_that("expr_rate default prefix is unchanged", {
  df <- data.frame(pkg = c("x", "y"), stringsAsFactors = FALSE)
  expect_equal(expr_rate(df, "pkg")$params, c("k_pkg_x", "k_pkg_y"))
  expect_equal(expr_rate(df, character())$params, "k_0")
  expect_equal(expr_rate(df, NULL)$params, "k_0")
  expect_equal(expr_rate(df, character(), prefix = "kB")$params, "kB_0")
})

test_that("build_parallel expression and params stay in sync", {
  # rate_cols and levels chosen to start with the "k_" the old gsub matched
  for (cs in list(list("k_lot", c("x", "y")),
                  list("grp",   c("k_1", "k_2")),
                  list("pkg",   c("x", "y")))) {
    d <- sim_parallel(cs[[1]], cs[[2]], c(0.03, 0.01))
    b <- build_parallel(d, rate_cols = cs[[1]],
                        temp_fml = expr_arrhenius()$formula)

    rp_B <- get("rp_B", envir = environment(b$fn))
    rp_C <- get("rp_C", envir = environment(b$fn))

    # No parameter may be referenced by the expression unless it is declared.
    # Symbols are taken from the parsed expression, not by splitting the text on
    # non-alphanumerics: a text split reads the contents of string literals (the
    # level values) as symbols. That is the same misclassification the tokeniser
    # fixes removed from the package itself (V10-05); all.vars() distinguishes
    # symbols from calls and literals.
    for (rp in list(rp_B, rp_C)) {
      toks  <- unique(all.vars(parse(text = rp$expr)))
      stray <- setdiff(toks, c(rp$params, cs[[1]], cs[[2]]))
      expect_equal(stray, character(0))
    }
    expect_true(all(rp_B$params %in% b$param_names))
    expect_true(all(rp_C$params %in% b$param_names))
  }
})

test_that("fit_parallel recovers truth for k_-prefixed columns and levels", {
  for (cs in list(list("grp",   c("L1", "L2")),
                  list("grp",   c("k_1", "k_2")),
                  list("k_lot", c("x", "y")),
                  list("pkg",   c("x", "y")))) {
    d <- sim_parallel(cs[[1]], cs[[2]], c(0.03, 0.01))
    m <- fit_parallel(d, kinetic_order = 1,
                      rate_cols = cs[[1]])
    cc <- coef(m)

    expect_equal(unname(cc[["EaB"]]), 80000, tolerance = 1e-3)
    expect_equal(unname(cc[["EaC"]]), 70000, tolerance = 1e-3)
    expect_equal(unname(cc[grep("^kB_", names(cc))]), c(0.03, 0.01),
                 tolerance = 1e-3)
    expect_equal(unname(cc[grep("^kC_", names(cc))]), c(0.015, 0.005),
                 tolerance = 1e-3)
    expect_lt(max(abs(fitted(m)[seq_len(nrow(d))] - d$yA)), 1e-6)
  }
})
