# expr_rate() cell enumeration.
#
# expr_rate() used to derive each cell by splitting interaction() labels on ".",
# assuming "." only ever appeared as interaction()'s separator.  A level value
# containing a "." produced more parts than there were rate_cols, so with one
# column the name was recycled into a product that is FALSE on every row --
# ((batch == 'A') * (batch == 1)) -- leaving that group's rate identically zero
# and its parameter unidentifiable; with two or more columns it errored on
# recycling.  Cells are now enumerated as rows of the data.
#
# The fit_*() entry points were shielded because sanitize_group_values() strips
# "." before the builders run, so these tests exercise expr_rate() directly --
# the exported, documented composition path where the bug was reachable.

# Evaluate a built rate expression against data with given parameter values.
eval_rate <- function(r, d, values) {
  p <- stats::setNames(values, r$params)
  rlang::eval_tidy(rlang::parse_expr(r$expr), data = d,
                   env = list2env(as.list(p), parent = baseenv()))
}

# ── indicator correctness ──────────────────────────────────────────────────

test_that("a dotted level yields an indicator that matches its own rows", {
  d <- data.frame(batch = c("A.1", "B"), stringsAsFactors = FALSE)
  r <- expr_rate(d, "batch")

  expect_equal(r$params, c("k_batch_A1", "k_batch_B"))
  # The level must be compared whole, not split into a bogus product. The
  # comparison target is the level as it appears in the data.
  expect_match(r$expr, "== 'A.1'", fixed = TRUE)
  expect_match(r$expr, "== 'B'", fixed = TRUE)
  # the split product compared the column against the fragment after the dot
  expect_false(grepl("== 1)", r$expr, fixed = TRUE))

  # with every parameter at 7, every row's multiplier must be 7 -- not 0
  expect_equal(eval_rate(r, d, rep(7, 2)), c(7, 7))
})

test_that("dotted levels work with two or more rate_cols", {
  d <- expand.grid(batch = c("A.1", "B"), orient = c("inv", "upr"),
                   stringsAsFactors = FALSE)
  # previously: "Can't recycle `.x` (size 2) to match `.y` (size 3)"
  expect_error(expr_rate(d, c("batch", "orient")), NA)

  r <- expr_rate(d, c("batch", "orient"))
  expect_equal(eval_rate(r, d, rep(7, length(r$params))), rep(7, nrow(d)))
})

test_that("multiple dots in a level are handled", {
  d <- data.frame(b = c("x.y.z", "p"), stringsAsFactors = FALSE)
  r <- expr_rate(d, "b")
  # again the level as in the data, compared whole rather than split on the dots
  expect_match(r$expr, "== 'x.y.z'", fixed = TRUE)
  expect_match(r$expr, "== 'p'", fixed = TRUE)
  expect_equal(eval_rate(r, d, rep(7, 2)), c(7, 7))
})

test_that("each indicator selects exactly its own rows", {
  cases <- list(
    list(data.frame(batch = c("A.1", "B"), stringsAsFactors = FALSE), "batch"),
    list(expand.grid(batch = c("A.1", "B"), orient = c("inv", "upr"),
                     stringsAsFactors = FALSE), c("batch", "orient")),
    list(data.frame(g = c("p", "q", "p"), stringsAsFactors = FALSE), "g"),
    list(data.frame(n = c(10, 2, 2)), "n")
  )

  for (cs in cases) {
    d <- cs[[1]]; cols <- cs[[2]]
    r <- expr_rate(d, cols)
    label <- paste(cols, collapse = ",")

    # every parameter must be reachable: setting it alone to 1 selects >=1 row
    for (i in seq_along(r$params)) {
      v <- rep(0, length(r$params)); v[i] <- 1
      expect_gt(sum(eval_rate(r, d, v)), 0,
                label = paste0(label, " / ", r$params[i]))
    }
    # and the indicators must partition the rows: exactly one matches each row
    expect_equal(eval_rate(r, d, rep(1, length(r$params))), rep(1, nrow(d)),
                 info = label)
  }
})

# ── cell set and ordering preserved ────────────────────────────────────────

test_that("observed_cells orders cells by parameter-name label", {
  # Cells are character labels ordered by parameter-name label, first column
  # varying fastest, regardless of column type or factor levels.
  cases <- list(
    list(data.frame(f = c("b", "a", "b", "a", "c"),
                    g = c("y", "x", "x", "x", "y"),
                    stringsAsFactors = FALSE), c("f", "g"),
         c("a.x", "b.x", "b.y", "c.y")),
    # numbers sort as text ("10" before "2"), factors ignore declared levels
    list(data.frame(n = c(10, 2, 2),
                    f = factor(c("q", "p", "q"), levels = c("q", "p"))),
         c("n", "f"), c("2.p", "10.q", "2.q")),
    list(data.frame(n = c(10, 2, 33, 2)), "n", c("10", "2", "33")),
    list(data.frame(f = factor(c("z", "y"), levels = c("z", "y"))), "f",
         c("y", "z")),
    list(expand.grid(a = c("A", "B"), b = c("1", "2"),
                     stringsAsFactors = FALSE), c("a", "b"),
         c("A.1", "B.1", "A.2", "B.2")),
    list(data.frame(l = c(TRUE, FALSE, TRUE)), "l", c("FALSE", "TRUE"))
  )

  for (cs in cases) {
    d <- cs[[1]]; cols <- cs[[2]]
    ce  <- observed_cells(d, cols)
    got <- do.call(paste, c(lapply(cols, function(v) ce[[v]]), sep = "."))
    expect_equal(got, cs[[3]], info = paste(cols, collapse = ","))
    expect_true(all(vapply(ce, is.character, logical(1))))
  }
})

# ── parameter order and grouping of numeric and factor levels ──────────────

sim_groups <- function(groups, k) {
  do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
    g  <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15))
    kk <- k[i] * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    cbind(groups[rep(i, nrow(g)), , drop = FALSE], g,
          y = 100 * exp(-kk * g$time), row.names = NULL)
  }))
}

test_that("numeric rate levels keep their parameter order", {
  d <- sim_groups(data.frame(n = c(2, 10)), c(0.01, 0.03))
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing", rate_cols = "n")
  expect_equal(grep("^k_", names(coef(fit)), value = TRUE), c("k_n_10", "k_n_2"))
  expect_equal(unname(coef(fit)["k_n_2"]), 0.01, tolerance = 1e-4)
})

test_that("suffixed rate names map to the same cells whatever the factor levels", {
  # Both cells are named k_f_g_a_b_b; the .1/.2 suffixes follow label order, not
  # factor level order.
  d <- sim_groups(data.frame(f = c("a_b", "a"), g = c("b", "b_b")), c(0.01, 0.03))
  d$g <- factor(d$g, levels = c("b_b", "b"))
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))
  expect_equal(unname(coef(fit)["k_f_g_a_b_b.1"]), 0.01, tolerance = 1e-4)  # (a_b, b)
  expect_equal(unname(coef(fit)["k_f_g_a_b_b.2"]), 0.03, tolerance = 1e-4)  # (a, b_b)
})

test_that("numeric values that print alike are one group", {
  # 0.3 and 0.1 * 3 differ but share the label "0.3", so they form one group.
  d <- sim_groups(data.frame(x = c(0.3, 0.1 * 3)), c(0.02, 0.02))
  expect_length(unique(d$x), 2L)
  fit <- fit_kinetic(d, kinetic_order = 1, trend = "decreasing",
                     intercept_cols = "x", rate_cols = "x")
  expect_equal(sort(names(coef(fit))), sort(c("Ea", "y0_x_03", "k_x_03")))
})

test_that("only observed cells get a parameter (incomplete designs)", {
  d <- data.frame(batch  = c("A", "B"),
                  orient = c("inv", "upr"), stringsAsFactors = FALSE)
  r <- expr_rate(d, c("batch", "orient"))
  # A/upr and B/inv were never observed
  expect_equal(r$params, c("k_batch_orient_A_inv", "k_batch_orient_B_upr"))
})

test_that("parameter names are unchanged for clean data", {
  d <- expand.grid(form = c("F1", "F2"), pkg = c("x", "y"),
                   stringsAsFactors = FALSE)
  expect_equal(expr_rate(d, c("form", "pkg"))$params,
               c("k_form_pkg_F1_x", "k_form_pkg_F2_x",
                 "k_form_pkg_F1_y", "k_form_pkg_F2_y"))

  p <- data.frame(pkg = c("x", "y"), stringsAsFactors = FALSE)
  expect_equal(expr_rate(p, "pkg")$params, c("k_pkg_x", "k_pkg_y"))
  expect_equal(expr_rate(p, "pkg", prefix = "kB")$params,
               c("kB_pkg_x", "kB_pkg_y"))
  expect_equal(expr_rate(p, character())$params, "k_0")
})

test_that("colliding cell names are disambiguated, not rejected", {
  # The "_"-joined recipe is not injective, but the design itself is perfectly
  # ordinary, so the duplicates are suffixed rather than refused -- see
  # test-expr-param-names.R.
  d <- expand.grid(f = c("a_b", "a"), f_a = c("b", "b_b"),
                   stringsAsFactors = FALSE)
  p <- expr_rate(d, c("f", "f_a"))$params
  expect_equal(length(p), 4L)
  expect_equal(length(p), length(unique(p)))
})

# ── end to end ─────────────────────────────────────────────────────────────

test_that("a manually composed model with dotted levels recovers truth", {
  kref <- c(A.1 = 0.03, B = 0.01)
  g <- do.call(rbind, lapply(names(kref), function(b) {
    gg <- expand.grid(time = c(0, 1, 3, 6, 9, 12), tempk = c(298.15, 313.15))
    k  <- kref[[b]] * exp(80000 / 8.314 * (1 / 298.15 - 1 / gg$tempk))
    data.frame(batch = b, time = gg$time, tempk = gg$tempk,
               y = 100 * exp(-k * gg$time))
  }))

  r   <- expr_rate(g, "batch")
  fml <- stats::as.formula(paste(
    "y ~ y0 * exp(-", r$expr, "* exp(Ea/8.314*(1/298.15 - 1/tempk)) * time)"))
  st  <- c(stats::setNames(rep(0.02, length(r$params)), r$params),
           y0 = 100, Ea = 80000)

  m  <- gslnls::gsl_nls(fml, data = g, start = st)
  cc <- coef(m)
  expect_equal(unname(cc[["k_batch_A1"]]),  0.03,  tolerance = 1e-4)
  expect_equal(unname(cc[["k_batch_B"]]),   0.01,  tolerance = 1e-4)
  expect_equal(unname(cc[["y0"]]),          100,   tolerance = 1e-4)
  expect_equal(unname(cc[["Ea"]]),          80000, tolerance = 1e-4)
})

test_that("fit_kinetic remains correct (sanitization strips the dot)", {
  g <- do.call(rbind, lapply(c("A.1", "B"), function(b) {
    gg <- expand.grid(time = c(0, 3, 6, 12), tempk = c(298.15, 313.15))
    k  <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / gg$tempk))
    data.frame(batch = b, time = gg$time, tempk = gg$tempk,
               y = 100 * exp(-k * gg$time))
  }))
  m <- fit_kinetic(g, kinetic_order = 1, trend = "decreasing",
                   rate_cols = "batch")
  expect_equal(unname(coef(m)[["Ea"]]), 80000, tolerance = 1e-4)
})
