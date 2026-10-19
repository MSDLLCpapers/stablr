# Data with only the cells listed in `cells` observed.  A rate parameter exists
# for each observed cell and for no others.
mk_cells <- function(cells) {
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    g  <- expand.grid(time = c(0, 3, 6, 9, 12), tempk = c(298.15, 313.15))
    k  <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    yA <- 100 * exp(-k * g$time)
    data.frame(f = cells$f[i], g = cells$g[i], time = g$time, tempk = g$tempk,
               y = yA, yA = yA, yB = 2 + (100 - yA), yC = 1 + (100 - yA) / 2,
               stringsAsFactors = FALSE)
  }))
}

# Observed: (A,x) and (B,y).  Unobserved: (A,y) and (B,x) -- each built from
# levels that ARE present in their own column, so the per-column check passes.
diag_data <- function() {
  mk_cells(data.frame(f = c("A", "B"), g = c("x", "y"), stringsAsFactors = FALSE))
}

newdata_cell <- function(f, g) {
  data.frame(f = f, g = g, time = c(0, 12, 24, 36), tempk = 298.15,
             stringsAsFactors = FALSE)
}

test_that("predict rejects an unobserved rate_cols combination", {
  fit <- fit_kinetic(diag_data(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  # Both levels are individually known, so the old per-column check passed and
  # the rate summed to exactly 0 -> exp(0) = 1 -> a flat, non-degrading line.
  testthat::expect_error(
    predict(fit, newdata = newdata_cell("A", "y")),
    "never observed when fitting"
  )
  testthat::expect_error(
    predict(fit, newdata = newdata_cell("B", "x")),
    "never observed when fitting"
  )
})

test_that("the error names the offending combination, not just the columns", {
  fit <- fit_kinetic(diag_data(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  err <- tryCatch(predict(fit, newdata = newdata_cell("A", "y")),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "f=")
  testthat::expect_match(err, "g=")
  testthat::expect_match(err, "A")
  testthat::expect_match(err, "y")
})

test_that("interval predictions reject it too, rather than returning zero width", {
  fit <- fit_kinetic(diag_data(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  # These were the worst symptom: fit == lwr == upr, i.e. certainty about a
  # product that supposedly never degrades.
  for (int in c("confidence", "prediction", "tolerance")) {
    testthat::expect_error(
      predict(fit, newdata = newdata_cell("A", "y"), interval = int),
      "never observed when fitting"
    )
  }
})

test_that("observed combinations still predict a real decay curve", {
  fit <- fit_kinetic(diag_data(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  for (cell in list(c("A", "x"), c("B", "y"))) {
    pred <- predict(fit, newdata = newdata_cell(cell[1], cell[2]))
    testthat::expect_length(pred, 4)
    testthat::expect_true(all(diff(pred) < 0))      # actually decreasing
    testthat::expect_false(all(pred == pred[1]))    # not the flat-line bug
  }
})

test_that("a full factorial design is never rejected", {
  df  <- mk_cells(expand.grid(f = c("A", "B"), g = c("x", "y"),
                              stringsAsFactors = FALSE))
  fit <- fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  grid <- expand.grid(f = c("A", "B"), g = c("x", "y"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(grid))) {
    testthat::expect_no_error(
      predict(fit, newdata = newdata_cell(grid$f[i], grid$g[i]))
    )
  }
})

test_that("a single rate column is unaffected: the marginal check still owns it", {
  # With one column, cells and levels coincide, so check_newdata_levels()'s
  # per-column check already covered this and must remain the one that reports.
  df  <- mk_cells(data.frame(f = c("A", "B"), g = c("x", "x"),
                             stringsAsFactors = FALSE))
  fit <- fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                     rate_cols = "f")

  testthat::expect_error(
    predict(fit, newdata = data.frame(f = "Z", time = 0, tempk = 298.15)),
    "unseen levels in column"
  )
  testthat::expect_no_error(
    predict(fit, newdata = newdata_cell("A", "x"))
  )
})

test_that("rate_cells is recorded only when there is a combination to check", {
  df <- mk_cells(expand.grid(f = c("A", "B"), g = c("x", "y"),
                             stringsAsFactors = FALSE))

  b0 <- build_kinetic(df, kinetic_order = 1, trend = "decreasing",
                      temp_fml = expr_arrhenius()$formula)
  testthat::expect_null(b0$settings$rate_cells)

  b1 <- build_kinetic(df, kinetic_order = 1, trend = "decreasing",
                      rate_cols = "f", temp_fml = expr_arrhenius()$formula)
  testthat::expect_null(b1$settings$rate_cells)

  b2 <- build_kinetic(df, kinetic_order = 1, trend = "decreasing",
                      rate_cols = c("f", "g"),
                      temp_fml = expr_arrhenius()$formula)
  testthat::expect_equal(nrow(b2$settings$rate_cells), 4L)
})

test_that("recorded cells match the parameters expr_rate actually built", {
  # If these two ever drift apart, the check would reject fitted cells or admit
  # unfitted ones, so tie them together explicitly.
  df <- diag_data()
  b  <- build_kinetic(df, kinetic_order = 1, trend = "decreasing",
                      rate_cols = c("f", "g"),
                      temp_fml = expr_arrhenius()$formula)

  testthat::expect_equal(nrow(b$settings$rate_cells), 2L)
  testthat::expect_equal(
    nrow(b$settings$rate_cells),
    length(expr_rate(df, c("f", "g"))$params)
  )
})

test_that("fit_singledeg rejects an unobserved combination", {
  fit <- fit_singledeg(diag_data(), kinetic_order = 1, trend = "decreasing",
                       Yinf = 0, rate_cols = c("f", "g"))

  testthat::expect_error(
    predict(fit, newdata = newdata_cell("A", "y")),
    "never observed when fitting"
  )
  testthat::expect_no_error(predict(fit, newdata = newdata_cell("A", "x")))
})

test_that("fit_parallel rejects an unobserved combination", {
  fit <- fit_parallel(diag_data(), kinetic_order = 1,
                      rate_cols = c("f", "g"))

  testthat::expect_error(
    predict(fit, newdata = newdata_cell("A", "y")),
    "never observed when fitting"
  )
  testthat::expect_no_error(predict(fit, newdata = newdata_cell("A", "x")))
})

test_that("predict_tos rejects an unobserved combination in the schedule", {
  sched <- function(f, g) {
    data.frame(t_start = c(0, 3), t_end = c(3, 12), temp = c(25, 40),
               f = f, g = g, stringsAsFactors = FALSE)
  }
  fit <- fit_kinetic(diag_data(), kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  testthat::expect_error(
    predict_tos(fit, sched("A", "y")),
    "never observed when fitting"
  )
  # The message must say "schedule", not "newdata", for this entry point.
  err <- tryCatch(predict_tos(fit, sched("A", "y")),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "^schedule")

  testthat::expect_no_error(predict_tos(fit, sched("A", "x")))
  testthat::expect_no_error(predict_tos(fit, sched("B", "y")))
})

test_that("levels with punctuation are matched as they appear in the data", {
  # Cells are recorded and checked using the levels as they appear in the data.
  df <- mk_cells(data.frame(f = c("A-1", "B-2"), g = c("x.1", "y.2"),
                            stringsAsFactors = FALSE))
  fit <- fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                     rate_cols = c("f", "g"))

  testthat::expect_no_error(predict(fit, newdata = newdata_cell("A-1", "x.1")))
  testthat::expect_error(
    predict(fit, newdata = newdata_cell("A-1", "y.2")),
    "never observed when fitting"
  )
})

test_that("cell keys cannot be confused by values that look like joined keys", {
  # A separator a value could contain would make two distinct cells compare
  # equal.  Grouping values are sanitized to alphanumerics and "_", so an
  # underscore is the adversarial case: ("a_b","c") vs ("a","b_c").
  cells <- data.frame(f = c("a_b", "a"), g = c("c", "b_c"),
                      stringsAsFactors = FALSE)
  keys  <- cell_keys(cells, c("f", "g"))

  testthat::expect_length(unique(keys), 2L)
})
