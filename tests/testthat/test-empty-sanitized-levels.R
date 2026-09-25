test_that("a level that sanitizes to the empty string is rejected", {
  # "##" contains no alphanumerics, so sanitize() leaves nothing behind.  There
  # is only one such level, so it is not a duplicate of anything and the
  # collision check cannot see it -- it needs its own rejection.
  df <- data.frame(batch = c("##", "B"), stringsAsFactors = FALSE)

  testthat::expect_error(
    sanitize_group_values(df, "batch"),
    "no letters, digits or underscores"
  )
  testthat::expect_error(
    sanitize_group_values(df, "batch"),
    "'##'"
  )
})

test_that("every punctuation-only level is named in the error", {
  df <- data.frame(batch = c("##", "(1)", "--", "B"), stringsAsFactors = FALSE)

  err <- tryCatch(sanitize_group_values(df, "batch"),
                  error = function(e) conditionMessage(e))

  # "##" and "--" both empty out; "(1)" survives as "1" and must not be listed.
  testthat::expect_match(err, "'##'")
  testthat::expect_match(err, "'--'")
  testthat::expect_false(grepl("\\(1\\)", err))
})

test_that("the empty-level check fires before the duplicate-pooling check", {
  # Two punctuation-only levels are *also* a sanitize collision.  Either message
  # is defensible, but the empty-level one names the actual defect, so it wins.
  df <- data.frame(batch = c("##", "--", "B"), stringsAsFactors = FALSE)

  testthat::expect_error(
    sanitize_group_values(df, "batch"),
    "no letters, digits or underscores"
  )
})

test_that("levels with punctuation but some alphanumerics still pass", {
  df <- data.frame(batch = c("Lot-1", "Lot#2", "(3)"), stringsAsFactors = FALSE)

  testthat::expect_no_error(out <- sanitize_group_values(df, "batch"))
  testthat::expect_equal(out$batch, c("Lot1", "Lot2", "3"))
})

test_that("NA levels do not trigger the empty-level check", {
  # sanitize(NA) is the string "NA", not "", and NAs are excluded from the
  # (orig, clean) pairs anyway.  NA handling is validated elsewhere; here we only
  # assert this check does not claim an NA is punctuation-only.
  df <- data.frame(batch = c(NA, "B"), stringsAsFactors = FALSE)

  err <- tryCatch({sanitize_group_values(df, "batch"); NA_character_},
                  error = function(e) conditionMessage(e))
  if (!is.na(err)) {
    testthat::expect_false(grepl("no letters, digits or underscores", err))
  } else {
    testthat::expect_true(TRUE)
  }
})

test_that("fit_kinetic rejects a punctuation-only group level", {
  g <- expand.grid(time = c(0, 1, 3, 6, 9, 12), tempk = c(298.15, 313.15))
  mk <- function(b, y0) {
    k <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    data.frame(batch = b, time = g$time, tempk = g$tempk,
               y = y0 * exp(-k * g$time), stringsAsFactors = FALSE)
  }
  df <- rbind(mk("##", 100), mk("B", 80))

  testthat::expect_error(
    fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                intercept_cols = "batch"),
    "no letters, digits or underscores"
  )
  testthat::expect_error(
    fit_kinetic(df, kinetic_order = 1, trend = "decreasing",
                rate_cols = "batch"),
    "no letters, digits or underscores"
  )
})

test_that("fit_sb rejects a punctuation-only level instead of aborting on subscripts", {
  # Previously: make_group() produced the level "", so the parameter was named
  # "k_"; get_param() stripped the prefix, leaving a zero-length name, and
  # k_vec[[""]] aborted with the opaque "subscript out of bounds".
  g <- expand.grid(time = c(0, 1, 3, 6, 9, 12), tempk = c(298.15, 313.15))
  mk <- function(b, y0) {
    k <- 0.02 * exp(80000 / 8.314 * (1 / 298.15 - 1 / g$tempk))
    data.frame(batch = b, time = g$time, tempk = g$tempk,
               y = y0 * exp(-k * g$time), stringsAsFactors = FALSE)
  }
  df <- rbind(mk("##", 1.0), mk("B", 0.8))

  err <- tryCatch(fit_sb(df, trend = "decreasing", Yinf = 0, rate_cols = "batch"),
                  error = function(e) conditionMessage(e))
  testthat::expect_match(err, "no letters, digits or underscores")
  testthat::expect_false(grepl("subscript out of bounds", err))
})

test_that("get_param cannot be handed a zero-length label any more", {
  # Direct statement of the mechanism the check now prevents.
  p <- c(k_ = 0.02, k_B = 0.03)
  k_vec <- get_param(p, "k_")

  testthat::expect_true("" %in% names(k_vec))
  testthat::expect_error(k_vec[[""]], "subscript out of bounds")

  # With the check in place, no fitter can construct such a name from data.
  testthat::expect_error(
    sanitize_group_values(data.frame(b = c("##", "B"), stringsAsFactors = FALSE), "b"),
    "no letters, digits or underscores"
  )
})
