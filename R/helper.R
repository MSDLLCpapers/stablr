# Tiny utilities                                                           
sanitize <- function(x) gsub("[^[:alnum:]_]", "", x)

# Yinf is the scalar asymptote of the degradation model. 
check_scalar_Yinf <- function(Yinf) {
  if (length(Yinf) != 1L || !is.numeric(Yinf) || !is.finite(Yinf)) {
    stop(
      sprintf(
        "`Yinf` must be a single finite numeric value, not length %d.",
        length(Yinf)
      ),
      call. = FALSE
    )
  }
  invisible(Yinf)
}

# check order
check_kinetic_order <- function(kinetic_order) {
  if (length(kinetic_order) != 1L || !is.numeric(kinetic_order) ||
      is.na(kinetic_order) || kinetic_order < 0 ||
      round(kinetic_order) != kinetic_order) {
    stop("`kinetic_order` must be a single non-negative integer.", call. = FALSE)
  }
  invisible(kinetic_order)
}

format_level <- function(lvl) {
  if (is.numeric(lvl)) return(as.character(lvl))
  if (is.character(lvl) && grepl("^[+-]?[0-9]*\\.?[0-9]+$", lvl)) return(lvl)
  shQuote(as.character(lvl))
}
# Quoted, escaped string literal for a grouping level as it appears in the data.
# Numeric and factor columns compare via as.character(), so decimals and leading
# zeros match.
level_literal <- function(lvl) {
  encodeString(as.character(lvl), quote = "'")
}

indicator_term <- function(var, lvl, par) {
  glue::glue("(({var} == {level_literal(lvl)}) * {par})")
}

# Check that grouping levels give valid, distinct parameter names. Levels are not
# trimmed, so 'Lot-1 ' and 'Lot-1' collide.
check_group_levels <- function(df, cols) {
  if (length(cols) == 0) return(invisible(NULL))
  cols <- unique(cols)
  check_syntactic_cols(cols)
  for (col in cols) {
    orig <- as.character(df[[col]])
    check_sanitize_collisions(col, orig, sanitize(orig))
  }
  invisible(NULL)
}

# Fitted levels per grouping column, as they appear in the data. predict() and
# predict_tos() check new data against these.
group_levels <- function(df, cols) {
  cols <- unique(cols[nzchar(cols)])
  stats::setNames(lapply(cols, function(v) unique(as.character(df[[v]]))), cols)
}

# Used by fit_sb(), which keys its groups on sanitized labels.
sanitize_group_values <- function(df, cols) {
  if (length(cols) == 0) return(df)
  cols <- unique(cols)
  check_group_levels(df, cols)
  for (col in cols) df[[col]] <- sanitize(as.character(df[[col]]))
  df
}

# Grouping column names are spliced into the expression as bare symbols (see
# indicator_term()), and also form parameter names, so a non-syntactic name ("test
# site", "y (%)", "2023 lot") breaks both the parser and gsl_nls. Backticking does
# not help; reject it here with a remedy.
check_syntactic_cols <- function(cols) {
  cols <- cols[nzchar(cols)]
  bad  <- unique(cols[make.names(cols) != cols])
  if (!length(bad)) return(invisible(NULL))

  stop(
    sprintf("Grouping column name(s) %s are not valid R names. ",
            paste(sQuote(bad), collapse = ", ")),
    call. = FALSE
  )
}

# Sanitizing strips punctuation, so distinct levels can collapse onto one label
# ("Lot-1" and "Lot 1" both become "Lot1") and share a parameter, silently
# pooling groups. Refuse to proceed instead.
check_sanitize_collisions <- function(col, orig, clean) {
  keep  <- !is.na(orig)
  pairs <- unique(data.frame(orig = orig[keep], clean = clean[keep]))

  check_empty_sanitized_levels(col, pairs)

  bad   <- unique(pairs$clean[duplicated(pairs$clean)])
  if (!length(bad)) return(invisible(NULL))

  groups <- vapply(bad, function(lbl) {
    sprintf("%s -> '%s'",
            paste(sQuote(pairs$orig[pairs$clean == lbl]), collapse = ", "),
            lbl)
  }, character(1))

  # Suggest trimws() when levels differ only by surrounding whitespace.
  hint <- if (anyDuplicated(trimws(pairs$orig))) {
    "Levels differing only by surrounding whitespace are distinct; remove it with trimws() before fitting."
  } else ""

  stop(
    sprintf("Grouping column '%s' has distinct levels that collide after removing punctuation: %s. %s",
            col, paste(groups, collapse = "; "), hint),
    call. = FALSE
  )
}

# Reject an intercept design whose indicator matrix is rank deficient.
#
# expr_intercept() uses reference coding for 2+ columns (first column keeps every
# level, later columns drop their first). Correct for a crossed design, but in a
# nested design (batches A,B only at site S1; C,D only at S2) site is a function
# of batch, so the matrix is collinear: only per-cell sums are identified and the
# individual coefficients depend on column order. Build the implied indicator
# matrix and compare its rank to the parameter count.
#
# @param df Data frame providing the grouping columns.
# @param specs Data frame of (col, lvl) pairs, one row per emitted parameter, in
#   parameter order.
# @param params Parameter names, aligned with `specs`.
# @param what Noun used in the error message.
# @return Invisibly NULL; errors when the design is rank deficient.
# @keywords internal
check_intercept_rank <- function(df, specs, params, what = "intercept") {
  if (nrow(specs) < 2L) return(invisible(NULL))

  # One indicator column per emitted parameter, exactly as the expression sums
  # them; rows are the distinct grouping cells (repeated rows add no rank).
  cols  <- unique(specs$col)
  cells <- unique(df[, cols, drop = FALSE])
  X <- vapply(seq_len(nrow(specs)), function(i) {
    as.numeric(as.character(cells[[specs$col[i]]]) == as.character(specs$lvl[i]))
  }, numeric(nrow(cells)))
  X <- matrix(X, nrow = nrow(cells))

  rank <- tryCatch(qr(X)$rank, error = function(e) NA_integer_)
  if (is.na(rank) || rank >= ncol(X)) return(invisible(NULL))

  stop(
    sprintf(
      "The %s design is rank deficient: %d parameters but only rank %d.",
      what, ncol(X), rank
    ),
    call. = FALSE
  )
}

# A level made up entirely of punctuation ("##", "--") sanitizes to "", so its
# parameter name is just the stem with a trailing "_" ("k_batch_", "y0_"), which
# cannot be looked up by group label. Reject it at the source.
check_empty_sanitized_levels <- function(col, pairs) {
  bad <- unique(pairs$orig[!nzchar(pairs$clean)])
  if (!length(bad)) return(invisible(NULL))

  stop(
    sprintf("Grouping column '%s' has level(s) %s with no letters, digits or underscores, so they cannot name a parameter.",
            col, paste(sQuote(bad), collapse = ", ")),
    call. = FALSE
  )
}

# Parameter names paste a prefix, column name(s) and level value(s) with "_".
# Because "_" is legal inside both, that recipe is not injective: column "f" level
# "a_b" and column "f_a" level "b" both collapse to "y0_f_a_b", so two distinct
# indicators would share one coefficient. Refuse to build such a model.
#
# `sources` is a character vector parallel to `params`, naming where each name
# came from, for an actionable message.
check_param_name_collisions <- function(params, sources, what = "parameter") {
  params <- as.character(params)
  if (!anyDuplicated(params)) return(invisible(NULL))

  dup <- unique(params[duplicated(params)])
  groups <- vapply(dup, function(nm) {
    sprintf("'%s' <- %s", nm,
            paste(unique(sources[params == nm]), collapse = " and "))
  }, character(1))

  stop(
    sprintf("Ambiguous %s names: %s. ", what, paste(groups, collapse = "; ")),
    "Distinct groups would be pooled into one coefficient; rename the columns ",
    "or levels so they do not collide when joined with '_'.",
    call. = FALSE
  )
}

# Make generated parameter names unique, leaving already-unique names untouched.
#
# The naming recipe is not injective (see check_param_name_collisions): distinct
# columns/levels can collapse to one name. Such designs are legitimate (fit_sb()
# fits them, joining with "||"), so disambiguate rather than error.
# Disambiguating keeps every non-colliding name byte-identical, so no existing
# fit is renamed. The suffix uses "." because sanitize() strips it, so ".1"
# cannot itself collide.
dedupe_param_names <- function(params) {
  params <- as.character(params)
  dups <- unique(params[duplicated(params)])
  for (nm in dups) {
    i <- which(params == nm)
    params[i] <- paste0(nm, ".", seq_along(i))
  }
  params
}

# A covariate column may not be named after a parameter of the same expression.
#
# temp_expr_params() classifies by subtraction, so a covariate named after a
# parameter would delete the parameter instead of clashing (the symbol is read as
# data). The check belongs here, the only place that knows the emitted names. A
# hand-built temp_expr is on its own.
#
# `covars` is a named character vector (names = argument names, for the message);
# `params` are the parameter names the expression will contain.
check_covar_names <- function(covars, params) {
  covars <- covars[!vapply(covars, is.null, logical(1))]
  clash  <- covars[as.character(covars) %in% params]
  if (!length(clash)) return(invisible(NULL))

  stop(
    sprintf("%s cannot be %s: that names a parameter in this expression (%s). ",
            paste(sprintf("`%s`", names(clash)), collapse = " and "),
            paste(sQuote(as.character(clash)), collapse = " and "),
            paste(params, collapse = ", ")),
    call. = FALSE
  )
}

# Check for duplicate parameter names across blocks, not within one.
#
# expr_intercept() and expr_rate() each check their own output but cannot see the
# other or the temperature expression, so a name duplicated across blocks (e.g. a
# custom temp_expr carrying its own `k_0` alongside expr_rate()'s) would reach the
# fit as two coefficients of the same name.
#
# `blocks` is a named list of character vectors; the names appear verbatim in the
# message, so they should read as phrases ("the rate term").
check_param_union <- function(blocks) {
  blocks <- blocks[lengths(blocks) > 0L]
  check_param_name_collisions(
    unlist(blocks, use.names = FALSE),
    rep(names(blocks), lengths(blocks)),
    what = "parameter"
  )
}

# Observed combinations of `cols` as character labels, one row per cell, ordered
# by parameter-name label.
observed_cells <- function(df, cols) {
  labels <- lapply(df[cols], as.character)
  cells  <- unique(data.frame(labels, stringsAsFactors = FALSE, check.names = FALSE))

  # First column varies fastest (order() treats its first argument as most
  # significant, so the keys are reversed).
  keys  <- lapply(cols, function(v) sanitize(cells[[v]]))
  cells <- cells[do.call(order, rev(keys)), , drop = FALSE]
  rownames(cells) <- NULL
  cells
}

# The rate cells a builder created parameters for, recorded in settings$rate_cells
# so predict can check against them. Deliberately the same observed_cells() call
# expr_rate() makes, so the record cannot drift from the fitted parameters. NULL
# when there is nothing to check (fewer than two rate columns).
rate_cells <- function(df, rate_cols) {
  if (is.null(rate_cols) || length(rate_cols) < 2L) return(NULL)
  if (any(!nzchar(rate_cols))) return(NULL)
  observed_cells(df, rate_cols)
}

# Row keys for comparing observed cells against prediction cells. encodeString()
# escapes control characters, so "\r" cannot occur in an encoded value and is a
# safe separator (unlike the lossy interaction()-label round-trip; see
# observed_cells()).
cell_keys <- function(df, cols) {
  do.call(paste, c(lapply(cols, function(v) encodeString(as.character(df[[v]]))),
                   sep = "\r"))
}

# Render cell keys as human-readable labels, e.g. "(batch='A', pkg='x')"
cell_labels <- function(keys, cols) {
  vapply(strsplit(keys, "\r", fixed = TRUE), function(parts) {
    sprintf("(%s)", paste(sprintf("%s=%s", cols, sQuote(parts)),
                          collapse = ", "))
  }, character(1))
}

# Reject prediction rows whose rate_cols combination was never observed.
#
# expr_rate() encodes rates as a sum of per-observed-cell indicators, so an
# unfitted combination matches none, giving exp(-0 * t) == 1: a flat prediction
# with zero-width intervals. check_newdata_levels()'s per-column check cannot catch
# this, because each marginal level is known while the pair is unfitted; it is a
# no-op for a single rate column.
check_observed_cells <- function(newdata, cells, rate_cols, what = "newdata",
                                 time_var = NULL) {
  if (is.null(cells) || length(rate_cols) < 2L) return(invisible(NULL))
  if (!all(rate_cols %in% names(newdata))) return(invisible(NULL))

  # A prediction at time 0 is answered entirely by the intercept (the rate term
  # is multiplied by t), so an unfitted rate cell is harmless there; the sibling
  # guard check_rate_support_newdata() exempts time 0 too. time_var is NULL where
  # there is no single time column (predict_tos() passes a schedule), and then
  # nothing is exempt.
  if (!is.null(time_var) && time_var %in% names(newdata)) {
    tm <- suppressWarnings(as.numeric(newdata[[time_var]]))
    newdata <- newdata[!(is.finite(tm) & tm == 0), , drop = FALSE]
    if (!nrow(newdata)) return(invisible(NULL))
  }

  new_keys <- cell_keys(newdata, rate_cols)
  unseen   <- unique(new_keys[!new_keys %in% cell_keys(cells, rate_cols)])
  if (!length(unseen)) return(invisible(NULL))

  labels <- cell_labels(unseen, rate_cols)

  stop(
    sprintf("%s contains rate_cols combination(s) never observed when fitting: %s. ",
            what, paste(labels, collapse = "; ")),
    "An unfitted combination has no rate, so it would be predicted as a flat ",
    "line with zero-width intervals. Restrict predictions to the fitted combinations.",
    call. = FALSE
  )
}

# Warn when the optimiser stopped without converging.
#
# gsl_nls() reports its exit status in convInfo. A fit that hit the iteration
# limit or stalled is otherwise indistinguishable from a converged one, so its
# last-visited parameter vector would be treated as an estimate. This matters most
# for models with many rate parameters, where crossing grouping columns leaves
# cells poorly determined.
#
# isConv alone is not sufficient: gsl_nls also reports FALSE when the step search
# stalls at a stationary point (a correctly solved but misspecified model), and
# warning there would cry wolf on ordinary lack of fit. So confirm with the scaled
# gradient J'r, which is ~0 near a solution regardless of residual size; only warn
# when the optimiser both stopped early and is still on a slope.
#
# @param mod Fitted gsl_nls object.
# @param n_par Number of parameters in the model.
# @param tol Relative-gradient threshold below which the stopping point is
#   treated as a genuine stationary point.
# @return Invisibly the model; warns when convInfo reports non-convergence.
# @keywords internal
check_convergence <- function(mod, n_par = NA_integer_, tol = 1e-5) {
  ci <- mod[["convInfo"]]
  if (is.null(ci) || is.null(ci$isConv) || isTRUE(ci$isConv)) {
    return(invisible(mod))
  }

  # Stalled at a stationary point: the fit may be poor, but it is the solution,
  # and lack of fit is not this function's business to report.
  rel_grad <- tryCatch({
    g <- mod$m$gradient()
    r <- stats::residuals(mod)
    if (is.null(g) || !length(r)) NA_real_
    else max(abs(crossprod(g, r))) / max(1, sum(r^2))
  }, error = function(e) NA_real_)
  if (is.finite(rel_grad) && rel_grad < tol) return(invisible(mod))

  detail <- if (!is.null(ci$stopMessage) && nzchar(ci$stopMessage)) {
    sprintf(" (%s", ci$stopMessage)
  } else {
    " (unknown reason"
  }
  detail <- if (!is.null(ci$finIter)) {
    sprintf("%s after %d iteration%s)", detail, ci$finIter,
            ifelse(ci$finIter == 1L, "", "s"))
  } else {
    paste0(detail, ")")
  }

  warning(sprintf("The optimiser did not converge%s. ", detail))

  invisible(mod)
}

# Warn when a fitted activation energy comes back negative.
#
# @param coefs Named numeric coefficient vector from the fitted model.
# @return Invisibly NULL; called for its side effect (a warning).
# @keywords internal
check_activation_energy_sign <- function(coefs) {
  if (is.null(coefs) || !length(coefs) || is.null(names(coefs))) {
    return(invisible(NULL))
  }
  nm <- names(coefs)
  # Match plain Ea (including branch-suffixed EaB/EaC) but exclude both spellings
  # of the log form: Ea_log/Ea_logB (kinetic family) and Ealog (fit_sb).
  is_log <- grepl("^Ea_?log", nm)
  is_Ea  <- grepl("^Ea", nm) & !is_log
  neg <- is_Ea & is.finite(coefs) & coefs < 0
  if (!any(neg)) return(invisible(NULL))

  warning(
    "Fitted activation energy is negative (",
    paste(sprintf("%s = %.4g", nm[neg], coefs[neg]), collapse = ", "), "). ",
    "Consider revising the model arguments or reparametrizing the model.",
    call. = FALSE
  )
  invisible(NULL)
}

# Verify a composite reconstruction evaluates every fitted parameter.
#
# predict_tos() does not reuse the fitted model function: each builder in
# predict_temp_deviation.R reassembles the model from parts to step through a
# schedule, and a reconstructed part can come back missing a term:
# .build_composite_fn_kinetic() rebuilt the intercept from the schedule, and
# because expr_intercept() derives reference coding from the levels it is given, a
# sub-schedule carrying one level of a second intercept column made that level its
# own reference and dropped its coefficient. Nothing in the output reveals it: the
# prediction is off by exactly the omitted coefficient and its delta-method
# interval is centred on the same wrong value. Comparing parameter sets catches
# it regardless of which part went missing.
#
# @param exprs Character vector of the expressions the composite evaluates.
# @param param_names Fitted parameter names, from build$param_names.
# @param what Noun naming the caller, used in the message.
# @return Invisibly NULL; errors when a fitted parameter appears in no
#   expression.
# @keywords internal
check_composite_params <- function(exprs, param_names,
                                   what = "predict_tos()") {
  exprs <- as.character(exprs)
  exprs <- exprs[!is.na(exprs) & nzchar(exprs)]
  if (!length(exprs) || !length(param_names)) return(invisible(NULL))

  seen <- unique(unlist(lapply(exprs, function(e) {
    tryCatch(all.vars(parse(text = e)), error = function(err) character())
  })))

  missing <- setdiff(param_names, seen)
  if (!length(missing)) return(invisible(NULL))

  stop(
    sprintf(
      paste0(
        "Fitted parameter(s) %s are missing from the %s expressions, so its ",
        "predictions and intervals would be silently wrong; please report it."
      ),
      paste(sprintf("'%s'", missing), collapse = ", "), what
    ),
    call. = FALSE
  )
}

# Reject non-numeric values in the columns that are used as numbers.
#
# Temperature and time are arithmetic inputs, so a factor or character column
# there is a type error, not a grouping. This is common from read.csv() with
# stringsAsFactors, a group_by()/summarise() table, or a factor-coded condition
# column.
#
# Checked rather than silently coerced: as.numeric() on a factor returns integer
# level codes, so factor(c(25, 40, 50)) would become 1, 2, 3, and coercing via
# as.character() quietly would also rescue genuinely mislabelled data.
#
# @param df Data frame supplied by the user.
# @param cols Columns that are used as numbers. Missing or empty names ignored.
# @param what Noun used in the error message.
# @return Invisibly NULL; errors when any named column is not numeric.
# @keywords internal
check_numeric_cols <- function(df, cols, what = "df") {
  cols <- unique(cols[nzchar(cols) & cols %in% names(df)])
  if (!length(cols)) return(invisible(NULL))

  bad <- cols[!vapply(cols, function(col) is.numeric(df[[col]]), logical(1))]
  if (!length(bad)) return(invisible(NULL))

  items <- vapply(bad, function(col) {
    sprintf("'%s' (%s)", col, class(df[[col]])[1])
  }, character(1))

  stop(
    sprintf(
      "%s column(s) %s must be numeric, check the column type",
      what, paste(items, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Reject negative times in the alpha-space ODE path (fit_sb only).
#
# fit_sb integrates dalpha/dt from an initial condition at time 0, so a negative
# time would re-reference that initial condition and change what y0 means.
# Integrating backward is not an option: deSolve needs a monotonic time vector,
# and more fundamentally alpha is the extent of degradation on [0, 1] with
# alpha(0) ~ 0, so a time before the reference implies alpha < 0, which the model
# cannot represent. Re-referencing the time axis is the caller's decision.
# Scoped to fit_sb: the closed-form fitters evaluate an algebraic expression in t
# where a negative time is well defined and disturbs no other row.
#
# @param df Data frame supplied by the user.
# @param time_var Name of the time column.
# @param what Noun used in the error message.
# @return Invisibly NULL; errors when the time column contains a negative value.
# @keywords internal
check_nonnegative_times <- function(df, time_var, what = "df") {
  if (!length(time_var) || !nzchar(time_var[1]) ||
      !time_var[1] %in% names(df)) {
    return(invisible(NULL))
  }

  tv  <- df[[time_var[1]]]
  if (!is.numeric(tv)) return(invisible(NULL))
  bad <- which(!is.na(tv) & tv < 0)
  if (!length(bad)) return(invisible(NULL))

  stop(
    sprintf(
      "The alpha-space ODE model (fit_sb) requires non-negative times, but %s column '%s' has %d negative value%s.",
      what, time_var[1], length(bad),
      if (length(bad) == 1L) "" else "s"
    ),
    call. = FALSE
  )
}

# Reject an empty data frame up front, so downstream code does not fail with a
# message that never mentions the empty input. Common from subset()/filter() with
# a mistyped condition.
#
# @param df Data frame supplied by the user.
# @param what Noun used in the error message.
# @return Invisibly NULL; errors when df has zero rows.
# @keywords internal
# @noRd
check_nonempty_df <- function(df, what = "df") {
  if (nrow(df) == 0L) {
    stop(what, " has no rows.", call. = FALSE)
  }
  invisible(NULL)
}

# Reject NA in the columns that define the model structure.
#
# An NA in a grouping column or covariate is not an observation of a group but
# unknown group membership, which no fitted parameter describes; left unchecked
# each fitter mishandles it differently (silently dropped rows, phantom levels, or
# estimates dragged toward zero). NA grouping values are plausible from an
# upstream join, so check them at the door.
#
# @param df Data frame supplied by the user.
# @param cols Columns that define the model structure (grouping columns,
#   temperature, time, humidity). Missing or empty names are ignored.
# @param what Noun used in the error message.
# @return Invisibly NULL; errors when any named column contains NA.
# @keywords internal
check_structural_nas <- function(df, cols, what = "df") {
  cols <- unique(cols[nzchar(cols) & cols %in% names(df)])
  if (!length(cols)) return(invisible(NULL))

  counts <- vapply(cols, function(col) sum(is.na(df[[col]])), integer(1))
  bad    <- counts[counts > 0L]
  if (!length(bad)) return(invisible(NULL))

  items <- sprintf("'%s' (%d row%s)", names(bad), bad,
                   ifelse(bad == 1L, "", "s"))

  consequence <- if (identical(what, "newdata")) {
    "An unknown group usually means a failed upstream join. "
  } else {
    "An unknown group matches no fitted parameter. "
  }

  stop(
    sprintf("NA values in %s column(s) that define the model: %s. ",
            what, paste(items, collapse = ", ")),
    consequence,
    "Remove those rows, or label them 'unknown'.",
    call. = FALSE
  )
}

# Identify rate cells that carry too little information to estimate their rate.
#
# expr_rate()/make_group() enumerate whatever cells are observed, without checking
# each has enough distinct times to identify a rate. A cell whose rows all sit at
# time == 0 gives residuals independent of its k (exp(-k*0) == 1), so k is
# unconstrained and its singularity also corrupts the whole covariance matrix.
#
# Two tiers of support failure:
#   "no elapsed time":  every row at time == 0; k carries zero information.
#   "single timepoint": one distinct time; k trades off against the intercept.
#
# Tier two is conditional: one timepoint suffices when k gets an equation it
# need not share with y0: either the cell's intercept is pooled with cells that do
# carry elapsed time (y0 pinned elsewhere), or the cell was observed at two or
# more temperatures at that time (a second equation, since Ea is shared). Tier one
# stays unconditional: no pooling or temperature spread recovers information that
# is not there.
#
# @param df Model data frame.
# @param rate_cols Rate grouping columns (may be empty for a single global rate).
# @param time_var Name of the time column.
# @param intercept_cols Intercept grouping columns, used to tell a pooled
#   intercept from a per-cell one. Empty means a single pooled intercept.
# @param temp_var Name of the temperature column; NULL skips the
#   two-temperature test (treating the cell as single-temperature).
# @return Data frame with columns `key` (the raw cell key, for matching),
#   `label` (for display), `reason`, and `obs_time` (the one time the cell was
#   observed at, for tier two; NA otherwise), one row per under-supported cell.
# @keywords internal
unsupported_rate_cells <- function(df, rate_cols, time_var = "time",
                                   intercept_cols = character(),
                                   temp_var = NULL) {
  empty <- data.frame(key = character(), label = character(),
                      reason = character(), obs_time = numeric(),
                      stringsAsFactors = FALSE)
  if (!time_var %in% names(df)) return(empty)

  rate_cols <- rate_cols[nzchar(rate_cols)]
  if (length(rate_cols) && !all(rate_cols %in% names(df))) return(empty)

  tm <- suppressWarnings(as.numeric(df[[time_var]]))

  # With no rate columns the whole data set is one cell.
  if (!length(rate_cols)) {
    keys   <- rep("0", nrow(df))
    labels <- c("0" = "the single global rate")
  } else {
    keys <- cell_keys(df, rate_cols)
    uk   <- unique(keys)
    labels <- cell_labels(uk, rate_cols)
    names(labels) <- uk
  }

  # Intercept cell of every row; one pooled intercept when there are no
  # intercept_cols (or they are absent from df, as in a legacy call).
  intercept_cols <- intercept_cols[nzchar(intercept_cols)]
  int_keys <- if (length(intercept_cols) &&
                  all(intercept_cols %in% names(df))) {
    cell_keys(df, intercept_cols)
  } else {
    rep("0", nrow(df))
  }

  temps <- if (!is.null(temp_var) && temp_var %in% names(df)) {
    suppressWarnings(as.numeric(df[[temp_var]]))
  } else {
    rep(NA_real_, nrow(df))
  }

  out <- lapply(unique(keys), function(k) {
    in_cell <- keys == k
    t_cell  <- tm[in_cell]
    t_cell  <- t_cell[is.finite(t_cell)]
    n_dist  <- length(unique(t_cell))

    if (!length(t_cell) || all(t_cell == 0)) {
      return(data.frame(key = k, label = unname(labels[k]),
                        reason = "no rows with elapsed time > 0",
                        obs_time = NA_real_, stringsAsFactors = FALSE))
    }
    if (n_dist >= 2L) return(NULL)

    # One distinct timepoint: identified if the intercept is pinned elsewhere by
    # rows that carry elapsed time, or if this cell spans two temperatures.
    shared_intercept <- any(
      int_keys[!in_cell & is.finite(tm) & tm > 0] %in% int_keys[in_cell]
    )
    n_temp <- length(unique(temps[in_cell & is.finite(temps)]))
    if (shared_intercept || n_temp >= 2L) return(NULL)

    data.frame(key = k, label = unname(labels[k]),
               reason = "only one distinct timepoint",
               obs_time = unique(t_cell), stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, out)
  if (is.null(out)) empty else out
}

# Warn about rate cells with too little support, naming each one.
#
# A warning rather than an error: incomplete sampling is the norm in stability
# data and the rest of the model is still worth fitting, but the affected rate is
# not an estimate, so say which cell and why.
#
# @inheritParams unsupported_rate_cells
# @return Invisibly the under-supported cells, for recording in settings.
# @keywords internal
check_rate_support <- function(df, rate_cols, time_var = "time",
                               intercept_cols = character(),
                               temp_var = NULL) {
  bad <- unsupported_rate_cells(df, rate_cols, time_var, intercept_cols,
                                temp_var)
  if (!nrow(bad)) return(invisible(bad))

  items <- sprintf("%s: %s", bad$label, bad$reason)
  warning(
    sprintf(
      "Rate parameter(s) fitted from insufficient data: %s.",
      paste(items, collapse = "; ")
    ),
    call. = FALSE
  )

  invisible(bad)
}

# Warn when rate cells are confined to a single temperature.
#
# check_rate_support() covers too few timepoints; a cell can have plenty of
# timepoints and still be ill-determined if they all sit at one temperature. Ea
# is shared while k_ref is per-cell, so a single-temperature cell cannot separate
# them: any Ea is absorbed by rescaling that cell's k_ref, the cell fits its own
# data at an arbitrary ridge point, and because Ea is shared the distortion leaks
# into every other cell. This is the usual cost of crossing a second column into
# rate_cols, which shatters cells that spanned the full temperature range into
# fragments that do not.
#
# @param df Model data frame.
# @param rate_cols Rate grouping columns.
# @param temp_var Name of the temperature column.
# @return Invisibly the affected cell labels.
# @keywords internal
check_rate_temp_support <- function(df, rate_cols, temp_var = "tempk") {
  rate_cols <- rate_cols[nzchar(rate_cols)]
  # A single column is enough to carry the problem: if it is the temperature arm
  # itself (or aligned with it) every cell sits at one temperature and Ea is
  # unidentified. Only a global rate has nothing to confound, so that is the only
  # case skipped here.
  if (length(rate_cols) < 1L) return(invisible(character()))
  if (!temp_var %in% names(df) || !all(rate_cols %in% names(df))) {
    return(invisible(character()))
  }
  # Only meaningful when the data as a whole spans temperatures; an isothermal
  # study has no Ea information anywhere and is a different problem.
  if (length(unique(df[[temp_var]][!is.na(df[[temp_var]])])) < 2L) {
    return(invisible(character()))
  }

  keys <- cell_keys(df, rate_cols)
  uk   <- unique(keys)
  bad  <- uk[vapply(uk, function(k) {
    tk <- df[[temp_var]][keys == k]
    length(unique(tk[!is.na(tk)])) < 2L
  }, logical(1))]
  if (!length(bad)) return(invisible(character()))

  labels <- vapply(strsplit(bad, "\r", fixed = TRUE), function(parts) {
    sprintf("(%s)", paste(sprintf("%s=%s", rate_cols, sQuote(parts)),
                          collapse = ", "))
  }, character(1))

  remedy <- if (length(rate_cols) > 1L) {
    sprintf(paste0(
      "Crossing %d columns made %d thin cells; use one column and check the ",
      "residual sum of squares first."),
      length(rate_cols), length(uk))
  } else {
    sprintf(paste0(
      "`rate_cols = \"%s\"` is aligned with `%s`, so Ea is unidentified. Drop ",
      "it from `rate_cols`, or group by a temperature-spanning column."),
      rate_cols, temp_var)
  }

  warning(
    sprintf(
      paste0(
        "%d of %d rate cell(s) observed at only one `%s` value: %s. Because ",
        "the activation energy is shared, these rates turn unphysical. %s"
      ),
      length(bad), length(uk), temp_var, paste(labels, collapse = "; "), remedy
    ),
    call. = FALSE
  )

  invisible(labels)
}

# Warn when the whole design spans a single temperature: the Arrhenius factor is
# then constant, so only k * A(Ea) is identifiable and Ea is start-driven with a
# huge SE. Complements check_rate_temp_support(), which handles the per-cell
# version once the design does span temperatures.
# @keywords internal
#
# @param temp_params Free parameters of the temperature expression, when known. A
#   hand-built temp_expr need not have an Ea at all, and warning about a missing
#   Ea would name a cause that does not exist. NULL means "unknown", keeping the
#   Arrhenius wording.
check_temp_identifiable <- function(df, temp_var = "tempk", temp_params = NULL) {
  temps <- unique(df[[temp_var]][!is.na(df[[temp_var]])])
  if (length(temps) >= 2L) return(invisible(NULL))
  # Nothing in the temperature term is estimated, so nothing about it is
  # unidentifiable: a single temperature is simply all this model needs.
  if (!is.null(temp_params) && !length(temp_params)) return(invisible(NULL))

  pars <- if (length(temp_params)) paste(temp_params, collapse = ", ") else "Ea"
  warning(sprintf(
    "Only one `%s` value in the data: temperature parameter(s) %s are not identifiable; supply >= 2 temperatures to estimate them.",
    temp_var, pars), call. = FALSE)
  invisible(NULL)
}


# The closed-form fitters pick the model RHS purely from the `trend` string;
# nothing compares it to the data. Getting `trend` backwards drives the fit to a
# negative rate constant and an uninterpretable Ea. Compare the sign of the
# least-squares time slope against the declared trend and warn on a clear
# reversal. Correlation (with a small floor) keeps assay noise on a genuine trend
# from tripping it and exempts near-flat outcomes, which carry no direction.
#
# The verdict is unanimous, not a median over the temperature arms: a median
# flips on a single odd arm, and one arm disagreeing is ordinary when a
# slow-degrading arm carries little signal against assay noise. Only a mismatch
# at *every* temperature is strong enough evidence to raise.
# @param df Data frame with the outcome and time columns.
# @param outcome_var,time_var Column names for the outcome and time.
# @param trend Declared trend, "increasing" or "decreasing".
# @param temp_var Temperature column; the correlation is computed within each of
#   its values rather than pooled.
# @param advice Optional string appended as the remedy sentence. Detection is
#   generic but the fix differs by model: callers with a `trend` argument want
#   the default, while a model that hardwires its direction (fit_parallel)
#   supplies its own.
# @return Invisibly NULL; warns when every temperature arm clearly contradicts
#   `trend`.
# @keywords internal
check_trend_direction <- function(df, outcome_var, time_var, trend,
                                  temp_var = NULL, advice = NULL) {
  y <- df[[outcome_var]]; t <- df[[time_var]]
  ok <- which(is.finite(y) & is.finite(t))
  grp <- if (!is.null(temp_var) && temp_var %in% names(df)) {
    as.character(df[[temp_var]])
  } else {
    rep("all", length(y))
  }

  # One correlation per temperature; arms with too little data to say anything
  # drop out as NA.
  rs <- vapply(split(ok, grp[ok]), function(i) {
    if (length(i) < 3L || stats::sd(t[i]) == 0 || stats::sd(y[i]) == 0) {
      return(NA_real_)
    }
    suppressWarnings(stats::cor(t[i], y[i]))
  }, numeric(1))
  rs <- rs[is.finite(rs)]
  if (!length(rs)) return(invisible(NULL))

  # An arm that is near-flat (|r| < 0.1) carries no direction, so it is not
  # evidence against `trend` and blocks the warning along with any arm that
  # actually agrees.
  contradicts <- if (trend == "decreasing") rs > 0.1 else rs < -0.1
  if (!all(contradicts)) return(invisible(NULL))

  if (is.null(advice)) advice <- "Consider revising `trend`."
  warning(sprintf(
    "Declared trend = \"%s\" but `%s` trends %s at every temperature. %s",
    trend, outcome_var, if (trend == "decreasing") "upward" else "downward",
    advice), call. = FALSE)
  invisible(NULL)
}


#
# Predicting an under-supported cell at time == 0 is harmless (the rate does not
# enter), but any positive time produces a forecast driven entirely by an
# unidentified parameter: a flat, zero-degradation line presented as fitted.
# Mirrors check_observed_cells(), which errors on cells never observed at all.
#
# @param newdata Prediction data, already sanitized.
# @param unsupported Data frame from unsupported_rate_cells(), recorded at fit
#   time in settings$unsupported_rate_cells.
# @param rate_cols Rate grouping columns.
# @param time_var Name of the time column.
# @return Invisibly NULL; errors when an under-supported cell is extrapolated.
# @keywords internal
check_rate_support_newdata <- function(newdata, unsupported, rate_cols,
                                       time_var = "time") {
  if (is.null(unsupported) || !nrow(unsupported)) return(invisible(NULL))
  if (!time_var %in% names(newdata)) return(invisible(NULL))

  rate_cols <- rate_cols[nzchar(rate_cols)]
  if (length(rate_cols) && !all(rate_cols %in% names(newdata))) {
    return(invisible(NULL))
  }

  if (!length(rate_cols)) {
    keys <- rep("0", nrow(newdata))
  } else {
    keys <- cell_keys(newdata, rate_cols)
  }

  # Match on the raw cell key, never on the label: the label is built with
  # sQuote(), whose output depends on options(useFancyQuotes) and so is not a
  # stable matching key. `key` is absent on fits recorded before it was stored,
  # so fall back to the old label comparison there.
  if (!is.null(unsupported$key)) {
    ref_keys <- unsupported$key
    new_keys <- keys
  } else {
    ref_keys <- unsupported$label
    new_keys <- cell_labels(keys, rate_cols)
  }

  idx <- match(new_keys, ref_keys)
  tm  <- suppressWarnings(as.numeric(newdata[[time_var]]))

  # t == 0 is exempt because the rate does not enter there, and the exemption is
  # two-sided so a backward extrapolation is caught too. Also exempt the one time
  # a tier-two cell was actually observed at, since reproducing an in-sample
  # fitted value is not extrapolation.
  obs_t <- if (!is.null(unsupported$obs_time)) unsupported$obs_time[idx] else NA_real_
  in_sample <- is.finite(obs_t) & abs(tm - obs_t) <= 1e-8 * pmax(1, abs(obs_t))

  flagged <- !is.na(idx) & is.finite(tm) & tm != 0 & !in_sample
  if (!any(flagged)) return(invisible(NULL))

  hit     <- unique(idx[flagged])
  labels  <- unsupported$label[hit]
  reasons <- unsupported$reason[hit]
  stop(
    sprintf(
      "newdata extrapolates rate parameter(s) that were never estimated from data: %s.",
      paste(sprintf("%s: %s", labels, reasons), collapse = "; ")
    ),
    call. = FALSE
  )
}

# Helpers for SB
alpha_to_Y <- function(alpha, Y0, Yinf) Y0 + (Yinf - Y0) * alpha

# Invert a symmetric positive-(semi)definite Hessian robustly.
invert_hessian <- function(hess) {
  n <- nrow(hess)
  # Diagonal (Jacobi) preconditioning: invert the dimensionless correlation-form
  # Hessian S = D^{-1} H D^{-1}, D = diag(sqrt(|diag H|)), then rescale back via
  # H^{-1} = D^{-1} S^{-1} D^{-1}.
  d <- sqrt(abs(diag(hess)))
  d[d == 0 | !is.finite(d)] <- 1
  scale <- outer(d, d)
  S <- hess / scale

  # Accept the Cholesky inverse only if it reconstructs the original Hessian.
  inv <- tryCatch(chol2inv(chol(S)), error = function(e) NULL)
  if (!is.null(inv)) {
    full <- inv / scale
    if (max(abs(hess %*% full - diag(n))) < 1e-7) return(full)
  }

  warning("Hessian is numerically singular or rank deficient; using a ",
          "pseudo-inverse. Standard errors/intervals may be unreliable.", call. = FALSE)
  MASS::ginv(S) / scale
}

# Validate a user-supplied `start` vector: names exactly the model parameters (no
# missing, unknown or duplicate names) and holds numbers. The gslnls fitters
# match start to the formula by name.
validate_start_names <- function(start, param_names) {
  # A quoted number is not a starting value: gsl_nls discards a non-numeric start
  # and self-starts. An all-NA start is meaningful (it is default_start()'s
  # self-start request), so allow it.
  vals <- if (is.list(start)) unlist(start, use.names = FALSE) else start
  if (!is.numeric(vals) && !(is.logical(vals) && all(is.na(vals)))) {
    stop(
      "`start` must be numeric, not ", class(vals)[1L],
      ". Quoted numbers such as \"80000\" are discarded, so the fit would ",
      "ignore the values you supplied.",
      call. = FALSE
    )
  }

  # A duplicate name adds a second slot for one parameter, inflating the
  # parameter count (hence residual df and sigma).
  dup_start <- unique(names(start)[duplicated(names(start))])
  if (length(dup_start)) {
    stop(
      "Duplicated name(s) in `start`: ", paste(dup_start, collapse = ", "),
      ". Supply exactly one starting value per parameter.",
      call. = FALSE
    )
  }

  missing_start <- setdiff(param_names, names(start))
  if (length(missing_start)) {
    stop(
      "Missing starting values for: ",
      paste(missing_start, collapse = ", "),
      ". Check group labels or supply start values.",
      call. = FALSE
    )
  }

  extra_start <- setdiff(names(start), param_names)
  if (length(extra_start)) {
    stop(
      "Unknown starting value name(s): ",
      paste(extra_start, collapse = ", "),
      ". Model parameters are: ", paste(param_names, collapse = ", "), ".",
      call. = FALSE
    )
  }

  invisible(start)
}

# Default activation-energy seed, in the energy units implied by the user's gas
# constant. The humidity fitters need an explicit Ea seed (gsl_nls' multi-start
# fails on them). Seeding the dimensionless ratio Ea/R keeps the same physical
# temperature sensitivity for any R.
#
# @param R Gas constant on the user's energy scale.
# @return Numeric Ea seed on the same energy scale as `R`.
# @keywords internal
default_Ea_seed <- function(R = 8.314) {
  if (!is.numeric(R) || length(R) != 1L || !is.finite(R) || R <= 0) {
    stop("`R` must be a single positive finite number.", call. = FALSE)
  }
  # 80000 J/mol / 8.314 J/mol/K = 9622.3 K, the seeded Ea/R ratio.
  (80000 / 8.314) * R
}

# Seed starting values for every model parameter.
#
# gsl_nls' self-start defaults every parameter to 1, a poor guess for this model
# family that diverges once the rate is grouped. Seeding each parameter on the
# right order of magnitude keeps grouped fits stable.
#
# The seeds are deliberately crude: the point is the right order of magnitude,
# not a good estimate. Intercepts come from the data, which self-start could not
# know.
#
# @param param_names Character vector of every parameter in the model.
# @param df Model data frame, used to seed intercepts from the response.
# @param R Gas constant, so the Ea seed matches its units.
# @param outcome_var Name of the response column; NULL to skip intercept seeding.
# @param trend "increasing" or "decreasing"; only a fallback, see intercept_seed().
# @param time_var Name of the time column, used to locate the intercept.
# @param k_seed Starting value for rate parameters. NULL (the default) seeds it
#   from the data via rate_seed(); pass a number to override.
# @param kinetic_order Reaction order of the model. 
# @return Named numeric vector of starting values.
# @keywords internal
default_start <- function(param_names, df = NULL, R = 8.314,
                          outcome_var = NULL,
                          trend = c("increasing", "decreasing"),
                          time_var = "time",
                          k_seed = NULL,
                          kinetic_order = 1L) {
  trend <- match.arg(trend)
  start <- stats::setNames(rep(NA_real_, length(param_names)), param_names)

  Ea_seed <- default_Ea_seed(R)
  # build_parallel() appends a branch letter, so match prefixes rather than
  # exact names: Ea -> EaB/EaC, Ea_log -> Ea_logB/Ea_logC, b -> bB/bC.
  is_log_Ea <- grepl("^Ea_?log", param_names)
  is_Ea     <- grepl("^Ea", param_names) & !is_log_Ea
  start[is_Ea]     <- Ea_seed
  start[is_log_Ea] <- log(Ea_seed)

  # A humidity coefficient of 0 means "no moisture effect until the data says
  # otherwise", the neutral starting point.
  start[grepl("^b$|^b[BC]$", param_names)] <- 0

  # Seed the intercept first: the rate seed below is put on the right scale for
  # the kinetic order using this value.
  int_par   <- grepl("^y0|^A0", param_names)
  intercept <- NA_real_
  if (any(int_par)) {
    intercept <- intercept_seed(df, outcome_var, time_var, trend)
    if (is.finite(intercept)) start[int_par] <- intercept
  }

  if (is.null(k_seed)) {
    k_seed <- rate_seed(df, outcome_var, time_var)
    # rate_seed() is a fractional rate (1/time); rescale to the units of k for
    # order != 1.
    if (kinetic_order != 1 && is.finite(intercept) && intercept > 0) {
      k_seed <- k_seed / intercept^(kinetic_order - 1)
    }
  }

  start[grepl("^k", param_names)] <- k_seed

  start
}

# Seed the rate constant from the observed decay.
#
# The seed only needs the right order of magnitude, but a fixed 0.01 is blind to
# the data: when the response nearly vanishes over the window the true rate is far
# larger. Reading the fractional change per unit time off the endpoints recovers
# the magnitude (for first-order kinetics y(t)/y(0) ~ 1 - k*t, so
# |y_end - y0| / (scale * span) ~ k).
#
# @param df Model data frame.
# @param resp_var Name of the response column; NULL to use the fallback.
# @param time_var Name of the time column.
# @return A single positive finite numeric rate seed.
# @keywords internal
rate_seed <- function(df, resp_var, time_var = "time") {
  floor_k <- 0.01
  ceil_k  <- 10
  if (is.null(df) || is.null(resp_var) || !resp_var %in% names(df) ||
      is.null(time_var) || !time_var %in% names(df)) {
    return(floor_k)
  }
  y  <- suppressWarnings(as.numeric(df[[resp_var]]))
  tv <- suppressWarnings(as.numeric(df[[time_var]]))
  ok <- is.finite(y) & is.finite(tv)
  if (!any(ok)) return(floor_k)
  tmin <- min(tv[ok]); tmax <- max(tv[ok])
  if (tmax <= tmin) return(floor_k)
  y0 <- mean(y[ok & tv == tmin])
  yT <- mean(y[ok & tv == tmax])
  scale <- max(diff(range(y[ok])), abs(y0), 1)
  k <- abs(yT - y0) / (scale * (tmax - tmin))
  if (is.finite(k)) min(max(k, floor_k), ceil_k) else floor_k
}

# Seed one intercept from the response at the earliest observed time.
#
# The intercept is the response at t = 0, so read it off the data rather than
# inferring it from `trend`. The trend's extreme (max for decreasing, min for
# increasing) agrees with the first timepoint only when the data moves the way
# `trend` says, which is not guaranteed (e.g. fit_parallel() defaults to
# "increasing" on data that may decay). Averaging over the t = 0 rows also handles
# an intercept split by grouping columns.
#
# @return A single finite numeric, or NA_real_ when the data cannot supply one
#   (the caller then leaves the parameter to gsl_nls self-start).
# @keywords internal
intercept_seed <- function(df, resp_var, time_var = "time",
                           trend = c("increasing", "decreasing")) {
  trend <- match.arg(trend)
  if (is.null(df) || is.null(resp_var) || !resp_var %in% names(df)) {
    return(NA_real_)
  }
  y <- suppressWarnings(as.numeric(df[[resp_var]]))
  ok <- is.finite(y)
  if (!any(ok)) return(NA_real_)

  if (!is.null(time_var) && time_var %in% names(df)) {
    tv <- suppressWarnings(as.numeric(df[[time_var]]))
    ok_t <- ok & is.finite(tv)
    if (any(ok_t)) return(mean(y[ok_t & tv == min(tv[ok_t])]))
  }
  # No usable time column: fall back to the end of the range the trend starts at.
  if (trend == "decreasing") max(y[ok]) else min(y[ok])
}

# Warn when the humidity column carries no information about `b`.
#
# The humidity term is exp(b * (rh - rh_ref)). If rh never varies, that factor is
# a single constant multiplying the rate, so b is perfectly aliased with the rate
# parameters and not identifiable: the optimiser never moves b off its seed and
# reports it as an estimate. Only a warning, since a constant-humidity study is
# reasonable and the temperature part is still estimable, but b must not be read
# as an estimate.
#
# @param df Model data frame.
# @param rh_col Humidity column name, or NULL when no humidity term is used.
# @param temp_fml Temperature expression; used to confirm `b` is actually a
#   free parameter of the model (a user expression may omit it).
# @return Invisibly NULL; called for its side effect (a warning).
# @keywords internal
check_rh_identifiable <- function(df, rh_col, temp_fml = NULL) {
  if (is.null(rh_col) || !rh_col %in% names(df)) return(invisible(NULL))

  # Nothing to warn about if the expression has no humidity coefficient.
  if (!is.null(temp_fml) &&
      !grepl("\\bb\\b", paste(as.character(temp_fml), collapse = " "))) {
    return(invisible(NULL))
  }

  rh <- suppressWarnings(as.numeric(df[[rh_col]]))
  rh <- rh[is.finite(rh)]
  if (length(rh) < 2L) return(invisible(NULL))

  # Tolerance-based: near-constant humidity is as unidentifiable in practice as
  # exactly-constant, and floating-point equality would miss it. The tolerance is
  # a design threshold, not a floating-point one: what matters is whether the
  # spread is large enough to read b off the data, 1% of the humidity scale, one
  # percentage point of RH on either coding. Real designs span 0.35 (40/60/75% RH).
  spread <- diff(range(rh))
  scale  <- max(abs(rh))
  tol    <- 0.01 * max(1, scale)
  if (spread >= tol) return(invisible(NULL))

  # Exactly constant and merely too-narrow fail the same way but report
  # differently, and the distinction is useful: with no spread the reported `b` is
  # the seed, whereas a too-narrow spread lets the optimiser chase noise and
  # return something large and confident-looking.
  extent <- if (spread == 0) {
    sprintf("is constant (%g)", rh[1])
  } else {
    sprintf("spans only %g (%g to %g), below the %g needed to estimate `b`",
            spread, min(rh), max(rh), tol)
  }
  reported <- if (spread == 0) {
    "The reported `b` is the starting value, not an estimate"
  } else {
    "The reported `b` is fitted to noise, not humidity"
  }

  warning(
    sprintf(
      paste0(
        "Humidity column `%s` %s, so the coefficient `b` is not identifiable: ",
        "exp(b * (%s - rh_ref)) is a constant absorbed into the rate. %s. Do ",
        "not read it as moisture sensitivity or predict at other humidities. ",
        "Drop the humidity term, or supply >= 2 humidity levels."
      ),
      rh_col, extent, rh_col, reported
    ),
    call. = FALSE
  )

  invisible(NULL)
}

# Warn/error when the humidity column and its reference value are on different
# scales. The humidity term is exp(b * (rh - rh_ref)), so pairing percent rh
# (40/60/75) with a fraction rh_ref (0.60) is, at fit time, an exact
# reparameterisation: b absorbs the factor of 100 and k_0 becomes the rate at the
# wrong RH, with residuals/sigma/AIC identical to the correct fit, so it must be
# caught at input time. At predict time the same mismatch is a wrong number: `b`
# is fixed, so the humidity term is evaluated on the wrong scale and the forecast
# can underflow to exactly 0. Hence `action`: warn for a fit that is at least
# internally consistent, stop when the answer would be wrong.
#
# @param df Model data frame, newdata, or schedule.
# @param rh_col Humidity column name, or NULL when no humidity term is used.
# @param rh_ref Reference humidity supplied alongside `rh_col`.
# @param what Label for the offending object in the message.
# @param action "warning" (fit time) or "error" (predict time).
# @param rh_fitted The humidity range the model was fitted on (settings$rh_range).
#   When given it replaces `rh_ref` as the scale `df` must agree with.
# @return Invisibly NULL; called for its side effect.
# @keywords internal
# The humidity range a model was fitted on. Recorded in settings so predict() can
# judge a newdata rh against the coding actually used rather than against
# `rh_ref`; see check_rh_scale(). NULL when there is no humidity term.
#
# @param df Model frame.
# @param rh_col Humidity column name, or NULL.
# @return Length-2 numeric range, or NULL.
# @keywords internal
rh_range <- function(df, rh_col) {
  if (is.null(rh_col) || !rh_col %in% names(df)) return(NULL)
  if (!is.numeric(df[[rh_col]])) return(NULL)
  rh <- df[[rh_col]][is.finite(df[[rh_col]])]
  if (!length(rh)) return(NULL)
  range(rh)
}

# The temperature range a model was fitted on, recorded in the builders'
# settings so predict() can tell an interpolation from an extrapolation.
# Same shape as rh_range(); NULL when the column is unusable.
temp_range <- function(df, temp_var) {
  if (is.null(temp_var) || !temp_var %in% names(df)) return(NULL)
  if (!is.numeric(df[[temp_var]])) return(NULL)
  tk <- df[[temp_var]][is.finite(df[[temp_var]])]
  if (!length(tk)) return(NULL)
  range(tk)
}

# Arrhenius extrapolation is exponential in 1/T, so a modest step outside the
# fitted range moves the rate constant by a large factor. The activation energy
# carries the whole extrapolation and was estimated from the fitted span alone,
# so this warns rather than errors.
check_temp_extrapolation <- function(newdata, settings, what = "newdata") {
  rng <- settings$temp_range
  temp_var <- settings$temp_var
  # Fits recorded before temp_range existed have none; skip rather than guess.
  if (is.null(rng) || !is.numeric(rng) || length(rng) != 2L ||
      !all(is.finite(rng)) || is.null(temp_var) ||
      !temp_var %in% names(newdata) || !is.numeric(newdata[[temp_var]])) {
    return(invisible(NULL))
  }
  tk <- newdata[[temp_var]][is.finite(newdata[[temp_var]])]
  if (!length(tk)) return(invisible(NULL))

  below <- tk[tk < rng[1]]
  above <- tk[tk > rng[2]]
  if (!length(below) && !length(above)) return(invisible(NULL))

  parts <- c(
    if (length(below)) sprintf("as low as %g (%g below)", min(below),
                               rng[1] - min(below)),
    if (length(above)) sprintf("as high as %g (%g above)", max(above),
                               max(above) - rng[2])
  )
  warning(
    sprintf(
      paste0("%s `%s` extrapolates outside the fitted temperature range ",
             "%g to %g: %s. The rate here rests on an activation energy ",
             "estimated from the fitted range. "),
      what, temp_var, rng[1], rng[2], paste(parts, collapse = " and ")
    ),
    call. = FALSE
  )
  invisible(NULL)
}

# A newdata column named after a fitted parameter silently wins. The closed-form
# fitters evaluate through rlang::eval_tidy(data = newdata, ...) and the formula
# path through predict.nls(); in both, the data mask outranks the parameter
# values, so such a column replaces the estimate.
# check_param_data_collisions() enforces this at fit time; this is the
# predict-time half. Errors rather than warns: no reading of a `Ea` column
# intends to silently override the fit.
check_newdata_param_shadow <- function(newdata, param_names, what = "newdata") {
  if (is.null(newdata) || is.null(param_names)) return(invisible(NULL))
  clash <- intersect(names(newdata), param_names)
  if (!length(clash)) return(invisible(NULL))
  stop(
    sprintf(
      paste0("%s column(s) %s share a name with a fitted parameter and would ",
             "take precedence over the estimate during evaluation, replacing ",
             "the fitted value. Rename or drop the column."),
      what, paste(sQuote(clash), collapse = ", ")
    ),
    call. = FALSE
  )
}

check_rh_scale <- function(df, rh_col, rh_ref, what = "df",
                           action = c("warning", "error"),
                           rh_fitted = NULL) {
  action <- match.arg(action)
  if (is.null(rh_col) || is.null(rh_ref)) return(invisible(NULL))
  if (!rh_col %in% names(df)) return(invisible(NULL))
  if (!is.numeric(df[[rh_col]])) return(invisible(NULL))  # check_numeric_cols' job

  rh <- suppressWarnings(as.numeric(df[[rh_col]]))
  rh <- rh[is.finite(rh)]
  if (!length(rh) || !is.numeric(rh_ref) || length(rh_ref) != 1L ||
      !is.finite(rh_ref)) {
    return(invisible(NULL))
  }

  rng <- range(rh)
  data_is_fraction <- rng[2] <= 1

  # What the incoming coding is judged against. At fit time `rh_ref` is the only
  # reference. At predict time it is the wrong one: pairing percent rh with
  # rh_ref = 0.60 is an exact reparameterisation of the fit, so that model's
  # predictions must be given percent rh too, and comparing against rh_ref would
  # refuse the correct coding. newdata must match the scale the model was fitted
  # on (rh_fitted). Fits recorded before rh_range was stored have none, so fall
  # back to rh_ref.
  fitted_rng <- if (is.numeric(rh_fitted) && length(rh_fitted) &&
                    all(is.finite(rh_fitted))) range(rh_fitted) else NULL
  if (is.null(fitted_rng)) {
    ref_is_fraction <- rh_ref <= 1
    ref_desc <- sprintf("`rh_ref` is %g", rh_ref)
    remedy   <- sprintf(", or pass rh_ref = %g",
                        if (ref_is_fraction) rh_ref * 100 else rh_ref / 100)
  } else {
    ref_is_fraction <- fitted_rng[2] <= 1
    ref_desc <- sprintf("the model was fitted on humidity ranging %g to %g",
                        fitted_rng[1], fitted_rng[2])
    # No rh_ref advice here: rh_ref is fixed by the fit and cannot be changed at
    # predict time, so suggesting it would send the reader down a dead end.
    remedy   <- " to match the fitted scale"
  }
  if (ref_is_fraction == data_is_fraction) return(invisible(NULL))

  # The consequence half of the message: what the mismatch does to a fit versus
  # to a prediction from a fit that was itself scaled correctly.
  consequence <- if (action == "warning") {
    "`k_0` then becomes the rate at the wrong RH; the fit looks flawless, with identical residuals, sigma and AIC."
  } else {
    paste0("`b` is fixed by the fit, so the humidity term is evaluated on the ",
           "wrong scale and the prediction can underflow to exactly 0.")
  }

  (if (action == "warning") warning else stop)(
    sprintf(
      paste0("%s humidity column `%s` ranges %g to %g but %s. This looks like ",
             "%s-coded humidity paired with a %s-scale reference. %s %s `%s` by ",
             "100%s."),
      what, rh_col, rng[1], rng[2], ref_desc,
      if (data_is_fraction) "fraction" else "percent",
      if (ref_is_fraction) "fraction" else "percent",
      consequence,
      if (data_is_fraction) "Multiply" else "Divide", rh_col, remedy
    ),
    call. = FALSE
  )

  invisible(NULL)
}

# Type- and scale-check the predictor columns of a prediction request.
#
# The humidity and numeric-type guards also run at fit time; this applies the same
# contract at predict time, so a percent-coded rh or a factor tempk in newdata is
# rejected rather than silently mishandled.
#
# @param newdata Prediction data or schedule.
# @param settings The fit's build$settings.
# @param what Label for the offending object in messages.
# @param cols Extra columns to type-check (predict_tos() passes the schedule's
#   own arithmetic columns, which are not named in settings).
# @return Invisibly NULL; errors on a non-numeric or mis-scaled predictor.
# @keywords internal
check_predictor_cols <- function(newdata, settings, what = "newdata",
                                 cols = character()) {
  # covar_cols covers rh_col plus any second covariate a custom temp_expr
  # declares; rh_col stays in the union for fits recorded before it existed.
  check_numeric_cols(newdata,
                     unique(c(cols, settings$time_var, settings$temp_var,
                              settings$rh_col, settings$covar_cols)),
                     what = what)
  # Judged against the humidity the model was fitted on, not against rh_ref;
  # see check_rh_scale()'s rh_fitted argument.
  check_rh_scale(newdata, settings$rh_col, settings$rh_ref, what = what,
                 action = "error", rh_fitted = settings$rh_range)
  invisible(NULL)
}

get_param <- function(p, prefix) {
  v           <- p[startsWith(names(p), prefix)]
  names(v)    <- sub(prefix, "", names(v))   # names == lot labels
  v
}
make_group <- function(x, cols, default_label = "0") {
  if (length(cols) == 0) {
    return(factor(rep(default_label, nrow(x)), levels = default_label))
  }
  x2 <- sanitize_group_values(x, cols)
  interaction(x2[cols], drop = TRUE, sep = "||")
}

#' Resolve a temperature expression into its component parts
#'
#' Shared by the \code{fit_*()} functions.  When \code{temp_expr} is supplied it
#' is unpacked into its components; when it is \code{NULL} a standard Arrhenius
#' expression is built from \code{R}, \code{trefk}, \code{temp_var} and
#' \code{log_Ea} (in which case \code{rh_col} and \code{rh_ref} are \code{NULL}).
#'
#' @section Supplied elements override the arguments; absent ones do not:
#' Only the elements the list actually carries are taken from it, so a hand-built
#' \code{temp_expr} that omits one keeps the caller's argument rather than
#' overwriting it with \code{NULL}.
#'
#' @param temp_expr Optional temperature expression list, or \code{NULL}.
#' @param R Gas constant. Used when \code{temp_expr} is \code{NULL} or omits it.
#' @param trefk Reference temperature in Kelvin, same rule.
#' @param temp_var Name of the temperature column, same rule.
#' @param log_Ea Logical; log-scale activation energy, same rule.
#'
#' @return A named list with elements \code{temp_expr}, \code{R}, \code{trefk},
#'   \code{temp_var}, \code{log_Ea}, \code{rh_col}, \code{rh_ref},
#'   \code{covar_cols} and \code{temp_fml}. \code{covar_cols} is every data
#'   column the expression is declared to reference besides \code{temp_var},
#'   \code{rh_col} included.
#' @keywords internal
resolve_temp_expr <- function(temp_expr = NULL, R = 8.314, trefk = 298.15,
                              temp_var = "tempk", log_Ea = FALSE) {
  if (!is.null(temp_expr)) {
    if (!is.list(temp_expr)) {
      stop("`temp_expr` must be a list such as expr_arrhenius() returns, not ",
           class(temp_expr)[1], ". Pass the whole list, not just its formula.",
           call. = FALSE)
    }
    # Keep the caller's value when the list is silent about an element.
    pick <- function(nm, default) {
      v <- unlist(temp_expr[[nm]])
      if (is.null(v)) default else v
    }
    R        <- pick("R", R)
    trefk    <- pick("trefk", trefk)
    temp_var <- pick("temp_var", temp_var)
    log_Ea   <- pick("log_Ea", log_Ea)

    # These have no meaningful default: absent means "this model does not use
    # one", so NULL is the answer rather than a fallback.
    rh_col     <- unlist(temp_expr[["rh_col"]])
    rh_ref     <- unlist(temp_expr[["rh_ref"]])
    covar_cols <- unlist(temp_expr[["covar_cols"]])
    temp_fml   <- unlist(temp_expr[["formula"]])

    # The formula is the one element with no substitute.
    if (is.null(temp_fml) || !length(temp_fml) || !nzchar(temp_fml[1])) {
      stop("`temp_expr` has no `formula` element. A temperature expression must ",
           "carry the formula string itself; see expr_arrhenius().",
           call. = FALSE)
    }
    if (length(covar_cols) && !is.character(covar_cols)) {
      stop("`temp_expr$covar_cols` must be a character vector of column names.",
           call. = FALSE)
    }
  } else {
    temp_expr <- expr_arrhenius(R = R, trefk = trefk, temp_var = temp_var, log_Ea = log_Ea)
    temp_fml  <- unlist(temp_expr[["formula"]])
    rh_col     <- NULL
    rh_ref     <- NULL
    covar_cols <- NULL
  }

  list(
    temp_expr = temp_expr,
    R         = R,
    trefk     = trefk,
    temp_var  = temp_var,
    log_Ea    = log_Ea,
    rh_col    = rh_col,
    rh_ref    = rh_ref,
    # One vector for "every declared data reference other than temp_var", so the
    # fitters have a single thing to check for existence, type and NAs, and a
    # single thing to hand the tokeniser.
    covar_cols = unique(c(rh_col, covar_cols)),
    temp_fml  = temp_fml
  )
}


#' Identify the free parameters of a custom temperature expression
#'
#' Parses the expression and returns every symbol that is not a declared data
#' column, the temperature variable, or a language constant.
#'
#' This replaces per-builder pattern sniffing.  Fixed greps such as
#' \code{grepl("\\\\bb\\\\s*\\\\*", expr)} are position-dependent:
#' \code{exp(b * (rh - 0.6))} finds \code{b} while the identical
#' \code{exp((rh - 0.6) * b)} does not, and a renamed parameter is missed.
#'
#' @section Why the parser and not a tokeniser:
#' Classification uses \code{all.vars()} on the parsed expression rather than
#' splitting text on \code{[^[:alnum:]_]} and subtracting a whitelist.  The
#' whitelist cannot be made correct: it over-reports (\code{tanh}, \code{pmax},
#' \code{ifelse}, \code{TRUE} were all read as parameters) and mangles scientific
#' notation (\code{1e-9} becomes \code{1e}).  R's parser already distinguishes
#' symbols from calls and literals, so none of that guesswork is needed.
#'
#' Data references are identified from the \emph{declared} covariates
#' (\code{temp_var}, \code{covar_cols}) rather than by matching columns of
#' \code{df}.  Testing "is there a column of that name?" is circular: a stray
#' column named \code{b} would reclassify the humidity coefficient as data and
#' drop it.
#'
#' @param temp_expr Character string, or NULL / function (returns \code{character(0)}).
#' @param temp_var Name of the temperature column.
#' @param covar_cols Character vector of additional declared data columns the
#'   expression may reference (e.g. \code{rh_col}), as assembled by
#'   \code{resolve_temp_expr()}.
#'
#' @return Character vector of parameter names, in order of first appearance.
#' @keywords internal
temp_expr_params <- function(temp_expr, temp_var = "tempk", covar_cols = NULL) {
  if (is.null(temp_expr) || is.function(temp_expr)) return(character())

  parsed <- tryCatch(
    parse(text = as.character(temp_expr)),
    error = function(e) {
      stop("The temperature expression is not valid R: ", conditionMessage(e),
           "\nExpression: ", as.character(temp_expr)[1], call. = FALSE)
    }
  )

  # all.vars() returns symbols used as variables, not function names or literals.
  # `pi`, `T` and `F` are the constants R exposes as symbols, so they are the
  # only subtraction still needed.
  setdiff(all.vars(parsed), c(temp_var, covar_cols, "pi", "T", "F"))
}


# A parameter name that also names a column of the fitting data is ambiguous.
# The residual-function fitters evaluate with rlang::eval_tidy(expr, data = data,
# env = penv); `data` is a mask that sits between the expression and `env`, so
# the column silently wins and the parameter it shadows drops out.
check_param_data_collisions <- function(params, df, what = "parameter") {
  clash <- intersect(unique(as.character(params)), names(df))
  if (!length(clash)) return(invisible(NULL))

  stop(
    sprintf("Data column(s) %s collide with model %s name(s). ",
            paste(sQuote(clash), collapse = ", "), what),
    "The column would mask the parameter during evaluation, silently removing ",
    "it. Rename the column(s), or rename the parameter(s) in your temperature ",
    "expression.",
    call. = FALSE
  )
}


#' Compute tolerance interval bounds via Krishnamoorthy & Mathew (2009)
#'
#' Shared low-level helper for \code{compute_tolerance_intervals()} (gsl_nls
#' models) and the Delta Method path in \code{predict.fit_sb()} (nls.lm models).
#' Model-specific extraction happens in the caller; this operates purely on
#' pre-computed quantities.
#'
#' @param mu_hat Numeric vector of point predictions.
#' @param d_squared Numeric vector of standardised distances
#'   \eqn{d^2 = \nabla f_0^\top (F^\top F)^{-1} \nabla f_0}.
#' @param sigma Residual standard error \eqn{\hat\sigma}.
#' @param df Residual degrees of freedom \eqn{n - p}.
#' @param coverage Coverage proportion \eqn{\gamma} (0, 1).
#' @param level Confidence level \eqn{1 - \alpha} (0, 1).
#' @param side One of \code{"two"}, \code{"lower"}, \code{"upper"}.
#'
#' @return A data.frame with columns \code{fit}, \code{lwr}, \code{upr}.
#' @keywords internal
#' @noRd
tolerance_bounds_km <- function(mu_hat, d_squared, sigma, df, coverage, level,
                                side = c("two", "lower", "upper")) {
  
  side <- match.arg(side)
  
  if (coverage <= 0 || coverage >= 1) {
    stop("'coverage' must be between 0 and 1", call. = FALSE)
  }
  if (level <= 0 || level >= 1) {
    stop("'level' must be between 0 and 1", call. = FALSE)
  }
  if (df <= 0) {
    stop("Insufficient degrees of freedom for tolerance intervals", call. = FALSE)
  }
  if (any(!is.finite(sigma)) || any(sigma < 0)) {
    stop("'sigma' must be non-negative and finite", call. = FALSE)
  }
  if (any(!is.finite(d_squared)) || any(d_squared < 0)) {
    stop("'d_squared' must be non-negative and finite", call. = FALSE)
  }
  
  d       <- sqrt(d_squared)
  alpha   <- 1 - level
  z_gamma <- stats::qnorm(coverage)
  
  if (side == "two") {
    # Wallis (1951) approximation (K&M Eq. 3.3.6, page 87)
    chi2_ncp     <- stats::qchisq(coverage, df = 1, ncp = d_squared)
    chi2_central <- stats::qchisq(alpha, df = df)
    k <- sqrt(df * chi2_ncp / chi2_central)
    
    lwr <- mu_hat - k * sigma
    upr <- mu_hat + k * sigma
  } else {
    # One-sided: noncentral t-distribution (K&M Eq. 3.2.6, page 63)
    ncp <- ifelse(d > 1e-10, z_gamma / d, 0)
    k   <- d * stats::qt(1 - alpha, df = df, ncp = ncp)
    
    if (side == "lower") {
      lwr <- mu_hat - k * sigma
      upr <- rep(Inf, length(mu_hat))
    } else {
      lwr <- rep(-Inf, length(mu_hat))
      upr <- mu_hat + k * sigma
    }
  }
  
  data.frame(fit = mu_hat, lwr = lwr, upr = upr)
}

#' Check newdata grouping levels against the fitted levels
#'
#' Grouping values must match the fitted labels exactly.
#'
#' @param newdata Data frame of new observations.
#' @param settings The build settings list containing \code{intercept_cols},
#'   \code{rate_cols}, \code{group_levels} (the fitted levels per column) and
#'   \code{rate_cells} (the observed rate combinations).
#'
#' @return Invisibly \code{NULL}; errors on a level or rate cell the model was not
#'   fitted on.
#' @keywords internal
#' @noRd
check_newdata_levels <- function(newdata, settings) {

  # Rate checks run even without grouping columns: a single global rate can also
  # be unidentified.
  group_cols <- unique(c(settings$intercept_cols, settings$rate_cols))
  group_cols <- group_cols[nzchar(group_cols)]

  if (length(group_cols)) {

    # Older fits do not record group_levels.
    if (is.null(settings$group_levels)) {
      stop("This model was fitted with an earlier version of stablr. ",
           "Refit it to predict on new data.", call. = FALSE)
    }

    # Check that grouping columns exist in newdata
    missing_group <- setdiff(group_cols, names(newdata))
    if (length(missing_group)) {
      stop("newdata missing grouping variables: ",
           paste(missing_group, collapse = ", "), call. = FALSE)
    }

    # Without this, an NA group value is reported as an "unseen level" named NA,
    # sending the user after a level they never fitted rather than a failed join.
    check_structural_nas(newdata, group_cols, what = "newdata")

    for (col in group_cols) {
      unseen <- setdiff(unique(as.character(newdata[[col]])),
                        settings$group_levels[[col]])
      if (length(unseen)) {
        stop(
          sprintf(
            "newdata contains unseen levels in column '%s': %s. ",
            col, paste(sQuote(unseen, FALSE), collapse = ", ")
          ),
          "Prediction is only possible for levels present in the training data.",
          call. = FALSE
        )
      }
    }
  }

  rate_cols <- if (is.null(settings$rate_cols)) character() else settings$rate_cols
  time_var  <- if (is.null(settings$time_var)) "time" else settings$time_var

  # The per-column loop validates each level in isolation. Rate parameters are
  # keyed on whole combinations, so check the tuple too.
  check_observed_cells(newdata, settings$rate_cells, rate_cols,
                       time_var = time_var)

  # A cell that was observed can still have had too little support to identify
  # its rate; extrapolating it is the same flat-forecast hazard.
  check_rate_support_newdata(newdata, settings$unsupported_rate_cells,
                            rate_cols, time_var)

  invisible(NULL)
}

