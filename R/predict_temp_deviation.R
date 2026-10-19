# ═══════════════════════════════════════════════════════════════════════════════
# Predict Time Out of Storage (TOS)
# ═══════════════════════════════════════════════════════════════════════════════

#' Identify data columns a fitted temperature expression depends on
#'
#' A custom temperature expression such as the one produced by
#' \code{expr_arrhenius_humidity()} may reference data columns in addition to
#' the temperature variable (i.e., \code{temp_var}). Those columns are neither model
#' parameters nor grouping variables, so the composite propagation used by
#' \code{predict_tos()} must supply them per segment.
#'
#' Variables are taken from \code{all.vars()} on the parsed expression, which
#' excludes function names and numeric literals by construction. 
#'
#' @param temp_fml Character scalar; the fitted temperature expression.
#' @param param_names Character vector of fitted parameter names.
#' @param temp_var Name of the temperature column.
#'
#' @return Character vector of data-column names (possibly empty).
#' @keywords internal
#' @noRd
.temp_expr_covariates <- function(temp_fml, param_names, temp_var) {

  temp_fml <- as.character(temp_fml)
  if (!length(temp_fml) || is.na(temp_fml[1]) || !nzchar(temp_fml[1])) {
    return(character())
  }

  vars <- unique(unlist(lapply(as.list(parse(text = temp_fml[1])), all.vars)))

  setdiff(vars, c(temp_var, param_names, "pi"))
}

#' Extract per-segment covariate values from a schedule
#'
#' Validates that every column required by the temperature expression is
#' present in the schedule and returns them as a plain list of vectors.
#'
#' @param covar_cols Character vector of required column names.
#' @param schedule The (already ordered) schedule data frame.
#'
#' @return Named list of vectors, one per required column.
#' @keywords internal
#' @noRd
.schedule_covariates <- function(covar_cols, schedule) {

  if (!length(covar_cols)) return(list())

  missing_cols <- setdiff(covar_cols, names(schedule))
  if (length(missing_cols)) {
    stop(
      sprintf(
        paste0("The fitted temperature expression depends on %s, so schedule ",
               "must supply %s a value for every segment."),
        paste(sprintf("'%s'", missing_cols), collapse = ", "),
        if (length(missing_cols) > 1) "those columns with" else "that column with"
      ),
      call. = FALSE
    )
  }

  vals <- lapply(covar_cols, function(cc) schedule[[cc]])
  names(vals) <- covar_cols

  na_cols <- covar_cols[vapply(vals, function(v) anyNA(v), logical(1))]
  if (length(na_cols)) {
    stop(
      sprintf("schedule contains NA values in temperature covariate(s): %s",
              paste(na_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  vals
}

#' Check schedule grouping columns against their fitted levels
#'
#' Schedule grouping values must match the fitted labels exactly. For
#' \code{fit_sb()}, values are sanitized before matching.
#'
#' @param schedule The schedule data frame.
#' @param group_cols Character vector of grouping columns to check.
#' @param known_levels Named list of the levels each column was fitted on
#'   (\code{build$settings$group_levels}, or for \code{fit_sb} the sanitized
#'   levels in \code{build$data}). A \code{NULL} element skips that column.
#' @param sanitize_values Logical; sanitize the schedule values before matching
#'   and return them sanitized (\code{fit_sb}).
#' @param rate_cols,rate_cells The fitted rate grouping columns and the
#'   combinations of them that were observed, from \code{build$settings}. Used to
#'   reject a schedule naming a combination that has no fitted rate; \code{NULL}
#'   skips the check.
#' @param unsupported Under-supported rate cells from
#'   \code{build$settings$unsupported_rate_cells}. Used to reject a schedule that
#'   extrapolates a rate which was never identified; \code{NULL} skips the check.
#'
#' @return \code{schedule}, with its grouping columns sanitized when
#'   \code{sanitize_values = TRUE} and otherwise unchanged.
#' @keywords internal
#' @noRd
.match_schedule_groups <- function(schedule, group_cols, known_levels,
                                   sanitize_values = FALSE,
                                   rate_cols = NULL, rate_cells = NULL,
                                   unsupported = NULL) {

  # Older fits do not record group_levels.
  if (!sanitize_values && is.null(known_levels) &&
      length(group_cols[nzchar(group_cols)])) {
    stop("This model was fitted with an earlier version of stablr. ",
         "Refit it to use predict_tos().", call. = FALSE)
  }

  for (col in group_cols) {
    if (!col %in% names(schedule)) next

    raw  <- as.character(schedule[[col]])
    vals <- if (sanitize_values) sanitize(raw) else raw

    unseen <- setdiff(unique(vals), known_levels[[col]])
    if (!is.null(known_levels[[col]]) && length(unseen)) {
      stop(
        sprintf("schedule contains unseen levels in column '%s': %s",
                col, paste(sQuote(unique(raw[vals %in% unseen]), FALSE),
                           collapse = ", ")),
        call. = FALSE
      )
    }
    if (sanitize_values) schedule[[col]] <- vals
  }

  # Per-column checks above cannot see an unobserved combination, which has no
  # fitted rate and would integrate to a flat, non-degrading profile.
  check_observed_cells(schedule, rate_cells, rate_cols, what = "schedule")

  # Same hazard for a cell that was observed but too thinly to identify its rate.
  # The guard keys on elapsed time, which for a schedule is the cumulative t_end.
  if (!is.null(unsupported) && nrow(unsupported)) {
    sched_time <- schedule
    sched_time$time <- schedule$t_end
    check_rate_support_newdata(sched_time, unsupported, rate_cols, "time")
  }

  schedule
}

# ═══════════════════════════════════════════════════════════════════════════════
# build function for fit_kinetic
# ═══════════════════════════════════════════════════════════════════════════════
.build_composite_fn_kinetic <- function(fit, schedule) {
  build     <- attr(fit, "build")

  # Recover the kinetic model configuration from the formula environment,
  # where build_kinetic() captured its arguments. 
  fn_env        <- environment(formula(fit))
  kinetic_order <- get("kinetic_order", envir = fn_env)
  trend         <- get("trend", envir = fn_env)
  Yinf          <- if (exists("Yinf", envir = fn_env)) get("Yinf", envir = fn_env) else build$settings$Yinf

  # Extract grouping columns from settings
  intercept_cols <- build$settings$intercept_cols
  rate_cols      <- build$settings$rate_cols
  group_cols <- unique(c(intercept_cols, rate_cols))

  # Reconstruct intercept and rate expressions from schedule
  temp_var <- build$settings$temp_var
  time_var <- build$settings$time_var

  tempk_vec <- schedule$tempk
  dt_vec    <- schedule$duration
  n_seg     <- length(dt_vec)

  # Check the schedule's grouping levels against the fitted levels.
  schedule_checked <- .match_schedule_groups(schedule, group_cols,
                                              build$settings$group_levels,
                                              rate_cols = rate_cols,
                                              rate_cells = build$settings$rate_cells,
                                              unsupported = build$settings$unsupported_rate_cells)

  # Reuse the intercept and rate expressions that were actually fitted, taken
  # from the formula's closure like kinetic_order / trend above. Re-deriving them
  # from the schedule cannot recover reference coding: expr_intercept() takes its
  # reference from the levels it is given, so a sub-schedule carrying one level of
  # a second intercept column would make that level its own reference and drop its
  # coefficient.
  int_part  <- get("int_part", envir = fn_env)
  rate_part <- get("rate_part", envir = fn_env)

  arr_expr <- as.character(get("arr_expr", envir = fn_env))

  # Every fitted parameter must appear in something the composite evaluates; one
  # that does not is silently ignored, and neither the prediction nor its interval
  # shows it.
  check_composite_params(c(int_part$expr, rate_part$expr, arr_expr),
                        build$param_names)

  # Data columns the temperature expression needs (e.g., relative humidity).
  covar_cols <- .temp_expr_covariates(arr_expr, build$param_names, temp_var)
  covar_vals <- .schedule_covariates(covar_cols, schedule)

  composite_fn <- function(params, data) {

    penv <- list2env(as.list(params), parent = baseenv())

    # Build initial data with grouping columns
    init_data <- data.frame(1)
    names(init_data) <- time_var
    init_data[[time_var]] <- 0
    init_data[[temp_var]] <- tempk_vec[1]

    # Add grouping columns to init_data
    for (col in group_cols) {
      if (col %in% names(schedule_checked)) {
        init_data[[col]] <- schedule_checked[[col]][1]
      }
    }

    # Evaluate intercept at initial conditions
    y0 <- rlang::eval_tidy(rlang::parse_expr(int_part$expr), data = init_data, env = penv)

    # Propagate the state sequentially through each segment, applying the
    # closed-form kinetic solution for the fitted kinetic_order / trend over
    # each segment's duration. 
    y_current <- y0
    for (i in seq_len(n_seg)) {
      seg_data <- data.frame(1)
      names(seg_data) <- time_var
      seg_data[[time_var]] <- dt_vec[i]
      seg_data[[temp_var]] <- tempk_vec[i]

      # Add grouping columns to seg_data
      for (col in group_cols) {
        if (col %in% names(schedule_checked)) {
          seg_data[[col]] <- schedule_checked[[col]][i]
        }
      }

      # Supply this segment's temperature-expression covariates (e.g., rh)
      for (cc in covar_cols) {
        seg_data[[cc]] <- covar_vals[[cc]][i]
      }

      # Evaluate rate and temperature dependence, then k(T) * dt for this segment
      k_eval    <- rlang::eval_tidy(rlang::parse_expr(rate_part$expr), data = seg_data, env = penv)
      arr_eval  <- rlang::eval_tidy(rlang::parse_expr(arr_expr), data = seg_data, env = penv)
      rate_time <- k_eval * arr_eval * dt_vec[i]

      if (kinetic_order == 0) {
        if (trend == "increasing") {
          y_current <- y_current + rate_time
        } else {
          y_current <- y_current - rate_time
        }
      } else if (kinetic_order == 1) {
        if (trend == "increasing") {
          y_current <- Yinf - (Yinf - y_current) * exp(-rate_time)
        } else {
          y_current <- y_current * exp(-rate_time)
        }
      } else {
        pow <- 1 - kinetic_order
        if (trend == "increasing") {
          y_current <- Yinf - ((Yinf - y_current)^pow - rate_time * pow)^(1 / pow)
        } else {
          y_current <- (y_current^pow - rate_time * pow)^(1 / pow)
        }
      }
    }

    y_current
  }

  composite_fn
}


# ═══════════════════════════════════════════════════════════════════════════════
# build function for single_deg (A -> B)
# ═══════════════════════════════════════════════════════════════════════════════
.build_composite_fn_singledeg <- function(fit, schedule) {

  build <- attr(fit, "build")

  # Extract grouping columns from settings
  intercept_cols <- build$settings$intercept_cols
  rate_cols      <- build$settings$rate_cols
  group_cols <- unique(c(intercept_cols, rate_cols))

  fn <- formula(fit)
  fn_env <- environment(fn)

  # Extract model configuration from the closure environment
  kinetic_order <- get("kinetic_order", envir = fn_env)
  trend         <- get("trend", envir = fn_env)
  ie_A          <- get("ie_A", envir = fn_env)
  ie_B          <- get("ie_B", envir = fn_env)
  rate_part     <- get("rate_part", envir = fn_env)
  arr_expr      <- get("arr_expr", envir = fn_env)
  temp_var      <- get("temp_var", envir = fn_env)
  time_var      <- get("time_var", envir = fn_env)
  Yinf          <- if (exists("Yinf", envir = fn_env)) get("Yinf", envir = fn_env) else NULL

  tempk_vec <- schedule$tempk
  dt_vec    <- schedule$duration
  n_seg     <- length(dt_vec)

  # Data columns the temperature expression needs (e.g. relative humidity).
  arr_expr   <- as.character(arr_expr)
  covar_cols <- .temp_expr_covariates(arr_expr, build$param_names, temp_var)
  covar_vals <- .schedule_covariates(covar_cols, schedule)

  check_composite_params(c(ie_A$expr, ie_B$expr, rate_part$expr, arr_expr),
                        build$param_names)

  # Check the schedule's grouping levels against the fitted levels.
  schedule_checked <- .match_schedule_groups(schedule, group_cols,
                                              build$settings$group_levels,
                                              rate_cols = rate_cols,
                                              rate_cells = build$settings$rate_cells,
                                              unsupported = build$settings$unsupported_rate_cells)

  composite_fn <- function(params, data) {

    penv <- list2env(as.list(params), parent = baseenv())

    init_data <- data.frame(1)
    names(init_data) <- time_var
    init_data[[time_var]] <- 0
    init_data[[temp_var]] <- tempk_vec[1]

    # Add grouping columns to init_data
    for (col in group_cols) {
      if (col %in% names(schedule_checked)) {
        init_data[[col]] <- schedule_checked[[col]][1]
      }
    }

    A_int <- rlang::eval_tidy(rlang::parse_expr(ie_A$expr), data = init_data, env = penv)
    B_int <- rlang::eval_tidy(rlang::parse_expr(ie_B$expr), data = init_data, env = penv)

    yA_current <- A_int

    for (i in seq_len(n_seg)) {
      seg_data <- data.frame(1)
      names(seg_data) <- time_var
      seg_data[[time_var]] <- dt_vec[i]
      seg_data[[temp_var]] <- tempk_vec[i]

      # Add grouping columns to seg_data
      for (col in group_cols) {
        if (col %in% names(schedule_checked)) {
          seg_data[[col]] <- schedule_checked[[col]][i]
        }
      }

      # Supply this segment's temperature-expression covariates (e.g. rh)
      for (cc in covar_cols) {
        seg_data[[cc]] <- covar_vals[[cc]][i]
      }

      # Evaluate rate at this segment's temperature
      rate_part_eval <- rlang::eval_tidy(rlang::parse_expr(rate_part$expr),
                                         data = seg_data, env = penv)
      arr_expr_eval  <- rlang::eval_tidy(rlang::parse_expr(arr_expr),
                                         data = seg_data, env = penv)
      rate <- rate_part_eval * arr_expr_eval
      rate_time <- rate * dt_vec[i]

      # Apply kinetic model for this segment
      if (kinetic_order == 0) {
        if (trend == "increasing") {
          yA_current <- yA_current + rate_time
        } else {
          yA_current <- yA_current - rate_time
        }
      } else if (kinetic_order == 1) {
        if (trend == "increasing") {
          yA_current <- Yinf - (Yinf - yA_current) * exp(-rate_time)
        } else {
          yA_current <- yA_current * exp(-rate_time)
        }
      } else {
        pow <- 1 - kinetic_order
        if (trend == "increasing") {
          yA_current <- Yinf - ((Yinf - yA_current)^pow - rate_time * pow)^(1/pow)
        } else {
          yA_current <- (yA_current^pow - rate_time * pow)^(1/pow)
        }
      }
    }

    yA_hat <- yA_current
    yB_hat <- B_int + (A_int - yA_hat)

    c(yA_hat, yB_hat)
  }

  composite_fn
}


# ═══════════════════════════════════════════════════════════════════════════════
# build function for parallel (A -> B + C)
# ═══════════════════════════════════════════════════════════════════════════════
.build_composite_fn_parallel <- function(fit, schedule) {

  build <- attr(fit, "build")

  # Extract grouping columns from settings
  intercept_cols <- build$settings$intercept_cols
  rate_cols      <- build$settings$rate_cols
  group_cols <- unique(c(intercept_cols, rate_cols))

  fn <- formula(fit)
  fn_env <- environment(fn)

  kinetic_order <- get("kinetic_order", envir = fn_env)
  ie_A          <- get("ie_A", envir = fn_env)
  ie_B          <- get("ie_B", envir = fn_env)
  ie_C          <- get("ie_C", envir = fn_env)
  rp_B          <- get("rp_B", envir = fn_env)
  rp_C          <- get("rp_C", envir = fn_env)
  arr_B         <- get("arr_B", envir = fn_env)
  arr_C         <- get("arr_C", envir = fn_env)
  temp_var      <- get("temp_var", envir = fn_env)
  time_var      <- get("time_var", envir = fn_env)

  tempk_vec <- schedule$tempk
  dt_vec    <- schedule$duration
  n_seg     <- length(dt_vec)

  # Data columns the branch temperature expressions need (e.g. relative
  # humidity). 
  arr_B <- as.character(arr_B)
  arr_C <- as.character(arr_C)
  covar_cols <- union(
    .temp_expr_covariates(arr_B, build$param_names, temp_var),
    .temp_expr_covariates(arr_C, build$param_names, temp_var)
  )
  covar_vals <- .schedule_covariates(covar_cols, schedule)
  check_composite_params(c(ie_A$expr, ie_B$expr, ie_C$expr, rp_B$expr,
                          rp_C$expr, arr_B, arr_C), build$param_names)

  # Check the schedule's grouping levels against the fitted levels.
  schedule_checked <- .match_schedule_groups(schedule, group_cols,
                                              build$settings$group_levels,
                                              rate_cols = rate_cols,
                                              rate_cells = build$settings$rate_cells,
                                              unsupported = build$settings$unsupported_rate_cells)

  composite_fn <- function(params, data) {

    penv <- list2env(as.list(params), parent = baseenv())

    # Initial conditions
    init_data <- data.frame(1)
    names(init_data) <- time_var
    init_data[[time_var]] <- 0
    init_data[[temp_var]] <- tempk_vec[1]

    # Add grouping columns to init_data
    for (col in group_cols) {
      if (col %in% names(schedule_checked)) {
        init_data[[col]] <- schedule_checked[[col]][1]
      }
    }

    A_int <- rlang::eval_tidy(rlang::parse_expr(ie_A$expr), data = init_data, env = penv)
    B_int <- rlang::eval_tidy(rlang::parse_expr(ie_B$expr), data = init_data, env = penv)
    C_int <- rlang::eval_tidy(rlang::parse_expr(ie_C$expr), data = init_data, env = penv)

    # Sequential propagation
    yA_current <- A_int
    cumul_B_gain <- 0
    cumul_C_gain <- 0

    for (i in seq_len(n_seg)) {
      seg_data <- data.frame(1)
      names(seg_data) <- time_var
      seg_data[[time_var]] <- dt_vec[i]
      seg_data[[temp_var]] <- tempk_vec[i]

      # Add grouping columns to seg_data
      for (col in group_cols) {
        if (col %in% names(schedule_checked)) {
          seg_data[[col]] <- schedule_checked[[col]][i]
        }
      }

      # Supply this segment's temperature-expression covariates (e.g., rh)
      for (cc in covar_cols) {
        seg_data[[cc]] <- covar_vals[[cc]][i]
      }

      kB_i <- rlang::eval_tidy(rlang::parse_expr(rp_B$expr), data = seg_data, env = penv) *
        rlang::eval_tidy(rlang::parse_expr(arr_B), data = seg_data, env = penv)
      kC_i <- rlang::eval_tidy(rlang::parse_expr(rp_C$expr), data = seg_data, env = penv) *
        rlang::eval_tidy(rlang::parse_expr(arr_C), data = seg_data, env = penv)

      kApp_i <- kB_i + kC_i
      rate_time <- kApp_i * dt_vec[i]

      yA_before <- yA_current
      if (kinetic_order == 0) {
        yA_current <- yA_current - rate_time
      } else if (kinetic_order == 1) {
        yA_current <- yA_current * exp(-rate_time)
      } else {
        pow <- 1 - kinetic_order
        yA_current <- (yA_current^pow - rate_time * pow)^(1/pow)
      }

      # Amount of A lost in this segment
      dA_lost_i <- yA_before - yA_current

      # Distribute to B and C based on this segment's selectivity
      cumul_B_gain <- cumul_B_gain + (kB_i / kApp_i) * dA_lost_i
      cumul_C_gain <- cumul_C_gain + (kC_i / kApp_i) * dA_lost_i
    }

    yA_hat <- yA_current
    yB_hat <- B_int + cumul_B_gain
    yC_hat <- C_int + cumul_C_gain

    c(yA_hat, yB_hat, yC_hat)
  }

  composite_fn
}


# ═══════════════════════════════════════════════════════════════════════════════
# build function for fit_sb (alpha-space ODE)
# ═══════════════════════════════════════════════════════════════════════════════
.build_composite_fn_sb <- function(fit, schedule) {

  build <- attr(fit, "build")

  # Extract grouping columns from settings
  intercept_cols <- build$settings$intercept_cols
  rate_cols      <- build$settings$rate_cols
  group_cols <- unique(c(intercept_cols, rate_cols))

  # Extract settings
  temp_var <- build$settings$temp_var
  time_var <- build$settings$time_var
  Yinf     <- build$settings$Yinf
  temp_params <- build$settings$temp_params
  temp_fml <- build$settings$temp_fml
  R <- build$settings$R
  trefk <- build$settings$trefk

  tempk_vec <- schedule$tempk
  dt_vec    <- schedule$duration
  n_seg     <- length(dt_vec)

  # Sanitize grouping columns from schedule. 
  known_levels <- stats::setNames(
    lapply(group_cols, function(col) {
      if (col %in% names(build$data)) unique(as.character(build$data[[col]]))
    }),
    group_cols
  )
  schedule_sanitized <- .match_schedule_groups(schedule, group_cols,
                                                known_levels,
                                                sanitize_values = TRUE,
                                                rate_cols = rate_cols,
                                                unsupported = build$settings$unsupported_rate_cells)

  # Parse temperature expression once
  if (!is.null(temp_fml)) {
    temp_parsed <- parse(text = temp_fml)
    temp_fml_chr <- as.character(temp_fml)
  } else {
    # Default Arrhenius, referring to the temperature by the fitted `temp_var`
    temp_fml_chr <- sprintf("exp(Ea / %f * (1/%f - 1/%s))", R, trefk, temp_var)
    temp_parsed  <- parse(text = temp_fml_chr)
  }

  # Data columns the temperature expression needs (e.g. relative humidity).
  covar_cols <- .temp_expr_covariates(temp_fml_chr, build$param_names, temp_var)
  covar_vals <- .schedule_covariates(covar_cols, schedule)

  # Define ODE function for sequential integration
  growth_ode <- function(t, state, pars) {
    alpha <- state[1]
    k     <- pars["k"]
    n     <- as.numeric(pars["n"])
    m     <- as.numeric(pars["m"])

    # Create evaluation environment with all parameters.  The temperature is
    # bound under the fitted `temp_var` name, which is what temp_fml references.
    eval_env <- list2env(as.list(pars), parent = baseenv())
    eval_env[[temp_var]] <- as.numeric(pars[[temp_var]])

    # Expose the temperature-expression covariates (e.g. rh) carried in `pars`
    for (cc in covar_cols) {
      eval_env[[cc]] <- as.numeric(pars[[cc]])
    }

    # Evaluate temperature factor
    temp_factor <- eval(temp_parsed, envir = eval_env)

    dalpha <- as.numeric(k) * temp_factor * (1 - alpha)^n * alpha^m

    list(dalpha)
  }

  composite_fn <- function(params, data) {

    # Extract parameters
    k_vec  <- get_param(params, "k_")
    y0_vec <- get_param(params, "y0_")
    n_par  <- unname(params["n"])
    m_par  <- unname(params["m"])

    # Extract temperature parameters
    temp_pars <- params[temp_params]

    # Build grouping IDs for the first segment
    rate_id <- make_group(schedule_sanitized[1, , drop = FALSE], rate_cols, default_label = "0")
    int_id  <- make_group(schedule_sanitized[1, , drop = FALSE], intercept_cols, default_label = "0")

    k_g  <- unname(k_vec[as.character(rate_id)[1]])
    y0_g <- unname(y0_vec[as.character(int_id)[1]])

    # Initialize alpha at near-zero (consistent with fit_sb)
    alpha_current <- 1e-6

    # Propagate through each segment using ODE solver
    for (i in seq_len(n_seg)) {
      T_i  <- tempk_vec[i]
      dt_i <- dt_vec[i]

      # Skip zero-duration segments
      if (dt_i <= 0) next

      # Build parameter vector for ODE solver, carrying this segment's
      # temperature-expression covariates (e.g. rh) alongside the parameters.
      pars_full <- c(k = k_g, temp_pars, n = n_par, m = m_par,
                     setNames(T_i, temp_var))
      for (cc in covar_cols) {
        pars_full <- c(pars_full, setNames(covar_vals[[cc]][i], cc))
      }
      sol <- deSolve::ode(
        y = c(alpha = alpha_current),
        times = c(0, dt_i),
        func = growth_ode,
        parms = pars_full,
        method = "lsoda",
        rtol = 1e-8,
        atol = 1e-10
      )

      # Update alpha to end of this segment
      alpha_current <- sol[nrow(sol), "alpha"]

      # Clamp to [0, 1] for safety
      alpha_current <- min(max(alpha_current, 0), 1)
    }

    # Convert alpha to Y
    Y_pred <- alpha_to_Y(alpha_current, Y0 = y0_g, Yinf = Yinf)
    Y_pred
  }

  composite_fn
}


# ═══════════════════════════════════════════════════════════════════════════════
# dispatch to correct composite builder
# ═══════════════════════════════════════════════════════════════════════════════
.build_composite_fn <- function(fit, schedule, model_type) {
  switch(model_type,
         "fit_kinetic" = .build_composite_fn_kinetic(fit, schedule),
         "fit_singledeg" = .build_composite_fn_singledeg(fit, schedule),
         "fit_parallel" = .build_composite_fn_parallel(fit, schedule),
         "fit_sb" = .build_composite_fn_sb(fit, schedule),
         stop("Unknown model type: ", model_type)
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# newdata for composite function
# ═══════════════════════════════════════════════════════════════════════════════
.build_composite_newdata <- function(schedule, model_type, fit) {

  build <- attr(fit, "build")

  # Extract grouping columns from settings
  intercept_cols <- build$settings$intercept_cols
  rate_cols      <- build$settings$rate_cols
  group_cols <- unique(c(intercept_cols, rate_cols))

  if (model_type == "fit_kinetic") {
    time_var <- build$settings$time_var
    temp_var <- build$settings$temp_var

    newdata <- data.frame(
      placeholder_time = sum(schedule$duration),
      placeholder_temp = schedule$tempk[1]
    )
    names(newdata) <- c(time_var, temp_var)

    # Add grouping columns
    for (col in group_cols) {
      if (col %in% names(schedule)) {
        newdata[[col]] <- as.character(schedule[[col]][1])
      }
    }
  } else if (model_type == "fit_sb") {
    time_var <- build$settings$time_var
    temp_var <- build$settings$temp_var

    newdata <- data.frame(
      placeholder_time = 0,
      placeholder_temp = schedule$tempk[1]
    )
    names(newdata) <- c(time_var, temp_var)

    # Add sanitized grouping columns (fit_sb keys groups on sanitized labels)
    for (col in group_cols) {
      if (col %in% names(schedule)) {
        val <- sanitize(as.character(schedule[[col]][1]))
        newdata[[col]] <- val
      }
    }
  } else {
    fn <- formula(fit)
    fn_env <- environment(fn)
    time_var <- get("time_var", envir = fn_env)
    temp_var <- get("temp_var", envir = fn_env)

    newdata <- data.frame(
      placeholder_time = 0,
      placeholder_temp = schedule$tempk[1]
    )
    names(newdata) <- c(time_var, temp_var)

    # Add grouping columns
    for (col in group_cols) {
      if (col %in% names(schedule)) {
        newdata[[col]] <- as.character(schedule[[col]][1])
      }
    }
  }

  newdata
}


#' Times at which to evaluate a thermal schedule
#'
#' \code{predict_tos()} evaluates the composite model at a set of absolute times,
#' each paired with the schedule segment it falls in. With
#' \code{time_step = NULL} that set is the segment boundaries, which is the
#' one-row-per-segment output; a numeric \code{time_step} adds interior points
#' spaced that far apart.
#'
#' A point at exactly \code{t_start} is not attributed to the segment starting
#' there: it is the previous segment's end (or, for the first segment, the time-0
#' baseline row), and would otherwise be reported twice with two temperatures.
#' Each segment therefore owns \code{(t_start, t_end]}, and \code{t_end} is
#' always included so that no boundary is lost to a step that does not divide
#' the segment evenly. A zero-duration segment owns no interior point, so it
#' falls back to its boundary and keeps its single row.
#'
#' @param schedule The ordered, validated schedule, carrying \code{duration}.
#' @param time_step Positive numeric spacing, or \code{NULL} for boundaries only.
#'
#' @return A data frame of \code{time} and \code{segment}, in increasing time
#'   order, whose first row is the time-0 baseline (\code{segment = 0}).
#' @keywords internal
#' @noRd
.tos_eval_points <- function(schedule, time_step) {

  # The baseline row, ahead of every segment.
  times <- 0
  segs  <- 0L

  for (i in seq_len(nrow(schedule))) {
    t0 <- schedule$t_start[i]
    t1 <- schedule$t_end[i]

    if (is.null(time_step)) {
      pts <- t1
    } else {
      tol <- 1e-8 * max(1, abs(t0), abs(t1))
      pts <- seq(t0, t1, by = time_step)
      # Drop anything at or past t_end -- seq() can land a hair inside it -- then
      # append the boundary itself, so it appears exactly once.
      pts <- pts[pts < t1 - tol]
      pts <- c(pts[pts > t0 + tol], t1)
    }

    times <- c(times, pts)
    segs  <- c(segs, rep.int(i, length(pts)))
  }

  data.frame(time = times, segment = segs)
}


#' Manual delta method prediction intervals for composite temperature schedule
#'
#' @param fit A gsl_nls fitted model object
#' @param schedule A data.frame with columns "t_start", "t_end", "temp", plus a
#'   column for every grouping variable used in the fit. If the model was fitted
#'   with a temperature expression that depends on additional covariates (e.g.
#'   the relative humidity column from \code{\link{expr_arrhenius_humidity}})
#'   the schedule must also carry that column, giving its value for each segment.
#'   \code{t_start}, \code{t_end} and \code{temp} must be \strong{numeric}.
#'   The schedule must begin at \code{t_start = 0} and its segments must be
#'   contiguous, such that each segment's \code{t_start} must equal the previous
#'   segment's \code{t_end}, leaving no gaps or overlaps in the timeline.
#' @param level Confidence level for prediction intervals
#' @param var Which variable to return (NULL for fit_kinetic)
#' @param temp_unit "C" or "K"
#' @param time_step Optional spacing, in the model's time units, at which to
#'   predict \emph{within} each segment. By default (\code{NULL}) only the segment
#'   boundaries are reported, one row per segment. With \code{time_step = 1} on a
#'   segment running from 0 to 12, the returned table carries times 0, 1, 2, ...,
#'   12. 
#' @return A data frame with columns \code{segment}, \code{time}, \code{temp},
#'   \code{variable}, \code{fit}, \code{lwr} and \code{upr}, in increasing
#'   \code{time} order, with the interval a prediction interval at \code{level}.
#'   The first row is the time-0 baseline (\code{segment = 0}), where \code{fit}
#'   is the model's fitted intercept for this unit. 
#' @export
predict_tos <- function(fit,
                                   schedule,
                                   level = 0.95,
                                   var = NULL,
                                   temp_unit = "C",
                                   time_step = NULL) {

  # ───1. Validate schedule ─────────────────────────────────────────────────────────────
  # Accept gsl_nls (fit_kinetic, fit_singledeg, fit_parallel) or nls.lm (fit_sb)
  if (!inherits(fit, "gsl_nls") && !inherits(fit, "nls.lm")) {
    stop("fit must be a gsl_nls or nls.lm object (from fit_kinetic, fit_singledeg, fit_parallel, or fit_sb)")
  }
  required_cols <- c("t_start", "t_end", "temp")
  if (!all(required_cols %in% names(schedule))) {
    stop("schedule must contain columns: ", paste(required_cols, collapse = ", "))
  }
  # Arithmetic inputs, so type-check before any of them is used: a factor or
  # character temp otherwise yields a full row of NA fit/se/lwr/upr, or a bare
  # coercion error, depending on temp_unit.
  check_numeric_cols(schedule, required_cols, what = "schedule")

  # No rows means nothing to rbind, and the attributes below land on NULL.
  if (nrow(schedule) == 0L) {
    stop("schedule has no rows; supply at least one segment.", call. = FALSE)
  }

  # NA here propagates through the composite into a full row of NA fit/se/lwr/upr
  na_cols <- required_cols[vapply(required_cols,
                                  function(cc) anyNA(schedule[[cc]]),
                                  logical(1))]
  if (length(na_cols)) {
    stop(sprintf("schedule contains NA in column(s): %s. Every segment needs a %s.",
                 paste(na_cols, collapse = ", "),
                 paste(required_cols, collapse = ", ")),
         call. = FALSE)
  }

  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("level must be a single number strictly between 0 and 1 (e.g. 0.95).",
         call. = FALSE)
  }

  # A non-positive or non-finite step makes seq() error deep inside the grid
  # builder, or spin forever; check it here where the remedy can be named.
  if (!is.null(time_step) &&
      (!is.numeric(time_step) || length(time_step) != 1L ||
       !is.finite(time_step) || time_step <= 0)) {
    stop(paste0("time_step must be a single positive number (e.g. time_step = 1), ",
                "or NULL to report only the segment boundaries."),
         call. = FALSE)
  }

  fit_build <- attr(fit, "build")

  # Covariate types and humidity scale, on the scale the fit used.
  check_predictor_cols(schedule, fit_build$settings, what = "schedule")

  # Every grouping column the fit uses must be supplied; .match_schedule_groups()
  fit_group_cols <- unique(c(fit_build$settings$intercept_cols,
                             fit_build$settings$rate_cols))
  fit_group_cols <- fit_group_cols[!is.na(fit_group_cols) & nzchar(fit_group_cols)]
  missing_groups <- setdiff(fit_group_cols, names(schedule))
  if (length(missing_groups)) {
    stop(sprintf(paste0("schedule is missing grouping column(s) the model was ",
                        "fitted with: %s. Supply a value for each segment."),
                 paste(missing_groups, collapse = ", ")),
         call. = FALSE)
  }
  varying_groups <- fit_group_cols[vapply(
    fit_group_cols,
    function(cc) length(unique(as.character(schedule[[cc]]))) > 1L,
    logical(1)
  )]
  if (length(varying_groups)) {
    stop(sprintf(paste0("schedule grouping column(s) %s vary across segments: %s. ",
                        "A schedule is one unit, so each must be constant; call ",
                        "predict_tos() once per unit."),
                 paste(sprintf("'%s'", varying_groups), collapse = ", "),
                 paste(vapply(varying_groups, function(cc)
                   sprintf("%s = %s", cc,
                           paste(unique(as.character(schedule[[cc]])),
                                 collapse = "/")),
                   character(1)), collapse = "; ")),
         call. = FALSE)
  }

  # Order by t_start
  schedule <- schedule[order(schedule$t_start, schedule$t_end), ]

  temp_unit <- match.arg(temp_unit, c("C", "K"))
  if (temp_unit == "C") {
    schedule$tempk <- schedule$temp + 273.15
  } else {
    schedule$tempk <- schedule$temp
  }
  if (any(schedule$tempk <= 0)) {
    stop(sprintf(paste0("schedule contains a non-physical absolute temperature: ",
                        "%s K. Check `temp` and `temp_unit` (temp_unit = '%s')."),
                 format(min(schedule$tempk)), temp_unit),
         call. = FALSE)
  }
  # Catches Celsius passed with temp_unit = "K"
  implausible <- schedule$tempk < 223.15 | schedule$tempk > 333.15
  if (any(implausible)) {
    warning(
      sprintf(paste0("schedule contains temperature(s) outside the usual ",
                     "stability range: %s K. ",
                     "If you meant Celsius, pass temp_unit = 'C' (temp_unit is ",
                     "currently '%s')."),
              paste(format(unique(schedule$tempk[implausible])), collapse = ", "),
              temp_unit),
      call. = FALSE
    )
  }
  schedule$duration <- schedule$t_end - schedule$t_start
  if (any(schedule$duration < 0)) {
     stop("schedule contains segments with t_end < t_start.", call. = FALSE)
   }

  # The timeline must be anchored at t = 0. 
  if (nrow(schedule) > 0) {
    t0 <- schedule$t_start[1]
    # Relative tolerance: accepts floating-point noise, rejects any real t0.
    if (!is.na(t0) && abs(t0) > 1e-8 * max(1, abs(t0))) {
      stop(
        sprintf(
          paste0("schedule must begin at time 0, but its first segment starts at ",
                 "%s."),
          format(t0), format(schedule$t_end[1] - t0), format(t0)
        ),
        call. = FALSE
      )
    }
  }

  # Segments must tile the timeline
  if (nrow(schedule) > 1) {
    n     <- nrow(schedule)
    starts <- schedule$t_start[-1]
    ends   <- schedule$t_end[-n]
    # Relative tolerance, so fractional times are not rejected for FP noise.
    scale  <- pmax(1, abs(starts), abs(ends))
    breaks <- which(abs(starts - ends) > 1e-8 * scale)
    if (length(breaks)) {
      i <- breaks[1]
      stop(
        sprintf(
          paste0("schedule segments must be contiguous: segment %d ends at %s ",
                 "but segment %d starts at %s (no gaps or overlaps)."),
          i, format(schedule$t_end[i]), i + 1L,
          format(schedule$t_start[i + 1L])
        ),
        call. = FALSE
      )
    }
  }

  # ───2. Detect model type and validate var call────────────────────────────────────────────────
  model_type <- class(fit)[1]
  if (model_type == "fit_singledeg") {
    if (is.null(var) || !var %in% c("yA", "yB")) {
      stop("For fit_singledeg models, must specify var = 'yA' or var = 'yB'")
    }
  }
  if (model_type == "fit_parallel") {
    if (is.null(var) || !var %in% c("yA", "yB","yC")) {
      stop("For fit_parallel models, must specify var = 'yA' or var = 'yB' or var = 'yC'")
    }
  }
  # fit_sb and fit_kinetic have single response, no var needed
  var_names <- switch(model_type,
                      "fit_kinetic" = "y",
                      "fit_singledeg" = c("yA", "yB"),
                      "fit_parallel" = c("yA", "yB", "yC"),
                      "fit_sb" = "y"
  )

  # ─── 3. Extract fitted parameters and covariance ─────────────────────────────
  theta_hat <- stats::coef(fit)
  V_theta   <- stats::vcov(fit)

  # capture.output(): summary.fit_sb() 
  utils::capture.output(fit_summary <- summary(fit))
  sigma2    <- fit_summary$sigma^2
  df_resid  <- fit_summary$df[2]
  free_cols <- seq_along(theta_hat)
  if (ncol(V_theta) != length(theta_hat)) {
    free_cols <- match(colnames(V_theta), names(theta_hat))
    if (anyNA(free_cols)) {
      stop(sprintf(paste0("cannot align vcov(fit) (%d parameter(s)) with ",
                          "coef(fit) (%d parameter(s)): %s. The delta-method ",
                          "variance needs the two to refer to the same ",
                          "parameters."),
                   ncol(V_theta), length(theta_hat),
                   paste(setdiff(colnames(V_theta), names(theta_hat)),
                         collapse = ", ")),
           call. = FALSE)
    }
  }

  alpha  <- 1 - level
  t_crit <- stats::qt(1 - alpha / 2, df = df_resid)

  # ─── 4. Evaluate the composite over [0, t] for each requested time t ─────────
  eval_pts <- .tos_eval_points(schedule, time_step)

  # Each point is a fresh composite build plus a Jacobian, so a step far finer
  # than the schedule is not a slow call but an unfinishable one.
  max_pts <- 5000L
  if (nrow(eval_pts) > max_pts) {
    span <- max(schedule$t_end)
    stop(sprintf(paste0("time_step = %s over length %s asks for %d prediction ",
                        "times, above the limit of %d. Each is a full delta-method ",
                        "evaluation; use a coarser time_step (%s stays under)."),
                 format(time_step), format(span), nrow(eval_pts),
                 max_pts, format(signif(span / max_pts, 2))),
         call. = FALSE)
  }

  # t = 0 is the baseline: the same composite evaluated over zero elapsed time,
  # which reduces to the fitted intercept for every model type. 
  predict_at <- function(j) {

    t_eval <- eval_pts$time[j]
    k      <- eval_pts$segment[j]

    if (k == 0L) {
      sub_schedule <- schedule[1L, , drop = FALSE]
      sub_schedule$t_start  <- 0
      sub_schedule$t_end    <- 0
      sub_schedule$duration <- 0
    } else {
      sub_schedule <- schedule[seq_len(k), , drop = FALSE]
      sub_schedule$t_end[k]    <- t_eval
      sub_schedule$duration[k] <- t_eval - sub_schedule$t_start[k]
    }

    composite_fn_k <- .build_composite_fn(fit, sub_schedule, model_type)
    # Placeholder newdata; the composite has the schedule baked in.
    newdata_k <- .build_composite_newdata(sub_schedule, model_type, fit)
    y_hat_k <- composite_fn_k(theta_hat, newdata_k)

    J_k <- numDeriv::jacobian(
      func = function(theta) composite_fn_k(theta, newdata_k),
      x = theta_hat
    )
    J_k <- J_k[, free_cols, drop = FALSE]   # see free_cols above

    # Delta method: mean variance, then prediction variance.
    var_mean_k <- diag(J_k %*% V_theta %*% t(J_k))
    var_pred_k <- var_mean_k + sigma2

    se_mean_k <- sqrt(var_mean_k)
    se_pred_k <- sqrt(var_pred_k)

    if (any(!is.finite(y_hat_k)) || any(!is.finite(se_mean_k))) {
      stop(
        sprintf(paste0("prediction for %s is not finite (fit = %s)."),
                if (k == 0L) "the time-0 baseline"
                else sprintf("time %s (segment %d)", format(t_eval), k),
                paste(format(y_hat_k), collapse = ", ")),
        call. = FALSE
      )
    }

    data.frame(
      segment  = k,
      time     = t_eval,
      temp     = schedule$temp[max(k, 1L)],
      variable = var_names,
      fit      = y_hat_k,
      lwr      = y_hat_k - t_crit * se_pred_k,
      upr      = y_hat_k + t_crit * se_pred_k,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, lapply(seq_len(nrow(eval_pts)), predict_at))

  # ───5. Filter by var if specified ───────────────────────────────────────────
  if (!is.null(var)) {
    # Match against the labels actually in `result$variable`
    var_match <- match.arg(var, var_names)
    result <- result[result$variable == var_match, , drop = FALSE]
  }

  attr(result, "schedule") <- schedule[, c("t_start", "t_end", "temp", "tempk", "duration")]
  attr(result, "level") <- level
  attr(result, "model_type") <- model_type

  rownames(result) <- NULL
  result
}
