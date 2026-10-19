
#' Build alpha-space ODE model objects (residual fn, param names, groups)
#'
#' Creates a residual function suitable for \code{minpack.lm::nls.lm()},
#' plus the fully enumerated parameter names given grouping in
#' \code{rate_cols} (-> k_) and \code{intercept_cols} (-> y0_).
#'
#' The state is alpha is 0 < alpha < 1; predictions in Y-space use
#' Y = Y0 + (Yinf - Y0) * alpha. Temperature dependence is supplied ready-built
#' via \code{temp_fml}.
#'
#' @param df Data frame with outcome, time, temperature, and optionally covariates.
#' @param intercept_cols,rate_cols Character vectors of grouping columns for
#'   intercept and rate effects.
#' @param outcome_var,temp_var,time_var Column names for outcome, temperature
#'   (Kelvin) and time.
#' @param R Gas constant (default  8.314 J K-1 mol-1).
#' @param trefk Reference temperature (Kelvin) used in Arrhenius.
#' @param log_Ea Logical; selects the activation-energy parameter name used when
#'   enumerating parameters ("Ealog" when \code{TRUE}, otherwise "Ea"). The
#'   parameterization itself comes from \code{temp_fml}, so this must agree with
#'   the expression that was passed in.
#' @param trend \code{"increasing"} or \code{"decreasing"} outcome trajectory.
#' @param Yinf Asymptote the trajectory travels toward. \code{NULL} (default)
#'   takes 100 when \code{trend = "increasing"} and 0 when
#'   \code{trend = "decreasing"}. It should lie beyond every observation in the
#'   declared direction: above \code{max(y)} for an increasing trend, below
#'   \code{min(y)} for a decreasing one. If it does not, the fit still proceeds
#'   but a warning is emitted.
#' @param temp_fml Character string from \code{expr_arrhenius()$formula} or \code{expr_arrhenius_humidity()$formula}.
#' @param rh_col Column name for relative humidity values (only needed when expr_arrhenius_humidity is used).
#' @param rh_ref Reference relative humidity value (only needed when expr_arrhenius_humidity is used).
#' @param covar_cols Character vector of every data column the temperature
#'   expression is declared to reference besides \code{temp_var} (\code{rh_col}
#'   included). Assembled by \code{resolve_temp_expr()}; symbols not named here are
#'   treated as free parameters.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{residuals}: function(par, data) -> numeric residual vector.
#'   \item \code{param_names}: character vector of parameter names (in expected order).
#'   \item \code{rate_levels}, \code{intercept_levels}: character vectors of grouping labels.
#'   \item \code{settings}: list with R, trefk, log_Ea, Ea_name, temp_fml, temp_params, etc.
#'   \item \code{predict_run}: low-level ODE predictor (time grid, tempk, pars).
#'   \item \code{data}: the model frame with sanitized grouping values and the
#'     internal \code{..rate_id} / \code{..intercept_id} columns, reused as the
#'     default \code{newdata} by \code{\link{predict.fit_sb}}.
#' }
#' @export
build_alpha_ode <- function(
    df,
    trend = c("increasing", "decreasing"),
    Yinf = NULL,
    intercept_cols = character(),
    rate_cols = character(),
    log_Ea = FALSE,
    outcome_var = "y",
    temp_var = "tempk",
    time_var = "time",
    R = 8.314,
    trefk = 298.15,
    rh_col = NULL,
    rh_ref = NULL,
    covar_cols = NULL,
    temp_fml
) {


  # ── 1. input validation ───────────────────────────────
  trend <- match.arg(trend)

  if (is.null(Yinf)) {
    Yinf <- if (trend == "increasing") 100 else 0
  } else {
    check_scalar_Yinf(Yinf)
  }

  # Check required and optional columns exist in df
  missing_cols <- setdiff(
    c(outcome_var, temp_var, time_var, intercept_cols, rate_cols, covar_cols),
    names(df)
  )
  if (length(missing_cols)) {
    stop(
      sprintf("Missing columns in df: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  check_nonempty_df(df)

  # Check NAs
  if (anyNA(df[[outcome_var]])) stop("Remove or impute NAs in the outcome before fitting.", call. = FALSE)
  if (!all(is.finite(df[[outcome_var]]))) stop("Non-finite (Inf) values in the outcome; remove or impute them before fitting.", call. = FALSE)

  # Check NA in temp_var, time_var, or a grouping column. 
  check_structural_nas(df, c(intercept_cols, rate_cols, temp_var, time_var, covar_cols))
  check_temp_identifiable(df, temp_var)

  # Temperature, time and humidity must be numeric. 
  check_numeric_cols(df, c(temp_var, time_var, covar_cols))

  # The ODE initial condition is placed at the first requested time, so a
  # negative time relocates it off t = 0 and redefines y0.
  check_nonnegative_times(df, time_var)

  y_obs <- df[[outcome_var]]; t_obs <- df[[time_var]]

  if (trend == "increasing" && Yinf <= max(y_obs)) {
    warning(sprintf(
      "Yinf (%g) does not lie above every observation when trend = \"increasing\"; max(%s) is %.4g. Fitting anyway.",
      Yinf, outcome_var, max(y_obs)), call. = FALSE)
  }
  if (trend == "decreasing" && Yinf >= min(y_obs)) {
    warning(sprintf(
      "Yinf (%g) does not lie below every observation when trend = \"decreasing\"; min(%s) is %.4g. Fitting anyway.",
      Yinf, outcome_var, min(y_obs)), call. = FALSE)
  }

  # Y0 is read off the observed time-zero readings, so replicate noise there can
  # push a neighbouring point a hair outside [Y0, Yinf]. Keep a tolerance so the
  # geometry warning flags real violations rather than ordinary assay scatter.
  tol <- max(1e-6, 0.05 * diff(range(y_obs)))

  grp_cols <- unique(c(intercept_cols, rate_cols))
  int_grp <- if (length(grp_cols)) interaction(df[grp_cols], drop = TRUE)
             else factor(rep("0", nrow(df)))
  outside <- character(0)
  for (g in levels(int_grp)) {
    yg <- y_obs[int_grp == g]; tg <- t_obs[int_grp == g]
    y0_obs <- yg[tg == min(tg)]
    y0_g <- if (trend == "increasing") min(y0_obs) else max(y0_obs)
    if (any(yg < min(y0_g, Yinf) - tol | yg > max(y0_g, Yinf) + tol)) {
      outside <- c(outside, sprintf("%s (Y0 = %.4g)", g, y0_g))
    }
  }
  if (length(outside)) {
    warning(sprintf(
      "Some observations fall outside the range spanned by Y0 and Yinf (%g) in: %s. Fitting anyway.",
      Yinf, paste(outside, collapse = "; ")), call. = FALSE)
  }

  df_model <- sanitize_group_values(df, c(intercept_cols, rate_cols))

  # Canonical row order so the fit is reproducible from any row permutation. 
  sort_cols  <- c(intercept_cols, rate_cols, temp_var, covar_cols,
                  time_var, outcome_var)
  sort_order <- do.call(order, df_model[sort_cols])
  df_model   <- df_model[sort_order, , drop = FALSE]
  rownames(df_model) <- NULL

  # ── 1. Grouping ───────────────────────────────
  
  rate_id      <- make_group(df_model, rate_cols,      default_label = "0")
  intercept_id <- make_group(df_model, intercept_cols, default_label = "0")
  
  rate_levels      <- levels(rate_id)
  intercept_levels <- levels(intercept_id)
  
  # add IDs to data (these are used inside residuals_fn)
  df_model$..rate_id      <- rate_id
  df_model$..intercept_id <- intercept_id
  
  # ── 3. Parameter naming ───────────────────────────────
  # Read the temperature parameters straight from the expression
  temp_param_names <- temp_expr_params(temp_fml, temp_var,
                                       unique(c(rh_col, covar_cols)))
  Ea_name <- grep(if (log_Ea) "^Ea_?log|^Ealog" else "^Ea",
                  temp_param_names, value = TRUE)[1]
  if (is.na(Ea_name)) Ea_name <- if (log_Ea) "Ealog" else "Ea"

  # global/shared: temperature params + n, m
  shared <- c(temp_param_names, "n", "m")

  # per-group: k_<rate>, y0_<intercept>
  k_names  <- paste0("k_",  rate_levels)
  int_names <- paste0("y0_", intercept_levels)

  param_names <- c(shared, k_names, int_names)


  # growth_ode() reads its parameters out of one named vector that also carries
  # temperature (and RH) under the caller's column name, so a temp_var/rh_col
  # named after a parameter would silently overwrite it, so reject it here. The
  # other builders get this from check_param_data_collisions().
  reserved <- unique(c("k", "n", "m", temp_param_names, Ea_name))
  clash    <- intersect(c(temp_var, covar_cols), reserved)
  if (length(clash)) {
    stop(
      sprintf("Column(s) %s share a name with a model parameter (%s). ",
              paste(sQuote(clash), collapse = ", "),
              paste(reserved, collapse = ", ")),
      "parameters, so the column would silently overwrite one and the fit would ",
      "converge to meaningless estimates. Rename the column.",
      call. = FALSE
    )
  }

  check_param_data_collisions(param_names, df)

  # ── 4. Helpers ───────────────────────────────

  # Parse temperature expression once if provided
  temp_parsed <- parse(text = temp_fml)
 
  growth_ode <- function(t, state, pars) {
    # Clamp the state into [0, 1]
    alpha <- min(max(state[1], 0), 1)
    k     <- pars["k"]
    n     <- as.numeric(pars["n"])
    m     <- as.numeric(pars["m"])

    # Create evaluation environment with all parameters
    eval_env <- list2env(as.list(pars), parent = baseenv())

    # Bind temperature under the caller's temp_var name (what temp_fml refers to),
    # not a hardcoded "tempk".
    eval_env[[temp_var]] <- as.numeric(pars[[temp_var]])

    # Bind every declared covariate the expression may reference, not just rh_col.
    for (cc in covar_cols) {
      eval_env[[cc]] <- as.numeric(pars[[cc]])
    }
    
    # Evaluate temperature factor
    temp_factor <- eval(temp_parsed, envir = eval_env)
    
    dalpha <- as.numeric(k) * temp_factor * (1 - alpha)^n * alpha^m
    
    list(dalpha)
  }
  
  # `covars`: named numeric vector of declared covariate values for this group,
  # one element per covar_cols entry, named by column.
  predict_run <- function(time_vec, tempk, pars, covars = NULL) {

    # deSolve places the initial condition (alpha = 1e-6) at times[1], so t = 0
    # must come first or y0 stops meaning Y(0).
    times_unique <- sort(unique(time_vec))
    if (times_unique[1] > 0)
      times_unique <- c(0, times_unique)
    # An all-zero cell leaves a length-1 grid, which deSolve rejects. Pad to two
    # times; approx() below still returns the baseline at t = 0.
    if (length(times_unique) < 2L)
      times_unique <- c(times_unique, times_unique[length(times_unique)] + 1)

    # Backstop for the t = 0 invariant; callers are screened by
    # check_nonnegative_times().
    if (times_unique[1] < 0) {
      stop("internal error: alpha-space ODE solved from a negative time (",
           format(times_unique[1]), "); the initial condition would not be at ",
           "t = 0. Please report this.", call. = FALSE)
    }

    # Temperature travels in `pars` keyed by temp_var (the name temp_fml uses);
    # the `tempk` argument name is this closure's public interface.
    names(pars) <- c("k", temp_param_names, "n", "m")
    pars_full <- c(pars, setNames(tempk, temp_var))

    # Covariates travel alongside, keyed by column name, for growth_ode() to bind.
    if (length(covars)) {
      pars_full <- c(pars_full, covars)
    }

    sol <- deSolve::ode(
      y      = c(alpha = 1e-6),
      times  = times_unique,
      func   = growth_ode,
      parms  = pars_full,
      method = "lsoda",
      rtol   = 1e-8,
      atol   = 1e-10
    )
    sol_df <- as.data.frame(sol)

    # find the non-time column dynamically
    state_cols <- setdiff(names(sol_df), "time")
    if (length(state_cols) != 1) stop("Expected one state column, got: ", paste(state_cols, collapse = ", "))

    # A truncated or blown-up solve returns too few rows and/or non-finite states.
    # Return NA so the residual guard steers the optimiser away, rather than let
    # approx(rule = 2) crash or fabricate a flat clamped tail.
    if (nrow(sol_df) < length(times_unique) ||
        sum(is.finite(sol_df[[state_cols]])) < 2L) {
      return(rep(NA_real_, length(time_vec)))
    }

    stats::approx(
      x    = sol_df$time,
      y    = sol_df[[state_cols]],
      xout = time_vec,
      method = "linear",
      rule   = 2
    )$y

  }
  
  # Residual function for nls.lm ----------------------------------------------
  residuals_fn <- function(p, data) {

    ## global (shared) kinetic parameters -------------------------
    # Extract temperature parameters (could be Ea, or Ea+b for humidity, etc.)
    temp_pars <- p[temp_param_names]
    n     <- p["n"]
    m     <- p["m"]

    ## vector of per‑lot y0 values
    ## per‑covariate parameters -----------------------------------------
    k_vec      <- get_param(p, "k_")         #  Aₗₒₜ
    y0_vec <- get_param(p, "y0_")
    
    ## placeholder for predictions
    pred_Y <- numeric(nrow(data))
    
    # Working copies of IDs (as characters to index the named vectors)
    rate_chr <- as.character(data$..rate_id)
    int_chr  <- as.character(data$..intercept_id)
    
    # Split by (rate_id, intercept_id, tempk, every covariate) to avoid repeated
    # ODE solves; grouping on rh_col alone would merge rows differing in a second
    # covariate.
    split_idx <- split(
      seq_len(nrow(data)),
      interaction(c(list(rate_chr, int_chr, data[[temp_var]]),
                    lapply(covar_cols, function(cc) data[[cc]])),
                  drop = TRUE, sep = "||")
    )
    
    pred_Y <- numeric(nrow(data))
    
    for (idx in split_idx) {
      
      if (length(idx) == 0) next
      
      rate_id  <- rate_chr[idx][1]
      int_id  <- int_chr[idx][1]
      k_g   <- k_vec[[rate_id]]
      y0_g  <- y0_vec[[int_id]]
      
      Tval  <- data[[temp_var]][idx][1]
      times <- data[[time_var]][idx]

      # Build parameter vector with appropriate temperature parameters
      if (!is.null(temp_fml)) {
        pars_run <- c(k = k_g, temp_pars, n = n, m = m)
      } else {
        pars_run <- c(k = k_g, setNames(temp_pars, Ea_name), n = n, m = m)
      }

      # One value per declared covariate for this group; constant within it
      # because the split above includes every covariate.
      covar_val <- if (length(covar_cols)) {
        stats::setNames(vapply(covar_cols,
                               function(cc) as.numeric(data[[cc]][idx][1]),
                               numeric(1)), covar_cols)
      } else NULL

      alpha_pred <- predict_run(times, Tval, pars_run, covars = covar_val)
      
      # guard alpha in [0,1] and finiteness
      if (any(!is.finite(alpha_pred)) || any(alpha_pred < -1e-6) || any(alpha_pred > 1 + 1e-6)) {
        return(Inf)  # forces nls.lm to move away
      }
      
      pred_Y[idx] <- alpha_to_Y(alpha_pred, Y0 = y0_g, Yinf = Yinf)
    }
    
    pred_Y - data[[outcome_var]]
  }
  
  list(
    residuals    = residuals_fn,
    param_names  = param_names,
    rate_levels  = rate_levels,
    intercept_levels = intercept_levels,
    settings     = list(
      R = R, trefk = trefk, log_Ea = log_Ea, Ea_name = Ea_name,
      temp_fml = temp_fml, temp_params = temp_param_names,
      outcome_var = outcome_var, time_var = time_var, temp_var = temp_var,
      trend = trend,
      Yinf = Yinf,
      rh_col = rh_col,
      covar_cols = covar_cols,
      # rh_ref lets predict.fit_sb() check a newdata rh against the fitted scale;
      # rh_range is the fitted humidity extent; see rh_range().
      rh_ref = rh_ref,
      rh_range = rh_range(df, rh_col),
      temp_range = temp_range(df, temp_var),
      rate_cols = rate_cols,
      # Recorded so predict.fit_sb() can enforce the rate-support guard.
      rate_cells = rate_cells(df, rate_cols),
      unsupported_rate_cells = unsupported_rate_cells(df, rate_cols, time_var,
                                                      intercept_cols, temp_var),
      intercept_cols = intercept_cols),
    predict_run  = predict_run,
    data         = df_model,
    sort_order   = sort_order
  )
}


#' Fit the alpha-space ODE kinetic model (LM algorithm with bounds)
#'
#' Wrapper that builds the ODE residual function, prepares starts/bounds,
#' and fits via \code{minpack.lm::nls.lm()}.
#'
#' @section Intercept coding differs from the closed-form fitters:
#' Grouping here is by \code{interaction()}, so both intercepts and rates are
#' **full cell means**: one parameter per observed combination of the columns,
#' named with a \code{"||"} separator (e.g. \code{y0_B1||inv}). The closed-form
#' fitters (\code{\link{fit_kinetic}}, \code{\link{fit_singledeg}},
#' \code{\link{fit_parallel}}) instead reference-code their intercepts via
#' \code{\link{expr_intercept}}.
#'
#' With two or more \code{intercept_cols}, an identical call therefore fits a
#' different model here than there. On a 2x2 design over 32 rows this function
#' fits 11 parameters (21 residual df) against \code{fit_kinetic}'s 8 (24
#' residual df). Cell means can represent an interaction between the grouping
#' columns and reference coding cannot, so the two are genuinely different
#' specifications, not two namings of one model. **Do not put fits from the two
#' paths into the same AIC/BIC/RMSE comparison**: the parameter-count gap
#' favours the closed-form entry for that reason alone. The fitters warn when
#' two or more \code{intercept_cols} are supplied. To make the paths comparable,
#' cross the columns yourself and pass one column:
#' \code{df$cell <- interaction(df$batch, df$orient, drop = TRUE)}.
#'
#' @section Time must be non-negative:
#' This model is an initial-value problem: the state is integrated from
#' \code{alpha(0) = 0}, so \code{y0} is by construction the value at time 0 and
#' the time axis must start there. **Negative times are rejected**, in both
#' \code{df} and the \code{newdata} of \code{\link{predict.fit_sb}}. If your
#' times are referenced to something other than the first measurement (a
#' retained sample read before the release test, or an axis anchored to a
#' stability start date), shift the axis so the earliest reading is time 0,
#' bearing in mind that \code{y0} then refers to that reading rather than to the
#' original origin. The closed-form fitters (\code{\link{fit_kinetic}},
#' \code{\link{fit_singledeg}}, \code{\link{fit_parallel}}) evaluate an algebraic
#' expression in \code{t} and accept negative times as they are.
#'
#' @param df Data frame with outcome, time, temperature, and optionally covariates.
#' @param intercept_cols,rate_cols Character vectors of grouping columns for
#'   intercept and rate effects. Both are coded as full cell means; see the
#'   section above, as this differs from the closed-form fitters.
#' @param outcome_var,temp_var,time_var Column names for outcome, temperature
#'   (Kelvin) and time. \code{temp_var} ignored when \code{temp_expr} is provided.
#' @param R Gas constant (default  8.314 J K-1 mol-1). Ignored when \code{temp_expr} is provided.
#' @param trefk Reference temperature (Kelvin) used in Arrhenius. Ignored when \code{temp_expr} is provided.
#' @param log_Ea Logical; use log(Ea) parameterization (parameter name "Ealog"). Ignored if \code{temp_expr} is provided.
#' @param trend \code{"increasing"} or \code{"decreasing"} outcome trajectory.
#' @param Yinf Asymptote the trajectory travels toward. \code{NULL} (default)
#'   takes 100 when \code{trend = "increasing"} and 0 when
#'   \code{trend = "decreasing"}, the asymptotes of a percentage-scale CQA. The
#'   value used, defaulted or supplied, must lie beyond every observation in the
#'   declared direction: above \code{max(y)} for an increasing trend, below
#'   \code{min(y)} for a decreasing one, since \code{Y = Y0 + (Yinf - Y0) * alpha}
#'   approaches \code{Yinf} without crossing it. Readings that graze the asymptote
#'   within a small tolerance are treated as noise rather than violations.
#' @param temp_expr Optional list output from \code{expr_arrhenius()} or \code{expr_arrhenius_humidity()}. If `NULL` (default) \code{expr_arrhenius()} is used.
#' @param start Optional named numeric vector of starting values. If NULL,
#'   defaults based on the temperature expression are used.
#' @param lower,upper Optional named numeric vectors for parameter bounds. If NULL,
#'   defaults based on the temperature expression are used; the y0 (intercept)
#'   bounds are derived from the observed outcome range.
#' @param control \code{nls.lm.control()} list (defaults to diagonal scaling across parameters).
#'
#' @return the nls.lm fit object
#' @examples
#' data(evers_supp)
#' dat <- subset(evers_supp, cqa == "Purity" & !is.na(y))
#' fit <- fit_sb(dat, trend = "decreasing")
#' summary(fit)
#'
#' \donttest{
#' fit_batch <- fit_sb(
#'   dat,
#'   rate_cols = "batch",
#'   intercept_cols = "batch",
#'   trend = "decreasing"
#' )
#' summary(fit_batch)
#' }
#'
#' @export
fit_sb <- function(
    df,
    trend          = c("increasing", "decreasing"),
    Yinf           = NULL,
    intercept_cols = character(),
    rate_cols      = character(),
    log_Ea         = FALSE,
    outcome_var    = "y",
    temp_var       = "tempk",
    time_var       = "time",
    R              = 8.314,
    trefk          = 298.15,
    temp_expr      = NULL,
    start          = NULL,
    lower          = NULL,
    upper          = NULL,
    control        = NULL
) {
  
  trend <- match.arg(trend)

  # ── 1. Expr params  ───────────────────────────────
  te <- resolve_temp_expr(temp_expr, R = R, trefk = trefk, temp_var = temp_var, log_Ea = log_Ea)
  temp_expr <- te$temp_expr
  R         <- te$R
  trefk     <- te$trefk
  temp_var  <- te$temp_var
  log_Ea    <- te$log_Ea
  rh_col    <- te$rh_col
  rh_ref    <- te$rh_ref
  # Every declared data reference besides temp_var, rh_col included. A custom
  # temp_expr may name a second covariate; without carrying it the tokeniser
  # reads it as a free parameter; see temp_expr_params().
  covar_cols <- te$covar_cols
  temp_fml  <- te$temp_fml

  if (isTRUE(log_Ea)) {
    log_Ea_names <- c("Ealog", "Ea_log", "log_Ea", "logEa")
    clash <- intersect(c(temp_var, rh_col, covar_cols), log_Ea_names)
    if (length(clash)) {
      stop(
        sprintf("Column(s) %s share a name with a model parameter (%s). ",
                paste(sQuote(clash), collapse = ", "),
                paste(log_Ea_names, collapse = ", ")),
        "parameters, so the column would silently overwrite one and the fit would ",
        "converge to meaningless estimates. Rename the column.",
        call. = FALSE
      )
    }
    # Replace all possible variants to the canonical name "Ealog"
    temp_fml <- gsub("\\b(Ea_log|log_Ea|logEa)\\b", "Ealog", temp_fml, perl = TRUE)
  }

  # ── 2. Build Equation ───────────────────────────────
  build <- build_alpha_ode(
    df = df,
    rate_cols = rate_cols,
    intercept_cols = intercept_cols,
    outcome_var = outcome_var,
    time_var = time_var,
    temp_var = temp_var,
    R = R,
    trefk = trefk,
    log_Ea = log_Ea,
    trend = trend,
    Yinf = Yinf,
    temp_fml = temp_fml,
    rh_col = rh_col,
    covar_cols = covar_cols,
    rh_ref = rh_ref
  )

  # After build_alpha_ode(), so its hard errors (missing columns, NAs,
  # parameter-name clashes) take precedence over these advisory warnings.
  check_rh_scale(df, rh_col, rh_ref, what = "Model")
  check_rh_identifiable(df, rh_col, temp_fml)
  check_rate_support(df, rate_cols, time_var, intercept_cols, temp_var)
  check_rate_temp_support(df, rate_cols, temp_var)
  check_trend_direction(df, outcome_var, time_var, trend, temp_var)

  # Seed/scale the rate from the observed time span (k*t ~ O(1))
  t_obs  <- build$data[[time_var]]
  tmax   <- max(t_obs, na.rm = TRUE)
  k_seed <- if (is.finite(tmax) && tmax > 0) 1 / tmax else 0.1

  # ── 3.  Generate default starts/bounds if missing ───────────────
  if (is.null(start) || is.null(lower) || is.null(upper)) {
    
    # group levels from build
    rate_lv  <- build$rate_levels
    int_lv   <- build$intercept_levels
    
    # per-group parameter names
    k_names     <- paste0("k_", rate_lv)
    y0_names <- paste0("y0_", int_lv)
    Ea_name <- build$settings$Ea_name

    # y0 defaults derived from the observed outcome range (a fixed [-20, 150]
    # box silently clips any CQA not on a percentage scale).
    y_obs <- build$data[[outcome_var]]
    y_rng <- range(y_obs, na.rm = TRUE)
    y_pad <- max(1, 0.5 * diff(y_rng))
    # Seed y0 from the outcome at the earliest observed time (trend-agnostic):
    # a fixed max(y) mis-seeds increasing-trend fits.
    y0_start <- mean(y_obs[t_obs == min(t_obs, na.rm = TRUE)], na.rm = TRUE)
    y0_lo    <- y_rng[1] - y_pad
    y0_hi    <- y_rng[2] + y_pad

    # Temperature parameter defaults, keyed on the parameters the expression

    Ea_seed      <- default_Ea_seed(R)
    covar_params <- setdiff(build$settings$temp_params, Ea_name)
    temp_start <- setNames(c(Ea_seed, rep(0.1,  length(covar_params))), c(Ea_name, covar_params))
    temp_lower <- setNames(c(0,       rep(-10,  length(covar_params))), c(Ea_name, covar_params))
    temp_upper <- setNames(c(1e6,     rep(10,   length(covar_params))), c(Ea_name, covar_params))

    if (log_Ea) {
      temp_start[Ea_name] <- log(Ea_seed)
      temp_lower[Ea_name] <- -10
      temp_upper[Ea_name] <- log(1e6)
    }

    start_def <- c(
      temp_start, n  = 1, m  = 0.1,
      setNames(rep(k_seed, length(k_names)), k_names),
      setNames(rep(y0_start, length(y0_names)), y0_names)
    )
    
    lower_def <- c(
      temp_lower, n  = 0, m  = 0,
      setNames(rep(0,   length(k_names)), k_names),
      setNames(rep(y0_lo, length(y0_names)), y0_names)
    )
    
    upper_def <- c(
      temp_upper, n  = 20, m  = 10,
      setNames(rep(100, length(k_names)), k_names),
      setNames(rep(y0_hi, length(y0_names)), y0_names)
    )
    
    # only fill the ones the user didn't provide
    if (is.null(start)) start <- start_def
    if (is.null(lower)) lower <- lower_def
    if (is.null(upper)) upper <- upper_def
  }

  # ── 3b.  Align start/lower/upper to the model parameters by name ──────
  pnames <- build$param_names
  align_par <- function(v, what) {
    if (is.null(names(v)) || anyNA(names(v)) || any(names(v) == "")) {
      stop(sprintf("`%s` must be a fully named vector; parameters are: %s",
                   what, paste(pnames, collapse = ", ")), call. = FALSE)
    }
    if (anyDuplicated(names(v))) {
      dup <- unique(names(v)[duplicated(names(v))])
      stop(sprintf("`%s` has duplicated names: %s", what,
                   paste(dup, collapse = ", ")), call. = FALSE)
    }
    if (anyNA(v)) {
      stop(sprintf("`%s` contains NA/NaN for: %s", what,
                   paste(names(v)[is.na(v)], collapse = ", ")), call. = FALSE)
    }
    extra   <- setdiff(names(v), pnames)
    missing <- setdiff(pnames, names(v))
    if (length(extra) || length(missing)) {
      stop(sprintf("`%s` does not match model parameters (%s).%s%s", what,
                   paste(pnames, collapse = ", "),
                   if (length(extra))   paste0(" Unknown: ",  paste(extra,   collapse = ", ")) else "",
                   if (length(missing)) paste0(" Missing: ", paste(missing, collapse = ", ")) else ""),
           call. = FALSE)
    }
    v[pnames]
  }
  start <- align_par(start, "start")
  lower <- align_par(lower, "lower")
  upper <- align_par(upper, "upper")
  if (any(lower > upper)) {
    stop(sprintf("lower > upper for: %s",
                 paste(pnames[lower > upper], collapse = ", ")), call. = FALSE)
  }

  # ── 4.  Column scaling + optimizer control ───────────────────────────
  # Column scaling (diag) conditions the problem: linear Ea spans ~1e4 while other
  # parameters are O(1). On the log scale Ealog is O(10) and must be scaled ~1.
  # Always computed so it survives a user control.
  Ea_name  <- build$settings$Ea_name
  Ea_scale <- if (log_Ea) 1 else 1/1e4
  # Scale every covariate coefficient the expression has, not just one named `b`.
  covar_params <- setdiff(build$settings$temp_params, Ea_name)
  temp_scale <- setNames(c(Ea_scale, rep(1/0.1, length(covar_params))),
                         c(Ea_name, covar_params))
  scale_vec <- c(temp_scale, n = 1, m = 1)
  k_idx   <- startsWith(names(start), "k_")
  scale_k <- 1 / pmax(abs(start[k_idx]), k_seed)
  scale_y0 <- rep(1/10, sum(startsWith(names(start), "y0_")))
  names(scale_y0) <- names(start)[startsWith(names(start), "y0_")]
  scale_vec <- c(scale_vec, scale_k, scale_y0)[names(start)]

  if (is.null(control)) {
    control <- minpack.lm::nls.lm.control(
      maxiter = 500,     # more steps
      maxfev  = 5000,    # more func evals
      ftol    = 1e-12,   # tighter SSR tolerance
      ptol    = 1e-12,   # tighter param tolerance
      gtol    = 1e-12,   # tighter gradient tolerance
      factor  = 5,       # smaller initial damping step
      nprint  = 0,
      diag    = scale_vec
    )
  } else {
    # Keep the package scaling unless the user passed a usable numeric diag.
    # nls.lm.control()$diag defaults to list(), which must not be subset by
    # name (that yields all-NA scaling and aborts inside the ODE).
    if (is.numeric(control$diag) && length(control$diag)) {
      if (is.null(names(control$diag)) || !setequal(names(control$diag), names(start))) {
        stop("`control$diag` must be named to match the model parameters: ",
             paste(names(start), collapse = ", "), call. = FALSE)
      }
      control$diag <- control$diag[names(start)]
    } else {
      control$diag <- scale_vec
    }
  }
  
  # ── 5.  Detect fixed parameters ───────────────────────────────
  # A parameter is "fixed" when lower == upper. Removed from the optimisation so
  # the Hessian is non-singular and summary() / intervals work.
  fixed <- (lower == upper) & is.finite(lower)
  fixed_vals <- lower[fixed]           # pin at the bound, not at start
  start[fixed] <- fixed_vals

  oob <- (start < lower) | (start > upper)
  if (any(oob)) {
    stop(sprintf("start value(s) outside [lower, upper] for: %s",
                 paste(pnames[oob], collapse = ", ")), call. = FALSE)
  }

  free_idx   <- which(!fixed)

  # Guard the degenerate case (lower == upper for all)
  if (length(free_idx) == 0L) {
    stop("No free parameters to estimate: every parameter is pinned via ",
         "lower == upper. Free at least one parameter to fit.", call. = FALSE)
  }

  start_free <- start[free_idx]
  lower_free <- lower[free_idx]
  upper_free <- upper[free_idx]

  # Guard against n < p 
  n_obs <- nrow(build$data)
  if (length(start_free) > n_obs) {
    stop(sprintf(
      "Negative residual degrees of freedom: %d free parameter(s) but only %d observation(s).",
      length(start_free), n_obs), call. = FALSE)
  }
  
  # Adjust control$diag if present (it must match free params only)
  if (!is.null(control$diag)) {
    control$diag <- control$diag[names(start_free)]
  }
  
  # Wrapper: inject fixed values, call full residual function
  resid_free <- function(p_free, data) {
    p_full <- numeric(length(start))
    names(p_full) <- names(start)
    p_full[free_idx] <- p_free
    p_full[fixed]    <- fixed_vals
    build$residuals(p_full, data)
  }
  
  # ── 6.  Estimate model ───────────────────────────────
  fit_raw <- minpack.lm::nls.lm(
    par     = start_free,
    fn      = resid_free,
    data    = build$data,
    lower   = lower_free,
    upper   = upper_free,
    control = control
  )
  
  # Reconstruct full parameter vector
  par_full <- numeric(length(start))
  names(par_full) <- names(start)
  par_full[free_idx] <- fit_raw$par
  par_full[fixed]    <- fixed_vals
  fit_raw$par <- par_full

  # Warn when a free estimate rests on a finite bound (silent clipping). n or m
  # resting at a lower bound of exactly 0 is a meaningful result (first-order /
  # no autocatalysis), not clipping, so it is exempt; any other bound still warns.
  fp  <- par_full[free_idx]
  lo  <- lower[free_idx]
  up  <- upper[free_idx]
  tol <- 1e-6 * pmax(1, abs(up - lo))
  hit_lo <- is.finite(lo) & abs(fp - lo) <= tol
  hit_up <- is.finite(up) & abs(fp - up) <= tol
  exempt <- names(fp) %in% c("n", "m") & lo == 0
  hit <- names(fp)[(hit_lo & !exempt) | hit_up]
  if (length(hit)) {
    warning("Estimate(s) rest on a bound: ", paste(hit, collapse = ", "),
            ". Consider widening `lower`/`upper`.", call. = FALSE)
  }
  
  # Store fixed-parameter metadata
  fit_raw$fixed      <- fixed
  fit_raw$fixed_vals <- fixed_vals
  fit_raw$n_free     <- length(free_idx)

  # A converged fit can still land on a negative Ea; warn
  check_activation_energy_sign(par_full)

  # Attach metadata + class, but return only the model object
  attr(fit_raw, "build") <- build
  class(fit_raw) <- c("fit_sb", class(fit_raw))
  fit_raw
  
}


#' Predict method for fit_sb models
#'
#' Computes predictions, confidence intervals, prediction intervals, or
#' tolerance intervals for objects fitted with \code{\link{fit_sb}}.
#'
#' @param object An object of class "fit_sb".
#' @param newdata Optional data.frame; if NULL, uses training data stored in build.
#' @param interval Interval to calculate. Default is \code{"none"} and the function
#'   will return the predictions as a numeric vector. Other options are
#'   \code{"confidence"} for a confidence interval, \code{"prediction"} for a
#'   prediction interval, or \code{"tolerance"} for a tolerance interval, all
#'   returning a data.frame with \code{fit}, \code{lwr}, and \code{upr} columns.
#' @param method Method for interval computation. \code{"delta"} (default) uses the Delta
#'   Method (first-order Taylor linearization) with Student's t-distribution. \code{"montecarlo"}
#'   propagates parameter uncertainty via Monte Carlo simulation from the
#'   multivariate normal sampling distribution. Tolerance intervals
#'   always use the Delta Method regardless of this setting.
#' @param level Confidence level for interval calculation. Default is 0.95.
#' @param n_sim Number of Monte Carlo simulations for interval estimation
#'   (only used when \code{method = "montecarlo"}). Default is 1000.
#' @param coverage Coverage proportion (\eqn{\gamma}) for tolerance intervals.
#'   The proportion of the population to be covered. Default is 0.90.
#'   Only used when \code{interval = "tolerance"}.
#' @param side Side for tolerance interval: \code{"two"} (default), \code{"lower"},
#'   or \code{"upper"}. Only used when \code{interval = "tolerance"}. Two-sided
#'   intervals use the Wallis (1951) approximation; one-sided intervals use the
#'   noncentral t-distribution.
#' @param rate_override Optional named numeric vector to override k_ per rate group.
#' @param intercept_override Optional named numeric vector to override y0_ per intercept group.
#' @param return_alpha If TRUE, return alpha(t,T); otherwise return Y(t,T).
#' @param ... Additional arguments (currently ignored; present for S3
#'   method consistency with the \code{predict} generic).
#'
#' @details
#' \strong{Delta method} (\code{method = "delta"}): The prediction variance is
#' computed as \eqn{\mathrm{Var}(\hat{y}) = J \, V_\theta \, J^\top} where
#' \eqn{J} is the numerical Jacobian of predictions with respect to the (free)
#' model parameters and \eqn{V_\theta = \hat\sigma^2 H^{-1}} is the parameter
#' covariance matrix derived from the Hessian. For prediction intervals, the
#' residual variance \eqn{\hat\sigma^2} is added. Intervals use Student's t
#' quantiles with \eqn{n - p} degrees of freedom.
#'
#' \strong{Monte Carlo method} (\code{method = "montecarlo"}): Parameters are
#' sampled from their multivariate normal distribution and propagated through
#' the ODE solver. For prediction intervals, residual noise is added to each
#' simulation. Quantiles of the simulated predictions form the interval bounds.
#'
#' \strong{Tolerance intervals} (\code{interval = "tolerance"}): Computed via
#' the Delta Method following Krishnamoorthy & Mathew (2009). See
#' \code{vignette("intervals-calculation", package = "stablr")} for details.
#'
#' @return If interval = "none", a numeric vector of predictions. Otherwise, a
#'   data.frame with columns \code{fit}, \code{lwr}, \code{upr}.
#' @export
predict.fit_sb <- function(object,
                           newdata = NULL,
                           interval = c("none", "confidence", "prediction",
                                        "tolerance"),
                           method = c("delta", "montecarlo"),
                           level = 0.95,
                           n_sim = 1000,
                           coverage = 0.90,
                           side = c("two", "lower", "upper"),
                           rate_override = NULL,
                           intercept_override = NULL,
                           return_alpha = FALSE,
                           ...) {
  
  interval <- match.arg(interval)
  method   <- match.arg(method)
  side     <- match.arg(side)


  if (interval != "none") {
    if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
        level <= 0 || level >= 1) {
      stop("level must be a single number strictly between 0 and 1 (e.g. 0.95).",
           call. = FALSE)
    }
    if (interval == "tolerance" &&
        (!is.numeric(coverage) || length(coverage) != 1L || is.na(coverage) ||
         coverage <= 0 || coverage >= 1)) {
      stop("coverage must be a single number strictly between 0 and 1 (e.g. 0.90).",
           call. = FALSE)
    }
    if (method == "montecarlo" &&
        (!is.numeric(n_sim) || length(n_sim) != 1L || is.na(n_sim) ||
         n_sim < 2 || n_sim != round(n_sim))) {
      stop("n_sim must be a single integer >= 2 for method = \"montecarlo\".",
           call. = FALSE)
    }
  }

  # Residual-variance intervals are defined only on the Y scale so will not be used for alpha predictions
  if (return_alpha && interval %in% c("prediction", "tolerance")) {
    stop("interval = \"", interval, "\" is not available with return_alpha = TRUE: ",
         "residual-variance intervals are only defined on the Y scale. ",
         "Use interval = \"confidence\", or return_alpha = FALSE.", call. = FALSE)
  }

  if (interval == "tolerance" && method == "montecarlo") {
    stop(
      "Tolerance intervals are only available with method = \"delta\". ",
      call. = FALSE
    )
  }
  
  if (!inherits(object, "fit_sb")) stop("object is not of class 'fit_sb'")
  build <- attr(object, "build")
  if (is.null(build)) stop("Missing 'build' metadata on fit object.")
  
  settings <- build$settings
  time_var <- settings$time_var
  temp_var <- settings$temp_var
  Ea_name  <- settings$Ea_name
  temp_params <- settings$temp_params
  Yinf     <- settings$Yinf
  rate_cols      <- settings$rate_cols
  intercept_cols <- settings$intercept_cols
  temp_fml <- settings$temp_fml
  rh_col <- settings$rh_col
  # Every declared data reference besides temp_var. Fits recorded before
  # covar_cols existed have only rh_col, so fall back to it.
  covar_cols <- settings$covar_cols
  if (is.null(covar_cols)) covar_cols <- rh_col

  # coefficients
  coef_est <- if (!is.null(object$par)) object$par else stop("No coefficients found in model.")
  n_hat    <- unname(coef_est["n"])
  m_hat    <- unname(coef_est["m"])

  # Extract temperature parameters
  if (!is.null(temp_fml)) {
    temp_pars_hat <- coef_est[temp_params]
  } else {
    temp_pars_hat <- unname(coef_est[Ea_name])
    names(temp_pars_hat) <- Ea_name
  }
  
  # fitted per-group params
  k_hat     <- get_param(coef_est, "k_")
  y0_hat <- get_param(coef_est, "y0_")
  
  # overrides. 
  if (!is.null(rate_override)) {
    onm <- names(rate_override)
    key <- ifelse(onm %in% names(k_hat), onm, sub("^k_", "", onm))
    hit <- key %in% names(k_hat)
    k_hat[key[hit]] <- rate_override[hit]
  }
  # Accept "y0_<level>" (as coef()/summary()/start= display it) as well as the
  # bare "<level>" used internally; prefer an exact bare-label match.
  if (!is.null(intercept_override)) {
    onm <- names(intercept_override)
    key <- ifelse(onm %in% names(y0_hat), onm, sub("^y0_", "", onm))
    hit <- key %in% names(y0_hat)
    y0_hat[key[hit]] <- intercept_override[hit]
  }
  
  # choose data
  if (is.null(newdata)) {
    inv <- if (is.null(build$sort_order)) {
      seq_len(nrow(build$data))
    } else {
      order(build$sort_order)
    }
    dat <- build$data[inv, , drop = FALSE]
  } else {
    dat <- newdata
  }
  
  # checks
  need_cols <- unique(c(time_var, temp_var, rate_cols, intercept_cols, covar_cols))
  if (length(setdiff(need_cols, names(dat))) > 0) {
    miss <- setdiff(need_cols, names(dat))
    stop("newdata is missing columns: ", paste(miss, collapse = ", "))

  }

  # The ODE is solved per group over that group's requested times, so a negative
  # time here moves the initial condition off t = 0 for the whole group and
  # shifts every other prediction in it, not just its own row.
  check_nonnegative_times(dat, time_var,
                          what = if (is.null(newdata)) "df" else "newdata")

  # Predictor types and humidity scale, the same contract fit_sb() enforced.
  # Only for supplied newdata: build$data already passed at fit time, and
  # re-running would turn a humidity-coding warning into a predict() error.
  if (!is.null(newdata)) {
    check_predictor_cols(newdata, settings, what = "newdata")

    # Applied here too so the newdata contract is identical across all four
    # fitters, though this model is structurally immune to the shadowing.
    check_newdata_param_shadow(newdata, build$param_names)

    check_temp_extrapolation(newdata, settings)

    # A rate cell observed but never identified at fit time (e.g. only at time 0).
    check_rate_support_newdata(newdata, settings$unsupported_rate_cells,
                               rate_cols, time_var)
  }

  # rebuild group ids exactly as in training
  rate_id <- make_group(dat, rate_cols, default_label = "0")
  int_id  <- make_group(dat, intercept_cols, default_label = "0")
  
  # guard against unseen groups unless overridden
  unseen_rate <- setdiff(levels(rate_id), names(k_hat))
  if (length(unseen_rate)) {
    stop("Unseen rate group(s): ",
         paste(unseen_rate, collapse = ", "),
         ". Provide rate_override for these labels.")
  }
  unseen_int <- setdiff(levels(int_id), names(y0_hat))
  if (length(unseen_int)) {
    stop("Unseen intercept group(s): ",
         paste(unseen_int, collapse = ", "),
         ". Provide intercept_override for these labels.")
  }
  
  # efficient ODE solves by (rate, intercept, temp)
  rate_chr <- as.character(rate_id)
  int_chr  <- as.character(int_id)

  # Group on every covariate, not just humidity: two rows differing only in a
  # second covariate would otherwise share one solve.
  idx_split <- split(
    seq_len(nrow(dat)),
    interaction(c(list(rate_chr, int_chr, dat[[temp_var]]),
                  lapply(covar_cols, function(cc) dat[[cc]])),
                drop = TRUE, sep = "||")
  )
  
  # ── Point predictions ───────────────────────────────
  pred_alpha <- rep(NA_real_, nrow(dat))
  pred_Y <- rep(NA_real_, nrow(dat))
  
  for (idx in idx_split) {
    
    if (!length(idx)) next
    
    rate_lbl <- rate_chr[idx][1]
    int_lbl <- int_chr[idx][1]
    
    k_g   <- unname(k_hat[rate_lbl])
    y0_g  <- unname(y0_hat[int_lbl])
    
    T_i  <- dat[[temp_var]][idx][1]
    t_i  <- dat[[time_var]][idx]

    pars_i <- c(k = k_g, temp_pars_hat, n = n_hat, m = m_hat)

    covar_val <- if (length(covar_cols)) {
      stats::setNames(vapply(covar_cols,
                             function(cc) as.numeric(dat[[cc]][idx][1]),
                             numeric(1)), covar_cols)
    } else NULL

    a_i <- build$predict_run(time_vec = t_i, tempk = T_i, pars = pars_i, covars = covar_val)
    if (any(!is.finite(a_i)) || any(a_i < -1e-6) || any(a_i > 1 + 1e-6)) {
      stop("Non-finite or out-of-bounds alpha predicted; check parameters or data.")
    }
    pred_alpha[idx] <- a_i
    pred_Y[idx] <- alpha_to_Y(a_i, Y0 = y0_g, Yinf = Yinf)
    
  }
  
  fit_vals <- if (return_alpha) pred_alpha else pred_Y
  
  # ── Return point predictions if no interval requested ───────────────────────────────
  if (interval == "none") {
    return(fit_vals)
  }
  
  # ── Delta method interval estimation ───────────────────────────────
  if (method == "delta") {
    
    p <- length(coef_est)
    
    # Handle fixed parameters
    fixed      <- if (!is.null(object$fixed)) object$fixed else rep(FALSE, p)
    free_idx   <- which(!fixed)
    n_free     <- length(free_idx)
    
    # Residual variance
    resid    <- object$fvec
    n_obs    <- length(resid)
    df_resid <- n_obs - n_free

    # df_resid <= 0; 
    if (df_resid <= 0) {
      stop(sprintf(
        "Cannot compute intervals: residual degrees of freedom is %d (%d observation(s), %d free parameter(s)). Refit with more data or fewer parameters.",
        df_resid, n_obs, n_free), call. = FALSE)
    }
    sigma2   <- sum(resid^2) / df_resid
    
    # Parameter covariance matrix (free params only; Hessian is n_free x n_free)
    hess <- object$hessian
    if (is.null(hess)) {
      stop("Hessian not available in fit object. Cannot compute intervals.")
    }
    
    vcov_free <- sigma2 * invert_hessian(hess)
    
    # Numerical Jacobian of predictions w.r.t. free parameters
    n_pred <- length(fit_vals)
    grad_matrix <- matrix(NA_real_, nrow = n_pred, ncol = n_free)
    
    for (j in seq_along(free_idx)) {
      par_plus <- coef_est
      # Relative finite-difference step. The ODE is solved to rtol = 1e-8, so the
      # reliable forward-difference step is sqrt(1e-8) = 1e-4 relative.
      theta_j <- abs(par_plus[free_idx[j]])
      eps_j <- if (theta_j > 0) theta_j * 1e-4 else 1e-4
      par_plus[free_idx[j]] <- par_plus[free_idx[j]] + eps_j
      
      obj_temp <- object
      obj_temp$par <- par_plus
      
      pred_plus <- predict(obj_temp, newdata = newdata,
                           rate_override = rate_override,
                           intercept_override = intercept_override,
                           return_alpha = return_alpha)
      grad_matrix[, j] <- (pred_plus - fit_vals) / eps_j
    }
    
    # ── Tolerance intervals (Krishnamoorthy & Mathew, 2009) ──────────
    if (interval == "tolerance") {
      
      # d^2 = J * (F'F)^{-1} * J'  where (F'F)^{-1} = H^{-1}
      Hinv <- vcov_free / sigma2
      d_squared <- rowSums((grad_matrix %*% Hinv) * grad_matrix)
      
      return(tolerance_bounds_km(
        mu_hat    = fit_vals,
        d_squared = d_squared,
        sigma     = sqrt(sigma2),
        df        = df_resid,
        coverage  = coverage,
        level     = level,
        side      = side
      ))
    }
    
    # ── Confidence / Prediction intervals ─────────────────────────────
    # Variance of predictions: diag(J V J')
    var_pred <- rowSums((grad_matrix %*% vcov_free) * grad_matrix)
    
    # For prediction intervals, add residual variance
    if (interval == "prediction") {
      var_pred <- var_pred + sigma2
    }
    
    se_pred <- sqrt(pmax(var_pred, 0))
    
    # Student's t quantile with residual degrees of freedom
    q_mult <- stats::qt(1 - (1 - level) / 2, df = df_resid)
    
    result <- data.frame(
      fit = fit_vals,
      lwr = fit_vals - q_mult * se_pred,
      upr = fit_vals + q_mult * se_pred
    )
    
    return(result)
  }
  
  # ── Monte Carlo interval estimation ───────────────────────────────
  p <- length(coef_est)
  
  # Handle fixed parameters
  fixed      <- if (!is.null(object$fixed)) object$fixed else rep(FALSE, p)
  fixed_vals <- if (!is.null(object$fixed_vals)) object$fixed_vals else numeric(0)
  free_idx   <- which(!fixed)
  n_free     <- length(free_idx)
  
  # Compute parameter covariance matrix from the Hessian (free params only)
  hess <- object$hessian
  if (is.null(hess)) {
    stop("Hessian not available in fit object. Cannot compute intervals.")
  }
  
  # Residual variance (df uses only free parameters)
  resid <- object$fvec
  n_obs <- length(resid)
  df_resid <- n_obs - n_free
  # sigma2 is undefined at df_resid <= 0; error rather than propagate NaN draws.
  if (df_resid <= 0) {
    stop(sprintf(
      "Cannot compute intervals: residual degrees of freedom is %d (%d observation(s), %d free parameter(s)). Refit with more data or fewer parameters.",
      df_resid, n_obs, n_free), call. = FALSE)
  }
  sigma2 <- sum(resid^2) / df_resid
  sigma <- sqrt(sigma2)
  
  # Unit-sigma parameter covariance H^{-1} (free params only). The per-draw
  # covariance is sigma_s^2 * H^{-1}, with sigma_s^2 drawn below.
  vcov_unit <- invert_hessian(hess)
  
  # Ensure positive definiteness
  eig <- eigen(vcov_unit, symmetric = TRUE)
  if (any(eig$values <= 0)) {
    min_eig <- min(eig$values)
    vcov_unit <- vcov_unit + diag(abs(min_eig) + 1e-6, n_free)
  }
  
  # Propagate the uncertainty in sigma_hat
  # Draw sigma_s^2 ~ df * sigma2_hat / chisq(df) per simulation 
  coef_free <- coef_est[free_idx]
  base_draw <- MASS::mvrnorm(n = n_sim, mu = rep(0, n_free), Sigma = vcov_unit)
  if (n_sim == 1L) base_draw <- matrix(base_draw, nrow = 1L)
  sigma2_s  <- df_resid * sigma2 / stats::rchisq(n_sim, df = df_resid)
  sigma_s   <- sqrt(sigma2_s)
  par_samples_free <- sweep(base_draw * sigma_s, 2, coef_free, `+`)
  
  # Expand to full parameter matrix
  par_samples <- matrix(NA_real_, nrow = n_sim, ncol = p)
  colnames(par_samples) <- names(coef_est)
  for (s in seq_len(n_sim)) {
    par_samples[s, free_idx] <- par_samples_free[s, ]
    par_samples[s, fixed]    <- fixed_vals
  }
  
  n_pred <- nrow(dat)
  pred_matrix <- matrix(NA_real_, nrow = n_pred, ncol = n_sim)
  n_valid <- 0L    # count simulations that produced a usable trajectory
  
  for (s in seq_len(n_sim)) {
    
    par_s <- par_samples[s, ]
    
    # Extract parameters for this simulation
    n_s  <- unname(par_s["n"])
    m_s  <- unname(par_s["m"])

    # Extract temperature parameters
    if (!is.null(temp_fml)) {
      temp_pars_s <- par_s[temp_params]
    } else {
      temp_pars_s <- unname(par_s[Ea_name])
      names(temp_pars_s) <- Ea_name
    }
    
    k_s     <- get_param(par_s, "k_")
    y0_s <- get_param(par_s, "y0_")
    
    # Apply overrides (accept "k_<label>" as well as the bare "<label>")
    if (!is.null(rate_override)) {
      onm_ov <- names(rate_override)
      key_ov <- ifelse(onm_ov %in% names(k_s), onm_ov, sub("^k_", "", onm_ov))
      hit_ov <- key_ov %in% names(k_s)
      k_s[key_ov[hit_ov]] <- rate_override[hit_ov]
    }
    if (!is.null(intercept_override)) {
      onm_ov <- names(intercept_override)
      key_ov <- ifelse(onm_ov %in% names(y0_s), onm_ov, sub("^y0_", "", onm_ov))
      hit_ov <- key_ov %in% names(y0_s)
      y0_s[key_ov[hit_ov]] <- intercept_override[hit_ov]
    }
    
    pred_s <- rep(NA_real_, n_pred)
    valid_sim <- TRUE
    
    for (idx in idx_split) {
      
      if (!length(idx)) next
      
      rate_lbl <- rate_chr[idx][1]
      int_lbl  <- int_chr[idx][1]
      
      k_g     <- unname(k_s[rate_lbl])
      y0_g <- unname(y0_s[int_lbl])
      
      T_i <- dat[[temp_var]][idx][1]
      t_i <- dat[[time_var]][idx]

      # Build parameter vector with appropriate temperature parameters
      pars_i <- c(k = k_g, temp_pars_s, n = n_s, m = m_s)

      covar_val <- if (length(covar_cols)) {
        stats::setNames(vapply(covar_cols,
                               function(cc) as.numeric(dat[[cc]][idx][1]),
                               numeric(1)), covar_cols)
      } else NULL

      a_i <- tryCatch({
        build$predict_run(time_vec = t_i, tempk = T_i, pars = pars_i, covars = covar_val)
      }, error = function(e) {
        rep(NA_real_, length(t_i))
      })
      
      if (any(!is.finite(a_i)) || any(a_i < -1e-6, na.rm = TRUE) || any(a_i > 1 + 1e-6, na.rm = TRUE)) {
        valid_sim <- FALSE
        break
      }
      
      if (return_alpha) {
        pred_s[idx] <- a_i
      } else {
        pred_s[idx] <- alpha_to_Y(a_i, Y0 = y0_g, Yinf = Yinf)
      }
    }
    
    if (!valid_sim) next
    
    # Add residual variance for prediction intervals
    if (interval == "prediction") {
      pred_s <- pred_s + stats::rnorm(n_pred, mean = 0, sd = sigma_s[s])
    }
    
    pred_matrix[, s] <- pred_s
    n_valid <- n_valid + 1L
  }
  
  # A discarded simulation contributes no column, and the quantiles below use
  # na.rm = TRUE, so a mostly-failed run would return a too-narrow interval.
  # Surface the discard rate: error if none remain, warn if a material fraction
  # was dropped.
  n_fail <- n_sim - n_valid
  if (n_valid == 0L) {
    stop("All ", n_sim, " Monte Carlo simulations were discarded (out-of-range ",
         "trajectories or ODE failures). Cannot form an interval; the fit is ",
         "too poorly constrained for method = \"montecarlo\". Try ",
         "method = \"delta\".", call. = FALSE)
  }
  if (n_fail > 0.05 * n_sim) {
    warning(sprintf(
      "%d of %d Monte Carlo simulations (%.0f%%) were discarded as invalid; the interval is based on %d draw(s) and may be unreliable. Consider more data, tighter bounds, or method = \"delta\".",
      n_fail, n_sim, 100 * n_fail / n_sim, n_valid), call. = FALSE)
  }
  
  # Compute quantiles
  alpha_q <- (1 - level) / 2
  probs <- c(alpha_q, 0.5, 1 - alpha_q)
  
  quantiles <- t(apply(pred_matrix, 1, stats::quantile, probs = probs, na.rm = TRUE))
  
  result <- data.frame(
    fit = fit_vals,
    lwr = quantiles[, 1],
    upr = quantiles[, 3]
  )
  
  result
}


#' Summary method for fit_sb models
#'
#' Temporarily restores the free-parameter-only \code{par} vector so that
#' the inherited \code{summary.nls.lm()} method sees consistent dimensions
#' between \code{par} and the Hessian.
#' Fixed parameters are printed separately.
#'
#' @param object An object of class "fit_sb".
#' @param ... Additional arguments passed to \code{summary.nls.lm()}.
#' @return The result of \code{summary.nls.lm()} (invisibly).
#' @export
summary.fit_sb <- function(object, ...) {
  
  fixed <- object$fixed
  
  if (!is.null(fixed) && any(fixed)) {
    # Print fixed params header
    cat("Fixed parameters:\n")
    for (nm in names(object$fixed_vals)) {
      cat(sprintf("  %s = %g\n", nm, object$fixed_vals[nm]))
    }
    cat("\n")
    
    # Temporarily set par to free-only so summary.nls.lm dimensions match
    free_idx <- which(!fixed)
    object$par <- object$par[free_idx]
  }
  
  # Delegate to summary.nls.lm
  class(object) <- setdiff(class(object), "fit_sb")
  tryCatch(
    summary(object, ...),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("leading minor of order 1 is not positive", msg, fixed = TRUE)) {
        stop(
          paste0(
            msg,
            "\nThe Hessian is not positive definite, so the uncertainty ",
            "estimates in summary() are unstable.",
            "Try widening `lower`/`upper`, supplying more informative ",
            "`start` values, or simplifying the grouping."
          ),
          call. = FALSE
        )
      }
      stop(e)
    }
  )
}


#' Variance-covariance matrix for fit_sb models
#'
#' Restores the free-parameter-only \code{par} vector before delegating to the
#' inherited \code{vcov.nls.lm()}, mirroring \code{\link{summary.fit_sb}}.
#'
#' @details
#' \code{fit_sb} refills \code{$par} to full length (including the fixed
#' parameters) while the Hessian stays \code{n_free x n_free}. The inherited
#' \code{vcov.nls.lm()} scales by \code{deviance / (length(fvec) - length(par))},
#' so the full-length \code{par} undercounts the residual degrees of freedom and
#' inflates every standard error by \code{sqrt((n - n_free) / (n - p_full))},
#' silently disagreeing with \code{summary()}, which already corrects for
#' this. Subsetting \code{par} to the free entries makes the covariance match
#' \code{summary()} and \code{predict(interval=)}, all using \code{n - n_free}.
#'
#' @param object An object of class "fit_sb".
#' @param ... Additional arguments passed to \code{vcov.nls.lm()}.
#' @return The free-parameter variance-covariance matrix.
#' @export
vcov.fit_sb <- function(object, ...) {
  fixed <- object$fixed
  if (!is.null(fixed) && any(fixed)) {
    object$par <- object$par[which(!fixed)]
  }
  class(object) <- setdiff(class(object), "fit_sb")
  stats::vcov(object, ...)
}


#' Residual degrees of freedom for fit_sb models
#'
#' Reports \code{n - n_free}, the residual degrees of freedom \code{fit_sb}
#' actually estimated with, rather than the inherited \code{df.residual.nls.lm}
#' value of \code{length(resid) - length(coef)}, which counts the fixed
#' parameters as estimated and so disagrees with \code{summary()} and
#' \code{predict.fit_sb()}.
#'
#' @param object An object of class "fit_sb".
#' @param ... Ignored.
#' @return Integer residual degrees of freedom.
#' @export
df.residual.fit_sb <- function(object, ...) {
  length(object$fvec) - object$n_free
}


#' Residuals for fit_sb models
#'
#' @param object An object of class "fit_sb".
#' @param ... Ignored.
#' @return Numeric vector of residuals (\code{observed - fitted}).
#' @export
residuals.fit_sb <- function(object, ...) {
  build    <- attr(object, "build")
  if (is.null(build)) stop("Missing 'build' metadata on fit object.")
  inv <- if (is.null(build$sort_order)) {
    seq_len(nrow(build$data))
  } else {
    order(build$sort_order)
  }
  observed <- build$data[[build$settings$outcome_var]][inv]
  observed - as.numeric(predict(object))
}


#' Confidence intervals for fit_sb parameters
#'
#' Wald confidence intervals that use the same residual degrees of freedom the
#' fit estimated with (\code{n - n_free}).
#'
#' @details
#' Without a \code{fit_sb} method, \code{confint()} dispatches to
#' \code{stats::confint.default}, which uses a normal quantile
#' (\code{qnorm}) regardless of sample size and returns \code{NA} for the pinned
#' parameters. On the small-n designs typical of predictive stability the normal
#' quantile understates the interval (e.g. \code{qnorm(0.975) / qt(0.975, 3)}
#' is about 0.62, a 38\% too-narrow interval). This method uses Student's
#' \code{t} on \code{n - n_free} degrees of freedom (matching
#' \code{summary()}, \code{vcov()} and \code{predict(interval=)}), and returns
#' each pinned parameter as a degenerate interval at its held value (zero
#' sampling variance) rather than \code{NA}.
#'
#' @param object An object of class "fit_sb".
#' @param parm Optional vector of parameter names or indices to report.
#' @param level Confidence level. Default 0.95.
#' @param ... Ignored.
#' @return A matrix with a row per parameter and lower/upper columns.
#' @export
confint.fit_sb <- function(object, parm, level = 0.95, ...) {
  cf     <- object$par
  pnames <- names(cf)
  fixed  <- if (!is.null(object$fixed)) object$fixed else rep(FALSE, length(cf))

  # Standard errors from the free-parameter covariance (vcov.fit_sb); pinned
  # parameters are held constants and so carry zero sampling variance.
  V  <- stats::vcov(object)
  se <- setNames(rep(0, length(cf)), pnames)
  se_free <- sqrt(diag(V))
  se[rownames(V)] <- se_free

  # Student's t on the residual df the fit actually used.
  a   <- (1 - level) / 2
  fac <- stats::qt(c(a, 1 - a), df = df.residual(object))

  ci <- cf + outer(se, fac)
  colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE), "%")
  rownames(ci) <- pnames

  if (!missing(parm)) {
    if (is.numeric(parm)) parm <- pnames[parm]
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}


#' Fitted values for fit_sb models
#'
#' @param object An object of class "fit_sb".
#' @param ... Ignored.
#' @return Numeric vector of fitted values, one per observation.
#' @export
fitted.fit_sb <- function(object, ...) {
  as.numeric(predict(object))
}


#' Number of observations for fit_sb models
#'
#' @param object An object of class "fit_sb".
#' @param ... Ignored.
#' @return Integer number of observations.
#' @export
nobs.fit_sb <- function(object, ...) {
  length(object$fvec)
}


#' Log-likelihood for fit_sb models
#'
#' @details
#' Gaussian least-squares log-likelihood, identical in form to
#' \code{stats:::logLik.nls} (which the three \code{gsl_nls} fitters inherit):
#' \eqn{\ell = -\frac{n}{2}\left(\log 2\pi + 1 + \log(\mathrm{RSS}/n)\right)}. 
#'
#' @param object An object of class "fit_sb".
#' @param ... Ignored.
#' @return An object of class "logLik" with \code{"df"} and \code{"nobs"} attributes.
#' @export
logLik.fit_sb <- function(object, ...) {
  n   <- length(object$fvec)
  rss <- sum(object$fvec^2)
  val <- -0.5 * n * (log(2 * pi) + 1 - log(n) + log(rss))
  attr(val, "df")   <- object$n_free + 1L        # +1 for the variance
  attr(val, "nobs") <- n
  class(val) <- "logLik"
  val
}


#' Residual standard error for fit_sb models
#'
#'
#' @param object An object of class "fit_sb".
#' @param ... Ignored.
#' @return Numeric residual standard error.
#' @export
sigma.fit_sb <- function(object, ...) {
  sqrt(sum(object$fvec^2) / df.residual(object))
}

