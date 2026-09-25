#' Assemble a Vector‑Valued Residual Function for an A → B Kinetic System
#'
#' One‑stop builder that combines intercept, rate, Arrhenius, and kinetic‑order
#' pieces for a two‑response single‑degradation model (A → B).
#' The result is ready for \code{\link[gslnls]{gsl_nls}}.
#'
#' @param df Data frame containing every variable referenced in the formula.
#' @param kinetic_order Integer reaction order (0, 1, 2, …).
#' @param trend Character string, `"increasing"` or `"decreasing"` outcome
#'   trajectory for \code{yA}.
#' @param Yinf Upper asymptote when \code{trend == "increasing"}.
#' @param intercept_cols,rate_cols Character vectors of grouping columns
#'   governing intercept and rate effects.
#' @param log_Ea Logical; recorded in \code{settings} only. The activation
#'   energy parameter name is taken from \code{temp_fml}, so build
#'   \code{temp_fml} with \code{expr_arrhenius(log_Ea = TRUE)} to obtain
#'   \code{Ea_log}.
#' @param yA_var,yB_var Column names for the primary response A and its
#'   companion response B (mass balance).
#' @param temp_var,time_var Column names for temperature (Kelvin) and time. 
#' @param R Gas constant. Defaults to 8.314 J/(mol * K).
#' @param trefk Reference temperature in Kelvin. 
#' @param temp_fml Character string from \code{expr_arrhenius()$formula} or \code{expr_arrhenius_humidity()$formula}.
#' @param rh_col Column name for relative humidity values (only needed when expr_arrhenius_humidity is used).
#' @param rh_ref Reference relative humidity value (only needed when expr_arrhenius_humidity is used).
#' @param covar_cols Character vector of every data column the temperature
#'   expression is declared to reference besides \code{temp_var} (\code{rh_col}
#'   included). Assembled by \code{resolve_temp_expr()}; symbols not named here are
#'   treated as free parameters.
#' @return A list with
#' * \code{fn} – function \code{f(params, data)} returning the stacked
#'              residual vector \eqn{(y_A, y_B)}.
#' * \code{param_names} – character vector of all parameter names.
#' * \code{response} – numeric vector \eqn{(y_A, y_B)} in the order expected by
#'                     \code{gslnls}.
#' * \code{settings} - list containing build settings.
#' 
#' @examples
#' data("nrcesds")
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' library(tidyverse)
#' nrcesds_wide <- nrcesds %>% tidyr::pivot_wider(names_from = "cqa", values_from = "y")
#' build <- build_singledeg(
#'           df = nrcesds_wide,
#'           kinetic_order = 2,
#'           trend = "decreasing",
#'           yA_var = "NR Major Peak",
#'           yB_var = "NR Minor Peaks",
#'           temp_fml = expr_arrhenius()$formula
#' )
#' str(build$param_names)
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' # This example uses the |> pipe and requires R version 4.1.0 or higher.
#' if (getRversion() >= "4.1.0") {
#'   library(tidyverse)
#'   build <- nrcesds |>
#'     tidyr::pivot_wider(names_from = "cqa", values_from = "y") |>
#'     dplyr::rename(
#'       yA = "NR Major Peak", yB = "NR Minor Peaks"
#'     ) |>
#'     build_singledeg(
#'       kinetic_order  = 2,
#'       trend          = "decreasing",
#'       intercept_cols = "batch",
#'       rate_cols      = character(),
#'       temp_fml = expr_arrhenius()$formula
#'     )
#'   str(build$param_names)
#' }
#' @export
build_singledeg <- function(
    df,
    kinetic_order = 1L,
    trend = c("increasing", "decreasing"),
    Yinf = 100,
    intercept_cols = character(),
    rate_cols     = character(),
    log_Ea = F,
    yA_var   = "yA",
    yB_var   = "yB",
    temp_var       = "tempk",
    time_var       = "time",
    R              = 8.314,
    trefk       = 298.15,
    rh_col = NULL,
    rh_ref = NULL,
    covar_cols = NULL,
    temp_fml
) {
  
  check_kinetic_order(kinetic_order)
  check_scalar_Yinf(Yinf)

  # ── 1. Build the two intercept expressions ───────────────────────────────
  ## Prefix each family up front (not renamed afterwards) so parameter names
  ## don't clash: a rename would also rewrite column names and level values
  ## containing "y0".
  ie_A <- expr_intercept(df, intercept_cols, prefix = "A0")
  ie_B <- expr_intercept(df, intercept_cols, prefix = "B0")

  ie_A$expr <- quote_numeric_levels(ie_A$expr)
  ie_B$expr <- quote_numeric_levels(ie_B$expr)
  
  
  ## ── 2. Rate‑constant & Arrhenius Equation ─────────────────────────────────
  rate_part <- expr_rate(df, rate_cols)                      # k(T, …)
  rate_part$expr <- quote_numeric_levels(rate_part$expr)                   # k(T, …)
  arr_expr <- temp_fml

  ## ── 3. Vector-valued residual function ─────────────────────────────────
  singledeg_fn <- function(params, data) {
    
    penv <- list2env(as.list(params), parent = baseenv())
    
    A_int  <- rlang::eval_tidy(rlang::parse_expr(ie_A$expr), data = data, env = penv)
    B_int  <- rlang::eval_tidy(rlang::parse_expr(ie_B$expr), data = data, env = penv)
    
    rate_part_eval <- rlang::eval_tidy(rlang::parse_expr(rate_part$expr),  data = data, env = penv)
    arr_expr_eval  <- rlang::eval_tidy(rlang::parse_expr(arr_expr), data = data, env = penv)
    rate <- rate_part_eval * arr_expr_eval
    
    tempk <- data[[temp_var]]
    time  <- data[[time_var]]
    rate_time <- rate * time
    
    ## ── 3a. General solution for yA ─────────────────────────────────
    if (kinetic_order == 0) {    
      if (trend == "increasing") {
        yA_hat <- A_int + rate_time
      }
      else{
        yA_hat <- A_int - rate_time 
      }
    } else if (kinetic_order == 1) {                              # first‑order
      if (trend == "increasing") {
        yA_hat <- Yinf - (Yinf - A_int) * exp(-rate_time)
      } else {
        yA_hat <- (A_int) * exp(-rate_time)
      }
    } else {                                          # higher‑order
      pow  <- 1 - kinetic_order
      if (trend == "increasing"){
        yA_hat <- Yinf - ( (Yinf - A_int)^(pow) - (rate_time) * (pow) )^(1 / pow)
      } else{
        yA_hat <- ( (A_int) ^ (pow) - (rate_time) * (pow) )^(1 / pow)
      }
    }
    
    ## ── 3b. Companion response (mass balance) ─────────────────────────────
    yB_hat <- B_int + (A_int - yA_hat)    
    ## gsl_nls wants a *vector* of predictions in the same order as y
    c(yA_hat, yB_hat)
  }
  
  ## ── 4. Parameter names ─────────────────────────────────
  # Tokenise rather than grep, so renamed or extra parameters still get seeded.
  temp_params <- temp_expr_params(arr_expr, temp_var, unique(c(rh_col, covar_cols)))
  par_names <- c(rate_part$params, temp_params, ie_A$params, ie_B$params)

  # Cross-block duplicates: the per-block checks cannot see across blocks.
  check_param_union(list("the rate term"              = rate_part$params,
                         "the temperature expression" = temp_params,
                         "the A intercept term"       = ie_A$params,
                         "the B intercept term"       = ie_B$params))

  # eval_tidy() puts the data mask first, so a column named after a parameter
  # would shadow it outright.
  check_param_data_collisions(par_names, df)

  ## ── 4. Build object ─────────────────────────────────
  list(fn          = singledeg_fn,
       param_names = par_names,
       response    = c(df[[yA_var]], df[[yB_var]]),
       settings     = list(
         R = R, trefk = trefk, log_Ea = log_Ea,
         yA_var = yA_var, yB_var=yB_var, time_var = time_var, temp_var = temp_var,
         Yinf = Yinf,
         rate_cols = rate_cols,
         rate_cells = rate_cells(df, rate_cols),
         unsupported_rate_cells = unsupported_rate_cells(df, rate_cols, time_var,
                                                         intercept_cols, temp_var),
         intercept_cols = intercept_cols,
         rh_col = rh_col,
         covar_cols = covar_cols,
         rh_ref = rh_ref,
         rh_range = rh_range(df, rh_col),
         temp_range = temp_range(df, temp_var)))
}

#' Fit a Two‑Response Single‑Degradation Kinetic Model with \code{gsl_nls}
#'
#' Convenience wrapper that
#' 1. builds the vector‑valued residual function for a coupled
#'    single‑degradation system (A → B) using
#'    \code{\link{build_singledeg}},
#' 2. creates default starting values when none are supplied, and
#' 3. fits the model with \link[gslnls]{gsl_nls}.
#'
#' @param df Data frame containing every variable referenced in the formula.
#' @param kinetic_order Integer reaction order (0, 1, 2, …).
#' @param trend Character string, `"increasing"` or `"decreasing"` outcome
#'   trajectory for \code{yA}.
#' @param Yinf Upper asymptote when \code{trend == "increasing"}.
#' @param intercept_cols,rate_cols Character vectors of grouping columns
#'   governing intercept and rate effects. Intercepts are reference coded
#'   (additive) when two or more columns are given (see
#'   \code{\link{expr_intercept}}); \code{\link{fit_sb}} codes the same call as
#'   full cell means, so the two fit different intercept models with different
#'   parameter counts. Rates are always cell means.
#' @param log_Ea Logical; if \code{TRUE} the activation energy parameter is
#'   named \code{Ea_log} and stored on the log‑scale. Ignored when \code{temp_expr} is provided.
#' @param yA_var,yB_var Column names for the primary response A and its
#'   companion response B (mass balance).
#' @param temp_var,time_var Column names for temperature (Kelvin) and time. \code{temp_var} ignored when \code{temp_expr} is provided
#' @param R Gas constant. Ignored when \code{temp_expr} is provided. Defaults to 8.314 J/(mol * K).
#' @param trefk Reference temperature in Kelvin. Ignored when \code{temp_expr} is provided
#' @param temp_expr Optional list output from \code{expr_arrhenius()} or \code{expr_arrhenius_humidity()}. If `NULL` (default) \code{expr_arrhenius()} is used.
#' @param start Optional **named** numeric vector of starting values.
#'   If `NULL` (default) every parameter is seeded explicitly: the activation
#'   energy from the gas constant you supply, any humidity coefficient at 0, each
#'   rate parameter at 0.01, and the intercepts \code{A0}/\code{B0} at the mean
#'   of \code{yA_var}/\code{yB_var} over the earliest observed \code{time_var}.
#'   Do not rely on \code{gsl_nls()} self-start here: it sets every unseeded
#'   parameter to 1, which for a percent-scale response with a split rate
#'   produces a converged-looking fit with an \code{NaN} residual sum of squares.
#' @param ... Additional arguments passed straight to
#'   \code{\link[gslnls]{gsl_nls}} (e.g. \code{algorithm =
#'   c("lm","lmaccel","dogleg","ddogleg","subspace2D")}).
#'
#' @return The fitted \code{gsl_nls} object.
#'
#' @examples
#' data("nrcesds")
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' library(tidyverse)
#' nrcesds_wide <- nrcesds %>% tidyr::pivot_wider(names_from = "cqa", values_from = "y")
#' fit <- fit_singledeg(
#'           df = nrcesds_wide,
#'           kinetic_order = 2,
#'           trend = "decreasing",
#'           yA_var = "NR Major Peak",
#'           yB_var = "NR Minor Peaks")
# summary(fit)
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' # This example uses the |> pipe and requires R version 4.1.0 or higher.
#' if (getRversion() >= "4.1.0") {
#'   library(tidyverse)
#'   fit <- nrcesds |>
#'     tidyr::pivot_wider(names_from = "cqa", values_from = "y") |>
#'     dplyr::rename(
#'       yA = "NR Major Peak", yB = "NR Minor Peaks"
#'     ) |>
#'     fit_singledeg(
#'       kinetic_order  = 2,
#'       trend          = "decreasing",
#'       intercept_cols = "batch",
#'       rate_cols      = character(),
#'       log_Ea         = TRUE
#'     )
#'   summary(fit)
#' }
#' @export
fit_singledeg <- function(
    df,
    kinetic_order  = 1L,
    trend          = c("increasing", "decreasing"),
    Yinf      = 100,
    intercept_cols = character(),
    rate_cols      = character(),
    log_Ea         = FALSE,
    yA_var         = "yA",
    yB_var         = "yB",
    temp_var       = "tempk",
    time_var       = "time",
    R              = 8.314,
    trefk          = 298.15,
    temp_expr      = NULL,
    start          = NULL,
    ...) {
  
  # ── 1. Expr params ───────────────────────────────
  te <- resolve_temp_expr(temp_expr, R = R, trefk = trefk, temp_var = temp_var, log_Ea = log_Ea)
  temp_expr <- te$temp_expr
  R         <- te$R
  trefk     <- te$trefk
  temp_var  <- te$temp_var
  log_Ea    <- te$log_Ea
  rh_col    <- te$rh_col
  rh_ref    <- te$rh_ref
  # Declared data references besides temp_var (rh_col included); anything not
  # listed here is read as a free parameter; see temp_expr_params().
  covar_cols <- te$covar_cols
  temp_fml  <- te$temp_fml
  
  # ── 2. Input validation ───────────────────────────────
  trend <- match.arg(trend)
  
  # Check required and optional columns exist in df
  missing_cols <- setdiff(
    c(yA_var, yB_var, temp_var, time_var, intercept_cols, rate_cols, covar_cols),
    names(df)
  )
  if (length(missing_cols)) {
    stop(
      sprintf("Missing columns in df: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  check_nonempty_df(df)

  # Before every structure-dependent check below: an NA grouping value otherwise
  # reaches the residual function and surfaces as the misleading "missing values
  # returned by 'fn' at starting values".
  check_structural_nas(df, c(intercept_cols, rate_cols, temp_var, time_var, covar_cols))

  # Arithmetic inputs; a non-numeric column produces the same misleading message.
  check_numeric_cols(df, c(temp_var, time_var, covar_cols))

  check_rh_scale(df, rh_col, rh_ref, what = "Model")
  check_rh_identifiable(df, rh_col, temp_fml)
  check_rate_support(df, rate_cols, time_var, intercept_cols, temp_var)
  check_rate_temp_support(df, rate_cols, temp_var)
  check_temp_identifiable(df, temp_var,
                          temp_expr_params(temp_fml, temp_var, covar_cols))
  check_trend_direction(df, yA_var, time_var, trend, temp_var)

  df_model <- sanitize_group_values(df, c(intercept_cols, rate_cols))

  # ──3. Record the mapping from original -> sanitized values ──────
  group_cols <- unique(c(intercept_cols, rate_cols))
  sanitize_map <- list()
  for (col in group_cols) {
    # names = original, values = sanitized; relies on sanitize_group_values()
    # preserving row order.
    lookup <- stats::setNames(
      as.character(df_model[[col]]),
      as.character(df[[col]])
    )
    lookup <- lookup[!duplicated(names(lookup))]
    sanitize_map[[col]] <- lookup
  }

  # !is.finite, not anyNA(): Inf would otherwise slip into a degenerate fit.
  resp_vars <- c(yA_var, yB_var)
  resp_with_na <- resp_vars[vapply(resp_vars, function(v) any(!is.finite(df[[v]])), logical(1))]
  if (length(resp_with_na)) {
    stop(
      sprintf("NA or non-finite values detected in response(s): %s", paste(resp_with_na, collapse = ", ")),
      call. = FALSE
    )
  } 
  
  # ── 4. Build model formula ───────────────────────────────
  build <- build_singledeg(
    df_model,
    kinetic_order = kinetic_order,
    trend  = match.arg(trend),
    Yinf = Yinf,
    intercept_cols = intercept_cols,
    rate_cols = rate_cols,
    log_Ea = log_Ea,
    yA_var = yA_var,
    yB_var = yB_var,
    temp_var = temp_var,
    time_var = time_var,
    R = R,
    trefk  = trefk,
    rh_col = rh_col,
    covar_cols = covar_cols,
    rh_ref = rh_ref,
    temp_fml = temp_fml
  )

  # ── 5. Estimate the model ───────────────────────────────
  # Seed every parameter (see default_start()); A0 from yA, B0 from its own
  # response yB, which mirrors yA and so trends the other way.
  if (is.null(start)) {
    start <- default_start(build$param_names, df_model, R, yA_var, trend,
                           time_var, kinetic_order = kinetic_order)
    b0 <- grepl("^B0", build$param_names)
    seed_B <- intercept_seed(df_model, yB_var, time_var,
                             if (trend == "decreasing") "increasing" else "decreasing")
    if (any(b0) && is.finite(seed_B)) start[b0] <- seed_B
  }

  validate_start_names(start, build$param_names)
  
  mod <- gslnls::gsl_nls(
    fn    = build$fn,
    y     = build$response,   # concatenated yA & yB
    data  = df_model,
    start = start,
    ...
  )
  # gsl_nls reports its exit status but does not warn, so a fit that hit the
  # iteration limit would look converged.
  check_convergence(mod, length(build$param_names))

  # A converged fit can still land on a negative Ea; warn
  check_activation_energy_sign(coef(mod))

  class(mod) <- c("fit_singledeg", class(mod))
  attr(mod, "build") <- build
  attr(mod, "sanitize_map") <- sanitize_map
  return(mod)
  
}
