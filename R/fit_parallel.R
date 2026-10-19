#' Assemble Vector-Valued Residuals for a Parallel A -> \{B, C\} Model
#'
#' One-stop builder that combines intercept, rate, Arrhenius, and kinetic-order
#' pieces for a three-response parallel-reaction system \eqn{A \to \{B, C\}}.
#' The result is ready for \code{\link[gslnls]{gsl_nls}}.
#'
#' @param df Data frame containing every variable referenced by the model.
#' @param kinetic_order Integer reaction order (0, 1, 2, …) applied to A, B and C.
#' @param Yinf Not used by the parallel model. \code{Yinf} sets the asymptote
#'   only for increasing trajectories; the parallel source always decreases
#'   (toward 0 for first- and higher-order kinetics, linearly for zero-order).
#'   Retained for signature symmetry with the other models; a non-zero value is
#'   ignored with a warning.
#' @param intercept_cols,rate_cols Character vectors of grouping columns
#'   governing intercept and rate effects (copied per branch to avoid clashes).
#' @param temp_fml Character string from \code{expr_arrhenius()$formula} or \code{expr_arrhenius_humidity()$formula}.
#' @param log_Ea Logical; recorded in \code{settings} only. Branch activation
#'   energy names are derived from \code{temp_fml}, so pass a \code{temp_fml}
#'   built with \code{expr_arrhenius(log_Ea = TRUE)} to obtain
#'   \code{Ea_logB} / \code{Ea_logC}.
#' @param yA_var,yB_var,yC_var Column names for the three responses. \code{yA} is
#'   the source species (the main attribute that degrades), and its lost mass
#'   is partitioned into the two products \code{yB} and \code{yC} in proportion to
#'   the branch rate constants.
#' @param temp_var,time_var Column names for temperature (Kelvin) and time.
#' @param R Gas constant. Defaults to 8.314 J/(mol * K).
#' @param trefk Reference temperature (Kelvin).
#' @param rh_col Column name for relative humidity values (only needed when expr_arrhenius_humidity is used).
#' @param rh_ref Reference relative humidity value (only needed when expr_arrhenius_humidity is used).
#' @param covar_cols Character vector of every data column the temperature
#'   expression is declared to reference besides \code{temp_var} (\code{rh_col}
#'   included). Assembled by \code{resolve_temp_expr()}; symbols not named here are
#'   treated as free parameters.
#' @param temp_expr Optional temperature-expression list, recorded in
#'   \code{settings} for the predict methods. The expression actually evaluated
#'   is \code{temp_fml}; \code{\link{fit_parallel}} passes both.
#'
#' @return A list with
#' * \code{fn} – function \code{f(params, data)} returning the stacked vector
#'               \eqn{(y_A, y_B, y_C)} in the order expected by
#'               \code{gslnls}.  
#' * \code{param_names} – character vector of all parameter names (branch-
#'               specific intercepts, rates, and activation energies).  
#' * \code{response} – numeric vector \eqn{(y_A, y_B, y_C)} to pass as \code{y}.
#' * \code{settings} - list containing build settings.
#'
#' @details The branch rate constants are built as
#' \eqn{k_B(T) = k_B(\cdot)\,\exp(\cdot)} and \eqn{k_C(T) = k_C(\cdot)\,\exp(\cdot)}.
#' Mass balance across branches is allocated in proportion to \eqn{k_B} and \eqn{k_C}.
#'
#' @examples
#' data("upsec")
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' upsec_wide <- tidyr::pivot_wider(upsec, names_from = "cqa", values_from = "y")
#' build <- build_parallel(
#'           df = upsec_wide,
#'           kinetic_order = 1,
#'           yA_var = "Total Main",
#'           yB_var = "HMWS", 
#'           yC_var = "LMWS",
#'           temp_fml = expr_arrhenius()$formula
#' )
#' str(build$param_names)
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' # This example uses the |> pipe and requires R version 4.1.0 or higher.
#' if (getRversion() >= "4.1.0") {
#'   build <- upsec |>
#'     tidyr::pivot_wider(names_from = "cqa", values_from = "y") |>
#'     dplyr::rename(
#'       yA = "Total Main", yB = "HMWS", yC = "LMWS"
#'     ) |>
#'     build_parallel(
#'       kinetic_order  = 1,
#'       intercept_cols = "batch",
#'       rate_cols      = character(),
#'       log_Ea         = TRUE,
#'       temp_fml = expr_arrhenius()$formula
#'     )
#'   str(build$param_names)
#' }
#' @seealso \code{\link{expr_intercept}}, \code{\link{expr_rate}},
#'   \code{\link{expr_arrhenius}}, \code{\link[gslnls]{gsl_nls}}
#' @export
build_parallel <- function(
    df,
    kinetic_order = 1L,
    Yinf = 0,
    intercept_cols = character(),
    rate_cols      = character(),
    log_Ea         = FALSE,
    yA_var = "yA",
    yB_var = "yB",
    yC_var = "yC",
    temp_var = "tempk",
    time_var = "time",
    R      = 8.314,
    trefk  = 298.15,
    temp_expr      = NULL,
    rh_col = NULL,
    rh_ref = NULL,
    covar_cols = NULL,
    temp_fml
) {
  
  check_kinetic_order(kinetic_order)
  check_scalar_Yinf(Yinf)
  # Yinf is an asymptote for increasing trajectories only, and this model always
  # decreases, so it is unused; warn rather than silently ignore it.
  if (Yinf != 0) {
    warning(
      "`Yinf` is not used by the parallel model, whose source decreases ",
      "toward 0. Remove `Yinf` or set it to 0.",
      call. = FALSE
    )
  }

  # ── 1. Build the three intercept expressions ───────────────────────────────
  ## Prefix each family up front (not renamed afterwards) so parameter names
  ## don't clash: a rename would also rewrite column names and level values
  ## containing "y0".
  ie_A <- expr_intercept(df, intercept_cols, prefix = "A0")
  ie_B <- expr_intercept(df, intercept_cols, prefix = "B0")
  ie_C <- expr_intercept(df, intercept_cols, prefix = "C0")
  
  ## ── 2. Rate-constant ─────────────────────────────────
  ## Prefix up front as with the intercepts: a gsub over the built expression
  ## cannot tell a "k_" parameter stem from a column name or level value.
  rp_B <- expr_rate(df, rate_cols, prefix = "kB")
  rp_C <- expr_rate(df, rate_cols, prefix = "kC")
  
  ## ── 3. Temperature Dependence Expressions ─────────────────────────────────
  arr_base <- temp_fml
  arr_B <- arr_base
  arr_C <- arr_base

  # Declared data columns (e.g. "rh") must not be read as parameters.
  param_names <- temp_expr_params(arr_base, temp_var, unique(c(rh_col, covar_cols)))

  # Must run before the branch suffixes below, the only place it can. Once each
  # parameter is renamed to `<p>B`/`<p>C` the original name is gone and a data
  # column shadowing it is invisible to the later check.
  check_param_data_collisions(param_names, df,
                              what = "temperature-expression parameter")

  # Word boundaries avoid partial replacements.
  for (p in param_names) {
    pattern <- paste0("\\b", p, "\\b")
    arr_B <- gsub(pattern, paste0(p, "B"), arr_B, perl = TRUE)
    arr_C <- gsub(pattern, paste0(p, "C"), arr_C, perl = TRUE)
  }

  temp_params <- c(paste0(param_names, "B"), paste0(param_names, "C"))

  ## ── 4. Vector-valued residual function ─────────────────────────────────
  parallel_fn <- function(params, data) {
    
    penv <- list2env(as.list(params), parent = baseenv())
    
    A_int <- rlang::eval_tidy(rlang::parse_expr(ie_A$expr), data = data, env = penv)
    B_int <- rlang::eval_tidy(rlang::parse_expr(ie_B$expr), data = data, env = penv)
    C_int <- rlang::eval_tidy(rlang::parse_expr(ie_C$expr), data = data, env = penv)
    
    kB <- ( rlang::eval_tidy(rlang::parse_expr(rp_B$expr), data = data, env = penv) ) *
      rlang::eval_tidy(rlang::parse_expr(arr_B),    data = data, env = penv)
    kC <- ( rlang::eval_tidy(rlang::parse_expr(rp_C$expr), data = data, env = penv) ) *
      rlang::eval_tidy(rlang::parse_expr(arr_C),    data = data, env = penv)
    
    kApp <- kB + kC
    time <- data[[time_var]]
    rate_time <- (kApp) * time
    
    if (kinetic_order == 0) {
      yA_hat <- A_int - rate_time
    } else if (kinetic_order == 1) {                              # first-order
      yA_hat <- (A_int) * exp(-rate_time)
    } else {                                          # higher-order
      pow  <- 1 - kinetic_order
      yA_hat <- ( (A_int) ^ (pow) - (rate_time) * (pow) )^(1 / pow)
    }
    
    yB_hat <- B_int + (kB / (kApp) ) * (A_int - yA_hat)
    yC_hat <- C_int + (kC / (kApp) ) * (A_int - yA_hat)
    
    ## gsl_nls wants a *vector* of predictions in the same order as y
    c(yA_hat, yB_hat, yC_hat)
    
  }
  
  
  ## Parameter names -------------------------------------------------------
  par_names <- c(rp_B$params, rp_C$params,
                 temp_params,
                 ie_A$params, ie_B$params, ie_C$params)

  # Cross-block duplicates: the per-block checks cannot see across blocks.
  check_param_union(list("the B rate term"            = rp_B$params,
                         "the C rate term"            = rp_C$params,
                         "the temperature expression" = temp_params,
                         "the A intercept term"       = ie_A$params,
                         "the B intercept term"       = ie_B$params,
                         "the C intercept term"       = ie_C$params))

  # eval_tidy() puts the data mask first, so a column named after a parameter
  # would shadow it outright.
  check_param_data_collisions(par_names, df)

  list(
    fn          = parallel_fn,
    param_names = par_names,
    response    = c(df[[yA_var]], df[[yB_var]], df[[yC_var]]),
    settings     = list(
      R = R, trefk = trefk, log_Ea = log_Ea, temp_expr = temp_expr,
      yA_var = yA_var, yB_var=yB_var, yC_var=yC_var, time_var = time_var, temp_var = temp_var,
      Yinf = Yinf,
      rate_cols = rate_cols,
      rate_cells = rate_cells(df, rate_cols),
      group_levels = group_levels(df, c(intercept_cols, rate_cols)),
      unsupported_rate_cells = unsupported_rate_cells(df, rate_cols, time_var,
                                                      intercept_cols, temp_var),
      intercept_cols = intercept_cols,
      rh_col = rh_col,
      covar_cols = covar_cols,
      rh_ref = rh_ref,
      rh_range = rh_range(df, rh_col),
      temp_range = temp_range(df, temp_var))
  )
}


#' Fit the Parallel A -> \{B, C\} Model with \code{gsl_nls}
#'
#' Convenience wrapper that (1) builds the vector-valued residual function via
#' \code{\link{build_parallel}}, (2) generates default starting values when none
#' are provided, and (3) fits the model with \link[gslnls]{gsl_nls}.
#'
#' @param df Data frame containing every variable referenced by the model.
#' @param kinetic_order Integer reaction order (0, 1, 2, …) applied to A, B and C.
#' @param Yinf Not used by the parallel model. \code{Yinf} sets the asymptote
#'   only for increasing trajectories; the parallel source always decreases
#'   (toward 0 for first- and higher-order kinetics, linearly for zero-order).
#'   Retained for signature symmetry with the other models; a non-zero value is
#'   ignored with a warning.
#' @param intercept_cols,rate_cols Character vectors of grouping columns
#'   governing intercept and rate effects (copied per branch to avoid clashes).
#'   Intercepts are reference coded (additive) when two or more columns are given
#'   (see \code{\link{expr_intercept}}); \code{\link{fit_sb}} codes the same call
#'   as full cell means, so the two fit different intercept models with different
#'   parameter counts. Rates are always cell means.
#' @param temp_expr Optional list output from \code{expr_arrhenius()} or \code{expr_arrhenius_humidity()}. If `NULL` (default) \code{expr_arrhenius()} is used.
#' @param log_Ea Logical; if \code{TRUE}, branch activation energies are named
#'   \code{Ea_logB} and \code{Ea_logC} and are on the log-scale. Ignored when
#'   \code{temp_expr} is provided.
#' @param yA_var,yB_var,yC_var Column names for the three responses. \code{yA} is
#'   the source species (the main attribute that degrades), and its lost mass
#'   is partitioned into the two products \code{yB} and \code{yC} in proportion to
#'   the branch rate constants.
#' @param temp_var,time_var Column names for temperature (Kelvin) and time. \code{temp_var} ignored when \code{temp_expr} is provided.
#' @param R Gas constant. Ignored when \code{temp_expr} is provided. Defaults to 8.314 J/(mol * K).
#' @param trefk Reference temperature (Kelvin). Ignored when \code{temp_expr} is provided.
#' @param start Optional **named** numeric vector of initial values. If
#'   \code{NULL} (default) every parameter is seeded explicitly: \code{EaB} /
#'   \code{EaC} from the gas constant you supply, \code{bB} / \code{bC} at 0, each
#'   branch rate at 0.01, and the intercepts \code{A0} / \code{B0} / \code{C0} at
#'   the mean of \code{yA_var} / \code{yB_var} / \code{yC_var} over the earliest
#'   observed \code{time_var}. Do not rely on \code{gsl_nls()} self-start here:
#'   it sets every unseeded parameter to 1, which on a percent-scale response is
#'   far enough from the data that the fit can diverge instead of converging.
#' @param ... Additional arguments passed to \code{\link[gslnls]{gsl_nls}}
#'   (e.g. \code{algorithm = c("lm","lmaccel","dogleg","ddogleg","subspace2D")}).
#'
#' @return The fitted \code{gsl_nls} object.  
#'
#' @examples 
#' data("upsec")
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' upsec_wide <- tidyr::pivot_wider(upsec, names_from = "cqa", values_from = "y")
#' fit <- fit_parallel(
#'           df = upsec_wide,
#'           kinetic_order = 1,
#'           yA_var = "Total Main",
#'           yB_var = "HMWS",
#'           yC_var = "LMWS")
# summary(fit)
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' # This example uses the |> pipe and requires R version 4.1.0 or higher.
#' if (getRversion() >= "4.1.0") {
#'   fit <- upsec |>
#'     tidyr::pivot_wider(names_from = "cqa", values_from = "y") |>
#'     dplyr::rename(
#'       yA = "Total Main", yB = "HMWS", yC = "LMWS"
#'     ) |>
#'     fit_parallel(
#'       kinetic_order  = 1,
#'       intercept_cols = "batch",
#'       rate_cols      = character(),
#'       log_Ea         = FALSE
#'     )
#'   summary(fit)
#' }
#' @seealso \code{\link{build_parallel}}, \code{\link[gslnls]{gsl_nls}}
#' @export
fit_parallel <- function(
    df,
    kinetic_order = 1L,
    Yinf = 0,
    intercept_cols = character(),
    rate_cols      = character(),
    log_Ea         = FALSE,
    yA_var = "yA",
    yB_var = "yB",
    yC_var = "yC",
    temp_var = "tempk",
    time_var = "time",
    R      = 8.314,
    trefk  = 298.15,
    temp_expr      = NULL,
    start  = NULL,
    ...
) {
  
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
  missing_cols <- setdiff(
    c(yA_var, yB_var, yC_var, temp_var, time_var, intercept_cols, rate_cols, covar_cols),
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
  # No `trend` argument to point at here, so name the likely cause instead.
  check_trend_direction(
    df, yA_var, time_var, "decreasing", temp_var,
    advice = "Check that `yA` is the source species, not a degradant."
  )

  # Check that grouping levels give valid, distinct parameter names.
  check_group_levels(df, c(intercept_cols, rate_cols))
  df_model <- df

  # !is.finite, not anyNA(): Inf would otherwise slip into a degenerate fit.
  resp_vars <- c(yA_var, yB_var, yC_var)
  resp_with_na <- resp_vars[vapply(resp_vars, function(v) any(!is.finite(df[[v]])), logical(1))]
  if (length(resp_with_na)) {
    stop(
      sprintf("NA or non-finite values detected in response(s): %s", paste(resp_with_na, collapse = ", ")),
      call. = FALSE
    )
  } 
  
  
  # ── 4. Build model formula ───────────────────────────────
  build <- build_parallel(
    df_model,
    kinetic_order = kinetic_order,
    Yinf = Yinf,
    intercept_cols = intercept_cols,
    rate_cols      = rate_cols,
    temp_expr      = temp_expr,
    log_Ea         = log_Ea,
    yA_var = yA_var,
    yB_var = yB_var,
    yC_var = yC_var,
    temp_var = temp_var,
    time_var = time_var,
    R = R,
    trefk = trefk,
    rh_col = rh_col,
    covar_cols = covar_cols,
    rh_ref = rh_ref,
    temp_fml = temp_fml
  )

  # ── 5. Estimate the model ───────────────────────────────
  # default_start() handles the branch-suffixed names (EaB/EaC, bB/bC, kB_*/kC_*)
  # and seeds the intercepts, which gsl_nls would otherwise leave at 1. B and C
  # accumulate what A loses, so each seeds from its own baseline and increases.
  if (is.null(start)) {
    start <- default_start(build$param_names, df_model, R, yA_var, "decreasing",
                           time_var, kinetic_order = kinetic_order)
    for (spec in list(list("^B0", yB_var), list("^C0", yC_var))) {
      hit  <- grepl(spec[[1]], build$param_names)
      seed <- intercept_seed(df_model, spec[[2]], time_var, "increasing")
      if (any(hit) && is.finite(seed)) start[hit] <- seed
    }
  }

  validate_start_names(start, build$param_names)
  
  mod <- gslnls::gsl_nls(
    fn    = build$fn,
    y     = build$response,
    data  = df_model,
    start = start,
    ...
  )
  
  # gsl_nls reports its exit status but does not warn, so a fit that hit the
  # iteration limit would look converged.
  check_convergence(mod, length(build$param_names))

  # A converged fit can still land on a negative Ea;
  check_activation_energy_sign(coef(mod))

  class(mod) <- c("fit_parallel", class(mod))
  attr(mod, "build") <- build
  
  return(mod)
}

