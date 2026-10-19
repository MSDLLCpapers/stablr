#' Assemble an NLS Kinetic Formula
#'
#' One‑stop builder that pieces together intercept, rate, Arrhenius and kinetic
#'‑order parts to create a self‑contained non‑linear least‑squares formula.
#'
#' @param df Data frame containing all variables referenced in the formula.
#' @param kinetic_order Integer reaction order (0, 1, 2, …).
#' @param trend \code{"increasing"} or \code{"decreasing"} outcome trajectory.
#' @param Yinf Asymptotic maximum for increasing trajectories.
#' @param intercept_cols Character vectors of grouping columns for intercept effects.
#' @param rate_cols Character vectors of grouping columns for rate effects.
#' @param log_Ea Logical; use `Ea_log` instead of `Ea`.
#' @param outcome_var,temp_var,time_var Column names for outcome, temperature (Kelvin) and time.
#' @param R Gas constant. Defaults to 8.314 J/(mol * K).
#' @param trefk Reference temperature (Kelvin).
#' @param temp_fml Character string from \code{expr_arrhenius()$formula} or \code{expr_arrhenius_humidity()$formula}.
#' @param rh_col Column name for relative humidity values (only needed when expr_arrhenius_humidity is used).
#' @param rh_ref Reference relative humidity value (only needed when expr_arrhenius_humidity is used).
#' @param covar_cols Character vector of every data column the temperature
#'   expression is declared to reference besides \code{temp_var} (\code{rh_col}
#'   included). Assembled by \code{resolve_temp_expr()}; symbols not named here are
#'   treated as free parameters.
#'
#' @return A list with
#' * \code{formula} – an object of class `formula`.
#' * \code{param_names} – character vector of all parameter names.
#' * \code{settings} - list containing build settings.
#' @examples
#' data(evers_supp)
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' build_kinetic(
#'   evers_supp, 
#'   kinetic_order = 1, 
#'   trend = "decreasing", 
#'   temp_fml = expr_arrhenius()$formula
#'   )
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' # This example uses the |> pipe and requires R version 4.1.0 or higher.
#' if (getRversion() >= "4.1.0") {
#'   evers_supp |>
#'     build_kinetic(
#'       kinetic_order = 1,
#'       trend = "increasing",
#'       Yinf = 100,
#'       intercept_cols = character(),
#'       rate_cols = c("pkg"),
#'       log_Ea = TRUE,
#'       temp_fml = expr_arrhenius()$formula
#'     )
#' }
#' @export
build_kinetic <- function(
    df,
    kinetic_order = 1L,
    trend = c("increasing", "decreasing"),
    Yinf = 100,
    intercept_cols = character(),
    rate_cols = character(),
    log_Ea = F,
    outcome_var = "y",
    temp_var = "tempk",
    time_var = "time",
    R = 8.314,
    trefk = 298.15,
    rh_col = NULL,
    rh_ref = NULL,
    covar_cols = NULL,
    temp_fml) {
  
  
  if(! trend %in% c("increasing", "decreasing")) stop("trend can be either `increasing` or `decreasing`")
  check_kinetic_order(kinetic_order)
  check_scalar_Yinf(Yinf)
  
  ## Building the component expressions --------------------------------
  int_part   <- expr_intercept(df, intercept_cols)
  rate_part <- expr_rate(df, rate_cols)
  arr_expr <- temp_fml

  ## Common sub‑expression: k(T) × t -----------------------------------
  rate_time <- glue::glue("{rate_part$expr}*{arr_expr}*{time_var}")
  int_expr <- glue::glue("({int_part$expr})")

  ## Right‑hand side (RHS) generator -----------------------------------
  if(kinetic_order == 0 & trend == "increasing"){
    rhs <- glue::glue("{int_expr} + {rate_time}")
  }else if(kinetic_order == 0 & trend == "decreasing"){
    rhs <- glue::glue("{int_expr} - {rate_time}")
  }else if(kinetic_order == 1 & trend == "increasing"){
    rhs <- glue::glue("{Yinf} - ({Yinf} - {int_expr})*exp(- {rate_time})")
  }else if(kinetic_order == 1 & trend == "decreasing"){
    rhs <- glue::glue("{int_expr}*exp(- {rate_time})")
  }else if(kinetic_order > 1 & trend == "increasing"){
    rhs <- glue::glue("{Yinf} - (( {Yinf} - {int_expr} )^(1-{kinetic_order}) - {rate_time}*(1-{kinetic_order}))^(1/(1-{kinetic_order}))    ")
  }else if(kinetic_order > 1 & trend == "decreasing"){
    rhs <- glue::glue("({int_expr}^(1-{kinetic_order}) - {rate_time}*(1-{kinetic_order}))^(1/(1-{kinetic_order}))")
  }
  
  ## Assemble formula and return ---------------------------------------
  fml <- stats::as.formula(glue::glue("`{outcome_var}` ~ {rhs}"))
  
  # Tokenise rather than grep, so renamed or extra parameters still get seeded.
  temp_params <- temp_expr_params(arr_expr, temp_var, unique(c(rh_col, covar_cols)))

  param_names <- c(temp_params, int_part$params, rate_part$params)

  # Cross-block duplicates: expr_intercept()/expr_rate() only see their own names.
  check_param_union(list("the temperature expression" = temp_params,
                         "the intercept term"         = int_part$params,
                         "the rate term"              = rate_part$params))

  # A data column named after a parameter would mask it and silently drop it.
  check_param_data_collisions(param_names, df)

  list(formula = fml,
       param_names = param_names,
       settings     = list(
         R = R, trefk = trefk, log_Ea = log_Ea,outcome_var = outcome_var, time_var = time_var, temp_var = temp_var,
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
         temp_range = temp_range(df, temp_var)))
}


#' Fit a Kinetic Model with \code{gslnls}
#'
#' Convenience wrapper that builds the formula, generates starting values and
#' fits the model via \link[gslnls]{gsl_nls}.
#'
#' @section Intercept coding differs from fit_sb():
#' With two or more \code{intercept_cols} the intercepts are reference coded (see
#' \code{\link{expr_intercept}}), giving \code{|A| + |B| - 1} additive
#' parameters, whereas \code{\link{fit_sb}} codes the same call as full cell
#' means (\code{|A| x |B|}). The two therefore fit different intercept models
#' with different parameter counts on identical input. Note also that
#' \code{rate_cols} is always cell means, so one call can mix cell-means rates
#' with reference-coded intercepts.
#'
#' @param df Data frame containing all variables referenced in the formula.
#' @param kinetic_order Integer reaction order (0, 1, 2, …).
#' @param trend \code{"increasing"} or \code{"decreasing"} outcome trajectory.
#' @param Yinf Asymptotic maximum for increasing trajectories.
#' @param intercept_cols Character vectors of grouping columns for intercept
#'   effects. Reference coded when two or more are given; see the section above.
#' @param rate_cols Character vectors of grouping columns for rate effects.
#' @param log_Ea Logical; use `Ea_log` instead of `Ea`. Ignored when \code{temp_expr} is provided.
#' @param outcome_var,temp_var,time_var Column names for outcome, temperature (Kelvin) and time. \code{temp_var} ignored when \code{temp_expr} is provided.
#' @param R Gas constant. Ignored when \code{temp_expr} is provided. Defaults to 8.314 J/(mol * K).
#' @param trefk Reference temperature (Kelvin). Ignored when \code{temp_expr} is provided.
#' @param temp_expr Optional list output from \code{expr_arrhenius()} or \code{expr_arrhenius_humidity()}. If `NULL` (default) \code{expr_arrhenius()} is used.
#' @param start Optional named numeric vector of starting values. If `NULL`
#'   (default) every parameter is seeded explicitly: the activation energy from
#'   the gas constant you supply (equivalent to 80000 J/mol at
#'   \code{R = 8.314}, on the log-scale when \code{log_Ea = TRUE}), any humidity
#'   coefficient at 0, each rate parameter at 0.01, and each intercept at the
#'   mean outcome over the earliest observed \code{time_var} value. Do not rely
#'   on \code{gsl_nls()} self-start here: it sets every unseeded parameter to 1,
#'   which for a percent-scale response with a split rate is far enough from the
#'   data that the fit diverges to an overflow instead of converging.
#' @param ... Additional arguments passed to \code{gsl_nls()}. For example
#'  algorithm = c("lm", "lmaccel", "dogleg", "ddogleg", "subspace2D"), or
#'  \code{weights} and \code{subset}, which may name columns of \code{df}.
#'
#' @return The fitted `gsl_nls` object.
#'
#' @examples
#' data(evers_supp)
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' dat <- subset(evers_supp, cqa == "Purity" & !is.na(y))
#' fit <- fit_kinetic(dat, kinetic_order = 1, trend = "decreasing")
#' summary(fit)
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' if (getRversion() >= "4.1.0") {
#'   library(dplyr)
#'   fit <- evers_supp |>
#'     dplyr::filter(cqa == "HMW", !is.na(y)) |>  # filter to one CQA
#'     dplyr::mutate(form = form |> as.factor() |> as.integer(),
#'                   pkg = pkg |> as.factor() |> as.integer()) |>
#'     fit_kinetic(kinetic_order = 1,
#'     trend = "increasing",
#'     intercept_cols = c("form"),
#'     rate_cols = c("form","pkg"),)
#'   summary(fit)
#' }
#' @export
fit_kinetic <- function(
    df,
    kinetic_order = 1L,
    trend = c("increasing", "decreasing"),
    Yinf = 100,
    intercept_cols = character(),
    rate_cols     = character(),
    log_Ea = FALSE,
    outcome_var    = "y",
    temp_var       = "tempk",
    time_var       = "time",
    R              = 8.314,
    trefk       = 298.15,
    temp_expr      = NULL,
    start          = NULL,
    ...) {
  
  # ── 1. Expr params  ───────────────────────────────
  te <- resolve_temp_expr(temp_expr, R = R, trefk = trefk, temp_var = temp_var, log_Ea = log_Ea)
  temp_expr <- te$temp_expr
  R         <- te$R
  trefk     <- te$trefk
  temp_var  <- te$temp_var
  log_Ea    <- te$log_Ea
  rh_col    <- te$rh_col
  rh_ref    <- te$rh_ref
  # Declared data references besides temp_var (rh_col included); anything not
  # listed here is read as a free parameter. See temp_expr_params().
  covar_cols <- te$covar_cols
  temp_fml  <- te$temp_fml
  
  # ── 2. Input validation ───────────────────────────────
  trend <- match.arg(trend)
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
  
  if (anyNA(df[[outcome_var]])) stop("Remove or impute NAs in the outcome before fitting.", call. = FALSE)
  # anyNA() misses Inf/-Inf, which yield a silent degenerate fit (sigma = Inf).
  if (!all(is.finite(df[[outcome_var]]))) stop("Non-finite (Inf) values in the outcome; remove or impute them before fitting.", call. = FALSE)
  check_structural_nas(df, c(intercept_cols, rate_cols, temp_var, time_var, covar_cols))

  # Arithmetic inputs
  check_numeric_cols(df, c(temp_var, time_var, covar_cols))

  check_rh_scale(df, rh_col, rh_ref, what = "Model")
  check_rh_identifiable(df, rh_col, temp_fml)
  check_rate_support(df, rate_cols, time_var, intercept_cols, temp_var)
  check_rate_temp_support(df, rate_cols, temp_var)
  check_temp_identifiable(df, temp_var,
                          temp_expr_params(temp_fml, temp_var, covar_cols))
  check_trend_direction(df, outcome_var, time_var, trend, temp_var)

  # Check that grouping levels give valid, distinct parameter names.
  check_group_levels(df, c(intercept_cols, rate_cols))
  df_model <- df
  
  # ── 4. Build model formula ───────────────────────────────
  build   <- build_kinetic(
    df_model,
    kinetic_order, trend, Yinf,
    intercept_cols, rate_cols,
    log_Ea,
    outcome_var, temp_var, time_var,
    R, trefk, rh_col, rh_ref, covar_cols, temp_fml)

  # ── 5. Estimate the model ───────────────────────────────
  # Seed every parameter: gsl_nls self-start uses 1 for each, which overflows as
  # soon as `rate_cols` splits the rate.
  if (is.null(start)){
    start <- default_start(build$param_names, df_model, R, outcome_var, trend,
                           time_var, kinetic_order = kinetic_order)
  }
  
  validate_start_names(start, build$param_names)
  
  # gsl_nls() rebuilds its model.frame() call via match.call()
  pf   <- parent.frame()
  dots <- lapply(as.list(substitute(list(...)))[-1L],
                 function(a) eval(a, df, pf))

  mod <- eval(as.call(c(list(quote(gslnls::gsl_nls),
                             fn    = quote(build$formula),
                             data  = quote(df_model),
                             start = quote(start)),
                       dots)))
  
  check_convergence(mod, length(build$param_names))
  check_activation_energy_sign(coef(mod))

  class(mod) <- c("fit_kinetic", class(mod))
  attr(mod, "build") <- build
  return(mod)
}