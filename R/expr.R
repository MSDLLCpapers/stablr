#' Construct the Intercept Expression
#'
#' Generates the fixed-intercept part of the kinetic model formula and the
#' corresponding parameter names.
#'
#' When a single grouping column is provided, all levels receive their own
#' intercept parameter (cell-means coding). When two or more grouping columns
#' are provided, reference-level coding is used: the first column retains all
#' levels while each subsequent column drops its first (sorted) level, avoiding
#' a linear dependency in the Jacobian that would inflate confidence intervals.
#'
#' @section Coding differs from fit_sb():
#' This is the intercept coding of the closed-form fitters
#' (\code{\link{fit_kinetic}}, \code{\link{fit_singledeg}},
#' \code{\link{fit_parallel}}). With two or more \code{intercept_cols} the two paths fit
#' different models and report different residual degrees of freedom. Reference
#' coding represents only additive intercept effects, so it cannot capture a
#' genuine interaction between grouping columns. To make both paths agree, cross
#' the columns yourself and pass the result as a single column:
#' \code{df$cell <- interaction(df$batch, df$orient, drop = TRUE)}.
#'
#' @section Crossed vs nested columns:
#' Reference coding removes the dependency only for crossed designs, where every
#' column varies independently. A nested design stays collinear: if batches `A`,
#' `B` only appear at site `S1` and `C`, `D` only at `S2`, then site is a
#' function of batch and the indicator matrix is rank deficient. This function
#' checks the rank and errors when it is deficient. For a nested layout, pass
#' the single finest-grained column (e.g. `intercept_cols = "batch"`).
#'
#' @param df A data frame providing the grouping variables.
#' @param intercept_cols Character vector; grouping columns whose levels get
#'   their own intercepts. Use `NULL` or `character(0)` for a single global
#'   intercept.
#' @param prefix Character scalar naming the intercept parameter stem. Multi-
#'   response fits (\code{build_singledeg()}, \code{build_parallel()}) pass
#'   `"A0"`, `"B0"`, `"C0"` so each response gets its own parameters.
#'
#' @return A list with two elements
#' * `expr`- character scalar containing the additive intercept expression.
#' * `params`- character vector of parameter names.
#'
#' @examples
#' data(evers_supp)
#' # -------------------------------------------------
#' # Example 1: Standard usage (single factor - all levels kept)
#' # -------------------------------------------------
#' expr_intercept(df = evers_supp, intercept_cols = "pkg")
#'
#' # -------------------------------------------------
#' # Example 2: Two factors - reference coding applied
#' # -------------------------------------------------
#' df <- data.frame(batch = rep(c("A","B"), each=2),
#'                  orient = rep(c("inv","upr"), 2))
#' expr_intercept(df, c("batch", "orient"))
#' # First factor (batch): A, B  - all levels
#' # Second factor (orient): upr - "inv" dropped as reference
#'
#' @export
expr_intercept <- function(df, intercept_cols = NULL, prefix = "y0") {

  if (is.null(intercept_cols) || length(intercept_cols) == 0L ||
      all(intercept_cols == ""))
    return(list(expr = prefix, params = prefix))

  multi_factor <- length(intercept_cols) >= 2

  params   <- c()
  lvl_vals <- c()
  cols     <- c()
  sources  <- c()
  specs    <- list()

  for (i in seq_along(intercept_cols)) {
    v    <- intercept_cols[i]
    # Levels in parameter-name order (see observed_cells()).
    lvls <- observed_cells(df, v)[[v]]

    # Reference coding: drop first level of 2nd+ factors to avoid rank deficiency
    if (multi_factor && i >= 2) {
      lvls <- lvls[-1]
    }

    # Names first, then terms, since each term embeds its parameter name.
    # Sanitize the level into the name so it is a syntactic R name.
    params   <- c(params, as.character(glue::glue("{prefix}_{v}_{sanitize(lvls)}")))
    lvl_vals <- c(lvl_vals, as.character(lvls))
    cols     <- c(cols, rep(v, length(lvls)))
    sources  <- c(sources, sprintf("column '%s' level '%s'", v, lvls))
    # Recycle `col = v` explicitly so a single-level later factor yields an
    # empty spec rather than an error.
    specs[[i]] <- data.frame(col = rep(v, length(lvls)),
                             lvl = as.character(lvls),
                             stringsAsFactors = FALSE)
  }

  # Suffix name collisions rather than reject the design.
  params <- dedupe_param_names(params)
  terms  <- purrr::pmap_chr(list(cols, lvl_vals, params), indicator_term)

  # Nested columns stay collinear; error on rank deficiency.
  if (multi_factor) {
    check_intercept_rank(df, do.call(rbind, specs), params)
  }

  list(expr = glue::glue_collapse(terms, sep = " + "), params = params)
}


#' Construct the Rate Expression
#'
#' Builds the temperature-independent part of the rate constant, allowing for
#' full interaction of one or more grouping variables.
#'
#' @inheritParams expr_intercept
#' @param rate_cols Character vector of grouping columns giving distinct rate
#'   parameters.
#' @param prefix Character scalar naming the rate parameter stem. Multi-branch
#'   fits (\code{build_parallel()}) pass `"kB"` and `"kC"` so each branch gets
#'   its own parameters.
#'
#' @return A list with components `expr` and `params` analogous to
#'   [expr_intercept()].
#'
#' @examples
#' data(evers_supp)
#' # -------------------------------------------------
#' # Example 1: Standard usage
#' # -------------------------------------------------
#' expr_rate(df = evers_supp, rate_cols = "pkg")
#'
#' # -------------------------------------------------
#' # Example 2: Usage with the native R pipe
#' # -------------------------------------------------
#' # This example uses the |> pipe and requires R version 4.1.0 or higher.
#' if (getRversion() >= "4.1.0") {
#'   library(dplyr)
#'   evers_supp |>
#'     dplyr::filter(cqa == "HMW") |>
#'     expr_rate(rate_cols = c("form","pkg"))
#' }
#' @export
expr_rate <- function(df, rate_cols = NULL, prefix = "k") {

  null_par <- glue::glue("{prefix}_0")
  if (is.null(rate_cols) | length(rate_cols) == 0 ) return(list(expr = null_par, params = null_par))
  if (all(rate_cols == "") ) return(list(expr = null_par, params = null_par))

  # Enumerate the observed cells as rows so that level values containing "."
  # are handled correctly.
  cells <- observed_cells(df, rate_cols)

  conds   <- character()
  params  <- character()
  sources <- character()

  for (j in seq_len(nrow(cells))) {
    # value of each grouping variable in this combination, taken from the data
    parts <- vapply(rate_cols, function(v) as.character(cells[[v]][j]),
                    character(1), USE.NAMES = FALSE)
    # Compared against the level as it appears in the data.
    cond  <- purrr::map2_chr(
      rate_cols, parts,
      ~ glue::glue("({.x} == {level_literal(.y)})"))
    cond  <- glue::glue_collapse(cond, sep = " * ")

    combo_label <- paste(sanitize(parts), collapse = "_")
    par_name <- glue::glue("{prefix}_{paste(rate_cols, collapse = '_')}_{combo_label}")

    # The term embeds the parameter name, so settle the names first.
    conds   <- c(conds, as.character(cond))
    params  <- c(params, as.character(par_name))
    sources <- c(sources,
                 sprintf("cell (%s)",
                         paste(sprintf("%s='%s'", rate_cols, parts), collapse = ", ")))
  }

  # Suffix name collisions rather than reject a design fit_sb() accepts.
  params <- dedupe_param_names(params)
  terms  <- sprintf("(%s) * %s", conds, params)

  list(expr = glue::glue("({glue::glue_collapse(terms, sep = ' + ')})"),
       params = params)
}


#' Arrhenius Expression
#'
#' Returns the exponential Arrhenius factor as a character string ready to be
#' inserted in a model formula.
#'
#' @param R Gas constant (default 8.314 J mol-1 K-1).
#' @param trefk Reference temperature in Kelvin.
#' @param temp_var Name of the temperature column in the data.
#' @param log_Ea Logical; if `TRUE`, the parameter is assumed to be
#'   `Ea_log` (log-scale) and is exponentiated inside the expression.
#'
#' @return A list with Arrhenius expression inputs including an Arrhenius formula.
#'
#' @examples
#' expr_arrhenius()
#' expr_arrhenius(log_Ea = TRUE)
#' @export
expr_arrhenius <- function(
    R      = 8.314,
    trefk  = 298.15,
    temp_var = "tempk",
    log_Ea = F) {


  Ea_par <- ifelse(log_Ea, "exp(Ea_log)", "Ea")
  check_covar_names(c(temp_var = temp_var), if (log_Ea) "Ea_log" else "Ea")
  return(list(formula = glue::glue("exp({Ea_par} / {R} * (1/{trefk} - 1/{temp_var}))"),
              R = R,
              trefk = trefk,
              temp_var = temp_var,
              log_Ea = log_Ea)
  )
}


#' Modified Arrhenius Expression with Humidity
#'
#' Returns the Arrhenius factor with humidity covariate as a character string
#' ready to be inserted in a model formula.
#' The equation is: k(T, RH) = k_ref * exp(Ea/R * (1/T_ref - 1/T)) * exp(b * (RH - RH_ref))
#'
#' @section Humidity scale:
#' The \code{rh_col} column and \code{rh_ref} must be on the same scale. With
#' the default \code{rh_ref = 0.60} the humidity column must be a fraction
#' (e.g. `0.40`, `0.60`, `0.75`), not a percent. A scale mismatch is an exact
#' reparameterisation, not a lack of fit, so it is invisible to goodness-of-fit
#' diagnostics. The \code{fit_*} functions warn when the humidity range is
#' inconsistent with \code{rh_ref}, but check the scale yourself. To work in
#' percent, pass \code{rh_ref = 60}.
#'
#' @section Identifiability of b:
#' The data must span at least two humidity levels. With humidity held constant,
#' \code{exp(b * (rh - rh_ref))} is a constant absorbed into the rate
#' parameters, so `b` is aliased with `k_ref` and cannot be estimated. The
#' \code{fit_*} functions warn in this case. For a single-humidity study, fit
#' without a humidity term using \code{\link{expr_arrhenius}}.
#'
#' @param R Gas constant (default 8.314 J mol-1 K-1).
#' @param trefk Reference temperature in Kelvin.
#' @param temp_var Name of the temperature column in the data.
#' @param rh_col Name of the relative humidity column in the data (default
#'   "rh"). Must be on the same scale as \code{rh_ref}, a fraction in
#'   \code{[0, 1]} for the default \code{rh_ref}. See the Humidity scale section.
#' @param rh_ref Reference relative humidity, on the same scale as the
#'   \code{rh_col} data column (default 0.60, i.e. 60\% RH expressed as a
#'   fraction). \code{k_0} is the rate constant at this humidity.
#' @param log_Ea Logical; if `TRUE`, the parameter is assumed to be
#'   `Ea_log` (log-scale) and is exponentiated inside the expression.
#'
#' @return A list with humidity-modified Arrhenius expression inputs including a humidity-modified Arrhenius formula.
#'
#' @examples
#' expr_arrhenius_humidity()
#' expr_arrhenius_humidity(log_Ea = TRUE)
#' expr_arrhenius_humidity(rh_col = "humidity", rh_ref = 0.75)
#' # Percent-coded humidity data: put rh_ref on the percent scale too.
#' expr_arrhenius_humidity(rh_ref = 60)
#' @export
expr_arrhenius_humidity <- function(
    R        = 8.314,
    trefk    = 298.15,
    temp_var = "tempk",
    rh_col   = "rh",
    rh_ref   = 0.60,
    log_Ea   = FALSE) {

  Ea_par <- ifelse(log_Ea, "exp(Ea_log)", "Ea")
  check_covar_names(c(temp_var = temp_var, rh_col = rh_col),
                    c(if (log_Ea) "Ea_log" else "Ea", "b"))
  return(list(formula = glue::glue("exp({Ea_par} / {R} * (1/{trefk} - 1/{temp_var})) * exp(b * ({rh_col} - {rh_ref}))"),
              R = R,
              trefk = trefk,
              temp_var = temp_var,
              rh_col = rh_col,
              rh_ref = rh_ref,
              log_Ea = log_Ea)
  )
}