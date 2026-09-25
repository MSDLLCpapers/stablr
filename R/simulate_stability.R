#### Simulate Kinetic Data Function ####

#' Simulate stability data for various kinetic reaction models
#'
#' This function generates simulated stability data for different kinetic
#' reaction models incorporating the modified Arrhenius equation. It supports
#' parallel reactions, single degradation reactions, and generic kinetic models
#' with zero-order, first-order, or second-order kinetics. Both increasing and
#' decreasing CQA trajectories are supported, consistent with the models fitted
#' by \code{\link{fit_kinetic}}, \code{\link{fit_singledeg}}, and
#' \code{\link{fit_parallel}}.
#'
#' @param factors Named vector or list specifying factors and their levels, or NULL.
#'   Names are factor names, values are the number of levels.
#'   Example: \code{c("conc" = 3, "pkg" = 2)} creates a design with 3
#'   concentration levels and 2 package types. If NULL (default), no factors
#'   are included and all observations share the same parameter values.
#' @param intercept_cols Character vector specifying which factors from \code{factors}
#'   should vary the intercept terms (initial CQA values). If NULL (default),
#'   no factors affect intercepts. Ignored if \code{factors} is NULL.
#' @param rate_cols Character vector specifying which factors from \code{factors}
#'   should vary the rate constant terms. If NULL (default), no factors affect
#'   rate constants. Ignored if \code{factors} is NULL.
#' @param time_points Numeric vector. Time points for sampling (e.g., in months).
#' @param temperatures_C Numeric vector. Storage temperatures in Celsius.
#' @param trefk Numeric. Reference temperature in Kelvin. Default is 298.15 K
#'   (25°C), matching the \code{trefk} default in \code{\link{fit_kinetic}},
#'   \code{\link{fit_singledeg}}, and \code{\link{fit_parallel}}. A value below
#'   100 warns of a likely Celsius slip.
#' @param reaction_type Character. Type of reaction model: "parallel" (A -> B and A -> C),
#'   "single_deg" (A -> B), or "kinetic" (generic kinetic model tracking only A).
#' @param kinetic_order Numeric. Order of the kinetic reaction: 0 for zero-order,
#'   1 for first-order, 2 for second-order.
#' @param trend Character. Direction of the CQA trajectory for species A:
#'   \code{"increasing"} (default) for accumulation/growth or \code{"decreasing"}
#'   for degradation. This matches the \code{trend} default in
#'   \code{\link{fit_kinetic}} and \code{\link{fit_singledeg}}. Not supported for
#'   \code{reaction_type = "parallel"}: an explicit \code{trend = "increasing"}
#'   with \code{"parallel"} is an error (it produces physically-impossible
#'   negative concentrations), and a bare parallel call is forced to
#'   \code{"decreasing"}.
#' @param Yinf Numeric. Asymptotic upper bound for increasing trajectories
#'   (first-order and higher-order only). Ignored when \code{trend = "decreasing"}
#'   and for \code{reaction_type = "parallel"} (whose source always decreases
#'   toward 0 for first- and higher-order kinetics, linearly for zero-order).
#'   Default is \code{A_max}. Corresponds to the \code{Yinf} parameter in the
#'   fitting functions.
#' @param params List. Model parameters. Each parameter can be specified as:
#'   \itemize{
#'     \item Single numeric value: Same value for all factor combinations
#'     \item Named list with one factor: \code{list(conc = c(94, 96, 98))} for
#'       factor-specific values (only when \code{factors} is not NULL and the
#'       factor is specified in \code{intercept_cols} or \code{rate_cols})
#'     \item Function: for one or more factors, e.g. \code{function(conc, pkg) 95 + 2*conc}
#'       (only when \code{factors} is not NULL and the factors are specified in
#'       \code{intercept_cols} or \code{rate_cols}). The function receives the
#'       factor level indices (1, 2, ...), not the level labels.
#'   }
#'   \code{sigma_err} (the marginal measurement-error SD) is the one exception:
#'   it must be a single non-negative numeric value, because the \code{fit_*}
#'   estimators assume a single pooled, homoscedastic error SD. Factor-varying
#'   or vector \code{sigma_err} is not yet supported.
#'   Each generating parameter is drawn from a distribution matched to its
#'   physical support, so all constraints hold by construction (no clamping,
#'   rescaling, or flooring). Initial values are compositional and drawn from a
#'   Dirichlet (Beta for two components; scaled Beta for a single component);
#'   rate constants are positive and drawn from a lognormal. Required parameters
#'   depend on \code{reaction_type}:
#'   \itemize{
#'     \item For "parallel": \code{A0_target}, \code{B0_target}, \code{C0_target}
#'       (composition centers that must sum to \code{A_max} per group),
#'       \code{A0_target_sd} (target SD on the \code{A0} anchor), \code{EaB},
#'       \code{EaC}, \code{krefB_median}, \code{krefB_gsd}, \code{krefC_median},
#'       \code{krefC_gsd}, \code{sigma_err}
#'     \item For "single_deg": \code{A0_target}, \code{B0_target} (must sum to
#'       \code{A_max} per group), \code{A0_target_sd}, \code{Ea},
#'       \code{kref_median}, \code{kref_gsd}, \code{sigma_err}
#'     \item For "kinetic": \code{A0_target} (in \code{[A_min, A_max]}),
#'       \code{A0_target_sd}, \code{Ea}, \code{kref_median}, \code{kref_gsd},
#'       \code{sigma_err}
#'   }
#'   \code{A0_target_sd} is the desired batch-to-batch SD for the \code{A0}
#'   anchor; it is converted to a Dirichlet concentration
#'   \eqn{\kappa = \mathrm{span}^2 \, p (1 - p) / \mathrm{sd}^2 - 1}
#'   (\eqn{\mathrm{span} = A_{max} - A_{min}}, \eqn{p} the \code{A0} fraction). A
#'   single \eqn{\kappa} governs the whole composition, so the minor components'
#'   SDs are then \emph{determined} by \eqn{\kappa} and their targets
#'   (\eqn{\mathrm{SD}_i = \mathrm{span}\sqrt{p_i(1 - p_i)/(\kappa + 1)}}), not
#'   independently specifiable. \code{A0_target_sd = 0} requests no batch-to-batch spread
#'   (the \code{kref_gsd = 1} analogue): every group is pinned to its exact
#'   target (\code{A0}, and for compositions \code{B0}/\code{C0} too), the
#'   closure-consistent zero-variance limit. Requesting an \code{A0_target_sd}
#'   at or above \eqn{\mathrm{span}\sqrt{p(1 - p)}} (the maximum feasible SD),
#'   or a negative value, is an error.
#'   Rate constants use a median and geometric SD (\code{gsd}, \eqn{\ge 1};
#'   \code{gsd = 1} means no batch-to-batch spread):
#'   \eqn{k_{ref} = \exp(\mathcal{N}(\log(\mathrm{median}), \log(\mathrm{gsd})^2))}.
#'   Activation energies (\code{Ea}, or \code{EaB}/\code{EaC}) accept the same
#'   factor-varying specifications and vary by the factors in \code{rate_cols},
#'   but are deterministic (there is no \code{sd_Ea}): an activation energy is
#'   fixed within a factor level rather than randomly perturbed.
#' @param component_names Named list. Optional custom names for output components.
#'   For "parallel": list(comp_A, comp_B, comp_C).
#'   For "single_deg" or "kinetic": list(comp_A, comp_B) or list(comp_A).
#'   If NULL, default names are used.
#' @param noise Character. Measurement-error model, one of
#'   \code{"mass_balance"} (default) or \code{"independent"}.
#'   \code{"mass_balance"} draws row-centered errors so the species sum is
#'   conserved exactly (a closed system, e.g. \code{A + B + C = A_max}); this
#'   reproduces realistic composition data that sums to \code{A_max}, at the cost
#'   of a small negative correlation between the species' errors.
#'   \code{"independent"} draws an independent \eqn{N(0, \sigma^2)} error per
#'   species, matching the unweighted, independent-error assumption of the
#'   \code{fit_*} least-squares estimators (so the fitted residual SD recovers
#'   \code{sigma_err} unbiasedly); the species sum is then only approximately
#'   \code{A_max}. Both models use \code{sigma_err} as the marginal error SD. For
#'   a single measured species (\code{reaction_type = "kinetic"}) the two models
#'   coincide.
#' @param clamp Logical. If \code{TRUE}, observed values are clamped to the
#'   physical bounds \code{[A_min, A_max]}. If \code{FALSE} (default), values may
#'   fall outside the bounds due to measurement noise, consistent with the
#'   unconstrained models fitted by the \code{fit_*} functions. Clamping is a
#'   simple bound (no post-clamp rescale), so it always keeps values within
#'   \code{[A_min, A_max]}; under \code{noise = "mass_balance"} it leaves the sum
#'   exactly \code{A_max} except on the rows where a value hit a bound. Censoring
#'   is applied silently to individual observations and is not recorded in the
#'   returned data (there is no indicator column identifying censored values); a
#'   \code{\link{warning}} reports how many observations were censored when any
#'   are.
#' @param A_max Numeric. The upper bound for the CQAs. Default is 100.
#'   Used to calculate initial CQA values of degradants and to constrain simulated
#'   values when \code{clamp = TRUE}.
#' @param A_min Numeric. The lower bound for CQAs. Default is 0.
#' @param R Numeric. Universal gas constant in J/(mol-K). Default is 8.314.
#' @param sdpt Integer. Seed for random number generation to ensure reproducibility.
#'   If NULL, no seed is set.
#' @param return_params Logical. If TRUE, returns a list containing both the
#'   simulated data and the parameters used to generate it. If FALSE (default),
#'   returns only the data frame.
#'
#' @return If \code{return_params = FALSE} (default): A data frame (tibble) in
#'   long format with columns:
#'   \describe{
#'     \item{factors}{Factor columns as specified in factors (if not NULL)}
#'     \item{time}{Numeric. Time point}
#'     \item{temp}{Numeric. Temperature in Celsius}
#'     \item{tempk}{Numeric. Temperature in Kelvin}
#'     \item{cqa}{Character. Component/quality attribute name}
#'     \item{y}{Numeric. Observed value (percentage)}
#'   }
#'
#'   If \code{return_params = TRUE}: A list with two elements:
#'   \describe{
#'     \item{data}{The data frame described above}
#'     \item{parameters}{A list containing:
#'       \itemize{
#'         \item \code{intercept_params}: Data frame of initial CQA
#'           parameters (A0, and B0/C0 for parallel reactions) for each factor
#'           combination affecting intercepts. If no factors affect intercepts,
#'           contains a single row with the parameter values used for all observations.
#'         \item \code{rate_params}: Data frame of rate constant (kref or krefB/krefC)
#'           and activation energy (Ea or EaB/EaC) parameters for each factor
#'           combination affecting rates. If no factors affect rates, contains a
#'           single row with the parameter values used for all observations.
#'         \item \code{fixed_params}: List of fixed parameters including trefk
#'           (reference temperature in K), R (gas constant), kinetic_order
#'           (reaction order), trend, Yinf, and sigma_err (measurement error
#'           standard deviation).
#'       }
#'     }
#'   }
#'
#' @details
#' The function uses the modified Arrhenius equation to calculate temperature-dependent
#' rate constants:
#' \deqn{k(T) = k_{ref} \exp\left(\frac{-E_a}{R}\left(\frac{1}{T} - \frac{1}{T_{ref}}\right)\right)}
#'
#' For decreasing CQA trajectories:
#' \itemize{
#'   \item Zero-order: \eqn{Y(t) = Y_0 - k(T) \cdot t}
#'   \item First-order: \eqn{Y(t) = Y_0 \exp(-k(T) \cdot t)}
#'   \item Second-order: \eqn{Y(t) = (Y_0^{1-n} - (1-n) k(T) t)^{1/(1-n)}}
#' }
#'
#' For increasing CQA trajectories:
#' \itemize{
#'   \item Zero-order: \eqn{Y(t) = Y_0 + k(T) \cdot t}
#'   \item First-order: \eqn{Y(t) = Y_\infty - (Y_\infty - Y_0) \exp(-k(T) \cdot t)}
#'   \item Second-order: \eqn{Y(t) = Y_\infty - ((Y_\infty - Y_0)^{1-n} - (1-n) k(T) t)^{1/(1-n)}}
#' }
#'
#' For parallel reactions, degradation products are apportioned based on the
#' relative magnitudes of their rate constants.
#'
#' Measurement error follows the \code{noise} argument. The default
#' \code{"mass_balance"} conserves the species sum exactly (row-centered errors)
#' for the parallel and single-degradation models, while keeping each species'
#' marginal SD equal to \code{sigma_err}; \code{"independent"} draws an
#' independent \eqn{N(0, \sigma^2)} error per species, matching the \code{fit_*}
#' least-squares assumptions. See the \code{noise} parameter for details.
#'
#' Activation energies (Ea, or EaB/EaC for parallel reactions) may vary by the
#' factors in \code{rate_cols}, using the same specification forms as the rate
#' medians (single value, named list, or function). Unlike rate constants, they
#' are deterministic: there is no \code{sd_Ea}, so an activation energy is fixed
#' within a factor level rather than randomly perturbed. Initial CQA values and
#' rate constants vary by factor levels as specified. Rate constants are drawn
#' from a lognormal and so are strictly positive by construction; initial values
#' are drawn from a Dirichlet/Beta distribution.
#'
#' When \code{return_params = TRUE}, the returned parameters show the actual
#' values drawn for each factor combination.
#'
#' @examples
#' \dontrun{
#' # Example 1: No factors; generic kinetic; 1st order; decreasing
#' params1 <- list(
#'   A0_target = 98,
#'   A0_target_sd = 0.5,
#'   Ea = 85000,
#'   kref_median = 2.5e-4,
#'   kref_gsd = 1.05,
#'   sigma_err = 0.3
#' )
#'
#' data1 <- simulate_kinetic_data(
#'   time_points = c(0, 1, 2, 3, 6, 9, 12, 18, 24),
#'   temperatures_C = c(5, 25, 40),
#'   reaction_type = "kinetic",
#'   kinetic_order = 1,
#'   trend = "decreasing",
#'   params = params1,
#'   sdpt = 123
#' )
#'
#' # Example 2: Generic kinetic; 1st order; increasing (e.g., HMW growth)
#' params2 <- list(
#'   A0_target = 2,
#'   A0_target_sd = 0.2,
#'   Ea = 85000,
#'   kref_median = 2.5e-4,
#'   kref_gsd = 1.05,
#'   sigma_err = 0.3
#' )
#'
#' data2 <- simulate_kinetic_data(
#'   time_points = c(0, 1, 2, 3, 6, 9, 12, 18, 24),
#'   temperatures_C = c(5, 25, 40),
#'   reaction_type = "kinetic",
#'   kinetic_order = 1,
#'   trend = "increasing",
#'   Yinf = 100,
#'   params = params2,
#'   sdpt = 234
#' )
#'
#' # Example 3: NR CE-SDS simulated data (single degradation, 2nd order)
#' # This generates the nrcesds dataset included with the package.
#' nrcesds <- simulate_kinetic_data(
#'   factors = c("batch" = 6),
#'   intercept_cols = "batch",
#'   time_points = c(0, 1, 2, 3, 6, 9, 12, 18, 24, 30, 36),
#'   temperatures_C = c(5, 25, 30, 37, 40),
#'   trefk = 278.15,
#'   reaction_type = "single_deg",
#'   kinetic_order = 2,
#'   trend = "decreasing",
#'   params = list(
#'     A0_target = 98,
#'     B0_target = 2,
#'     A0_target_sd = 0.4,
#'     Ea = 80000,
#'     kref_median = 5e-4,
#'     kref_gsd = 1.05,
#'     sigma_err = 0.3
#'   ),
#'   component_names = list(comp_A = "NR Major Peak", comp_B = "NR Minor Peaks"),
#'   clamp = TRUE,
#'   sdpt = 42
#' )
#'
#' # Example 4: UPSEC simulated data (parallel reaction, 1st order)
#' # This generates the upsec dataset included with the package.
#' upsec <- simulate_kinetic_data(
#'   factors = c("batch" = 6),
#'   intercept_cols = "batch",
#'   time_points = c(0, 1, 2, 3, 6, 9, 12, 18, 24, 30, 36),
#'   temperatures_C = c(5, 25, 30, 37, 40),
#'   trefk = 278.15,
#'   reaction_type = "parallel",
#'   kinetic_order = 1,
#'   trend = "decreasing",
#'   params = list(
#'     A0_target = 98,
#'     B0_target = 0.8,
#'     C0_target = 1.2,
#'     A0_target_sd = 0.5,
#'     EaB = 85000,
#'     EaC = 80000,
#'     krefB_median = 4e-4,
#'     krefB_gsd = 1.05,
#'     krefC_median = 2e-4,
#'     krefC_gsd = 1.05,
#'     sigma_err = 0.4
#'   ),
#'   component_names = list(comp_A = "Total Main", comp_B = "HMWS", comp_C = "LMWS"),
#'   clamp = TRUE,
#'   sdpt = 123
#' )
#' }
#'
#' @export
simulate_kinetic_data <- function(factors = NULL,
                                  intercept_cols = NULL,
                                  rate_cols = NULL,
                                  time_points,
                                  temperatures_C,
                                  trefk = 298.15,
                                  reaction_type = c("parallel", "single_deg", "kinetic"),
                                  kinetic_order = 1,
                                  trend = c("increasing", "decreasing"),
                                  Yinf = NULL,
                                  params,
                                  component_names = NULL,
                                  noise = c("mass_balance", "independent"),
                                  clamp = FALSE,
                                  A_max = 100,
                                  A_min = 0,
                                  R = 8.314,
                                  sdpt = NULL,
                                  return_params = FALSE) {

  # Match and validate reaction_type, trend, and noise
  reaction_type <- match.arg(reaction_type)
  trend_supplied <- !missing(trend)
  trend <- match.arg(trend)
  noise <- match.arg(noise)

  # trefk is in Kelvin; a value this low is almost certainly a Celsius slip.
  if (is.numeric(trefk) && length(trefk) == 1 && is.finite(trefk) && trefk < 100) {
    warning(sprintf(
      "trefk (%g) is in Kelvin, not Celsius. Did you mean %g?", trefk, trefk + 273.15
    ), call. = FALSE)
  }

  # The parallel model has no valid `trend = "increasing"` mechanism.
  if (reaction_type == "parallel") {
    if (trend_supplied && trend == "increasing") {
      stop(
        "trend = \"increasing\" is not supported for reaction_type = ",
        "\"parallel\" (it produces negative concentrations). Use ",
        "reaction_type = \"single_deg\" with trend = \"increasing\" for a single ",
        "growing species, or reaction_type = \"parallel\" with trend = ",
        "\"decreasing\" for two products forming from a source.",
        call. = FALSE
      )
    }
    trend <- "decreasing"
  }

  # Reserved column names owned by the design/output
  reserved_cols <- c("time", "temp", "tempk", "cqa", "y")

  # Validate params: must be a list containing exactly the parameters required
  # by reaction_type.
  if (!is.list(params)) {
    stop("params must be a list of model parameters.")
  }
  required_params <- switch(reaction_type,
    parallel   = c("A0_target", "B0_target", "C0_target", "A0_target_sd",
                   "EaB", "EaC", "krefB_median", "krefB_gsd",
                   "krefC_median", "krefC_gsd", "sigma_err"),
    single_deg = c("A0_target", "B0_target", "A0_target_sd",
                   "Ea", "kref_median", "kref_gsd", "sigma_err"),
    kinetic    = c("A0_target", "A0_target_sd",
                   "Ea", "kref_median", "kref_gsd", "sigma_err")
  )
  missing_params <- setdiff(required_params, names(params))
  if (length(missing_params) > 0) {
    stop(sprintf(
      "params is missing required parameter(s) for reaction_type '%s': %s",
      reaction_type, paste(missing_params, collapse = ", ")
    ))
  }

  # Validate sigma_err
  if (!is.numeric(params$sigma_err) || length(params$sigma_err) != 1 ||
      !is.finite(params$sigma_err) || params$sigma_err < 0) {
    stop(sprintf(
      paste0("Parameter 'sigma_err' must be a single non-negative numeric value ",
             "(the marginal measurement-error SD); per-factor measurement error ",
             "is not yet supported. Got: %s."),
      paste(utils::head(format(params$sigma_err), 3), collapse = ", ")
    ))
  }

  extra_params <- setdiff(names(params), required_params)
  if (length(extra_params) > 0) {
    warning(sprintf(
      "params contains unrecognized name(s) for reaction_type '%s' (ignored): %s",
      reaction_type, paste(extra_params, collapse = ", ")
    ))
  }

  # Validate component_names, if supplied: must have exactly the elements
  # required by reaction_type, with unique values that do not collide with the
  # reserved output columns.
  if (!is.null(component_names)) {
    needed_components <- switch(reaction_type,
      parallel   = c("comp_A", "comp_B", "comp_C"),
      single_deg = c("comp_A", "comp_B"),
      kinetic    = "comp_A"
    )
    missing_cn <- setdiff(needed_components, names(component_names))
    extra_cn   <- setdiff(names(component_names), needed_components)
    if (length(missing_cn) > 0 || length(extra_cn) > 0) {
      stop(sprintf(
        "component_names for reaction_type '%s' must have exactly: %s%s%s",
        reaction_type, paste(needed_components, collapse = ", "),
        if (length(missing_cn) > 0)
          sprintf(" (missing: %s)", paste(missing_cn, collapse = ", ")) else "",
        if (length(extra_cn) > 0)
          sprintf(" (unexpected: %s)", paste(extra_cn, collapse = ", ")) else ""
      ))
    }
    cn_values <- unlist(component_names[needed_components])
    if (anyDuplicated(cn_values) > 0) {
      stop(sprintf(
        "component_names values must be unique; got: %s",
        paste(cn_values, collapse = ", ")
      ))
    }
    cn_clash <- intersect(cn_values, reserved_cols)
    if (length(cn_clash) > 0) {
      stop(sprintf(
        "component_names collide with reserved output columns: %s",
        paste(cn_clash, collapse = ", ")
      ))
    }
  }

  # Set Yinf default if not provided
  if (is.null(Yinf)) {
    Yinf <- A_max
  }

  # Validate kinetic_order
  if (!kinetic_order %in% c(0, 1, 2)) {
    stop("kinetic_order must be 0 (zero-order), 1 (first-order), or 2 (second-order)")
  }

  # Handle NULL factors case
  if (is.null(factors)) {
    # No factors - set intercept_cols and rate_cols to empty
    intercept_cols <- character(0)
    rate_cols <- character(0)
    factors <- list()  # Empty list for consistency
  } else {
    # Convert factors to named list if vector
    if (!is.list(factors)) factors <- as.list(factors)

    # Validate that factors is named
    if (is.null(names(factors)) || any(names(factors) == "")) {
      stop("factors must be a named vector or list")
    }

    # Validate that each level count is a single positive whole number.
    bad_counts <- !vapply(factors, function(n) {
      length(n) == 1 && is.numeric(n) && is.finite(n) && n >= 1 && n == round(n)
    }, logical(1))
    if (any(bad_counts)) {
      offenders <- names(factors)[bad_counts]
      stop(sprintf(
        "factor level counts must each be a single positive integer; got %s.",
        paste(sprintf("%s = %s", offenders,
                      vapply(factors[bad_counts],
                             function(n) paste(format(n), collapse = ", "),
                             character(1))),
              collapse = ", ")
      ))
    }

    # Reject factor names that collide with reserved design/output columns
    factor_clash <- intersect(names(factors), reserved_cols)
    if (length(factor_clash) > 0) {
      stop(sprintf(
        "factor name(s) collide with reserved design/output columns (%s): %s. Please rename these factors.",
        paste(reserved_cols, collapse = ", "), paste(factor_clash, collapse = ", ")
      ))
    }

    # Handle NULL intercept_cols and rate_cols
    if (is.null(intercept_cols)) {
      intercept_cols <- character(0)
    }
    if (is.null(rate_cols)) {
      rate_cols <- character(0)
    }

    # Validate that intercept_cols and rate_cols are in factors
    invalid_intercept <- setdiff(intercept_cols, names(factors))
    if (length(invalid_intercept) > 0) {
      stop(sprintf(
        "intercept_cols contains factor names not in factors: %s",
        paste(invalid_intercept, collapse = ", ")
      ))
    }

    invalid_rate <- setdiff(rate_cols, names(factors))
    if (length(invalid_rate) > 0) {
      stop(sprintf(
        "rate_cols contains factor names not in factors: %s",
        paste(invalid_rate, collapse = ", ")
      ))
    }

    # Warning if factors specified but not used
    if (length(intercept_cols) == 0 && length(rate_cols) == 0) {
      warning(
        "factors were specified but neither intercept_cols nor rate_cols were provided. ",
        "All observations will share the same parameter values. ",
        "The factors will only create different combinations in the experimental design."
      )
    }
  }

  # Set seed if provided
  if (!is.null(sdpt)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(sdpt)
  }

  # Create study design with factors
  if (length(factors) > 0) {
    factor_grid <- lapply(factors, function(n_levels) factor(1:n_levels))
    names(factor_grid) <- names(factors)
    factor_grid$time <- time_points
    factor_grid$temp <- temperatures_C
    factor_grid$stringsAsFactors <- FALSE
    design <- do.call(expand.grid, factor_grid)
  } else {
    # No factors - just time and temperature
    design <- expand.grid(
      time = time_points,
      temp = temperatures_C,
      stringsAsFactors = FALSE
    )
  }

  # Simulate based on reaction type
  if (reaction_type == "parallel") {
    sim_result <- simulate_parallel_kinetics(
      design, factors, intercept_cols, rate_cols,
      trefk, kinetic_order, trend, Yinf, params, component_names, R,
      A_max, A_min, noise, clamp,
      return_params = return_params
    )
  } else if (reaction_type == "single_deg") {
    sim_result <- simulate_single_deg_kinetics(
      design, factors, intercept_cols, rate_cols,
      trefk, kinetic_order, trend, Yinf, params, component_names, R,
      A_max, A_min, noise, clamp,
      return_params = return_params
    )
  } else if (reaction_type == "kinetic") {
    sim_result <- simulate_generic_kinetics(
      design, factors, intercept_cols, rate_cols,
      trefk, kinetic_order, trend, Yinf, params, component_names, R,
      A_max, A_min, noise, clamp,
      return_params = return_params
    )
  }

  if (return_params) {
    # Return list with data and parameters
    return(list(
      data = sim_result$data,
      parameters = sim_result$parameters
    ))
  } else {
    # Return only data
    if (is.list(sim_result) && "data" %in% names(sim_result)) {
      return(sim_result$data)
    } else {
      return(sim_result)
    }
  }
}


#### Helper function to parse parameters ####

#' Parse parameter specification and generate values for factor combinations
#'
#' @param param_spec Parameter specification (numeric, list, or function)
#' @param param_name Name of the parameter (for error messages)
#' @param factor_combos Data frame of factor combinations
#' @param relevant_factors Character vector of factor names that affect this parameter
#' @param factors List of all factors with their levels
#'
#' @return Numeric vector of parameter values for each factor combination
#' @keywords internal
#' @noRd
parse_parameter_spec <- function(param_spec, param_name, factor_combos,
                                 relevant_factors, factors) {

  # Case 1: Single numeric value - use for all combinations
  if (is.numeric(param_spec) && length(param_spec) == 1) {
    return(rep(param_spec, nrow(factor_combos)))
  }

  # Case 2: Named list with single factor
  if (is.list(param_spec) && !is.function(param_spec)) {
    # Check if factors are being used
    if (length(relevant_factors) == 0) {
      stop(sprintf(
        "Parameter '%s': Cannot use named list specification when no factors affect this parameter. Use a single numeric value.",
        param_name
      ))
    }

    if (length(param_spec) != 1) {
      stop(sprintf(
        "Parameter '%s': Named list specification must contain exactly one factor. For multiple factors, use a function.",
        param_name
      ))
    }

    factor_name <- names(param_spec)[1]
    factor_values <- param_spec[[1]]

    # Validate factor name
    if (!factor_name %in% relevant_factors) {
      stop(sprintf(
        "Parameter '%s': Factor '%s' is not in the relevant factors for this parameter (%s)",
        param_name, factor_name, paste(relevant_factors, collapse = ", ")
      ))
    }

    # Validate number of values matches factor levels
    expected_levels <- factors[[factor_name]]
    if (length(factor_values) != expected_levels) {
      stop(sprintf(
        "Parameter '%s': Factor '%s' has %d levels but %d values were provided",
        param_name, factor_name, expected_levels, length(factor_values)
      ))
    }

    # Map values to factor combinations
    factor_indices <- as.numeric(factor_combos[[factor_name]])
    return(factor_values[factor_indices])
  }

  # Case 3: Function - evaluate for each combination
  if (is.function(param_spec)) {
    # Check if factors are being used
    if (length(relevant_factors) == 0) {
      stop(sprintf(
        "Parameter '%s': Cannot use function specification when no factors affect this parameter. Use a single numeric value.",
        param_name
      ))
    }

    # Get function arguments
    func_args <- names(formals(param_spec))

    # Check that all function arguments are in relevant_factors
    invalid_args <- setdiff(func_args, relevant_factors)
    if (length(invalid_args) > 0) {
      stop(sprintf(
        "Parameter '%s': Function arguments (%s) include factors not relevant to this parameter (%s)",
        param_name, paste(invalid_args, collapse = ", "),
        paste(relevant_factors, collapse = ", ")
      ))
    }

    # Evaluate for each combination, passing factor level indices (not labels).
    param_values <- vapply(seq_len(nrow(factor_combos)), function(i) {
      args <- lapply(factor_combos[i, func_args, drop = FALSE],
                     function(v) if (is.factor(v)) as.numeric(v) else v)
      do.call(param_spec, args)
    }, numeric(1))

    return(param_values)
  }

  # Case 4: Numeric vector (not named list) - error
  if (is.numeric(param_spec) && length(param_spec) > 1) {
    stop(sprintf(
      "Parameter '%s': Numeric vectors must be wrapped in a named list. Use list(factorname = c(...)) or provide a single value.",
      param_name
    ))
  }

  # Unknown specification type
  stop(sprintf(
    "Parameter '%s': Invalid specification. Must be a single numeric value, named list with one factor, or function.",
    param_name
  ))
}


#### Unified parameter-draw layer (compositions and rates) ####

#' Dirichlet/Beta concentration from a target SD on the anchor component
#'
#' Inverts the marginal variance of a Beta component of a closed composition,
#' \eqn{\mathrm{Var} = \mathrm{span}^2 \, p (1 - p) / (\kappa + 1)}, to recover the
#' concentration \eqn{\kappa} that yields a requested standard deviation
#' \code{sd_anchor} on the anchor component (whose target is fraction \code{p} of
#' the range \code{span}). Vectorised over groups.
#'
#' A single \eqn{\kappa} governs the whole composition, so honouring the anchor
#' SD determines (does not free) the minor components' SDs
#' (\eqn{\mathrm{SD}_i = \mathrm{span}\sqrt{p_i(1 - p_i)/(\kappa + 1)}}).
#'
#' @param sd_anchor Numeric vector. Target SD on the anchor component, per group.
#' @param p Numeric vector. Anchor target as a fraction of \code{span}, per group.
#' @param span Numeric scalar. Range width the composition spans.
#'
#' @return Numeric vector of concentrations \eqn{\kappa}, one per group.
#' @keywords internal
#' @noRd
kappa_from_sd <- function(sd_anchor, p, span) {
  # Feasibility: kappa > 0 requires sd_anchor < span * sqrt(p(1-p)), the largest
  # SD a closed composition with this anchor fraction can have.
  max_sd <- span * sqrt(p * (1 - p))
  bad <- which(sd_anchor >= max_sd | !is.finite(max_sd) | max_sd <= 0)
  if (length(bad) > 0) {
    i <- bad[1]
    stop(sprintf(
      paste0("A0_target_sd (%g) is too large for the anchor target: it must be ",
             "< %g (= span * sqrt(p*(1-p)), the maximum feasible SD for a closed ",
             "composition with this A0 fraction). Use a smaller A0_target_sd or ",
             "move A0_target away from the bounds."),
      sd_anchor[i], max_sd[i]
    ))
  }
  span^2 * p * (1 - p) / sd_anchor^2 - 1
}

#' Draw compositions/intercepts from their natural support (one draw per group)
#'
#' Draws initial values on the support:
#' \itemize{
#'   \item \code{m == 1} (generic kinetic): a scalar \code{A0} on
#'     \code{[A_min, A_max]}, drawn as a scaled Beta of
#'     \code{(A0 - A_min, A_max - A0)}, so a non-zero \code{A_min} is respected.
#'   \item \code{m > 1} (single_deg, parallel): a closed percent-composition on
#'     the simplex \code{{x >= 0, sum(x) = A_max}}, drawn as a Dirichlet.
#'     \code{A_min} plays no role for a composition anchored at 0.
#' }
#' The anchor is always component 1 (\code{A0}); \code{sd_anchor} sets its target
#' SD via \code{\link{kappa_from_sd}}.
#'
#' @param targets Numeric matrix (ng x m) of per-group component targets. For
#'   \code{m > 1} each row sums to \code{A_max}; for \code{m == 1} the single
#'   column lies in \code{[A_min, A_max]}.
#' @param sd_anchor Numeric vector (length ng). Target SD on the anchor, per group.
#' @param A_min,A_max Numeric scalars. Physical bounds for the CQA.
#'
#' @return Numeric matrix (ng x m) of drawn values.
#' @keywords internal
#' @noRd
draw_simplex <- function(targets, sd_anchor, A_min, A_max) {
  ng <- nrow(targets)
  m  <- ncol(targets)
  out <- matrix(NA_real_, ng, m)

  # A negative SD is meaningless.
  neg <- which(sd_anchor < 0)
  if (length(neg) > 0) {
    stop(sprintf(
      "A0_target_sd (%g) must be >= 0 (0 = no batch-to-batch spread).",
      sd_anchor[neg[1]]
    ))
  }
  # sd_anchor == 0 requests no batch-to-batch spread (the kref_gsd = 1 analogue
  # for compositions). Handled per group so a mixed spec (e.g. sd = c(0, 0.5)) pins the
  # zero-SD groups and randomly draws the rest. kappa is computed only where
  # sd_anchor > 0.
  pos <- sd_anchor > 0

  if (m == 1L) {
    # Generic kinetic: scaled Beta of (A0 - A_min, A_max - A0) on [A_min, A_max].
    span <- A_max - A_min
    p <- (targets[, 1] - A_min) / span
    kappa <- rep(NA_real_, ng)
    if (any(pos)) kappa[pos] <- kappa_from_sd(sd_anchor[pos], p[pos], span)
    for (i in seq_len(ng)) {
      out[i, 1] <- if (sd_anchor[i] == 0) {
        targets[i, 1]
      } else {
        A_min + span * stats::rbeta(1, p[i] * kappa[i], (1 - p[i]) * kappa[i])
      }
    }
    return(out)
  }

  # Closed %-composition on {>= 0, sum = A_max}: Dirichlet anchored at 0.
  span <- A_max
  frac <- targets / span                        # rows sum to 1
  kappa <- rep(NA_real_, ng)
  if (any(pos)) {
    # anchor = component 1 (A0)
    kappa[pos] <- kappa_from_sd(sd_anchor[pos], frac[pos, 1], span)
  }
  for (i in seq_len(ng)) {
    if (sd_anchor[i] == 0) {
      out[i, ] <- targets[i, ]                  # zero spread: exact targets (sum = A_max)
    } else {
      g <- stats::rgamma(m, shape = frac[i, ] * kappa[i])
      out[i, ] <- g / sum(g) * span
    }
  }
  out
}

#' Draw rate constants from a lognormal
#'
#' Rate constants are positive and multiplicative, so the lognormal is their
#' natural law (and matches the Arrhenius \eqn{k_{ref} = \exp(\cdot)} form).
#' Parameterised by a median and geometric SD in the mechanism's own units:
#' \eqn{k_{ref} = \exp(\mathrm{rnorm}(\log(\mathrm{median}),\ \log(\mathrm{gsd})))}.
#' Strictly positive, so no floor is ever needed.
#'
#' @param median_vec Numeric vector. Per-group rate median (> 0).
#' @param gsd_vec Numeric vector. Per-group geometric SD (>= 1; 1 = no spread).
#'
#' @return Numeric vector of drawn rate constants.
#' @keywords internal
#' @noRd
draw_rate_lognormal <- function(median_vec, gsd_vec) {
  if (any(median_vec <= 0)) {
    stop("rate median(s) must be > 0.")
  }
  if (any(gsd_vec < 1)) {
    stop("rate gsd (geometric SD) must be >= 1 (1 = no batch-to-batch spread).")
  }
  stats::rlnorm(length(median_vec), meanlog = log(median_vec), sdlog = log(gsd_vec))
}

#' Validate that per-group composition targets sum to A_max
#'
#' @param targets Numeric matrix (ng x m) of component targets.
#' @param A_max Numeric scalar. Required per-group sum.
#' @param reaction_type Character. Used in the error message.
#'
#' @return Invisibly \code{TRUE}; stops on the first offending group.
#' @keywords internal
#' @noRd
validate_composition_targets <- function(targets, A_max, reaction_type) {
  s <- rowSums(targets)
  bad <- which(abs(s - A_max) > 1e-8)
  if (length(bad) > 0) {
    i <- bad[1]
    stop(sprintf(
      "composition targets for reaction_type '%s' must sum to A_max (%g) per group; group %d sums to %g.",
      reaction_type, A_max, i, s[i]
    ))
  }
  invisible(TRUE)
}


#### Helper function to draw measurement noise ####

#' Draw additive measurement noise for the observed species
#'
#' Generates an \code{n_obs} x \code{k} matrix of additive Gaussian measurement
#' errors under one of two models:
#' \describe{
#'   \item{\code{"independent"}}{Each species gets an independent
#'     \eqn{N(0, \sigma^2)} draw. Matches the unweighted, homoscedastic,
#'     independent-error assumption of the \code{fit_*} least-squares estimators,
#'     so the fitted residual SD recovers \code{sigma_err} unbiasedly. Row sums
#'     are only approximately conserved (off by ~\code{sigma_err}).}
#'   \item{\code{"mass_balance"}}{Row-centered noise, so the per-row errors sum
#'     to zero and the species sum is conserved exactly (a closed system,
#'     \eqn{A + B + \dots = A_{max}}). Centering shrinks each marginal variance to
#'     \eqn{(k-1)/k \sigma^2}, so the draws are rescaled by \eqn{\sqrt{k/(k-1)}}
#'     to keep the marginal SD equal to \code{sigma_err}. The trade-off is a
#'     \eqn{-1/(k-1)} correlation between species' errors.}
#' }
#' For \code{k = 1} (a single measured species) there is nothing to balance, so
#' both models reduce to a single independent draw.
#'
#' @param n_obs Integer. Number of observations (rows).
#' @param k Integer. Number of species (columns): 3 for parallel, 2 for
#'   single_deg, 1 for generic.
#' @param sigma Numeric. Marginal measurement-error SD (\code{sigma_err}).
#' @param noise Character. \code{"mass_balance"} or \code{"independent"}.
#'
#' @return Numeric matrix (\code{n_obs} x \code{k}) of measurement errors.
#' @keywords internal
#' @noRd
draw_noise <- function(n_obs, k, sigma, noise) {
  eps <- matrix(stats::rnorm(n_obs * k, mean = 0, sd = sigma), ncol = k)
  if (noise == "mass_balance" && k > 1) {
    # Center each row so errors sum to zero (exact mass balance), then rescale to
    # restore the marginal SD that row-centering would otherwise shrink.
    eps <- (eps - rowMeans(eps)) * sqrt(k / (k - 1))
  }
  eps
}


#### Parallel Reaction ####

#' Helper function for parallel reactions (A -> B and A -> C)
#'
#' This internal function simulates parallel degradation reactions where a main
#' component A degrades into two products B and C simultaneously. The source A
#' always decreases (there is no valid \code{trend = "increasing"} parallel
#' mechanism; the caller fixes \code{trend = "decreasing"}). The kinetics follow
#' the modified Arrhenius equation with temperature-dependent rate constants.
#'
#' @param design Data frame with experimental design including factor levels,
#'   time points, and temperatures.
#' @param factors List of all factors with their levels.
#' @param intercept_cols Character vector of factor names that affect intercepts
#'   (initial CQA values).
#' @param rate_cols Character vector of factor names that affect rate constants.
#' @param trefk Numeric. Reference temperature in Kelvin.
#' @param kinetic_order Numeric. Order of kinetics (0, 1, or 2).
#' @param trend Character. Always \code{"decreasing"} for the parallel model.
#' @param Yinf Numeric. Unused for the parallel model (the source decays to 0).
#' @param params List. Model parameters including A0_target, B0_target,
#'   C0_target, A0_target_sd, EaB, EaC, krefB_median, krefB_gsd, krefC_median,
#'   krefC_gsd, and sigma_err.
#' @param component_names Named list. Optional component names for output.
#' @param R Numeric. Universal gas constant in J/(mol-K).
#' @param A_max Numeric. The upper bound for the CQAs.
#' @param A_min Numeric. The lower bound for CQAs.
#' @param noise Character. Measurement-error model, \code{"mass_balance"} or
#'   \code{"independent"} (see \code{draw_noise()}).
#' @param clamp Logical. Whether to clamp observed values to \code{[A_min, A_max]}.
#' @param return_params Logical. If TRUE, returns both data and parameters.
#'   If FALSE, returns only data. Default is FALSE.
#'
#' @return If \code{return_params = FALSE}: A data frame in long format with
#'   simulated observations.
#'
#'   If \code{return_params = TRUE}: A list with two elements:
#'   \describe{
#'     \item{data}{Data frame in long format with columns for factors, time,
#'       temp, tempk, cqa (component name), and y (observed value)}
#'     \item{parameters}{List containing:
#'       \itemize{
#'         \item \code{intercept_params}: Data frame with columns for factor
#'           combinations (if applicable), A0, B0, and C0
#'         \item \code{rate_params}: Data frame with columns for factor
#'           combinations (if applicable), Ea_B, Ea_C, kref_B, and kref_C
#'         \item \code{fixed_params}: List with trefk, R, kinetic_order,
#'           trend, Yinf, and sigma_err
#'       }
#'     }
#'   }
#'
#' @details
#' The function implements parallel reactions where component A
#' converts into B and C with rate constants kB and kC. The total rate
#' is k_total = kB + kC, and products are apportioned based on the ratio
#' of individual rate constants to the total rate.
#'
#' Initial CQA values (A0, B0, C0) are drawn from a Dirichlet distribution so
#' they are non-negative and sum to \code{A_max}. Rate constants
#' are drawn from a lognormal, so they are strictly positive.
#'
#' Measurement error is added per the \code{noise} model (see
#' \code{draw_noise()}). Under \code{"mass_balance"} the three species' errors
#' are row-centered so they sum to zero, preserving mass balance
#' (\eqn{A + B + C = A_0 + B_0 + C_0}) while keeping each species' marginal SD
#' equal to \code{sigma_err}. Under \code{"independent"} each species gets an
#' independent \eqn{N(0, \sigma^2)} error.
#'
#' @keywords internal
#' @noRd
simulate_parallel_kinetics <- function(design, factors, intercept_cols,
                                       rate_cols, trefk, kinetic_order,
                                       trend, Yinf,
                                       params, component_names, R,
                                       A_max, A_min, noise, clamp,
                                       return_params = FALSE) {

  # Create unique combinations for intercept parameters
  if (length(intercept_cols) > 0) {
    intercept_grid <- lapply(factors[intercept_cols], function(n) factor(1:n))
    names(intercept_grid) <- intercept_cols
    intercept_combos <- do.call(expand.grid, c(intercept_grid, stringsAsFactors = FALSE))
  } else {
    intercept_combos <- data.frame(dummy = 1, stringsAsFactors = FALSE)
  }

  # Create unique combinations for rate parameters
  if (length(rate_cols) > 0) {
    rate_grid <- lapply(factors[rate_cols], function(n) factor(1:n))
    names(rate_grid) <- rate_cols
    rate_combos <- do.call(expand.grid, c(rate_grid, stringsAsFactors = FALSE))
  } else {
    rate_combos <- data.frame(dummy = 1, stringsAsFactors = FALSE)
  }

  # Parse per-group composition targets (each sums to A_max) and the anchor SD,
  # then draw one Dirichlet composition per group.
  A0_target <- parse_parameter_spec(param_spec=params$A0_target, param_name="A0_target", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)
  B0_target <- parse_parameter_spec(param_spec=params$B0_target, param_name="B0_target", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)
  C0_target <- parse_parameter_spec(param_spec=params$C0_target, param_name="C0_target", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)
  A0_sd     <- parse_parameter_spec(param_spec=params$A0_target_sd, param_name="A0_target_sd", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)

  targets <- cbind(A0_target, B0_target, C0_target)
  validate_composition_targets(targets, A_max, "parallel")
  drawn <- draw_simplex(targets, A0_sd, A_min, A_max)
  A0_values <- drawn[, 1]
  B0_values <- drawn[, 2]
  C0_values <- drawn[, 3]

  intercept_params <- data.frame(
    A0 = A0_values,
    B0 = B0_values,
    C0 = C0_values,
    stringsAsFactors = FALSE
  )

  if (length(intercept_cols) > 0) {
    intercept_params <- cbind(intercept_combos, intercept_params)
  }

  # Parse per-group rate medians and geometric SDs, then draw lognormal rate
  # constants.
  krefB_median <- parse_parameter_spec(param_spec=params$krefB_median, param_name="krefB_median", factor_combos=rate_combos,
                                       relevant_factors=rate_cols, factors=factors)
  krefB_gsd    <- parse_parameter_spec(param_spec=params$krefB_gsd, param_name="krefB_gsd", factor_combos=rate_combos,
                                       relevant_factors=rate_cols, factors=factors)
  krefC_median <- parse_parameter_spec(param_spec=params$krefC_median, param_name="krefC_median", factor_combos=rate_combos,
                                       relevant_factors=rate_cols, factors=factors)
  krefC_gsd    <- parse_parameter_spec(param_spec=params$krefC_gsd, param_name="krefC_gsd", factor_combos=rate_combos,
                                       relevant_factors=rate_cols, factors=factors)

  krefB_values <- draw_rate_lognormal(krefB_median, krefB_gsd)
  krefC_values <- draw_rate_lognormal(krefC_median, krefC_gsd)

  # Activation energies vary by rate_cols like the rate constants, but are
  # deterministic.
  EaB_values <- parse_parameter_spec(param_spec=params$EaB, param_name="EaB", factor_combos=rate_combos,
                                     relevant_factors=rate_cols, factors=factors)
  EaC_values <- parse_parameter_spec(param_spec=params$EaC, param_name="EaC", factor_combos=rate_combos,
                                     relevant_factors=rate_cols, factors=factors)

  rate_params <- data.frame(
    Ea_B = EaB_values,
    Ea_C = EaC_values,
    kref_B = krefB_values,
    kref_C = krefC_values,
    stringsAsFactors = FALSE
  )

  if (length(rate_cols) > 0) {
    rate_params <- cbind(rate_combos, rate_params)
  }

  # Join parameters with design
  sim_data <- design

  if (length(intercept_cols) > 0) {
    sim_data <- dplyr::left_join(sim_data, intercept_params, by = intercept_cols)
  } else {
    sim_data$A0 <- intercept_params$A0[1]
    sim_data$B0 <- intercept_params$B0[1]
    sim_data$C0 <- intercept_params$C0[1]
  }

  if (length(rate_cols) > 0) {
    sim_data <- dplyr::left_join(sim_data, rate_params, by = rate_cols)
  } else {
    sim_data$Ea_B <- rate_params$Ea_B[1]
    sim_data$Ea_C <- rate_params$Ea_C[1]
    sim_data$kref_B <- rate_params$kref_B[1]
    sim_data$kref_C <- rate_params$kref_C[1]
  }

  # Calculate kinetics based on order and trend
  sim_data <- sim_data %>%
    dplyr::mutate(
      tempk = temp + 273.15,
      k_B = kref_B * exp((-Ea_B / R) * (1/tempk - 1/trefk)),
      k_C = kref_C * exp((-Ea_C / R) * (1/tempk - 1/trefk)),
      k_total = k_B + k_C,
      rate_time = k_total * time,
      True_A = if (kinetic_order == 0) {
        A0 - rate_time
      } else if (kinetic_order == 1) {
        A0 * exp(-rate_time)
      } else if (kinetic_order == 2) {
        (A0^(1 - kinetic_order) - (1 - kinetic_order) * rate_time)^(1 / (1 - kinetic_order))
      } else {
        NA_real_
      },
      Degraded_A = A0 - True_A,
      # Handle division by zero if k_total is 0
      True_B = B0 + ifelse(k_total != 0, (k_B / k_total) * Degraded_A, 0),
      True_C = C0 + ifelse(k_total != 0, (k_C / k_total) * Degraded_A, 0)
    )

  # Add measurement error.
  n_obs <- nrow(sim_data)
  eps <- draw_noise(n_obs, k = 3, sigma = params$sigma_err, noise = noise)

  sim_data <- sim_data %>%
    dplyr::mutate(
      Observed_A = True_A + eps[, 1],
      Observed_B = True_B + eps[, 2],
      Observed_C = True_C + eps[, 3]
    )

  # Optionally clamp values to the physical bounds [A_min, A_max].
  if (clamp) {
    n_clamped <- sum(
      sim_data$Observed_A < A_min | sim_data$Observed_A > A_max,
      sim_data$Observed_B < A_min | sim_data$Observed_B > A_max,
      sim_data$Observed_C < A_min | sim_data$Observed_C > A_max
    )
    sim_data <- sim_data %>%
      dplyr::mutate(
        Observed_A = pmin(pmax(Observed_A, A_min), A_max),
        Observed_B = pmin(pmax(Observed_B, A_min), A_max),
        Observed_C = pmin(pmax(Observed_C, A_min), A_max)
      )
    if (n_clamped > 0) {
      warning(sprintf(
        "clamp = TRUE censored %d observation(s) that fell outside [%g, %g].",
        n_clamped, A_min, A_max
      ))
    }
  }

  # Set default component names if not provided
  if (is.null(component_names)) {
    component_names <- list(
      comp_A = "Total Main",
      comp_B = "HMWS",
      comp_C = "LMWS"
    )
  }

  # Build a named vector for renaming: names = new names, values = old names
  rename_map <- stats::setNames(
    c("Observed_A", "Observed_B", "Observed_C"),
    c(component_names$comp_A, component_names$comp_B, component_names$comp_C)
  )

  # Select and format output
  select_cols <- c(names(factors), "time", "temp", "tempk")

  output_data <- sim_data %>%
    # 1) Select fixed columns + the original observed columns
    dplyr::select(
      dplyr::all_of(select_cols),
      "Observed_A", "Observed_B", "Observed_C"
    ) %>%
    # 2) Rename to user-provided names using the named vector
    dplyr::rename(!!!rename_map) %>%
    # 3) Pivot the renamed columns
    tidyr::pivot_longer(
      cols = c(
        component_names$comp_A,
        component_names$comp_B,
        component_names$comp_C
      ),
      names_to = "cqa",
      values_to = "y"
    )

  if (return_params) {
    param_list <- list(
      intercept_params = intercept_params,
      rate_params = rate_params,
      fixed_params = list(
        trefk = trefk,
        R = R,
        kinetic_order = kinetic_order,
        trend = trend,
        Yinf = Yinf,
        sigma_err = params$sigma_err
      )
    )

    return(list(
      data = output_data,
      parameters = param_list
    ))
  } else {
    return(output_data)
  }
}



#### Single Degradation Reaction ####

#' Helper function for single degradation reactions (A -> B)
#'
#' This internal function simulates a single degradation reaction where a main
#' component A degrades into a single product B (or, with
#' \code{trend = "increasing"}, A accumulates). The kinetics follow the
#' modified Arrhenius equation with temperature-dependent rate constants.
#'
#' @param design Data frame with experimental design including factor levels,
#'   time points, and temperatures.
#' @param factors List of all factors with their levels.
#' @param intercept_cols Character vector of factor names that affect intercepts
#'   (initial CQA values).
#' @param rate_cols Character vector of factor names that affect rate constants.
#' @param trefk Numeric. Reference temperature in Kelvin.
#' @param kinetic_order Numeric. Order of kinetics (0, 1, or 2).
#' @param trend Character. \code{"decreasing"} or \code{"increasing"}.
#' @param Yinf Numeric. Asymptotic bound for increasing trajectories.
#' @param params List. Model parameters including A0_target, B0_target,
#'   A0_target_sd, Ea, kref_median, kref_gsd, and sigma_err.
#' @param component_names Named list. Optional component names for output.
#' @param R Numeric. Universal gas constant in J/(mol-K).
#' @param A_max Numeric. The upper bound for the CQAs.
#' @param A_min Numeric. The lower bound for CQAs.
#' @param noise Character. Measurement-error model, \code{"mass_balance"} or
#'   \code{"independent"} (see \code{draw_noise()}).
#' @param clamp Logical. Whether to clamp observed values to \code{[A_min, A_max]}.
#' @param return_params Logical. If TRUE, returns both data and parameters.
#'   If FALSE, returns only data. Default is FALSE.
#'
#' @return If \code{return_params = FALSE}: A data frame in long format with
#'   simulated observations.
#'
#'   If \code{return_params = TRUE}: A list with two elements:
#'   \describe{
#'     \item{data}{Data frame in long format with columns for factors, time,
#'       temp, tempk, cqa (component name), and y (observed value)}
#'     \item{parameters}{List containing:
#'       \itemize{
#'         \item \code{intercept_params}: Data frame with columns for factor
#'           combinations (if applicable), A0, and B0
#'         \item \code{rate_params}: Data frame with columns for factor
#'           combinations (if applicable), Ea, and kref
#'         \item \code{fixed_params}: List with trefk, R, kinetic_order,
#'           trend, Yinf, and sigma_err
#'       }
#'     }
#'   }
#'
#' @details
#' The function implements a single degradation reaction where component A
#' degrades into B with rate constant k. Initial values (A0, B0) are drawn from
#' a two-part Dirichlet (Beta) distribution so they are non-negative and sum to
#' \code{A_max}. The rate constant is drawn from a lognormal, so
#' it is strictly positive.
#'
#' \code{True_B} is derived by mass balance as \eqn{B_0 + (A_0 - True_A)},
#' matching the equation used in \code{\link{fit_singledeg}}
#' (\eqn{Y_B = B_0 + (A_0 - Y_A)}). Measurement error is then added per the
#' \code{noise} model (see \code{draw_noise()}). Under \code{"mass_balance"} the
#' two species' errors are row-centered so they sum to zero, preserving
#' \eqn{A + B = A_{max}} exactly while keeping each species' marginal SD equal to
#' \code{sigma_err}. Under \code{"independent"} each species gets an independent
#' \eqn{N(0, \sigma^2)} error, so the sum is only approximately conserved.
#'
#' @keywords internal
#' @noRd
simulate_single_deg_kinetics <- function(design, factors, intercept_cols,
                                         rate_cols, trefk, kinetic_order,
                                         trend, Yinf,
                                         params, component_names, R,
                                         A_max, A_min, noise, clamp,
                                         return_params = FALSE) {

  # Create unique combinations for intercept parameters
  if (length(intercept_cols) > 0) {
    intercept_grid <- lapply(factors[intercept_cols], function(n) factor(1:n))
    names(intercept_grid) <- intercept_cols
    intercept_combos <- do.call(expand.grid, c(intercept_grid, stringsAsFactors = FALSE))
  } else {
    intercept_combos <- data.frame(dummy = 1, stringsAsFactors = FALSE)
  }

  # Create unique combinations for rate parameters
  if (length(rate_cols) > 0) {
    rate_grid <- lapply(factors[rate_cols], function(n) factor(1:n))
    names(rate_grid) <- rate_cols
    rate_combos <- do.call(expand.grid, c(rate_grid, stringsAsFactors = FALSE))
  } else {
    rate_combos <- data.frame(dummy = 1, stringsAsFactors = FALSE)
  }

  # Parse per-group composition targets (A0+B0 = A_max) and the anchor SD, then
  # draw one 2-part Dirichlet (Beta) composition per group.
  A0_target <- parse_parameter_spec(param_spec=params$A0_target, param_name="A0_target", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)
  B0_target <- parse_parameter_spec(param_spec=params$B0_target, param_name="B0_target", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)
  A0_sd     <- parse_parameter_spec(param_spec=params$A0_target_sd, param_name="A0_target_sd", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)

  targets <- cbind(A0_target, B0_target)
  validate_composition_targets(targets, A_max, "single_deg")
  drawn <- draw_simplex(targets, A0_sd, A_min, A_max)
  A0_values <- drawn[, 1]
  B0_values <- drawn[, 2]

  intercept_params <- data.frame(
    A0 = A0_values,
    B0 = B0_values,
    stringsAsFactors = FALSE
  )

  if (length(intercept_cols) > 0) {
    intercept_params <- cbind(intercept_combos, intercept_params)
  }

  # Parse per-group rate median and geometric SD, then draw lognormal rate
  # constants.
  kref_median <- parse_parameter_spec(param_spec=params$kref_median, param_name="kref_median", factor_combos=rate_combos,
                                      relevant_factors=rate_cols, factors=factors)
  kref_gsd    <- parse_parameter_spec(param_spec=params$kref_gsd, param_name="kref_gsd", factor_combos=rate_combos,
                                      relevant_factors=rate_cols, factors=factors)
  kref_values <- draw_rate_lognormal(kref_median, kref_gsd)

  # Activation energy varies by rate_cols like the rate constant, but is
  # deterministic.
  Ea_values <- parse_parameter_spec(param_spec=params$Ea, param_name="Ea", factor_combos=rate_combos,
                                    relevant_factors=rate_cols, factors=factors)

  rate_params <- data.frame(
    Ea = Ea_values,
    kref = kref_values,
    stringsAsFactors = FALSE
  )

  if (length(rate_cols) > 0) {
    rate_params <- cbind(rate_combos, rate_params)
  }

  # Join parameters with design
  sim_data <- design

  if (length(intercept_cols) > 0) {
    sim_data <- dplyr::left_join(sim_data, intercept_params, by = intercept_cols)
  } else {
    sim_data$A0 <- intercept_params$A0[1]
    sim_data$B0 <- intercept_params$B0[1]
  }

  if (length(rate_cols) > 0) {
    sim_data <- dplyr::left_join(sim_data, rate_params, by = rate_cols)
  } else {
    sim_data$Ea <- rate_params$Ea[1]
    sim_data$kref <- rate_params$kref[1]
  }

  # Calculate kinetics based on order and trend
  sim_data <- sim_data %>%
    dplyr::mutate(
      tempk = temp + 273.15,
      k = kref * exp((-Ea / R) * (1/tempk - 1/trefk)),
      rate_time = k * time,
      # kinetic_order and trend are scalars, so pick the (vectorised) closed form
      # with if/else. case_when() here would trigger the dplyr >= 1.2.0 warning
      # about size-1 LHS with size >1 RHS.
      True_A = if (kinetic_order == 0) {
        if (trend == "decreasing") A0 - rate_time else A0 + rate_time
      } else if (kinetic_order == 1) {
        if (trend == "decreasing") A0 * exp(-rate_time)
        else Yinf - (Yinf - A0) * exp(-rate_time)
      } else if (kinetic_order == 2) {
        if (trend == "decreasing")
          (A0^(1 - kinetic_order) - (1 - kinetic_order) * rate_time)^(1 / (1 - kinetic_order))
        else
          Yinf - ((Yinf - A0)^(1 - kinetic_order) - (1 - kinetic_order) * rate_time)^(1 / (1 - kinetic_order))
      } else {
        NA_real_
      },
      # Mass balance: B absorbs whatever A loses (or gives up)
      True_B = B0 + (A0 - True_A)
    )

  # Add measurement error.
  n_obs <- nrow(sim_data)
  eps <- draw_noise(n_obs, k = 2, sigma = params$sigma_err, noise = noise)

  sim_data <- sim_data %>%
    dplyr::mutate(
      Observed_A = True_A + eps[, 1],
      Observed_B = True_B + eps[, 2]
    )

  # Optionally clamp values to the physical bounds [A_min, A_max].
  if (clamp) {
    n_clamped <- sum(
      sim_data$Observed_A < A_min | sim_data$Observed_A > A_max,
      sim_data$Observed_B < A_min | sim_data$Observed_B > A_max
    )
    sim_data <- sim_data %>%
      dplyr::mutate(
        Observed_A = pmin(pmax(Observed_A, A_min), A_max),
        Observed_B = pmin(pmax(Observed_B, A_min), A_max)
      )
    if (n_clamped > 0) {
      warning(sprintf(
        "clamp = TRUE censored %d observation(s) that fell outside [%g, %g].",
        n_clamped, A_min, A_max
      ))
    }
  }

  # Set default component names if not provided
  if (is.null(component_names)) {
    component_names <- list(
      comp_A = "NR Major Peak",
      comp_B = "NR Minor Peaks"
    )
  }

  # Build a named vector for renaming: names = new names, values = old names
  rename_map <- stats::setNames(
    c("Observed_A", "Observed_B"),
    c(component_names$comp_A, component_names$comp_B)
  )

  # Select and format output
  select_cols <- c(names(factors), "time", "temp", "tempk")

  output_data <- sim_data %>%
    # 1) Select fixed columns + the original observed columns
    dplyr::select(
      dplyr::all_of(select_cols),
      "Observed_A", "Observed_B"
    ) %>%
    # 2) Rename to user-provided names using the named vector
    dplyr::rename(!!!rename_map) %>%
    # 3) Pivot the renamed columns
    tidyr::pivot_longer(
      cols = c(
        component_names$comp_A,
        component_names$comp_B
      ),
      names_to = "cqa",
      values_to = "y"
    )

  if (return_params) {
    param_list <- list(
      intercept_params = intercept_params,
      rate_params = rate_params,
      fixed_params = list(
        trefk = trefk,
        R = R,
        kinetic_order = kinetic_order,
        trend = trend,
        Yinf = Yinf,
        sigma_err = params$sigma_err
      )
    )

    return(list(
      data = output_data,
      parameters = param_list
    ))
  } else {
    return(output_data)
  }
}



#### Generic Reaction ####

#' Helper function for generic kinetic reactions (no specific degradant naming)
#'
#' This internal function simulates a generic kinetic reaction tracking only the
#' CQA value of component A over time. Unlike the other reaction types, this
#' does not explicitly model degradation products. The reaction follows the
#' modified Arrhenius equation with temperature-dependent rate constants.
#' Both decreasing and increasing trends are supported.
#'
#' @param design Data frame with experimental design including factor levels,
#'   time points, and temperatures.
#' @param factors List of all factors with their levels.
#' @param intercept_cols Character vector of factor names that affect intercepts
#'   (initial CQA values).
#' @param rate_cols Character vector of factor names that affect rate constants.
#' @param trefk Numeric. Reference temperature in Kelvin.
#' @param kinetic_order Numeric. Order of kinetics (0, 1, or 2).
#' @param trend Character. \code{"decreasing"} or \code{"increasing"}.
#' @param Yinf Numeric. Asymptotic bound for increasing trajectories.
#' @param params List. Model parameters including A0_target, A0_target_sd, Ea,
#'   kref_median, kref_gsd, and sigma_err.
#' @param component_names Named list. Optional component names for output.
#' @param R Numeric. Universal gas constant in J/(mol-K).
#' @param A_max Numeric. The upper bound for the CQAs.
#' @param A_min Numeric. The lower bound for CQAs.
#' @param noise Character. Measurement-error model, \code{"mass_balance"} or
#'   \code{"independent"} (see \code{draw_noise()}).
#' @param clamp Logical. Whether to clamp observed values to \code{[A_min, A_max]}.
#' @param return_params Logical. If TRUE, returns both data and parameters.
#'   If FALSE, returns only data. Default is FALSE.
#'
#' @return If \code{return_params = FALSE}: A data frame in long format with
#'   simulated observations.
#'
#'   If \code{return_params = TRUE}: A list with two elements:
#'   \describe{
#'     \item{data}{Data frame in long format with columns for factors, time,
#'       temp, tempk, cqa (component name), and y (observed value)}
#'     \item{parameters}{List containing:
#'       \itemize{
#'         \item \code{intercept_params}: Data frame with columns for factor
#'           combinations (if applicable) and A0
#'         \item \code{rate_params}: Data frame with columns for factor
#'           combinations (if applicable), Ea, and kref
#'         \item \code{fixed_params}: List with trefk, R, kinetic_order,
#'           trend, Yinf, and sigma_err
#'       }
#'     }
#'   }
#'
#' @details
#' The function implements a generic kinetic model that tracks only component A
#' without explicitly modeling degradation products. This is useful for modeling
#' CQA changes over time when the specific degradation pathways are
#' not of interest or when modeling non-degradation processes.
#'
#' Initial CQA value A0 is drawn from a scaled Beta distribution on
#' \code{[A_min, A_max]}, so it stays within the physical bounds by construction.
#' The rate constant is drawn from a lognormal, so it is strictly positive by
#' construction. Measurement error is added to the observed CQA value.
#'
#' @keywords internal
#' @noRd
simulate_generic_kinetics <- function(design, factors, intercept_cols,
                                      rate_cols, trefk, kinetic_order,
                                      trend, Yinf,
                                      params, component_names, R,
                                      A_max, A_min, noise, clamp,
                                      return_params = FALSE) {

  # Create unique combinations for intercept parameters
  if (length(intercept_cols) > 0) {
    intercept_grid <- lapply(factors[intercept_cols], function(n) factor(1:n))
    names(intercept_grid) <- intercept_cols
    intercept_combos <- do.call(expand.grid, c(intercept_grid, stringsAsFactors = FALSE))
  } else {
    intercept_combos <- data.frame(dummy = 1, stringsAsFactors = FALSE)
  }

  # Create unique combinations for rate parameters
  if (length(rate_cols) > 0) {
    rate_grid <- lapply(factors[rate_cols], function(n) factor(1:n))
    names(rate_grid) <- rate_cols
    rate_combos <- do.call(expand.grid, c(rate_grid, stringsAsFactors = FALSE))
  } else {
    rate_combos <- data.frame(dummy = 1, stringsAsFactors = FALSE)
  }

  # Parse per-group A0 target and anchor SD, then draw one scaled-Beta A0 per
  # group on [A_min, A_max].
  A0_target <- parse_parameter_spec(param_spec=params$A0_target, param_name="A0_target", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)
  A0_sd     <- parse_parameter_spec(param_spec=params$A0_target_sd, param_name="A0_target_sd", factor_combos=intercept_combos,
                                    relevant_factors=intercept_cols, factors=factors)

  A0_values <- draw_simplex(cbind(A0_target), A0_sd, A_min, A_max)[, 1]

  intercept_params <- data.frame(
    A0 = A0_values,
    stringsAsFactors = FALSE
  )

  if (length(intercept_cols) > 0) {
    intercept_params <- cbind(intercept_combos, intercept_params)
  }

  # Parse per-group rate median and geometric SD, then draw lognormal rate
  # constants.
  kref_median <- parse_parameter_spec(param_spec=params$kref_median, param_name="kref_median", factor_combos=rate_combos,
                                      relevant_factors=rate_cols, factors=factors)
  kref_gsd    <- parse_parameter_spec(param_spec=params$kref_gsd, param_name="kref_gsd", factor_combos=rate_combos,
                                      relevant_factors=rate_cols, factors=factors)
  kref_values <- draw_rate_lognormal(kref_median, kref_gsd)

  # Activation energy varies by rate_cols like the rate constant, but is
  # deterministic.
  Ea_values <- parse_parameter_spec(param_spec=params$Ea, param_name="Ea", factor_combos=rate_combos,
                                    relevant_factors=rate_cols, factors=factors)

  rate_params <- data.frame(
    Ea = Ea_values,
    kref = kref_values,
    stringsAsFactors = FALSE
  )

  if (length(rate_cols) > 0) {
    rate_params <- cbind(rate_combos, rate_params)
  }

  # Join parameters with design
  sim_data <- design

  if (length(intercept_cols) > 0) {
    sim_data <- dplyr::left_join(sim_data, intercept_params, by = intercept_cols)
  } else {
    sim_data$A0 <- intercept_params$A0[1]
  }

  if (length(rate_cols) > 0) {
    sim_data <- dplyr::left_join(sim_data, rate_params, by = rate_cols)
  } else {
    sim_data$Ea <- rate_params$Ea[1]
    sim_data$kref <- rate_params$kref[1]
  }

  # Calculate kinetics based on order and trend
  sim_data <- sim_data %>%
    dplyr::mutate(
      tempk = temp + 273.15,
      k = kref * exp((-Ea / R) * (1/tempk - 1/trefk)),
      rate_time = k * time,
      # kinetic_order and trend are scalars, so pick the (vectorised) closed form
      # with if/else. case_when() here would trigger the dplyr >= 1.2.0 warning
      # about size-1 LHS with size >1 RHS.
      True_A = if (kinetic_order == 0) {
        if (trend == "decreasing") A0 - rate_time else A0 + rate_time
      } else if (kinetic_order == 1) {
        if (trend == "decreasing") A0 * exp(-rate_time)
        else Yinf - (Yinf - A0) * exp(-rate_time)
      } else if (kinetic_order == 2) {
        if (trend == "decreasing")
          (A0^(1 - kinetic_order) - (1 - kinetic_order) * rate_time)^(1 / (1 - kinetic_order))
        else
          Yinf - ((Yinf - A0)^(1 - kinetic_order) - (1 - kinetic_order) * rate_time)^(1 / (1 - kinetic_order))
      } else {
        NA_real_
      }
    )

  # Add measurement error.
  n_obs <- nrow(sim_data)
  eps <- draw_noise(n_obs, k = 1, sigma = params$sigma_err, noise = noise)

  sim_data <- sim_data %>%
    dplyr::mutate(
      Observed_A = True_A + eps[, 1]
    )

  # Optionally clamp values to the physical bounds [A_min, A_max].
  if (clamp) {
    n_clamped <- sum(
      sim_data$Observed_A < A_min | sim_data$Observed_A > A_max
    )
    sim_data <- sim_data %>%
      dplyr::mutate(
        Observed_A = pmin(pmax(Observed_A, A_min), A_max)
      )
    if (n_clamped > 0) {
      warning(sprintf(
        "clamp = TRUE censored %d observation(s) that fell outside [%g, %g].",
        n_clamped, A_min, A_max
      ))
    }
  }

  # Set default component name if not provided
  if (is.null(component_names)) {
    component_names <- list(comp_A = "Purity")
  }

  # Build a named vector for renaming: names = new names, values = old names
  rename_map <- stats::setNames(
    c("Observed_A"),
    c(component_names$comp_A)
  )

  # Select and format output
  select_cols <- c(names(factors), "time", "temp", "tempk")

  output_data <- sim_data %>%
    # 1) Select fixed columns + the original observed columns
    dplyr::select(
      dplyr::all_of(select_cols),
      "Observed_A"
    ) %>%
    # 2) Rename to user-provided names using the named vector
    dplyr::rename(!!!rename_map) %>%
    # 3) Pivot the renamed columns
    tidyr::pivot_longer(
      cols = c(
        component_names$comp_A
      ),
      names_to = "cqa",
      values_to = "y"
    )

  if (return_params) {
    param_list <- list(
      intercept_params = intercept_params,
      rate_params = rate_params,
      fixed_params = list(
        trefk = trefk,
        R = R,
        kinetic_order = kinetic_order,
        trend = trend,
        Yinf = Yinf,
        sigma_err = params$sigma_err
      )
    )

    return(list(
      data = output_data,
      parameters = param_list
    ))
  } else {
    return(output_data)
  }
}
