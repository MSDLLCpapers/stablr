#' Calculate Tolerance Intervals for gsl_nls Models
#'
#' Internal helper function to compute tolerance intervals following
#' Krishnamoorthy & Mathew (2009) for nonlinear regression models fitted
#' with gslnls.
#'
#' @param model A gsl_nls model object
#' @param newdata_list Data for prediction (in list format for gradient1)
#' @param coverage Coverage proportion (gamma)
#' @param level Confidence level (1-alpha)
#' @param side "two", "lower", or "upper"
#'
#' @return Data frame with fit, lwr, upr columns
#' @keywords internal
#' @noRd
compute_tolerance_intervals <- function(model, newdata, coverage, level, side) {
  
  predict_gsl_nls <- getFromNamespace("predict.gsl_nls", "gslnls")
  
  # Get point predictions and Jacobian
  mu_hat <- as.numeric(predict_gsl_nls(model, newdata))
  
  # Get gradient (Jacobian) at new points
  Fdot <- model$m$gradient1(newdata)
  
  # Get variance-covariance matrix: (F^T F)^{-1}
  Finv <- tryCatch(
    chol2inv(model$m$Rmat()),
    error = function(e) {
      stop("Singular gradient matrix at parameter estimates. ",
           "Tolerance intervals cannot be computed.", call. = FALSE)
    }
  )
  
  # Compute d^2 for each prediction point
  d_squared <- rowSums(Fdot %*% Finv * Fdot)
  
  # K&M helper function
  tolerance_bounds_km(
    mu_hat    = mu_hat,
    d_squared = d_squared,
    sigma     = sigma(model),
    df        = df.residual(model),
    coverage  = coverage,
    level     = level,
    side      = side
  )
}

#' Predict function for the fit_parallel model
#'
#' Computes predictions, confidence intervals, prediction intervals, or tolerance
#' intervals for objects fitted with \code{\link{fit_parallel}}. For tolerance
#' intervals, implements the noncentral t and Wallis approximation methods from
#' Krishnamoorthy & Mathew (2009).
#'
#' @param object Model output from stablr::fit_parallel function
#' @param newdata New dataframe containing time and tempk columns to be used for prediction. If NULL, the data used for model training will be used. Default is NULL.
#' @param interval Interval to calculate. Default is 'none' and the function will return the predictions as 'fit'. Other options are 'confidence' for a confidence interval containing 'fit', 'lwr', and 'upr', 'prediction' for a prediction interval containing 'fit', 'lwr', and 'upr', or 'tolerance' for tolerance intervals.
#' @param level Confidence level for the interval (default 0.95).
#' @param coverage Coverage proportion (\eqn{\gamma}) for tolerance intervals.
#'   The proportion of the population to be covered. Default is 0.90.
#'   Only used when \code{interval = "tolerance"}.
#' @param side Side for tolerance interval: \code{"two"} (default), \code{"lower"},
#'   or \code{"upper"}. Only used when \code{interval = "tolerance"}. Two-sided
#'   intervals use the Wallis (1951) approximation (Equation 3.3.6 in K&M);
#'   one-sided intervals use the noncentral t-distribution (Equation 3.2.6 in K&M).
#' @param ... Additional arguments (currently ignored; present for S3
#'   method consistency with the \code{predict} generic).
#'
#' @details
#' For \code{interval = "none"}, \code{"confidence"}, or \code{"prediction"},
#' this function delegates to \code{gslnls:::predict.gsl_nls()}, which implements
#' the Delta Method for nonlinear regression.
#'
#' For \code{interval = "tolerance"}, this function computes tolerance intervals
#' following Krishnamoorthy & Mathew (2009), adapted for nonlinear regression
#' via the Delta Method. The tolerance interval covers at least proportion
#' \eqn{\gamma} (specified by \code{coverage}) of the population with confidence
#' level \eqn{1-\alpha} (specified by \code{level}).
#'
#' Two-sided tolerance intervals (\code{side = "two"}) use the Wallis (1951)
#' approximation with noncentral chi-square distributions. One-sided tolerance
#' intervals (\code{side = "lower"} or \code{side = "upper"}) use the noncentral
#' t-distribution method.
#'
#' See \code{vignette("intervals-calculation", package = "stablr")} for detailed
#' mathematical background.
#'
#' @return Dataframe containing predictions for yA, yB, and yC at tempk and time values specified in newdata
#'
#' @references
#' Krishnamoorthy, K., & Mathew, T. (2009). \emph{Statistical Tolerance Regions:
#' Theory, Applications, and Computation}. John Wiley & Sons.
#'
#' @export
predict.fit_parallel <- function(object, newdata = NULL,
                                  interval = c("none", "confidence", "prediction", "tolerance"),
                                  level = 0.95, coverage = 0.90, 
                                  side = c("two", "lower", "upper"),
                                 ...
                                ){
  
  interval <- match.arg(interval)
  side <- match.arg(side)
  if (!is.null(newdata)) {
    build <- attr(object, "build")
    if (is.null(build)) {
      stop("Model object missing 'build' attribute. Was it fitted with fit_parallel()?", call. = FALSE)
    }
    
    settings <- build$settings
    # covar_cols covers rh_col and any second covariate; rh_col stays in the
    # union for fits recorded before it existed.
    required_vars <- unique(c(settings$time_var, settings$temp_var,
                              settings$rh_col, settings$covar_cols))
    missing_vars <- setdiff(required_vars, names(newdata))
    if (length(missing_vars)) {
      stop("newdata missing required variables: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    # Check NA in a predictor column (time/temp/rh)
    check_structural_nas(newdata, required_vars, what = "newdata")

    # Type and humidity-scale checks, the contract the fitter enforced.
    check_predictor_cols(newdata, settings, what = "newdata")

    # A newdata column named after a parameter outranks the estimate during
    # evaluation; see check_newdata_param_shadow().
    check_newdata_param_shadow(newdata, build$param_names)

    # Warns rather than errors: extrapolation is legitimate, just not silent.
    check_temp_extrapolation(newdata, settings)

    # ── Apply the same sanitization to newdata ──────────────────
    newdata <- apply_sanitize_map(
      newdata,
      sanitize_map = attr(object, "sanitize_map"),
      settings     = settings
    )
  }
  
  predict_gsl_nls <- getFromNamespace("predict.gsl_nls", "gslnls")

  # gsl_nls's fn/data interface needs newdata wrapped as list(data = ...).
  newdata_wrapped <- if (!is.null(newdata)) list(data = newdata) else NULL

  # ── Tolerance Intervals ────────────────────────────────────────────
  if (interval == "tolerance") {
    preds <- compute_tolerance_intervals(object, newdata_wrapped, coverage, level, side)
  } else {
    # ── Standard Intervals (none, confidence, prediction) ──────────────
    
    if(interval == "none" & !is.null(newdata)){
      preds <- predict_gsl_nls(object, newdata_wrapped)
    }
    else if(interval=="none" & is.null(newdata)){
      preds <- predict_gsl_nls(object)
    }
    else if (!is.null(newdata) && interval %in% c("confidence", "prediction")){
      preds <- predict_gsl_nls(object, newdata = newdata_wrapped, interval = interval, level = level)
    }
    else if (is.null(newdata) && interval %in% c("confidence", "prediction")){
      preds <- predict_gsl_nls(object, interval = interval, level = level)
    }
    
    preds <- as.data.frame(preds)
  }
  
  # ── Separate A, B, and C predictions ────────────────────────────────
  # Calculate n from predictions (3 responses: yA, yB, yC)
  n <- ifelse(!is.null(newdata), nrow(newdata),(nrow(preds) / 3))
  preds_A <- preds[seq_len(n), , drop = FALSE]
  preds_B <- preds[n + seq_len(n), , drop = FALSE]
  preds_C <- preds[2 * n + seq_len(n), , drop = FALSE]
  
  if(interval == "none"){
    preds_A <- as.data.frame(preds_A)
    preds_B <- as.data.frame(preds_B)
    preds_C <- as.data.frame(preds_C)
    colnames(preds_A) <- c("fit_yA")
    colnames(preds_B) <- c("fit_yB")
    colnames(preds_C) <- c("fit_yC")
  }
  else if(interval == "confidence" | interval == "prediction" | interval == "tolerance"){
    colnames(preds_A) <- c("fit_yA","lwr_yA","upr_yA")
    colnames(preds_B) <- c("fit_yB","lwr_yB","upr_yB")
    colnames(preds_C) <- c("fit_yC","lwr_yC","upr_yC")
  }
  
  # preds dataframe
  preds <- cbind(preds_A, preds_B, preds_C)
  
  return(preds)
}

#' Predict function for the fit_singledeg model
#'
#' Computes predictions, confidence intervals, prediction intervals, or tolerance
#' intervals for objects fitted with \code{\link{fit_singledeg}}. For tolerance
#' intervals, implements the noncentral t and Wallis approximation methods from
#' Krishnamoorthy & Mathew (2009).
#'
#' @param object Model output from stablr::fit_singledeg function
#' @param newdata New dataframe containing time and tempk columns to be used for prediction. If NULL, data used for model training will be used for prediction. Default is NULL. 
#' @param interval Interval to calculate. Default is 'none' and the function will return the predictions as 'fit'. Other options are 'confidence' for a confidence interval containing 'fit', 'lwr', and 'upr', 'prediction' for a prediction interval containing 'fit', 'lwr', and 'upr', or 'tolerance' for tolerance intervals.
#' @param level Confidence level for the interval (default 0.95).
#' @param coverage Coverage proportion (\eqn{\gamma}) for tolerance intervals.
#'   The proportion of the population to be covered. Default is 0.90.
#'   Only used when \code{interval = "tolerance"}.
#' @param side Side for tolerance interval: \code{"two"} (default), \code{"lower"},
#'   or \code{"upper"}. Only used when \code{interval = "tolerance"}. Two-sided
#'   intervals use the Wallis (1951) approximation (Equation 3.3.6 in K&M);
#'   one-sided intervals use the noncentral t-distribution (Equation 3.2.6 in K&M).
#' @param ... Additional arguments (currently ignored; present for S3
#'   method consistency with the \code{predict} generic).
#'
#' @details
#' For \code{interval = "none"}, \code{"confidence"}, or \code{"prediction"},
#' this function delegates to \code{gslnls:::predict.gsl_nls()}, which implements
#' the Delta Method for nonlinear regression.
#'
#' For \code{interval = "tolerance"}, this function computes tolerance intervals
#' following Krishnamoorthy & Mathew (2009), adapted for nonlinear regression
#' via the Delta Method. The tolerance interval covers at least proportion
#' \eqn{\gamma} (specified by \code{coverage}) of the population with confidence
#' level \eqn{1-\alpha} (specified by \code{level}).
#'
#' Two-sided tolerance intervals (\code{side = "two"}) use the Wallis (1951)
#' approximation with noncentral chi-square distributions. One-sided tolerance
#' intervals (\code{side = "lower"} or \code{side = "upper"}) use the noncentral
#' t-distribution method.
#'
#' See \code{vignette("intervals-calculation", package = "stablr")} for detailed
#' mathematical background.
#'
#' @return Dataframe containing predictions for yA and yB at tempk and time values specified in newdata
#'
#' @references
#' Krishnamoorthy, K., & Mathew, T. (2009). \emph{Statistical Tolerance Regions:
#' Theory, Applications, and Computation}. John Wiley & Sons.
#'
#' @examples
#' data("nrcesds")
#' library(tidyr)
#' nrcesds_wide <- pivot_wider(nrcesds, names_from = "cqa", values_from = "y")
#' fit <- fit_singledeg(
#'   df = nrcesds_wide,
#'   kinetic_order = 2,
#'   trend = "decreasing",
#'   yA_var = "NR Major Peak",
#'   yB_var = "NR Minor Peaks"
#' )
#'
#' # Point predictions
#' predict(fit)
#'
#' # Confidence intervals
#' predict(fit, interval = "confidence", level = 0.95)
#'
#' # Prediction intervals
#' predict(fit, interval = "prediction", level = 0.95)
#'
#' # Tolerance intervals (two-sided, 90% coverage, 95% confidence)
#' predict(fit, interval = "tolerance", coverage = 0.90, level = 0.95)
#'
#' # One-sided lower tolerance limit (90% coverage, 95% confidence)
#' predict(fit, interval = "tolerance", coverage = 0.90, level = 0.95, side = "lower")
#'
#' @export
predict.fit_singledeg <- function(object, newdata = NULL,
                                   interval = c("none", "confidence", "prediction", "tolerance"),
                                   level = 0.95, coverage = 0.90, 
                                   side = c("two", "lower", "upper"),
                                  ...){
  
  interval <- match.arg(interval)
  side <- match.arg(side)
  if (!is.null(newdata)) {
    build <- attr(object, "build")
    if (is.null(build)) {
      stop("Model object missing 'build' attribute. Was it fitted with fit_kinetic()?",
           call. = FALSE)
    }
    
    settings <- build$settings
    # covar_cols covers rh_col and any second covariate; rh_col stays in the
    # union for fits recorded before it existed.
    required_vars <- unique(c(settings$time_var, settings$temp_var,
                              settings$rh_col, settings$covar_cols))
    missing_vars <- setdiff(required_vars, names(newdata))
    if (length(missing_vars)) {
      stop("newdata missing required variables: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    # Check NA in a predictor column (time/temp/rh)
    check_structural_nas(newdata, required_vars, what = "newdata")

    # Type and humidity-scale checks, the contract the fitter enforced.
    check_predictor_cols(newdata, settings, what = "newdata")

    # A newdata column named after a parameter outranks the estimate during
    # evaluation; see check_newdata_param_shadow().
    check_newdata_param_shadow(newdata, build$param_names)

    # Warns rather than errors: extrapolation is legitimate, just not silent.
    check_temp_extrapolation(newdata, settings)

    # ── Apply the same sanitization to newdata ──────────────────
    newdata <- apply_sanitize_map(
      newdata,
      sanitize_map = attr(object, "sanitize_map"),
      settings     = settings
    )
  }
  
  predict_gsl_nls <- getFromNamespace("predict.gsl_nls", "gslnls")
  
  

  # gsl_nls's fn/data interface needs newdata wrapped as list(data = ...).
  newdata_wrapped <- if (!is.null(newdata)) list(data = newdata) else NULL

  # ── Tolerance Intervals ────────────────────────────────────────────
  if (interval == "tolerance") {
    preds <- compute_tolerance_intervals(object, newdata_wrapped, coverage, level, side)
  } else {
    # ── Standard Intervals (none, confidence, prediction) ──────────────
    
    if(interval == "none" & !is.null(newdata)){
      preds <- predict_gsl_nls(object, newdata_wrapped)
    }
    else if(interval=="none" & is.null(newdata)){
      preds <- predict_gsl_nls(object)
    }
    else if (!is.null(newdata) && interval %in% c("confidence", "prediction")){
      preds <- predict_gsl_nls(object, newdata = newdata_wrapped, interval = interval, level = level)
    }
    else if (is.null(newdata) && interval %in% c("confidence", "prediction")){
      preds <- predict_gsl_nls(object, interval = interval, level = level)
    }
    
    preds <- as.data.frame(preds)
  }
  
  # ── Separate A and B predictions ────────────────────────────────────
  # Calculate n from predictions (2 responses: yA, yB)
  n <- ifelse(!is.null(newdata), nrow(newdata),(nrow(preds) / 2))
  preds_A <- preds[seq_len(n), , drop = FALSE]
  preds_B <- preds[n + seq_len(n), , drop = FALSE]
  
  if(interval == "none"){
    preds_A <- as.data.frame(preds_A)
    preds_B <- as.data.frame(preds_B)
    colnames(preds_A) <- c("fit_yA")
    colnames(preds_B) <- c("fit_yB")
  }
  else if(interval == "confidence" | interval == "prediction" | interval == "tolerance"){
    colnames(preds_A) <- c("fit_yA","lwr_yA","upr_yA")
    colnames(preds_B) <- c("fit_yB","lwr_yB","upr_yB")
  }
  
  # Combine a and b predictions
  preds <- cbind(preds_A, preds_B)
  
  return(preds)
}


#' Predict Method for fit_kinetic Objects
#'
#' Computes predictions, confidence intervals, prediction intervals, or tolerance
#' intervals for objects fitted with \code{\link{fit_kinetic}}. For tolerance
#' intervals, implements the noncentral t and Wallis approximation methods from
#' Krishnamoorthy & Mathew (2009).
#'
#' @param object A \code{fit_kinetic} object (output from \code{\link{fit_kinetic}}).
#' @param newdata Optional data frame of new observations. If \code{NULL}, uses
#'   the training data.
#' @param interval Type of interval to compute: \code{"none"} (default),
#'   \code{"confidence"}, \code{"prediction"}, or \code{"tolerance"}.
#' @param level Confidence level for the interval (default 0.95).
#' @param coverage Coverage proportion (\eqn{\gamma}) for tolerance intervals.
#'   The proportion of the population to be covered. Default is 0.90.
#'   Only used when \code{interval = "tolerance"}.
#' @param side Side for tolerance interval: \code{"two"} (default), \code{"lower"},
#'   or \code{"upper"}. Only used when \code{interval = "tolerance"}. Two-sided
#'   intervals use the Wallis (1951) approximation (Equation 3.3.6 in K&M);
#'   one-sided intervals use the noncentral t-distribution (Equation 3.2.6 in K&M).
#' @param ... Additional arguments (currently ignored; present for S3
#'   method consistency with the \code{predict} generic).
#'
#' @details
#' For \code{interval = "none"}, \code{"confidence"}, or \code{"prediction"},
#' this function delegates to \code{gslnls:::predict.gsl_nls()}, which implements
#' the Delta Method for nonlinear regression.
#'
#' For \code{interval = "tolerance"}, this function computes tolerance intervals
#' following Krishnamoorthy & Mathew (2009), adapted for nonlinear regression
#' via the Delta Method. The tolerance interval covers at least proportion
#' \eqn{\gamma} (specified by \code{coverage}) of the population with confidence
#' level \eqn{1-\alpha} (specified by \code{level}).
#'
#' Two-sided tolerance intervals (\code{side = "two"}) use the Wallis (1951)
#' approximation with noncentral chi-square distributions. One-sided tolerance
#' intervals (\code{side = "lower"} or \code{side = "upper"}) use the noncentral
#' t-distribution method.
#'
#' The implementation uses:
#' \itemize{
#'   \item The Jacobian matrix \eqn{F} computed via \code{object$m$gradient1(newdata)}
#'   \item The standardized distance \eqn{d^2 = \nabla f_0^T(F^TF)^{-1}\nabla f_0}
#'   \item Noncentral t or noncentral chi-square distributions for the tolerance factor
#' }
#'
#' See \code{vignette("intervals-calculation", package = "stablr")} for detailed
#' mathematical background.
#'
#' @return
#' \itemize{
#'   \item If \code{interval = "none"}: numeric vector of fitted values
#'   \item Otherwise: data frame with columns \code{fit}, \code{lwr}, \code{upr}
#' }
#'
#' @references
#' Krishnamoorthy, K., & Mathew, T. (2009). \emph{Statistical Tolerance Regions:
#' Theory, Applications, and Computation}. John Wiley & Sons.
#'
#' @examples
#' data(evers_supp)
#' dat <- subset(evers_supp, cqa == "Purity" & !is.na(y))
#' fit <- fit_kinetic(dat, kinetic_order = 1, trend = "decreasing")
#'
#' # Point predictions
#' predict(fit)
#'
#' # Confidence intervals
#' predict(fit, interval = "confidence", level = 0.95)
#'
#' # Prediction intervals
#' predict(fit, interval = "prediction", level = 0.95)
#'
#' # Tolerance intervals (two-sided, 90% coverage, 95% confidence)
#' predict(fit, interval = "tolerance", coverage = 0.90, level = 0.95)
#'
#' # One-sided lower tolerance limit (90% coverage, 95% confidence)
#' predict(fit, interval = "tolerance", coverage = 0.90, level = 0.95, side = "lower")
#'
#' @export
predict.fit_kinetic <- function(
    object,
    newdata = NULL,
    interval = c("none", "confidence", "prediction", "tolerance"),
    level = 0.95,
    coverage = 0.90,
    side = c("two", "lower", "upper"),
    ...) {
  
  interval <- match.arg(interval)
  side <- match.arg(side)
  
  # Validate and sanitize newdata
  if (!is.null(newdata)) {
    build <- attr(object, "build")
    if (is.null(build)) {
      stop("Model object missing 'build' attribute. Was it fitted with fit_kinetic()?",
           call. = FALSE)
    }
    
    settings <- build$settings
    # covar_cols covers rh_col and any second covariate; rh_col stays in the
    # union for fits recorded before it existed.
    required_vars <- unique(c(settings$time_var, settings$temp_var,
                              settings$rh_col, settings$covar_cols))
    missing_vars <- setdiff(required_vars, names(newdata))
    if (length(missing_vars)) {
      stop("newdata missing required variables: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    # Check NA in a predictor column (time/temp/rh) 
    check_structural_nas(newdata, required_vars, what = "newdata")

    # Type and humidity-scale checks, the contract the fitter enforced.
    check_predictor_cols(newdata, settings, what = "newdata")

    # A newdata column named after a parameter outranks the estimate during
    # evaluation; see check_newdata_param_shadow().
    check_newdata_param_shadow(newdata, build$param_names)

    # Warns rather than errors: extrapolation is legitimate, just not silent.
    check_temp_extrapolation(newdata, settings)

    # ── Apply the same sanitization to newdata ──────────────────
    newdata <- apply_sanitize_map(
      newdata,
      sanitize_map = attr(object, "sanitize_map"),
      settings     = settings
    )
  }
  
  # ── Delegate to gslnls for standard intervals ──────────────────────
  if (interval != "tolerance") {
    predict_gsl_nls <- getFromNamespace("predict.gsl_nls", "gslnls")
    
    if (interval == "none") {
      return(predict_gsl_nls(object = object, newdata = newdata))
    } else {
      return(predict_gsl_nls(
        object = object,
        newdata = newdata,
        interval = interval,
        level = level
      ))
    }
  }
  
  # ── Tolerance Intervals ────────────────────────────────────────────
  return(compute_tolerance_intervals(object, newdata, coverage, level, side))
}
