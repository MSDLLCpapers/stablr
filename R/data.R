#' Purity and HMWS data from Evers et al. 2022 publication in Pharmaceutics
#'
#' HPLC purity and SEC high‑molecular‑weight (HMW) stability data for the therapeutic peptide SAR441255, measured for different combinations of three formulations and two primary packaging materials over time at various temperatures.
#'
#' @format A data frame with 630 rows and 9 variables:
#' \itemize{
#'   \item \code{batch} - The batch ID; combination of peptide formulation and primary packaging materials
#'   \item \code{form} - Peptide formulation
#'   \item \code{pkg} - Primary packaging material
#'   \item \code{temp} - Storage temperature in Celsius
#'   \item \code{tempk} - Storage temperature in Kelvin
#'   \item \code{timed} - Timepoints CQA measurements were taken (days)
#'   \item \code{time} - Timepoints CQA measurements were taken (months)
#'   \item \code{y} - CQA value
#'   \item \code{cqa} - Critical Quality Attribute
#' }
#'
#' @source
#' Evers, A., Clénet, D., & Pfeiffer-Marek, S. (2022). Long-Term Stability 
#' Prediction for Developability Assessment of Biopharmaceutics Using 
#' Advanced Kinetic Modeling. *Pharmaceutics*, 14(2), 375. 
#' \doi{10.3390/pharmaceutics14020375}.
#'
#' The data is available under the Creative Commons Attribution (CC BY 4.0)
#' license. This dataset is a subset of the original supplementary material.
#' @docType data
#' @usage data("evers_supp")
#' @keywords datasets
"evers_supp"

#' Simulated UPSEC Example Data
#'
#' A dataset containing simulated values mimicking the UPSEC assay for demonstration purposes. This data was generated using a parallel first-order reactions model.
#'
#' @format A data frame with 990 rows and 6 variables:
#' \itemize{
#'   \item \code{batch} - The batch ID; an integer identifier from 1 to 6
#'   \item \code{time} - Timepoints CQA measurements were taken (months)
#'   \item \code{temp} - Storage temperature in Celsius
#'   \item \code{tempk} - Storage temperature in Kelvin
#'   \item \code{cqa} - Critical Quality Attribute
#'   \item \code{y} - CQA value
#' }
#' 
#' @docType data
#' @usage data("upsec")
#' @keywords datasets
"upsec"

#' Simulated NR CE-SDS Example Data
#'
#' A dataset containing simulated values mimicking the NR CE-SDS assay for demonstration purposes. This data was generated using a single degradation second-order reaction model.
#'
#' @format A data frame with 660 rows and 6 variables:
#' \itemize{
#'   \item \code{batch} - The batch ID; an integer identifier from 1 to 6
#'   \item \code{time} - Timepoints CQA measurements were taken (months)
#'   \item \code{temp} - Storage temperature in Celsius
#'   \item \code{tempk} - Storage temperature in Kelvin
#'   \item \code{cqa} - Critical Quality Attribute
#'   \item \code{y} - CQA value
#' }
#' 
#' @docType data
#' @usage data("nrcesds")
#' @keywords datasets
"nrcesds"

#' Colour-change stability data with a relative humidity covariate
#'
#' Colour difference (\code{delta_E*}) measurements taken over time at four
#' elevated temperatures and ten relative humidity levels. Because storage
#' conditions vary in both temperature and humidity, this dataset is the example
#' used to demonstrate the humidity-dependent Arrhenius model built by
#' [expr_arrhenius_humidity()].
#'
#' @format A data frame with 125 rows and 6 variables:
#' \itemize{
#'   \item \code{temp} - Storage temperature in Celsius (60, 70, 80, 90)
#'   \item \code{tempk} - Storage temperature in Kelvin
#'   \item \code{time} - Timepoints CQA measurements were taken
#'   \item \code{rh} - Relative humidity, as a proportion (0 to 0.5)
#'   \item \code{cqa} - Critical Quality Attribute (\code{"delta_E*"})
#'   \item \code{y} - CQA value
#' }
#'
#' @docType data
#' @usage data("flavier_supp")
#' @keywords datasets
"flavier_supp"