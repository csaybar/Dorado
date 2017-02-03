#' Altiplane Peru Bolivia Rainfall Information
#'
#' A dataset containing rainfall information and other spatial covariables.
#'  Dorado. The variables are as follows:
#' @format A list with 2 camps: rain and cov:
#' \itemize{
#'   \item rain: 72 rainfall station is a sp object.
#'   \item cov$long: longitude raster.
#'   \item cov$lat: latitude raster.
#'   \item cov$prec: rainfall raster.
#'   \item cov$dem: Digital Elevation Model raster.
#'   \item cov$dist: Distance to sea raster}
#' @examples
#' library(Dorado)
#' data("Dorado")
"Dorado"
