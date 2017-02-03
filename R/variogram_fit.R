#' Variogram Estimate based in RasterLayer and sp Objects
#' @author Cesar Aybar <aybar1994@gmail.com>
#' @param gauge Is a object of SpatialPointsDataFrame-class.
#' @param cov Is a object of RasterLayer.
#' @param formula that defines the dependent variable as a linear model
#' of independent variables; suppose the dependent variable has
#' name 'z', for Kriging with Extrenal Drift (KED) use the formula
#' 'z~x+y+....', you do not need define
#' @param ... parameters that are passed on to \link[gstat]{variogram} variogram when calculating the sample variogram
#' @importFrom raster extract
#' @seealso \link[gstat]{fit.variogram} \link[automap]{autofitVariogram}
#' @examples
#' library(raster)
#' data(Dorado)
#' variogram_fit(gauge = Dorado$rain,cov = stack(Dorado$cov),rain~prec+dem)
#' @export
#'
variogram_fit <- function(gauge, cov, formula, ...) {
  ext <- raster::extract(cov, gauge, cellnumber = F, sp = T)
  list(ftvariogram = FitVariogram(formula, ext, ...), ext = ext)
}

