#' KED (Kriging with Extrenal Drift or Ordinary Kriging) using Raster and sp Objects
#' @author Cesar Aybar <aybar1994@gmail.com>
#' @param gauge Is an object of SpatialPointsDataFrame-class.
#' @param cov Is an object of RasterLayer.
#' @param formula that defines the dependent variable as a linear model
#' @param crossval logical, Generate crosvalidation
#' of independent variables; suppose the dependent variable has
#' name 'z', for Kriging with Extrenal Drift (KED) use the formula
#' 'z~x+y+....', you do not need define
#' @param model variogram model of dependent variable (or its residuals),
#' defined by a call to \code{\link[Dorado]{FitVariogram}}
#' @param ... parameters that are passed on to \link[gstat]{variogram} variogram when calculating the sample variogram
#' @return a List that contains: \code{Interpol} is the KED result in Raster,
#'  \code{params} being \code{residual} that are spatial residual obtains in the cross-validation (residual),
#'  \code{MSE} is the Residual Mean squared error and finally \code{var} is  the variogram model.
#' @seealso \link[gstat]{krige} \link[automap]{autoKrige}
#' @examples
#' library(raster)
#' library(Dorado)
#' data('Dorado')
#' #k = KED(gauge = Dorado$rain,cov = stack(Dorado$cov),formula = rain~prec+dem,crossval=T)
#' @importFrom automap autofitVariogram
#' @importFrom raster extract projection writeRaster
#' @importFrom sp coordinates
#' @importFrom gstat krige.cv idw krige
#' @export

KED <- function(gauge, cov, formula, model,crossval=F,...) {
  if(missing(model)){
    vm.fit <- variogram_fit(gauge, cov, formula)
    ext <- vm.fit$ext
    model <- vm.fit$ftvariogram
  }

  # Define Grid
  # -------------------------------------------------------------
  point <- rasterToPoints(cov) %>% data.frame
  coordinates(point) <- ~x + y
  projection(point) <- projection(cov)

  # Interpolation
  # ----------------------------------------------------------
  Zs <- krige(formula, locations = ext, newdata = point, model = model$var_model)
  map <- as(Zs[1], "SpatialPixelsDataFrame")
  gridded(map) <- TRUE
  mapa <- raster(map)
  plot(mapa)
  if(crossval==T){
    Zs.cv <- krige.cv(formula, ext,model$var_model, nfold = nrow(ext), nmax = Inf)
    Zs.cvresidual <- Zs.cv["residual"]
    list(Interpol = mapa, params = list(residual = Zs.cvresidual, MSE = mean(Zs.cvresidual$residual^2),
                                      var = vm.fit))
  } else  mapa
}
