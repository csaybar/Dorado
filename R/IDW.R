#' IDW Optimizing inverse distance weighting power
#' @description This function use gstat packages for interpolate spatial point data (\link[sp]{SpatialPoints})
#' and RasterLayer data (see \link[raster]{raster}).
#' @author Cesar Aybar <aybar1994@gmail.com>
#' @seealso \link[gstat]{idw}
#' @param gauge Is an object of SpatialPointsDataFrame class.
#' @param idpR Is vector numeric of the power coeficient to evaluate.
#' @param cov Is An object of RasterLayer.
#' @param formula that defines the dependent variable as a linear model
#' of independent variables; suppose the dependent variable has
#' name 'z', for Regression Inverse Distance Weigthing (RIDW) use the formula
#' 'z~x+y+....', you do not need define.
#' @details R_IDW use crossvalidation Leave-p-out cross-validation (LpO CV) and force brute (optimize MSE)
#'  for estimate the best idp power coeficient.
#' @return a List that contains: \code{Interpol} is the RIDW result in Raster,
#'  \code{params} being \code{bestp} is the best distance weighting power,
#'  \code{MSE} is the Residual Mean squared error of the residuals and
#'   finally \code{linear_Model} is  the adjusted linear Model.
#' @examples
#'  library(raster)
#'  data(Titicaca)
#'  x <- IDW(gauge = Titicaca$rain,cov = Titicaca$cov$long,formula = rain~1)
#'  plot(x$Interpol)
#' @importFrom automap autofitVariogram
#' @importFrom dplyr %>% tbl_df mutate_all
#' @importFrom raster extract projection writeRaster stack nlayers
#' @importFrom sp coordinates gridded
#' @importFrom gstat krige.cv idw
#' @importFrom methods as is
#' @importFrom stats lm median na.omit
#' @export
#'
IDW <- function(gauge,cov,formula, idpR = seq(0.8, 3.5, 0.1)) {
  # Define Grid -------------------------------------------------------------
  point <- rasterToPoints(cov) %>% data.frame
  coordinates(point) <- ~x + y
  projection(point) <- projection(cov)

  # Estimate Best Parameter -------------------------------------------------

  idpRange <- idpR
  mse <- rep(NA, length(idpRange))
  for (i in 1:length(idpRange)) {
    mse[i] <- mean(krige.cv(as.formula(formula), gauge, nfold = nrow(gauge),
                            nmax = Inf, set = list(idp = idpRange[i]), verbose = F)$residual^2)
  }

  poss <- which(mse %in% min(mse))
  bestparam <- idpRange[poss]
  residual.best <- krige.cv(as.formula(formula), gauge, nfold = nrow(gauge),
                            nmax = Inf, set = list(idp = idpRange[poss]), verbose = F)$residual

  # Interpolation ----------------------------------------------------------
  idwError <- idw(as.formula(formula), gauge, point, idp = bestparam)
  idwError <- idwError["var1.pred"]
  gridded(idwError) <- TRUE
  mapa <- raster(idwError)
  # Save Data ---------------------------------------------------------------
  list(Interpol = mapa, params = list(bestp = bestparam, MSE = mean(residual.best^2)))
}
