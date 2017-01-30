#' Regression IDW Optimizing inverse distance weighting power
#' @author Cesar Aybar
#' @description This function use gstat packages for interpolate spatial point data (see \code{\link[sp]{SpatialPointsDataFrame}} )
#' and RasterLayer data  (see \code{\link[raster]})
#' @seealso \link[gstat]{idw}
#' @param gauge Is an object of SpatialPointsDataFrame class.
#' @param newdata Is An object of RasterLayer.
#' @param idpR Is vector numeric of the power coeficient to evaluate.
#' @param formula that defines the dependent variable as a linear model
#' of independent variables; suppose the dependent variable has
#' name 'z', for ordinary and simple kriging use the formula
#' 'z~1'; for simple kriging also define 'beta' (see below); for
#' universal kriging, suppose 'z' is linearly dependent on 'x'
#' and 'y', use the formula 'z~x+y'.
#' @details R_IDW use crossvalidation Leave-p-out cross-validation (LpO CV) and force brute (optimize MSE)
#'  for estimate the best idp power coeficient.
#' @importFrom automap autofitVariogram
#' @importFrom raster extract projection writeRaster
#' @importFrom sp coordinates
#' @importFrom gstat krige.cv idw
#' @export
IDW_R <- function(gauge,
                  newdata,
                  formula,
                  idpR = seq(0.8, 3.5,0.1)
                  ){
  ext <- raster::extract(cov, gauge, cellnumber = F, sp = T)
  names(ext) <- c("gauge", "cov")
  station <- gauge
  llm <- lm(ext$gauge ~ ext$cov)
  station$residuals <- llm$residuals

  # Define Grid
  # -------------------------------------------------------------

  point <- rasterToPoints(cov)
  point <- as.data.frame(point)
  colnames(point) <- c("x", "y", "cov")
  coordinates(point) <- ~x + y
  projection(point) <- projection(cov)
  names(ext) <- c("gauge", "cov")

  # Estimate Best Parameter
  # -------------------------------------------------

  idpRange <- idpR
  mse <- rep(NA, length(idpRange))
  for (i in 1:length(idpRange)) {
    mse[i] <- mean(krige.cv(residuals ~ 1, station, nfold = nrow(station),
                            nmax = Inf, set = list(idp = idpRange[i]), verbose = F)$residual^2)
  }
  poss <- which(mse %in% min(mse))
  bestparam <- idpRange[poss]
  residual.best <- krige.cv(residuals ~ 1, station, nfold = nrow(station),
                            nmax = Inf, set = list(idp = idpRange[poss]), verbose = F)$residual

  # Interpolation
  # ----------------------------------------------------------

  idwError <- idw(residuals ~ 1, station, point, idp = bestparam)
  idwError <- idwError["var1.pred"]
  gridded(idwError) <- TRUE
  mapa <- raster(idwError)

  OBSp <- cov * llm$coefficients[2] + llm$coefficients[1]
  Ridw <- OBSp + mapa
  Ridw[Ridw < 0] <- 0

  # Save Data
  # ---------------------------------------------------------------
  list(Interpol=Ridw,params=list(bestp=bestparam,MSE=sqrt(mean(residual.best^2))))
}
