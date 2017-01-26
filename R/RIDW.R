#' Regression-IDW optimizado
#' @author Cesar Aybar
#' @param gauge An object of SpatialPointsDataFrame-class.
#' @param sat An object of RasterLayer.
#' @param idpR Is the range of the power coeficient
#' @importFrom automap autofitVariogram
#' @importFrom raster extract projection writeRaster
#' @importFrom sp coordinates
#' @importFrom gstat krige.cv idw
#' @export
RegressionIDW <- function(gauge = dperu2, sat = rSAT,
                          idpR = seq(0.8, 3.5,0.1), saveparam, saverast) {
  ext <- raster::extract(sat, gauge, cellnumber = F, sp = T)
  names(ext) <- c("gauge", "sat")
  station <- gauge
  llm <- lm(ext$gauge ~ ext$sat)
  station$residuals <- llm$residuals

  # Define Grid
  # -------------------------------------------------------------

  point <- rasterToPoints(sat)
  point <- as.data.frame(point)
  colnames(point) <- c("x", "y", "sat")
  coordinates(point) <- ~x + y
  projection(point) <- projection(sat)
  names(ext) <- c("gauge", "sat")

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

  OBSp <- sat * llm$coefficients[2] + llm$coefficients[1]
  Ridw <- OBSp + mapa
  Ridw[Ridw < 0] <- 0

  # Save Data
  # ---------------------------------------------------------------
  writeRaster(Ridw, filename = saverast, overwrite = TRUE)
  estimatepowerweigth <- bestparam
  RMSD <- sqrt(mean(residual.best^2))
  save(estimatepowerweigth, residual.best, RMSD, file = saveparam)
}
