#'Estimación de la distancia Euclideana en Raster
#'@author Cesar Luis Aybar Camacho
#'@param shp: Es un objeto de clase SpatialPolygonsDataFrame, SpatialLineDataFrame, SpatialPointsDataFrame.
#'@param rasterModel: Es un objeto clase Raster
#'@param output: Es el carpeto de destino
#'@importFrom raster shapefile raster mask distance writeRaster
#'@export
EuclideanRaster <- function(shp, rasterModel, output) {
  rasterModel[rasterModel != 0] <- 0
  mk <- mask(rasterModel, shp)
  mk[mk == 0] <- 1
  euclidean <- distance(mk)
  writeRaster(euclidean, output)
}
