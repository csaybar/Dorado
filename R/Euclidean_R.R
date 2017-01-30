#'Estimate Euclidean Distance for sp Objects
#'@description This function implement a quick method for estimate the euclidean distance using
#'a RasterLayer object for defined the boundaries.
#'@author Cesar Aybar <aybar1994@gmail.com>
#'@param sp: Is a object sp (\code{SpatialPolygonsDataFrame}, \code{SpatialLineDataFrame}, \code{SpatialPointsDataFrame}).
#'@param raster: Is a RasterLayer object.
#'@details This function use the raster only for extract the coordinates referent system and geometries attributes,
#'the values are not used.
#'@importFrom raster shapefile raster mask distance writeRaster
#'@export
Euclidean_R <- function(sp, raster) {
  raster[raster != 0] <- 0
  mk <- mask(raster, sp)
  mk[mk == 0] <- 1
  euclidean <- distance(mk)
  return(euclidean)
}
