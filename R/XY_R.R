#' Estimate Longitude and Latitude RasterLayer object
#'@author Cesar Aybar <aybar1994@gmail.com>
#'@param x: RasterLayer object
#'@importFrom dplyr %>%
#'@importFrom raster rasterToPoints
#'@importFrom sp coordinates
#'@examples
#' data(Prec_Bolivia)
#' generateXYraster(Prec_Bolivia)
#'@export
XY_R <- function(x) {
  long <- rasterToPoints(x)[, c(1, 2, 1)] %>% data.frame
  coordinates(long) <- ~x + y
  gridded(long) = T
  long <- raster(long)
  lat <- rasterToPoints(x)[, c(1, 2, 2)] %>% data.frame
  coordinates(lat) <- ~x + y
  gridded(lat) = T
  lat <- raster(lat)
  return(list(lat = lat, long = long))
}




raster()
