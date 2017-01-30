#'Fast Resample RasterLayer object
#'@description This function allows manipulating the geometry of one RasterLayer object.
#'@author Cesar Aybar <aybar1994@gmail.com>
#'@param R1  Is a RasterLayer object to resample.
#'@param R2  Is a RasterLayer object that contains the CRS and desired geometry.
#'@param t_res: Numeric. Is the pixel resolution (c(xres,yres))
#'@param t_ext: Numeric.(c(xmin,ymin,xmax,ymax)). set georeferenced extents of output file to be created (in target SRS).
#'@param Routput: Character. The destination file name.
#'@details
#'This function use \code{\link[gdalUtils]{gdalwarp}} for resampling RasterLayer object and
#'implement the two next possible situations. The first is that you have two RasterLayer objects
#'and you want to change the geometries and SRC of the first using the second. The second
#'possibility is that you have one RasterLayer object and you know the pixel resolution
#'and the extents.
#'@return A RasterLayer object
#'@importFrom raster projection xmin ymin xmax ymax res
#'@importFrom gdalUtils gdalwarp
#'@examples
#' #First Case
#'resampleR(R1=Prec_Bolivia,R2=dem,Routput = "~/")
#' #Second Case
#'resampleR(R1=Prec_Bolivia,t_res=c(0.1,0.1),t_ext=c(-71,-65,-23,-13),Routput = "~/")
#'@export

resampleR <- function(R1 = NULL, R2 = NULL, t_res = NULL, t_ext = NULL,
                        Routput = "~/") {
  if (!missing(R2)) {
    if (projection(R1) != projection(R2)) {
      message("RASTERS CON DIFERENTES PROYECCIONES")
      #----------------------------------------------------------------
      t1 <- c(xmin(R2), ymin(R2), xmax(R2), ymax(R2))  #coordenadas de la imagen destino
      t2 <- c(res(R2)[1], res(R2)[2])  #resolucion de la imagen destino
      namess <- paste0(.Internal(dirname(R1@file@name)), "/R.", names(R1),
                       ".tif")  #NAME
      #----------------------------------------------------------------

      gdalwarp(R1@file@name, dstfile = namess, t_srs = R2@crs)

      gdalwarp(namess, dstfile = Routput, tr = t2, te = t1, output_Raster = T,
               overwrite = T, verbose = T)

      file.remove(namess)

    } else {
      #----------------------------------------------------------------
      t1 <- c(xmin(R2), ymin(R2), xmax(R2), ymax(R2))  #coordenadas de la imagen destino
      t2 <- c(res(R2)[1], res(R2)[2])  #resolucion de la imagen destino
      #----------------------------------------------------------------

      gdalwarp(R1@file@name, dstfile = Routput, tr = t2, te = t1,
               output_Raster = T, overwrite = T, verbose = T)
    }
  } else {

    gdalwarp(srcfile = R1@file@name, dstfile = Routput, tr = t_res, te = t_ext, output_Raster = T,
             overwrite = T, verbose = T)
    }
  }
