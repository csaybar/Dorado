#'Resamplear raster en funcion a otro
#'@author Cesar Luis Aybar Camacho
#'@param R1  Es un objeto raster que se desea resamplear.
#'@param R2  Es un objeto raster de la resolucion deseada
#'@param t_res: Opcional parametro necesario si no se ingresa raster con que comparar.
#'@param t_ext: Opcional parametro necesario si no se ingresa raster con que comparar.
#'@param Routput: Es el carpeto de destino del raster
#'@importFrom raster projection xmin ymin xmax ymax res
#'@importFrom gdalUtils gdalwarp
#'@examples
#' #Cuando tienes un raster con que comparar.
#'resampleR(R1=Prec_Bolivia,R2=dem,Routput = "~/")
#' #Cuando no tienes un raster con que comparar.
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
