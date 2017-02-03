#' Extract Dorado data
#' @author Cesar Aybar <aybar1994@gmail.com>
#' @param data Data to extract
#' @param output is the folder directory to save the data.
#' @examples
#' data(Dorado)
#' data(Titicaca)
#' #give_me_the_data(Titicaca,output=getwd())
#' @importFrom raster writeRaster
#' @importFrom rgdal writeOGR
#' @export
give_me_the_data<-function(data="Titicaca",output="~/"){
dir.create(output)
setwd(output)
shp_folder=paste0(output,"/shapefile")
raster_folder=paste0(output,"/raster")
dir.create(shp_folder)
dir.create(raster_folder)
setwd(shp_folder)
x=eval(parse(text = data))
writeOGR(x$rain, ".",names(x$rain), driver="ESRI Shapefile")
setwd(raster_folder)
mapply(function(i) writeRaster(x$cov[[i]],paste0(names(x$cov[[i]]),".tif")),1:length(x$cov))
}
