#'Fast Extract RasterLayer data
#'@params R1 RasterLayer Object
#'@params shpR es un shapefile
#'@Routput Routput es el Raster de Salida
#'@ImportFrom gdalUtils gdalwarp
#'@export
Extract_Fast <- function(R1, shpR,Routput) {
  shape_cut <- try(shapefile(shpR))
  if (class(shape_cut) == "try-error")
    stop("ShpR no es un shapefile")
  te <- c(xmin(shape_cut), ymin(shape_cut), xmax(shape_cut), ymax(shape_cut))
  #Raster_Output_name <- sprintf("%s/R-%s", dirname(R1), basename(R1))
  gdalwarp(srcfile = R1, te = te, output_Raster = T, overwrite = T, verbose = T,
           cutline = shpR, dstfile = Routput)
}

#'Fast Extract RasterLayer data through shapefile with one geometry
#'@params R1 RasterLayer Object
#'@params shpR es un shapefile
#'@Routput Routput es el Raster de Salida
#'@params fun is a function to apply
#'@ImportFrom gdalUtils gdalwarp
#'@export
Extract_Fast_One <- function(R1, shpR,Routput, fun = mean) {
  shape_cut <- try(shapefile(shpR))
  if (class(shape_cut) == "try-error")
    stop("ShpR no es un shapefile")
  te <- c(xmin(shape_cut), ymin(shape_cut), xmax(shape_cut), ymax(shape_cut))
  gdalwarp(srcfile = R1, te = te, output_Raster = T, overwrite = T, verbose = T,
           cutline = shpR, dstfile = Routput)
     raster_Data <- na.omit(getValues(raster(Raster_Output_name)))
     if (is.list(fun))
      mapply(function(i) fun[[i]](raster_Data), 1:length(fun)) else fun(raster_Data)
}


#'Fast Extract RasterLayer data through shapefile with many geometries
#'@params R1 RasterLayer Object
#'@params shpR es un shapefile
#'@params Routput es el Raster de Salida
#'@params fun is a function to apply
#'@ImportFrom gdalUtils gdalwarp
#'@ImportFrom parallel detectCores
#'@ImportFrom doMC registerDoMC
#'@ImportFrom foreach foreach
Extract_Fast_Multiple<-function(R1,shpR,dirtrash,fun=list(mean,max,min)){
  dir.create(dirtrash)
  setwd(dirtrash)
  message("Creando shapefile a utilizar espere....")
  dirname_shp<-sprintf("%s/%s",dirtrash,"shapefiles")
  dir.create(dirname_shp)
  setwd(dirname_shp)
  shp_ogr<-shapefile(shpR)
  zero_len<-nchar(length(shp_ogr))
  zero_lenF<-paste0("%0",zero_len,"d")
  foreach(i=1:length(shp_ogr))%dopar%{
  writeOGR(shp_ogr[i,1],".",sprintf(paste0("%s/",zero_lenF),dirname_shp,i),"ESRI Shapefile")
  }
  lt_shapefile<-list.files(dirname_shp,pattern = "\\.shp$",full.names = T)
  message("Se creo los shp sin problemas cortando el Raster ......")
  dirname_raster<-sprintf("%s/%s",dirtrash,"rasters")
  dir.create(dirname_raster)
  Raster_Imput_name <- R1

  no_cores <- detectCores() - 2
  registerDoMC(no_cores)
  fun_Extract=function(i){
    Raster_Output_name<- paste0(dirname_raster,"/",i,".tif")
    Extract_Fast(Raster_Imput_name,shpR = lt_shapefile[i],fun = fun,Routput = Raster_Output_name)
    x=na.omit(getValues(raster(Raster_Output_name)))
    mapply(function(i) fun[[i]](x), 1:length(fun))
  }
  foreach(i=1:length(shp_ogr),.combine = cbind) %dopar%{
    fun_Extract(i)
    }
  }


#Extract_Fast2(R1,shpR,dirtrash)



600/60

