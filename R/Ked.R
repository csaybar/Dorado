#' KED (Kriging with Extrenal Drift)
#' @author Cesar Aybar
#' @param gauge An object of SpatialPointsDataFrame-class.
#' @param cov An object of RasterLayer.
#' @param idpR Is the range of the power coeficient
#' @importFrom automap autofitVariogram
#' @importFrom raster extract projection writeRaster
#' @importFrom sp coordinates
#' @importFrom gstat krige.cv idw
#' @export

KED<-function(gauge,
              cov,
              boundariess,
              formula,
              fx = c(NA,NA,NA),
              saveparam="/home/senamhi-cesar/Escritorio/name.Rdata",
              saverast="/home/senamhi-cesar/Escritorio/name.tif"
){
  ext<-raster::extract(cov,gauge,cellnumber=F,sp=T)
  names(ext)<-c("gauge","cov")
  vm.fit<-Rvariogrm(gauge,cov,formula = formula)

  #------Define grid----------------------
  point<-readGDAL(cov@file@name)
  names(point)<-"cov"
  #---------------------------------------
  Zs<-krige(gauge~cov, locations = ext, newdata = point, model = vm.fit$var_model)
  Zs.cv<-krige.cv(gauge~cov, ext, nfold = nrow(ext),nmax = Inf)
  Zs.cvresidual<-Zs.cv["residual"]
  map 			<- as(Zs[1],"SpatialPixelsDataFrame")
  gridded(map) <- TRUE
  mapa <- raster(map); mapa[mapa<=0]=0
  save(vm.fit,Zs.cvresidual,file = saveparam)
  writeRaster(mapa,filename = saverast,overwrite=TRUE)
}

