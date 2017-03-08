#' Mean Doble Station
#' @author Cesar Aybar <aybar1994@gmail.com>
#' @param gauge Is an object of SpatialPointsDataFrame-class.
#' @param cov Is an object of RasterLayer.
#' @param longlat Using WGS84 system
#' @examples
#' library(raster)
#' library(Dorado)
#' data("Dorado")
#' gauge <- mean_doble_Station(gauge = Dorado$gauge,cov = Dorado$TRMM)
#' @importFrom automap autofitVariogram
#' @importFrom raster extract projection writeRaster  rasterToPoints
#' @importFrom sp coordinates spDists
#' @importFrom gstat krige.cv idw krige
#' @export


mean_doble_Station = function (gauge = NULL, cov = NULL, longlat = TRUE) {
  #----------------------------------------------------------------------
  projection(gauge) <- projection(cov)
  point <- data.frame(rasterToPoints(cov))
  colnames(point) <- c("x","y","cov")
  coordinates(point) <-  ~ x + y
  projection(point) <- projection(cov)

  distances <-
    function(x,ptsat = point)
      which.min(spDists(ptsat,gauge[x,],longlat = T))
  #----------------------------------------------------------------------
  loc <- do.call("c",lapply(1:length(gauge), distances))
  duplicates <- loc[which(duplicated(loc))]
  gauge2 <- cbind(coordinates(gauge),gauge@data)
  colnames(gauge2) <- c("x","y",colnames(gauge2)[3:length(gauge2)])
  gauge2p <- gauge2 %>% data.frame
  list<-list()
  for(i in 1:length(duplicates)){
    dupliStation <- which(loc==duplicates[i])
    gaugeD <- gauge2p[dupliStation,]
    PromStation <- colMeans(x = gaugeD,na.rm=T)
    list[[i]]<-list(Prom=PromStation,position=dupliStation)
  }
  stat<-do.call("rbind",lapply(1:length(duplicates),function(x) list[[x]]$position))
  stat2<-do.call("rbind",lapply(1:length(duplicates),function(x) list[[x]]$Prom))
  gauge2p[stat[,1],]<-stat2
  newG<-gauge2p[-stat[,2],]
  coordinates(newG)<-~x+y
  projection(newG) <- projection(cov)
  return(newG)
}



