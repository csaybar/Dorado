library(Dorado)
library(raster)
library(sp)
library(testthat)

data(meuse)
data(meuse.grid)
coordinates(meuse.grid)<-~x+y
dist<-meuse.grid[3]
gridded(dist)=T
dist<-raster(dist)
coordinates(meuse)<-~x+y
cadmium<-meuse["cadmium"]

x<-RIDW(gauge = cadmium,cov = dist,formula = cadmium ~ dist)
expect_is(mean(x$Interpol@data@values,na.rm=T),"numeric")
expect_is(x$params$bestp,"numeric")
expect_is(x$Interpol,"RasterLayer")
