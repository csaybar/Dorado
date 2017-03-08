#'Download CHIRPS
#'@author Cesar  Aybar
#'@param day Dia que se desea descargar
#'@param saveCHIRPS Lugar en donde se guardara el archivo
#'@param BBlonMin Crop-Box long min
#'@param BBlonMax Crop-Box long max
#'@param BBlatMin Crop-Box lat min
#'@param BBlatMax Crop-Box lat max
#'@param res CHIRP type 005 or 025.
#'@importFrom raster crop extent writeRaster res
#'@importFrom R.utils gunzip
#'@importFrom RCurl getURL
#'@importFrom R.utils gunzip
#'@importFrom utils download.file
#'@export
downloadCHIRPS <- function(res = "p05",
                           day = as.Date("2015-10-01"),
                           saveCHIRPS = "~/",
                           BBlonMin = -86,
                           BBlonMax = -66,
                           BBlatMin = -19.25,
                           BBlatMax = 1.25){
  url <- "ftp://ftp.chg.ucsb.edu/pub/org/chg/products/CHIRPS-2.0/global_daily/tifs/"
  year <- year(as.Date(day))
  month <- sprintf("%02d", month(as.Date(day)))
  num_month <- as.numeric(month)
  myURL <- sprintf("%s/%s/%s/", url, res, year)
  filenames <- getURL(myURL,
                             ftp.use.epsv = FALSE,
                             ftplistonly = TRUE,
                             crlf = TRUE)
  filePaths <-
    paste(myURL, strsplit(filenames, "\r*\n")[[1]], sep = "")
  selectedfilePaths <-
    filePaths[grep(filePaths, pattern = paste("\\.tif.gz$"))]
  select <-
    selectedfilePaths[grep(format(day, "%Y.%m.%d"), selectedfilePaths)]
  download.file(select, basename(select))
  rastecomprss <- paste0(getwd(), "/", basename(select))
  gunzip(rastecomprss)
  Rraster <-
    raster(paste0(getwd(), "/", gsub(".gz", "", basename(select))))
  crpPP <-
    crop(Rraster, extent(BBlonMin, BBlonMax, BBlatMin, BBlatMax))
  crpPP[crpPP < 0] = 0
  writeRaster(crpPP, paste0(saveCHIRPS,"/",names(crpPP),".tif"),overwrite=T)
  file.remove(gsub(".gz","",rastecomprss))
}
