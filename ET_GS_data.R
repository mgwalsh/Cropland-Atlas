# Ethiopia GeoSurvey 250m resolution data setup 
# M. Walsh, March 2018

# Required packages
# install.packages(c("downloader","rgdal","raster","caret")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("ET_GS250", showWarnings=F)
setwd("./ET_GS250")

# download current Ethiopia GeoSurvey data
download("https://www.dropbox.com/s/mfw02vsrrit674z/ET_geos_2017.csv.zip?raw=1", "ET_geos_2017.csv.zip", mode="wb")
unzip("ET_geos_2017.csv.zip", overwrite=T)
geos <- read.table("ET_geos_2017.csv", header=T, sep=",")

# download Ethiopia 250m Gtifs and stack in raster (note this is a big 950+ Mb download)
download("https://www.dropbox.com/s/iqix6sn66w04jo0/ET_250m_2017.zip?raw=1", "ET_250m_2017.zip", mode="wb")
unzip("ET_250m_2017.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)
gsdat <- as.data.frame(cbind(geos, geosgrid)) 
gsdat <- na.omit(gsdat) ## includes only complete cases
gsdat <- gsdat[!duplicated(gsdat), ] ## removes any duplicates 

# Write output file -------------------------------------------------------
dir.create("Results", showWarnings=F)
write.csv(gsdat, "./Results/ET_gsdat.csv", row.names = FALSE)
