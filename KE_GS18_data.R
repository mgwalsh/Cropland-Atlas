# Kenya GeoSurvey 2018 250m resolution GS data setup
# M. Walsh, April 2020

# Required packages
# install.packages(c("downloader","rgdal","raster","leaflet","htmlwidgets","wordcloud")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(leaflet)
  require(htmlwidgets)
  require(wordcloud)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("KE_GS18", showWarnings = F)
setwd("./KE_GS18")

# download GeoSurvey data
download("https://osf.io/rx3ge?raw=1", "KE_geos_2018.csv.zip", mode = "wb")
unzip("KE_geos_2018.csv.zip", overwrite = T)
geos <- read.table("KE_geos_2018.csv", header = T, sep = ",")

# download raster stack (note this is a big 1+ Gb download)
download("https://osf.io/b2xnr?raw=1", "KE_250m_2018.zip", mode = "wb")
unzip("KE_250m_2018.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)
gsdat <- as.data.frame(cbind(geos, geosgrid)) 
# gsdat <- gsdat[!duplicated(gsdat), ] ## removes any duplicates ... if needed
gsdat <- gsdat[complete.cases(gsdat[ ,c(1:56)]),] ## removes incomplete cases
gsdat$observer <- sub("@.*", "", as.character(gsdat$observer)) ## shortens observer ID's

# Write data frame --------------------------------------------------------
write.csv(gsdat, "./Results/KE_gsdat_2018.csv", row.names = F)
