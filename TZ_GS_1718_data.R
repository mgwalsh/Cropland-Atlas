# Tanzania combined GeoSurvey 2017/2018 250m resolution GS data setup
# M. Walsh, March 2019 (updated grids to 2019 version)

# Required packages
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(leaflet)
  require(htmlwidgets)
  require(wordcloud)
})

# Data downloads -----------------------------------------------------------
# download GeoSurvey data
download("https://osf.io/uhx9b?raw=1", "TZ_geos_1718.csv.zip", mode = "wb")
unzip("TZ_geos_1718.csv.zip", overwrite = T)
geos <- read.table("TZ_geos_1718.csv", header = T, sep = ",")

# download raster stack (note this is a big 1+ Gb download)
download("https://osf.io/ke5ya?raw=1", "TZ_250m_2019.zip", mode = "wb")
unzip("TZ_250m_2019.zip", overwrite = T)
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
gsdat <- gsdat[complete.cases(gsdat[ ,c(1:58)]),] ## removes incomplete cases
gsdat$observer <- sub("@.*", "", as.character(gsdat$observer)) ## shortens observer ID's

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/TZ_gsdat_1718.csv", row.names = F)

