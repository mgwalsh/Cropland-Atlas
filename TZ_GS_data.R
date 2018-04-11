# Tanzania GeoSurvey 250m resolution data setups
# M. Walsh, March 2018

# Required packages
# install.packages(c("downloader","rgdal","jsonlite","raster","leaflet","htmlwidgets","wordcloud")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(jsonlite)
  require(raster)
  require(leaflet)
  require(htmlwidgets)
  require(wordcloud)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("TZ_GS250", showWarnings = F)
setwd("./TZ_GS250")

# download GeoSurvey data
# GeoSurvey 2017 (baseline)
download("https://www.dropbox.com/s/94d68wrq93dj7te/TZ_geos_2017.csv.zip?raw=1", "TZ_geos_2017.csv.zip", mode = "wb")
unzip("TZ_geos_2017.csv.zip", overwrite = T)
geos <- read.table("TZ_geos_2017.csv", header = T, sep = ",")
geos$BIC <- as.factor(ifelse(geos$CP == "Y" & geos$BP == "Y", "Y", "N")) ## identifies croplands with buildings

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/bhefsc8u120uqwp/TZA_adm3.zip?raw=1", "TZA_adm3.zip", mode = "wb")
unzip("TZA_adm3.zip", overwrite = T)
shape <- shapefile("TZA_adm3.shp")

# download raster stack (note this is a big 860+ Mb download)
download("https://www.dropbox.com/s/pshrtvjf7navegu/TZ_250m_2017.zip?raw=1", "TZ_250m_2017.zip", mode = "wb")
unzip("TZ_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(5,7,9)], geos)
colnames(geos) <- c("region", "district", "ward", "survey", "observer", "lat", "lon", "BP", "CP", "WP", "BIC")

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
gsdat$observer <- sub("@.*", "", as.character(gsdat$observer)) ## shortens observer ID's

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/TZ_gsdat17.csv", row.names = F) ## 2017 baseline
write.csv(gsdat18, "./Results/TZ_gsdat18.csv", row.names = F) ## 2018 cropland GS

# GeoSurvey map widget ----------------------------------------------------
# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'TZ_GS.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
# 2017 baseline survey
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(3,0.1), random.order = T)

