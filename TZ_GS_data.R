# Tanzania GeoSurvey 250m resolution data setup 
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

# cropland-focused GeoSurvey 2018
# see sampling frame @ https://github.com/mgwalsh/Sampling/blob/master/TZ_GS_sample.R
download("https://www.dropbox.com/s/0x4y4j6ifqidmhh/TZ_geos_2018.csv.zip?raw=1", "TZ_geos_2018.csv.zip", mode = "wb")
unzip("TZ_geos_2018.csv.zip", overwrite = T)
geos18 <- read.table("TZ_geos_2018.csv", header = T, sep = ",")
geos18$BIC <- as.factor(ifelse(geos18$CP == "Y" & geos18$BP == "Y", "Y", "N")) ## identifies croplands with buildings

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/bhefsc8u120uqwp/TZA_adm3.zip?raw=1", "TZA_adm3.zip", mode = "wb")
unzip("TZA_adm3.zip", overwrite = T)
shape <- shapefile("TZA_adm3.shp")

# download raster stack (note this is a big 860+ Mb download)
download("https://www.dropbox.com/s/pshrtvjf7navegu/TZ_250m_2017.zip?raw=1", "TZ_250m_2017.zip", mode = "wb")
unzip("TZ_250m_2017.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Baseline data setup -----------------------------------------------------
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

# 2018 cropland survey data setup -----------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(geos18) <- ~lon+lat
projection(geos18) <- projection(shape)
gadm <- geos18 %over% shape
geos18 <- as.data.frame(geos18)
geos18 <- cbind(gadm[ ,c(5,7,9)], geos18)
colnames(geos18) <- c("region","district","ward","survey","time","id","observer","lat","lon","BP","CP","WP","rice","bloc","cgrid","BIC")

# project GeoSurvey coords to grid CRS
geos18.proj <- as.data.frame(project(cbind(geos18$lon, geos18$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos18.proj) <- c("x","y")
geos18 <- cbind(geos18, geos18.proj)
coordinates(geos18) <- ~x+y
projection(geos18) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos18)
gsdat18 <- as.data.frame(cbind(geos18, geosgrid)) 
gsdat18 <- gsdat18[!duplicated(gsdat18), ] ## removes any duplicates 
gsdat18$observer <- sub("@.*", "", as.character(gsdat18$observer)) ## shortens observer ID's

# Write output files ------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/TZ_gsdat.csv", row.names = F) ## baseline
write.csv(gsdat18, "./Results/TZ_gsdat18.csv", row.names = F) ## 2018 cropland GS

# GeoSurvey map widget ----------------------------------------------------
# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'TZ_GS.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
# Baseline survey
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(3,0.1), random.order = T)

# 2018 cropland survey
gscon18 <- as.data.frame(table(gsdat18$observer))
set.seed(1235813)
wordcloud(gscon18$Var1, freq = gscon18$Freq, scale = c(3,0.1), random.order = T)

