# Tanzania GeoSurvey 250m resolution cropland data setup
# M. Walsh, April 2018

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
dir.create("TZ_GS18", showWarnings = F)
setwd("./TZ_GS18")

# download GeoSurvey data
# see sampling frame @ https://github.com/mgwalsh/Sampling/blob/master/TZ_GS_sample.R
download("https://www.dropbox.com/s/0x4y4j6ifqidmhh/TZ_geos_2018.csv.zip?raw=1", "TZ_geos_2018.csv.zip", mode = "wb")
unzip("TZ_geos_2018.csv.zip", overwrite = T)
geos <- read.table("TZ_geos_2018.csv", header = T, sep = ",")
geos$BIC <- as.factor(ifelse(geos$CP == "Y" & geos$BP == "Y", "Y", "N")) ## identifies croplands with buildings

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/bhefsc8u120uqwp/TZA_adm3.zip?raw=1", "TZA_adm3.zip", mode = "wb")
unzip("TZA_adm3.zip", overwrite = T)
shape <- shapefile("TZA_adm3.shp")

# download raster stack (note this is a big 1+ Gb download)
download("https://www.dropbox.com/s/y1ht9tazxqvsggp/TZ_250m_2018.zip?raw=1", "TZ_250m_2018.zip", mode = "wb")
unzip("TZ_250m_2018.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(5,7,9)], geos)
colnames(geos) <- c("region","district","ward","survey","time","id","observer","lat","lon","BP","CP","WP","rice","bloc","cgrid","BIC")

# Count number of buildings per quadrat -----------------------------------
bp <- geos[which(geos$BP == "Y"), ] ## identify quadrats with buildings
bp$bloc <- as.character(bp$bloc)

# Counting tagged building locations from quadrats with buildings
n <- rep(NA, nrow(bp))
for(i in 1:nrow(bp)) {
  t <- fromJSON(bp$bloc[i])
  n[i] <- nrow(t$features)
}
n ## vector of number of buildings per quadrats with buildings
ba <- geos[which(geos$BP == "N"), ]
ba$n <- 0
bp <- cbind(bp, n)
geos <- rbind(ba, bp)
geos <- geos[order(geos$id),] ## sort in original sample order

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
gsdat<- gsdat[complete.cases(gsdat[ ,c(10:13,18:61)]),] ## removes incomplete cases
gsdat$observer <- sub("@.*", "", as.character(gsdat$observer)) ## shortens observer ID's

# Write data frames -------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/TZ_gsdat18.csv", row.names = F)

# GeoSurvey map widget ----------------------------------------------------
# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'TZ_GS18.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(3,0.1), random.order = T)

