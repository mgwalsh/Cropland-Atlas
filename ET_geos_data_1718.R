# Ethiopia GeoSurvey 2017/2018 250m resolution data setup
# M. Walsh, January 2019

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
dir.create("ET_geos_1718", showWarnings = F)
setwd("./ET_geos_1718")

# download GeoSurvey data
download("https://www.dropbox.com/s/3t9khr81sk801pr/ET_geos_1718.csv.zip?raw=1", "ET_geos_1718.csv.zip", mode = "wb")
unzip("ET_geos_1718.csv.zip", overwrite = T)
geos <- read.table("ET_geos_1718.csv", header = T, sep = ",")
geos$BIC <- as.factor(ifelse(geos$CP == "Y" & geos$BP == "Y", "Y", "N")) ## identifies croplands with buildings

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/25kw13359f8l5tr/ET_adm_shp.zip?raw=1", "ET_adm_shp.zip", mode = "wb")
unzip("ET_adm_shp.zip", overwrite = T)
shape <- shapefile("ETH_adm3.shp")

# download raster stack (note this is a big 1+ Gb download)
download("https://www.dropbox.com/s/iqix6sn66w04jo0/ET_250m_2018.zip?raw=1", "ET_250m_2018.zip", mode = "wb")
unzip("ET_250m_2018.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(5,7,9)], geos)
colnames(geos) <- c("region","zone","woreda","survey","observer","lat","lon","BP","CP","WP","BIC")

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
gsdat <- gsdat[complete.cases(gsdat[ ,c(8:55)]),] ## removes incomplete cases
# gsdat <- gsdat[ which(gsdat$CP=='Y'), ] ## selects croplands only
gsdat$observer <- sub("@.*", "", as.character(gsdat$observer)) ## shortens observer ID's

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/ET_gsdat_1718.csv", row.names = F)

# GeoSurvey map widget ----------------------------------------------------
w <- leaflet() %>%
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'ET_geos_1718.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(2,0.1), random.order = T)
