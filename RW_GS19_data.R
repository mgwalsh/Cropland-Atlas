# Rwanda GeoSurvey 2019 250m resolution GS data setup
# M. Walsh, April 2019

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
dir.create("RW_GS19", showWarnings = F)
setwd("./RW_GS19")
dir.create("Results", showWarnings = F)

# download GeoSurvey data
download("https://www.dropbox.com/s/oqao51hxxvc09ec/RW_geos_2019.csv.zip?raw=1", "RW_geos_2019.csv.zip", mode = "wb")
unzip("RW_geos_2019.csv.zip", overwrite = T)
geos <- read.table("RW_geos_2019.csv", header = T, sep = ",")

# download GADM-L5 shapefile (courtesy of: http://www.gadm.org)
download("https://www.dropbox.com/s/fhusrzswk599crn/RWA_level5.zip?raw=1", "RWA_level5.zip", mode = "wb")
unzip("RWA_level5.zip", overwrite = T)
shape <- shapefile("gadm36_RWA_5.shp")

# download raster stack
download("https://osf.io/xts2y?raw=1", "RW_250m_2020.zip", mode = "wb")
unzip("RW_250m_2020.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L5 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(4,6,8,10,12)], geos)
colnames(geos) <- c("region","district","sector","cell", "village","time","observer","id","lat","lon","BP","CP","TP","WP","bloc","cgrid")

# Coordinates and number of buildings per quadrat -------------------------
bp <- geos[which(geos$BP == "Y"), ] ## identify quadrats with buildings
bp$bloc <- as.character(bp$bloc)

# coordinates of tagged building locations from quadrats with buildings
c <- fromJSON(bp$bloc[1])
bcoord <- do.call("rbind", c$feature$geometry$coordinates)
for(i in 2:nrow(bp)) {
  c <- fromJSON(bp$bloc[i])
  bcoord_temp <- do.call("rbind", c$feature$geometry$coordinates)
  bcoord <- rbind(bcoord, bcoord_temp)
}
bcoord <- as.data.frame(bcoord) ## vector of coordinates per quadrats with buildings
colnames(bcoord) <- c("lon","lat")

# number of tagged building locations from quadrats with buildings
bcount <- rep(NA, nrow(bp))
for(i in 1:nrow(bp)) {
  t <- fromJSON(bp$bloc[i])
  bcount[i] <- nrow(t$features)
}
bcount ## vector of number of buildings per quadrats with buildings
ba <- geos[which(geos$BP == "N"), ]
ba$bcount <- 0
bp <- cbind(bp, bcount)
geos <- rbind(ba, bp)
geos <- geos[order(geos$time),] ## sort in original sample order

# cropland grid count
cp <- geos[which(geos$CP == "Y"), ] ## identify quadrats with cropland
cp$cgrid <- as.character(cp$cgrid)

# number of tagged grid locations from quadrats with cropland
ccount <- rep(NA, nrow(cp))
for(i in 1:nrow(cp)) {
  t <- fromJSON(cp$cgrid[i])
  ccount[i] <- nrow(t$features)
}
ccount ## cropland grid count
ca <- geos[which(geos$CP == "N"), ]
ca$ccount <- 0
cp <- cbind(cp, ccount)
geos <- rbind(ca, cp)
geos <- geos[order(geos$time),] ## sort in original sample order

# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)
gsdat <- as.data.frame(cbind(geos, geosgrid))
gsdat <- gsdat[ which(gsdat$ccount < 17), ]

# Write data frame --------------------------------------------------------
write.csv(bcoord, "./Results/RW_bcoord.csv", row.names = F)
write.csv(gsdat, "./Results/RW_gsdat_2019.csv", row.names = F)

# GeoSurvey map widgets ---------------------------------------------------
# number of GeoSurvey quadrats
w <- leaflet() %>%
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 9) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'RW_GS19.html', selfcontained = T) ## save widget

# number of building tags
b <- leaflet() %>%
  setView(lng = mean(bcoord$lon), lat = mean(bcoord$lat), zoom = 9) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(bcoord$lon, bcoord$lat, clusterOptions = markerClusterOptions())
b ## plot widget 
saveWidget(b, 'RW_GS19_buildings.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(4,0.1), random.order = T)
