# Mozambique GeoSurvey 2019 250m resolution Level-1 data setup
# M. Walsh, Jan 2020

# Required packages
# install.packages(c("downloader","rgdal","jsonlite","raster","leaflet","htmlwidgets","wordcloud"), dependencies=TRUE)
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
dir.create("MZ_GS19", showWarnings = F)
setwd("./MZ_GS19")

# download GeoSurvey data
download("https://osf.io/2evzk?raw=1", "MZ_geos_L1_2019.csv.zip", mode = "wb")
unzip("MZ_geos_L1_2019.csv.zip", overwrite = T)
geos <- read.table("MZ_geos_L1_2019.csv", header = T, sep = ",")

# download Mozambique GADM-L3 shapefile (courtesy of: http://www.gadm.org)
download("https://osf.io/yw4vr?raw=1", "MZ_GADM_L3.zip", mode = "wb")
unzip("MZ_GADM_L3.zip", overwrite = T)
shape <- shapefile("gadm36_MOZ_3.shp")

# download raster stack
download("https://osf.io/dp6ez$raw=1", "MZ_250m_2019.zip", mode = "wb")
unzip("MZ_250m_2019.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(4,7,10)], geos)
colnames(geos) <- c("province","district","locality","today","observer","id","lat","lon","BP","CP","WP","TP","bloc","cgrid")

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
geos <- geos[order(geos$today),] ## sort in original sample order

# Cropland grid count
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
geos <- geos[order(geos$today),] ## sort in original sample order

# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)
gsdat <- as.data.frame(cbind(geos, geosgrid))
gsdat$ccount[is.na(gsdat$ccount)] <- 1

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(bcoord, "./Results/MZ_bcoord.csv", row.names = F)
write.csv(gsdat, "./Results/MZ_gsdat_2019.csv", row.names = F)

# GeoSurvey map widgets ---------------------------------------------------
# number of GeoSurvey quadrats
w <- leaflet() %>%
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'MZ_GS19.html', selfcontained = T) ## save widget

# number of building tags
b <- leaflet() %>%
  setView(lng = mean(bcoord$lon), lat = mean(bcoord$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(bcoord$lon, bcoord$lat, clusterOptions = markerClusterOptions())
b ## plot widget 
saveWidget(b, 'MZ_GS19_buildings.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(4,0.1), random.order = T)

