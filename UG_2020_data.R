# Uganda GeoSurvey 2020 250m resolution L1-GS data setup
# M. Walsh, May 2020

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
dir.create("UG_GS20", showWarnings = F)
setwd("./UG_GS20")

# download GeoSurvey data
download("https://osf.io/zpbwm?raw=1", "UG_GS_2020.csv.zip", mode = "wb")
unzip("UG_GS_2020.csv.zip", overwrite = T)
geos <- read.table("UG_GS_2020.csv", header = T, sep = ",")

# download GADM-L4 shapefile (courtesy of: http://www.gadm.org)
download("https://osf.io/f9tr3?raw=1", "UG_GADM_L4.zip", mode = "wb")
unzip("UG_GADM_L4.zip", overwrite = T)
shape <- shapefile("gadm36_UGA_4.shp")

# download raster stack
download("https://osf.io/t9dvz?raw=1", "UG_250m_2020.zip", mode = "wb")
unzip("UG_250m_2020.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L4 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(4,6,8,10)], geos)
colnames(geos) <- c("district","county","subcounty","parish","today","observer","id","lat","lon","BP","CP","WP","CS","bloc","cgrid")

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
write.csv(bcoord, "./Results/UG_bcoord.csv", row.names = F)
write.csv(gsdat, "./Results/UG_gsdat_2020.csv", row.names = F)

# GeoSurvey map widgets ---------------------------------------------------
# number of GeoSurvey quadrats
w <- leaflet() %>%
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'UG_GS20.html', selfcontained = T) ## save widget

# number of building tags
b <- leaflet() %>%
  setView(lng = mean(bcoord$lon), lat = mean(bcoord$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(bcoord$lon, bcoord$lat, clusterOptions = markerClusterOptions())
b ## plot widget 
saveWidget(b, 'UG_GS20_buildings.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(4,0.1), random.order = T)

