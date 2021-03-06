# Tanzania GeoSurvey 2018 250m resolution GS-L2 data setup
# version 2, w. updated grid data
# M. Walsh, May 2019

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
download("https://osf.io/9gndf/?raw=1", "TZ_geos_2018.csv.zip", mode = "wb")
unzip("TZ_geos_2018.csv.zip", overwrite = T)
geos <- read.table("TZ_geos_2018.csv", header = T, sep = ",")

# download GADM-L3 shapefile (courtesy: http://www.gadm.org)
download("https://www.dropbox.com/s/bhefsc8u120uqwp/TZA_adm3.zip?raw=1", "TZA_adm3.zip", mode = "wb")
unzip("TZA_adm3.zip", overwrite = T)
shape <- shapefile("TZA_adm3.shp")

# download raster stack (note this is a big 1+ Gb download)
download("https://osf.io/ke5ya?raw=1", "TZ_250m_2019.zip", mode = "wb")
unzip("TZ_250m_2019.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# Data setup --------------------------------------------------------------
# attach GADM-L3 admin unit names from shape
coordinates(geos) <- ~lon+lat
projection(geos) <- projection(shape)
gadm <- geos %over% shape
geos <- as.data.frame(geos)
geos <- cbind(gadm[ ,c(5,7,9)], geos)
colnames(geos) <- c("region","district","ward","survey","time","id","observer","lat","lon","BP","CP","WP","rice","bloc","cgrid")

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
geos <- geos[order(geos$id),] ## sort in original sample order

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
gsdat <- gsdat[complete.cases(gsdat[ ,c(18:67)]),] ## removes incomplete cases
gsdat$observer <- sub("@.*", "", as.character(gsdat$observer)) ## shortens observer ID's

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(bcoord, "./Results/TZ_bcoord.csv", row.names = F)
write.csv(gsdat, "./Results/TZ_gsdat_2018.csv", row.names = F)

# GeoSurvey map widget ----------------------------------------------------
w <- leaflet() %>%
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(gsdat$lon, gsdat$lat, clusterOptions = markerClusterOptions())
w ## plot widget 
saveWidget(w, 'TZ_GS18.html', selfcontained = T) ## save widget

# GeoSurvey contributions -------------------------------------------------
gscon <- as.data.frame(table(gsdat$observer))
set.seed(1235813)
wordcloud(gscon$Var1, freq = gscon$Freq, scale = c(3,0.1), random.order = T)

