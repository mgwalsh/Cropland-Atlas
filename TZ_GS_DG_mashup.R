# Tanzania GeoSurvey building density / DigitalGlobe building footprint density mashup
# M. Walsh, January 2020

# Required packages
# install.packages(c("downloader","rgdal","raster","MASS","caret","dismo","doParallel")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(MASS)
  require(caret)
  require(dismo)
  require(doParallel)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("TZ_GSDG", showWarnings = F)
setwd("./TZ_GSDG")
dir.create("Results", showWarnings = F)

# download TZ GeoSurvey building data
download("https://osf.io/8jma4?raw=1", "TZ_buildings.csv.zip", mode = "wb")
unzip("TZ_buildings.csv.zip", overwrite = T)
geos <- read.table("TZ_buildings.csv", header = T, sep = ",")

# download raster stack
download("https://osf.io/r8xgf?raw=1", "TZ_GS_DG_comp.zip", mode = "wb")
unzip("TZ_GS_DG_comp.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)

# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)
gsdat <- as.data.frame(cbind(geos, geosgrid)) 

# Model setup -------------------------------------------------------------
# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$BP, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# start doParallel
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# ZIP model ---------------------------------------------------------------
# calibration labels
lcal <- gs_cal$bcount

# prediction features
fcal <- gs_cal[,7:9]

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
zip <- train(fcal, lcal, 
             method = "glm",
             family = "poisson",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="RMSE")

# model predictions
print(zip)
summary(zip)
zip.pred <- (predict(grids, zip)) ## spatial predictions of building densities/ha
stopCluster(mc)
saveRDS(zip, "./Results/zip_bdens.rds")
plot(zip.pred, axes=F)

# Receiver-operator characteristics ---------------------------------------
gs_val$zip_pred <- predict(zip, gs_val) ## predictions on validation set
p <- gs_val[ which(gs_val$BP=="Y"), ]
p <- p[,12]
a <- gs_val[ which(gs_val$BP=="N"), ]
a <- a[,12]
e <- evaluate(p=p, a=a) ## calculate ROC
plot(e, 'ROC') ## plot ROC curve

# Generate settlement mask
gsdat$zip_pred <- predict(zip, gsdat) ## predictions on all of the GeoSurvey data
p <- gsdat[ which(gsdat$BP=="Y"), ]
p <- p[,12]
a <- gsdat[ which(gsdat$BP=="N"), ]
a <- a[,12]
e <- evaluate(p=p, a=a) ## calculate ROC
plot(e, 'ROC') ## plot ROC curve
t <- threshold(e) ## calculate thresholds based on ROC
mk <- reclassify(zip.pred, c(-Inf, t[,1], 0, t[,1], Inf, 1)) ## reclassify map based on kappa
plot(mk, axes=F)
gsdat$zip_pa <- ifelse(gsdat$zip_pred > t[,1], "Y", "N")
confusionMatrix(data = gsdat$zip_pa, reference = gsdat$BP, positive = "Y")

# Write files -------------------------------------------------------------
write.csv(gsdat, "./Results/TZ_building_zip.csv", row.names = F)
gspred <- stack(zip.pred, mk)
writeRaster(gspred, filename="./Results/TZ_bcount_250m.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
