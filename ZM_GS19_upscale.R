# Zambia GeoSurvey 2019 L1/L2-GS upscale to 100m resolution
# M. Walsh, November 2019

# Required packages
# install.packages(c("downloader","rgdal","raster","dismo","doParallel"), dependencies=TRUE)
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
dir.create("ZM_GS19_100m", showWarnings = F)
setwd("./ZM_GS19_100m")
dir.create("Results", showWarnings = F)

# download GeoSurvey data
download("https://osf.io/6srce?raw=1", "ZM_buildings.csv.zip", mode = "wb")
unzip("ZM_buildings.csv.zip", overwrite = T)
geos <- read.table("ZM_buildings.csv", header = T, sep = ",")

# download raster stack (upscaled 250m building predictions)
download("https://osf.io/4m3x2?raw=1", "ZM_building_preds_2019.zip", mode = "wb")
unzip("ZM_building_preds_2019.zip", overwrite = T)
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

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$BP, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# calibration labels
lcal <- gs_cal$bcount

# prediction features
fcal <- gs_cal[,10:15]

# Upscaling model ---------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
up <- train(fcal, lcal, 
            method = "glm",
            family = "poisson",
            preProc = c("center","scale"), 
            trControl = tc,
            metric ="RMSE")

# model predictions
print(up)
summary(up)
up.pred <- (predict(grids, up))/6.25 ## spatial predictions of building densities/ha
stopCluster(mc)
saveRDS(up, "./Results/up_bdens.rds")

# Receiver-operator characteristics ---------------------------------------
gs_val$up_pred <- predict(up, gs_val) ## predictions on validation set
p <- gs_val[ which(gs_val$BP=="Y"), ]
p <- p[,18]
a <- gs_val[ which(gs_val$BP=="N"), ]
a <- a[,18]
e <- evaluate(p=p, a=a) ## calculate ROC
plot(e, 'ROC') ## plot ROC curve

# Generate settlement mask
gsdat$up_pred <- predict(up, gsdat) ## predictions on all of the GeoSurvey data
p <- gsdat[ which(gsdat$BP=="Y"), ]
p <- p[,18]
a <- gsdat[ which(gsdat$BP=="N"), ]
a <- a[,18]
e <- evaluate(p=p, a=a) ## calculate ROC
plot(e, 'ROC') ## plot ROC curve
t <- threshold(e) ## calculate thresholds based on ROC
mk <- reclassify(up.pred, c(-Inf, t[,1], 0, t[,1], Inf, 1)) ## reclassify map based on kappa
plot(mk, axes=F)
gsdat$up_pa <- ifelse(gsdat$up_pred > t[,1], "Y", "N")
confusionMatrix(data = gsdat$up_pa, reference = gsdat$BP, positive = "Y")

# Write files -------------------------------------------------------------
write.csv(gsdat, "./Results/ZM_building_upscale.csv", row.names = F)
gspred <- stack(up.pred, mk)
writeRaster(gspred, filename="./Results/ZM_bcount_100m.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

