# Stacked predictions of Tanzania 2018 GeoSurvey building density observations
# M. Walsh, January 2019

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","plyr","doParallel")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(plyr)
  require(doParallel)
})

# Data setup --------------------------------------------------------------
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/blob/master/TZ_GS18_data.R"
# source_url(SourceURL)
rm(list=setdiff(ls(), c("gsdat","grids"))) ## scrubs extraneous objects in memory
gsdat <- gsdat[ which(gsdat$CP=='Y'), ] ## selects observations in croplands only, could be changed to $BP, $BIC

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$bcount, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- log(((gs_cal$bcount)/6.25)+1) ## log transform of the building count data to buildings/ha

# raster calibration features
gf_cal <- gs_cal[,18:62]

# Central place theory model <glm> -----------------------------------------
# select central place covariates
gf_cpv <- gs_cal[,27:37]

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
gl1 <- train(gf_cpv, cp_cal, 
             method = "glmStepAIC",
             family = "gaussian",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="RMSE")

# model outputs & predictions
gl1
summary(gl1)
gl1.pred <- predict(grids, gl1) ## spatial predictions

stopCluster(mc)

# GLM with all covariates -------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
gl2 <- train(gf_cal, cp_cal, 
             method = "glmStepAIC",
             family = "gaussian",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="RMSE")

# model outputs & predictions
gl2
summary(gl2)
gl2.pred <- predict(grids, gl2) ## spatial predictions

stopCluster(mc)

# Random forest <randomForest> --------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(mtry = seq(1,8, by=1)) ## model tuning steps

# model training
rf <- train(gf_cal, cp_cal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## RMSEs accross tuning parameters
rf.pred <- predict(grids, rf) ## spatial predictions

stopCluster(mc)

# Generalized boosting <gbm> ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

## for initial <gbm> tuning guidelines see @ https://stats.stackexchange.com/questions/25748/what-are-some-useful-guidelines-for-gbm-parameters
tg <- expand.grid(interaction.depth = seq(6,12, by=2), shrinkage = seq(0.02,0.1, by=0.02), n.trees = 501,
                  n.minobsinnode = 25) ## model tuning steps

# model training
gb <- train(gf_cal, cp_cal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg)

# model outputs & predictions
print(gb) ## RMSEs accross tuning parameters
plot(varImp(gb)) ## relative variable importance
gb.pred <- predict(grids, gb) ## spatial predictions

stopCluster(mc)

# Neural network <nnet> ---------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(size = seq(2,10, by=2), decay = c(0.001, 0.01, 0.1)) ## model tuning steps

# model training
nn <- train(gf_cal, cp_cal, 
            method = "nnet",
            preProc = c("center","scale"), 
            tuneGrid = tg,
            trControl = tc,
            metric ="RMSE")

# model outputs & predictions
print(nn) ## RMSEs accross tuning parameters
plot(varImp(nn)) ## relative variable importance
nn.pred <- predict(grids, nn) ## spatial predictions

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(gl1.pred, gl2.pred, rf.pred, gb.pred, nn.pred)
names(preds) <- c("gl1","gl2","rf","gb","nn")
plot(preds, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
cp_val <- gspred$bcount
cp_val <- log(((gs_val$bcount)/6.25)+1)
gf_val <- gspred[,63:67] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# model training
st <- train(gf_val, cp_val,
            method = "glm",
            trControl = tc)

# model outputs & predictions
summary(st)
plot(varImp(st))
st.pred <- predict(preds, st) ## spatial predictions
plot(st.pred, axes = F)

stopCluster(mc)

# Write prediction grids --------------------------------------------------
gspreds <- stack(preds, 1-st.pred, mask)
names(gspreds) <- c("rr","rf","gb","nn","st","mk")
# change this to include other dependent variables e.g, $BP, $BIC
writeRaster(gspreds, filename="./Results/TZ_rice_preds_2018.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)## ... change feature names here

# Write output data frame -------------------------------------------------
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(grids)
gspre <- extract(gspreds, gsdat)
gsout <- as.data.frame(cbind(gsdat, gspre))
write.csv(gsout, "./Results/TZ_bcount_out.csv", row.names = F)

# Prediction map widget ---------------------------------------------------
pred <- st.pred ## GeoSurvey ensemble probability
pal <- colorBin("Reds", domain = 0:5) ## set color palette
w <- leaflet() %>% 
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 6) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.6, maxBytes=6000000) %>%
  addLegend(pal = pal, values = values(pred), title = "Probability")
w ## plot widget 
saveWidget(w, 'TZ_bcount.html', selfcontained = T) ## save html ... change feature names here

