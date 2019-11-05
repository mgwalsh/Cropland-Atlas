# Stacked predictions of Zambia GeoSurvey observations
# M. Walsh, July 2019

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","plyr","doParallel","dismo")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(plyr)
  require(doParallel)
  require(dismo)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Cropland-Atlas/blob/master/ZM_GS19_data.R
rm(list=setdiff(ls(), c("gsdat","grids","glist"))) ## scrub extraneous objects in memory
gsdat <- gsdat[complete.cases(gsdat[ ,c(9:11,17:55)]),] ## removes incomplete cases

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$BP, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
labs <- c("BP") ## insert other labels (CP,WP ...) here!
lcal <- as.vector(t(gs_cal[labs]))

# raster calibration features
fcal <- gs_cal[,17:55]

# Central place theory model <glm> -----------------------------------------
# select central place covariates
fcpv <- gs_cal[,22:31,45]

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl1 <- train(fcpv, lcal, 
             method = "glmStepAIC",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model outputs & predictions
summary(gl1)
print(gl1) ## ROC's accross cross-validation
gl1.pred <- predict(grids, gl1, type = "prob") ## spatial predictions
stopCluster(mc)
fname <- paste("./Results/", labs, "_gl1.rds", sep = "")
saveRDS(gl1, fname)

# GLM with all covariates -------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl2 <- train(fcal, lcal, 
             method = "glmStepAIC",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model outputs & predictions
summary(gl2)
print(gl2) ## ROC's accross cross-validation
gl2.pred <- predict(grids, gl2, type = "prob") ## spatial predictions
stopCluster(mc)
fname <- paste("./Results/", labs, "_gl2.rds", sep = "")
saveRDS(gl2, fname)

# Random forest <randomForest> --------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(mtry = seq(1,5, by=1)) ## model tuning steps

# model training
rf <- train(fcal, lcal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            metric = "ROC",
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## ROC's accross tuning parameters
rf.pred <- predict(grids, rf, type = "prob") ## spatial predictions
stopCluster(mc)
fname <- paste("./Results/", labs, "_rf.rds", sep = "")
saveRDS(rf, fname)

# Generalized boosting <gbm> ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, summaryFunction = twoClassSummary,
                   allowParallel = T)

## for initial <gbm> tuning guidelines see @ https://stats.stackexchange.com/questions/25748/what-are-some-useful-guidelines-for-gbm-parameters
tg <- expand.grid(interaction.depth = seq(2,5, by=1), shrinkage = 0.01, n.trees = seq(101,501, by=50),
                  n.minobsinnode = 50) ## model tuning steps

# model training
gb <- train(fcal, lcal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg,
            metric = "ROC")

# model outputs & predictions
print(gb) ## ROC's accross tuning parameters
gb.pred <- predict(grids, gb, type = "prob") ## spatial predictions
stopCluster(mc)
fname <- paste("./Results/", labs, "_gb.rds", sep = "")
saveRDS(gb, fname)

# Neural network <nnet> ---------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(size = seq(2,10, by=2), decay = c(0.001, 0.01, 0.1)) ## model tuning steps

# model training
nn <- train(fcal, lcal, 
            method = "nnet",
            preProc = c("center","scale"), 
            tuneGrid = tg,
            trControl = tc,
            metric ="ROC")

# model outputs & predictions
print(nn) ## ROC's accross tuning parameters
nn.pred <- predict(grids, nn, type = "prob") ## spatial predictions
stopCluster(mc)
fname <- paste("./Results/", labs, "_nn.rds", sep = "")
saveRDS(nn, fname)

# Model stacking setup ----------------------------------------------------
preds <- stack(1-gl1.pred, 1-gl2.pred, 1-rf.pred, 1-gb.pred, 1-nn.pred)
names(preds) <- c("gl1","gl2","rf","gb","nn")
plot(preds, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
gs_val <- as.data.frame(gs_val)
lval <- as.vector(t(gs_val[labs]))
fval <- gspred[,56:60] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
st <- train(fval, lval,
            method = "glm",
            family = "binomial",
            metric = "ROC",
            trControl = tc)

# model outputs & predictions
print(st)
st.pred <- predict(preds, st, type = "prob") ## spatial predictions
plot(1-st.pred, axes = F)
stopCluster(mc)
fname <- paste("./Results/", labs, "_st.rds", sep = "")
saveRDS(st, fname)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(st, fval, type="prob")
cp_val <- cbind(lval, cp_pre)
cpp <- subset(cp_val, cp_val=="Y", select=c(Y))
cpa <- subset(cp_val, cp_val=="N", select=c(Y))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC's on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# Generate feature mask ---------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,1], 0, t[,1], 1, 1), ncol=3, byrow = T) ## set threshold value <kappa>
mask <- reclassify(1-st.pred, r) ## reclassify stacked predictions
plot(mask, axes=F)

# Write prediction grids --------------------------------------------------
gspreds <- stack(preds, 1-st.pred, mask)
names(gspreds) <- c("gl1","gl2","rf","gb","nn","st","mk")
fname <- paste("./Results/","ZM_", labs, "_preds_2019.tif", sep = "")
writeRaster(gspreds, filename=fname, datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Write output data frame -------------------------------------------------
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(grids)
gspre <- extract(gspreds, gsdat)
gsout <- as.data.frame(cbind(gsdat, gspre))
gsout$mzone <- ifelse(gsout$mk == 1, "Y", "N")
# confusionMatrix(data = gsout$mzone, reference = gsout$BP, positive = "Y")
fname <- paste("./Results/","ZM_", labs, "_out.tif", sep = "")
write.csv(gsout, fname, row.names = F)

# Prediction map widget ---------------------------------------------------
pred <- 1-st.pred ## GeoSurvey ensemble probability
pal <- colorBin("Reds", domain = 0:1) ## set color palette
w <- leaflet() %>% 
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 7) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.4, maxBytes=6000000) %>%
  addLegend(pal = pal, values = values(pred), title = "Probability")
w ## plot widget 
saveWidget(w, 'ZW_BP_2019.html', selfcontained = T) ## save html ... change feature names here

