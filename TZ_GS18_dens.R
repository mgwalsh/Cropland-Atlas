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
# gsdat <- gsdat[ which(gsdat$CP=='Y'), ] ## selects cropland observations

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
gf_cal <- gs_cal[,19:67]

# Central place theory model <glm> -----------------------------------------
# select central place covariates
gf_cpv <- gs_cal[,26:37]

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
saveRDS(gl1, "./Results/gl1_bdens.rds")

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
saveRDS(gl2, "./Results/gl2_bdens.rds")

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
saveRDS(rf, "./Results/rf_bdens.rds")

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
gb.pred <- predict(grids, gb) ## spatial predictions
stopCluster(mc)
saveRDS(gb, "./Results/gb_bdens.rds")

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
nn.pred <- predict(grids, nn) ## spatial predictions
stopCluster(mc)
saveRDS(nn, "./Results/nn_bdens.rds")

# Model stacking setup ----------------------------------------------------
preds <- stack(gl1.pred, gl2.pred, rf.pred, gb.pred, nn.pred)
names(preds) <- c("gl1","gl2","rf","gb","nn")
plot(preds, axes=F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# Model stacking ----------------------------------------------------------
# poisson model
summary(st <- glm(bcount ~ gl1+gl2+rf+gb+nn, family=poisson, gspred))
(est <- cbind(Estimate = coef(st), confint(st))) ## standard 95% confidence intervals
st.pred <- predict(preds, st, type="response")
plot(st.pred, axes=F)

# Write prediction grids --------------------------------------------------
gspreds <- stack(gl1.pred, gl2.pred, rf.pred, gb.pred, nn.pred, st.pred)
names(gspreds) <- c("gl1","gl2","rf","gb","nn","st")
writeRaster(gspreds, filename="./Results/ZM_bcount_2019.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)## ... change feature names here

# Write output data frame -------------------------------------------------
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(grids)
gspre <- extract(gspreds, gsdat)
gsout <- as.data.frame(cbind(gsdat, gspre))
write.csv(gsout, "./Results/RW_bcount_out.csv", row.names = F)

# Prediction plot checks
require(devtools)
require(quantreg)

par(pty="s")
par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(bcount~st, xlab="Ensemble prediction", ylab="GeoSurvey building density", cex.lab=1.3, 
     xlim=c(-1,300), ylim=c(-1,300), gsout)
stQ <- rq(bcount~st, tau=c(0.05,0.5,0.95), data=gsout)
print(stQ)
curve(stQ$coefficients[2]*x+stQ$coefficients[1], add=T, from=0, to=300, col="blue", lwd=1)
curve(stQ$coefficients[4]*x+stQ$coefficients[3], add=T, from=0, to=300, col="red", lwd=1)
curve(stQ$coefficients[6]*x+stQ$coefficients[5], add=T, from=0, to=300, col="blue", lwd=1)
abline(c(0,1), col="grey", lwd=2)
