# Stacked predictions of Tanzania 2018 GeoSurvey building count observations
# M. Walsh, April 2018

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","glmnet","plyr","doParallel","quantreg")), dependencies=T)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(glmnet)
  require(plyr)
  require(doParallel)
  require(quantreg)
})

# Data setup --------------------------------------------------------------
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/blob/master/TZ_GS18_data.R"
# source_url(SourceURL)
rm(list=setdiff(ls(), c("gsdat","grids"))) ## scrub extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$bcount, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- gs_cal$bcount 

# raster calibration features
gf_cal <- gs_cal[,c(18:39,41:62)]

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
             family = "poisson",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="RMSE")

# model outputs & predictions
summary(gl1)
print(gl1) ## RMSE accross cross-validation
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
             family = "poisson",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="RMSE")

# model outputs & predictions
summary(gl2)
print(gl2) ## RMSE accross cross-validation
gl2.pred <- predict(grids, gl2) ## spatial predictions

stopCluster(mc)

# Random forest <randomForest> --------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
tg <- expand.grid(mtry = seq(5,25, by=5)) ## model tuning steps

# model training
rf <- train(gf_cal, cp_cal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 501,
            metric = "RMSE",
            tuneGrid = tg,
            trControl = tc)

# model outputs & predictions
print(rf) ## RMSE accross tuning parameters
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
tg <- expand.grid(interaction.depth = seq(5,25, by=5), shrinkage = 0.01, n.trees = 501,
                  n.minobsinnode = 25) ## model tuning steps

# model training
gb <- train(gf_cal, cp_cal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg,
            metric = "RMSE")

# model outputs & predictions
print(gb) ## ROC's accross tuning parameters
gb.pred <- predict(grids, gb) ## spatial predictions

stopCluster(mc)


# Quantile regression against DigitalGlobe counts ---------------------------
qrl <- rq(GDB ~ bcount, tau = 0.05, data = gsdat)
coef(qrl)
qrm <- rq(GDB ~ bcount, tau = 0.5, data = gsdat)
coef(qrm)
qrh <- rq(GDB ~ bcount, tau = 0.95, data = gsdat)
coef(qrh)

# Quantile plot
par(pty="s")
plot(GDB~bcount, gsdat, ylab="DigitalGlobe building count", xlab="GeoSurvey building count", xlim=c(0,350), ylim=c(0,350), cex.lab=1.3)
lines(gsdat$bcount, qrl$fitted.values, col = "blue")
lines(gsdat$bcount, qrh$fitted.values, col = "red")
abline(c(0,1))
