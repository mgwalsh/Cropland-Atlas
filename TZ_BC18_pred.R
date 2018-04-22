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
gf_cal <- gs_cal[,18:61]

# Quantile regression against DigitalGlobe counts ---------------------------
qrl <- rq(GDB ~ bcount, tau = 0.05, data = gsdat)
coef(qrl)
qrm <- rq(GDB ~ bcount, tau = 0.5, data = gsdat)
coef(qrm)
qrh <- rq(GDB ~ bcount, tau = 0.95, data = gsdat)
coef(qrh)

# Plot
par(pty="s")
plot(GDB~bcount, gsdat, ylab="DigitalGlobe building count", xlab="GeoSurvey building count", xlim=c(0,350), ylim=c(0,350), cex.lab=1.3)
lines(gsdat$bcount, qrl$fitted.values, col = "blue")
lines(gsdat$bcount, qrh$fitted.values, col = "red")
abline(c(0,1))
