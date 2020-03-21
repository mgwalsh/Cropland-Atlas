# Tanzania 2019 GeoSurvey 250m small area estimates (SAE)
# M. Walsh, December 2019

# Required packages
# install.packages(c("downloader","rgdal","raster","MASS","arm",leaflet","htmlwidgets")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(MASS)
  require(arm)
  require(leaflet)
  require(htmlwidgets)
})

# set working directory
dir.create("TZ_sae", showWarnings = F)
setwd("./TZ_sae")

# Data dow!nloads -----------------------------------------------------------
# download GeoSurvey data
download("https://osf.io/27avb?raw=1", "TZ_gsdat_2018.csv.zip", mode = "wb")
unzip("TZ_gsdat_2018.csv.zip", overwrite = T)
geos <- read.table("TZ_gsdat_2018.csv", header = T, sep = ",")
vars <- c("region","district","ward","lat","lon","BP","CP","bcount","ccount","PH")
geos <- geos[vars] ## removes extraneous variables
geos <- geos[ which(geos$ccount < 17), ] ## drops cropland miscounts

# download GeoSurvey prediction rasters
download("https://osf.io/fdkz8?raw=1", "TZ_GS_preds.zip", mode = "wb")
unzip("TZ_GS_preds.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)
# (gave <- cellStats(grids, mean)) ## calculates mean grids values

# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)
saedat <- as.data.frame(cbind(geos, geosgrid)) 

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(saedat, "./Results/TZ_sae.csv", row.names = F)

# Cropland area models ----------------------------------------------------
# binomial models of GeoSurvey cropland grid counts
# cp <-  saedat[which(saedat$CP=='Y'), ] ## actual cropland observations only
summary(m0 <- glm(cbind(ccount, 16-ccount) ~ 1, family=binomial, saedat)) ## mean model
(est0 <- cbind(Estimate = coef(m0), confint(m0))) ## standard 95% confidence intervals
# summary(mq <- glm(cbind(ccount, 16-ccount) ~ 1, family=quasibinomial(link="logit"), gsdat)) ## overdispersed model

# with cropland spatial presence prediction (CP18)
summary(m1 <- glm(cbind(ccount, 16-ccount) ~ CP18, family=binomial, saedat)) ## scaling model
(est1 <- cbind(Estimate = coef(m1), confint(m1))) ## standard 95% confidence intervals
m1.pred <- predict(grids, m1, type="response")
plot(m1.pred, axes=F)
# gsdat$m1 <- predict(m1, gsdat, type="response")

# with additional LCC covariates
summary(m2 <- glm(cbind(ccount, 16-ccount) ~ BP18*CP18*WP18, family=binomial, saedat))
(est2 <- cbind(Estimate = coef(m2), confint(m2))) ## standard 95% confidence intervals
anova(m1, m2) ## model comparison
m2.pred <- predict(grids, m2, type="response")
plot(m2.pred, axes=F)
# gsdat$m2 <- predict(m2, gsdat, type="response")

# Write prediction grids
gspreds <- stack(m1.pred, m2.pred)
names(gspreds) <- c("m1","m2")
writeRaster(gspreds, filename="./Results/TZ_cp_sae.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Small area estimates (SAE)
# post-stratified by admin units
summary(m3 <- glmer(cbind(ccount, 16-ccount) ~ 1 + (1|district), family=binomial, saedat))

#  with additional LCC covariates
summary(m4 <- glmer(cbind(ccount, 16-ccount) ~ BP19*CP19*WP19 + (1|district), family=binomial, saedat))
ran <- ranef(m4) ## extract regional random effects
ses <- se.coef(m4) ## extract regional standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## regional-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-0.5,0.5), CI=2, main="") ## district coefficient plot
write.csv(sae, "./Results/TZ_crop_area_sae.csv", row.names = F)

# Building count models ---------------------------------------------------
# Poisson models of GeoSurvey building counts
summary(m5 <- glm(bcount ~ 1, family=poisson, saedat)) ## country mean model
(est5 <- cbind(Estimate = coef(m5), confint(m5))) ## standard 95% confidence intervals
summary(mnb <- glm.nb(bcount ~ 1, saedat)) ## overdispersed negative binomial model
anova(m5, mnb)

# with building presence prediction (BP18, building probability)
summary(m6 <- glm(bcount ~ BP18, family=poisson, saedat)) ## scaling model
(est6 <- cbind(Estimate = coef(m6), confint(m6))) ## standard 95% confidence intervals
m6.pred <- predict(grids, m6, type="response")
dev.off()
plot(m6.pred, axes=F)
# gsdat$m6 <- predict(m6, gsdat, type="response")

# with additional LCC covariates
summary(m7 <- glm(bcount ~ BP18*CP18*WP18, family=poisson, saedat))
(est7 <- cbind(Estimate = coef(m7), confint(m7))) ## standard 95% confidence intervals
m7.pred <- predict(grids, m7, type="response")
plot(m7.pred, axes=F)
saedat$m7 <- predict(m7, gsdat, type="response")

# Write prediction grids
gspreds <- stack(m6.pred, m7.pred)
names(gspreds) <- c("m6","m7")
writeRaster(gspreds, filename="./Results/TZ_bcount_sae.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Small area estimates (SAE)
# post-stratified by admin units (72 districts)
summary(m8 <- glmer(bcount ~ 1 + (1|district), family=poisson, saedat))
ran <- ranef(m8) ## extract district random effects
ses <- se.coef(m8) ## extract district standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## district-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-2,2), CI=2, main="") ## district coefficient plot

# with additional LCC covariates
summary(m9 <- glmer(bcount ~ BP18*CP18*WP18 + (1|district), family=poisson, saedat))
anova(m8, m9) ## model comparison
ran <- ranef(m9) ## extract district random effects
ses <- se.coef(m9) ## extract district standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## district-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-1,1), CI=2, main="") ## district coefficient plot
write.csv(sae, "./Results/TZ_bcount_sae.csv", row.names = F)

