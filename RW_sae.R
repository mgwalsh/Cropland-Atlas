# Rwanda 2019 GeoSurvey 250m small area estimates (SAE)
# M. Walsh, July 2019

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
dir.create("RW_sae", showWarnings = F)
setwd("./RW_sae")

# Data dow!nloads -----------------------------------------------------------
# download GeoSurvey data
# see sampling frame @ https://github.com/mgwalsh/Sampling/blob/master/TZ_GS_sample.R
download("https://www.dropbox.com/s/e4k2vuc3homygfi/RW_gsdat_2019.csv.zip?raw=1", "RW_gsdat_2019.csv.zip", mode = "wb")
unzip("RW_gsdat_2019.csv.zip", overwrite = T)
geos <- read.table("RW_gsdat_2019.csv", header = T, sep = ",")

# download GeoSurvey prediction rasters
download("https://www.dropbox.com/s/fif34ecgagcfpej/RW_GS_preds_2019.zip?raw=1", "RW_GS_preds_2019.zip", mode = "wb")
unzip("RW_GS_preds_2019.zip", overwrite = T)
glist <- list.files(pattern="tif", full.names = T)
grids <- stack(glist)
(gave <- cellStats(grids, mean)) ## calculates mean grids values

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
gsdat <- gsdat[complete.cases(gsdat[ ,c(59:63)]),] ## removes incomplete cases
gsdat <- gsdat[c(1:63)] 

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/TZ_sae.csv", row.names = F)

# Cropland area models ----------------------------------------------------
# binomial models of GeoSurvey cropland grid counts
# cp <-  gsdat[which(gsdat$CP=='Y'), ] ## actual cropland observations only
summary(m0 <- glm(cbind(ccount, 16-ccount) ~ 1, family=binomial, gsdat)) ## mean model
(est0 <- cbind(Estimate = coef(m0), confint(m0))) ## standard 95% confidence intervals
# summary(mq <- glm(cbind(ccount, 16-ccount) ~ 1, family=quasibinomial(link="logit"), gsdat)) ## overdispersed model

# with cropland spatial presence prediction (CP19)
summary(m1 <- glm(cbind(ccount, 16-ccount) ~ CP19, family=binomial, gsdat)) ## scaling model
(est1 <- cbind(Estimate = coef(m1), confint(m1))) ## standard 95% confidence intervals
m1.pred <- predict(grids, m1, type="response")
(m1.area <- cellStats(m1.pred*6.25, sum)) ## calculates total cropland area (ha)
plot(m1.pred, axes=F)
gsdat$m1 <- predict(m1, gsdat, type="response")

# +additional LCC covariates
summary(m2 <- glm(cbind(ccount, 16-ccount) ~ BC19+BP19+CP19+TP19+WP19, family=binomial, gsdat))
(est2 <- cbind(Estimate = coef(m2), confint(m2))) ## standard 95% confidence intervals
anova(m1, m2) ## model comparison
m2.pred <- predict(grids, m2, type="response")
(m2.area <- cellStats(m2.pred*6.25, sum)) ## calculates total cropland area (ha)
plot(m2.pred, axes=F)
gsdat$m2 <- predict(m2, gsdat, type="response")

# Small area estimates (SAE)
# post-stratified by districts
summary(m3 <- glmer(cbind(ccount, 16-ccount) ~ 1 + (1|district), family=binomial, gsdat))
summary(m4 <- glmer(cbind(ccount, 16-ccount) ~ 1 + (1|district/sector), family=binomial, gsdat))

# +additional LCC covariates
summary(m5 <- glmer(cbind(ccount, 16-ccount) ~ BC19+BP19+CP19+TP19+WP19 + (1|district), family=binomial, gsdat))
anova(m3, m5) ## model comparison
ran <- ranef(m5) ## extract regional random effects
ses <- se.coef(m5) ## extract regional standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## regional-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-1,1), CI=2, main="") ## district coefficient plot
dev.off()
write.csv(sae, "./Results/RW_crop_area_sae.csv", row.names = F)

# Write prediction grids
gspreds <- stack(m1.pred, m2.pred)
names(gspreds) <- c("m1","m2")
writeRaster(gspreds, filename="./Results/TZ_cp_area.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Cropland area map widget
pred <- m2.pred*100 ## GeoSurvey cropland percentage
pal <- colorBin("Reds", domain = 0:100, na.color = "light grey") ## set color palette
w <- leaflet() %>% 
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 9) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.6, maxBytes=6000000) %>%
  addLegend(pal = pal, values = values(pred), title = "Cropland area (%)")
w ## plot widget 
saveWidget(w, 'RW_cp_area.html', selfcontained = T)

# Building count models ---------------------------------------------------
# Poisson models of GeoSurvey building counts
# bp <-  gsdat[which(gsdat$BP=='Y'), ] ## actual settlement observations only
summary(m6 <- glm(bcount ~ 1, family=poisson, gsdat)) ## country mean model
(est6 <- cbind(Estimate = coef(m6), confint(m6))) ## standard 95% confidence intervals
summary(mnb <- glm.nb(bcount ~ 1, gsdat)) ## overdispersed negative binomial model
anova(m6, mnb)

# with building count prediction (BC19)
summary(m7 <- glm(bcount ~ BC19, family=poisson, gsdat)) ## scaling model
(est7 <- cbind(Estimate = coef(m7), confint(m7))) ## standard 95% confidence intervals
m7.pred <- predict(grids, m7, type="response")
plot(m7.pred, axes=F)
gsdat$m7 <- predict(m7, gsdat, type="response")

# +additional LCC covariates
summary(m8 <- glm(bcount ~ BC19+BP19+CP19+TP19+WP19, family=poisson, gsdat))
(est2 <- cbind(Estimate = coef(m8), confint(m8))) ## standard 95% confidence intervals
anova(m7, m8) ## model comparison
m8.pred <- predict(grids, m8, type="response")
plot(m8.pred, axes=F)
gsdat$m8 <- predict(m8, gsdat, type="response")

# Small area estimates (SAE)
# post-stratified by districts
summary(m9 <- glmer(bcount ~ 1 + (1|district), family=poisson, gsdat))
summary(m10 <- glmer(bcount ~ 1 + (1|district/sector), family=poisson, gsdat))

# +additional LCC covariates
summary(m11 <- glmer(bcount ~ BC19+BP19+CP19+TP19+WP19 + (1|district), family=poisson, gsdat))
anova(m9, m11) ## model comparison
ran <- ranef(m11) ## extract district random effects
ses <- se.coef(m11) ## extract district standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## district-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-0.2,0.2), CI=2, main="") ## district coefficient plot
write.csv(sae, "./Results/RW_bcount_sae.csv", row.names = F)

# Write prediction grids
gspreds <- stack(m1.pred, m2.pred, m7.pred, m8.pred)
names(gspreds) <- c("m1","m2","m7","m8")
writeRaster(gspreds, filename="./Results/RW_sae.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

