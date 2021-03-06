# Zambia 2019 GeoSurvey 250m small area estimates (SAE)
# M. Walsh, August 2019

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
dir.create("ZM_sae", showWarnings = F)
setwd("./ZM_sae")

# Data dow!nloads -----------------------------------------------------------
# download GeoSurvey data
# see sampling frame @ https://github.com/mgwalsh/Sampling/blob/master/ZM_GS_sample.R
download("https://www.dropbox.com/s/ph52baolxmg7v4l/ZM_gsdat_2019.csv.zip?raw=1", "ZM_gsdat_2019.csv.zip", mode = "wb")
unzip("ZM_gsdat_2019.csv.zip", overwrite = T)
geos <- read.table("ZM_gsdat_2019.csv", header = T, sep = ",")

# download GeoSurvey prediction rasters
download("https://www.dropbox.com/s/ks5eg6grmodl17g/ZM_GS_preds.zip?raw=1", "ZM_GS_preds.zip", mode = "wb")
unzip("ZM_GS_preds.zip", overwrite = T)
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
gsdat <- as.data.frame(cbind(geos, geosgrid)) 
# gsdat <- gsdat[!duplicated(gsdat), ] ## removes any duplicates ... if needed
gsdat <- gsdat[complete.cases(gsdat[ ,c(57:63)]),] ## removes incomplete cases
gsdat <- gsdat[c(1:63)] 

# Write data frame --------------------------------------------------------
dir.create("Results", showWarnings = F)
write.csv(gsdat, "./Results/ZM_sae.csv", row.names = F)

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
plot(m1.pred, axes=F)
# gsdat$m1 <- predict(m1, gsdat, type="response")

# +additional LCC covariates
summary(m2 <- glm(cbind(ccount, 16-ccount) ~ BP19*CP19*WP19, family=binomial, gsdat))
(est2 <- cbind(Estimate = coef(m2), confint(m2))) ## standard 95% confidence intervals
anova(m1, m2) ## model comparison
m2.pred <- predict(grids, m2, type="response")
plot(m2.pred, axes=F)
# gsdat$m2 <- predict(m2, gsdat, type="response")

# Write prediction grids
gspreds <- stack(m1.pred, m2.pred)
names(gspreds) <- c("m1","m2")
writeRaster(gspreds, filename="./Results/ZM_cp_sae.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Small area estimates (SAE)
# post-stratified by admin units
summary(m3 <- glmer(cbind(ccount, 16-ccount) ~ 1 + (1|district), family=binomial, gsdat))

# +additional LCC covariates
summary(m4 <- glmer(cbind(ccount, 16-ccount) ~ BP19*CP19*WP19 + (1|district), family=binomial, gsdat))
ran <- ranef(m4) ## extract regional random effects
ses <- se.coef(m4) ## extract regional standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## regional-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-0.5,0.5), CI=2, main="") ## district coefficient plot
write.csv(sae, "./Results/ZM_crop_area_sae.csv", row.names = F)

# Building count models ---------------------------------------------------
# Poisson models of GeoSurvey building counts
# bp <-  gsdat[which(gsdat$BP=='Y'), ] ## actual settlement observations only
summary(m5 <- glm(bcount ~ 1, family=poisson, gsdat)) ## country mean model
(est5 <- cbind(Estimate = coef(m5), confint(m5))) ## standard 95% confidence intervals
summary(mnb <- glm.nb(bcount ~ 1, gsdat)) ## overdispersed negative binomial model
anova(m5, mnb)

# with building presence prediction (BM19, building mask)
summary(m6 <- glm(bcount ~ BM19, family=poisson, gsdat)) ## scaling model
(est6 <- cbind(Estimate = coef(m6), confint(m6))) ## standard 95% confidence intervals
m6.pred <- predict(grids, m6, type="response")
dev.off()
plot(m6.pred, axes=F)
# gsdat$m6 <- predict(m6, gsdat, type="response")

# +additional LCC covariates
summary(m7 <- glm(bcount ~ BP19*CP19*WP19, family=poisson, gsdat))
(est7 <- cbind(Estimate = coef(m7), confint(m7))) ## standard 95% confidence intervals
m7.pred <- predict(grids, m7, type="response")
plot(m7.pred, axes=F)
gsdat$m7 <- predict(m7, gsdat, type="response")

# Write prediction grids
gspreds <- stack(m6.pred, m7.pred)
names(gspreds) <- c("m6","m7")
writeRaster(gspreds, filename="./Results/ZM_bcount_sae.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Small area estimates (SAE)
# post-stratified by admin units (72 districts)
summary(m8 <- glmer(bcount ~ 1 + (1|district), family=poisson, gsdat))

# +additional LCC covariates
summary(m9 <- glmer(bcount ~ BP19*CP19*WP19 + (1|district), family=poisson, gsdat))
anova(m8, m9) ## model comparison
ran <- ranef(m9) ## extract district random effects
ses <- se.coef(m9) ## extract district standard errors
nam <- rownames(ran$district)
sae <- as.data.frame(cbind(ran$district, ses$district)) ## district-level small area estimates
colnames(sae) <- c("ran","se")
par(pty="s", mar=c(10,10,1,1))
coefplot(ran$district[,1], ses$district[,1], varnames=nam, xlim=c(-1,1), CI=2, main="") ## district coefficient plot
write.csv(sae, "./Results/ZM_bcount_sae.csv", row.names = F)

# Building density map widget
pred <- m7.pred/6.25 ## GeoSurvey building densities (per ha)
pal <- colorBin("Reds", domain = 0:maxValue(pred), na.color = "light grey") ## set color palette
w <- leaflet() %>% 
  setView(lng = mean(gsdat$lon), lat = mean(gsdat$lat), zoom = 7) %>%
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.6, maxBytes=6000000) %>%
  addLegend(pal = pal, values = values(pred), title = "Building density")
w ## plot widget 
saveWidget(w, 'ZM_bcount_sae.html', selfcontained = T)
