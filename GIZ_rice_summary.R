# Tanzania, Nigeria, Ghana, Ethiopia "paddy rice" 2018 GeoSurvey summary
# M. Walsh, June 2018

# Required packages
# install.packages(c("downloader","MASS")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(MASS)
})

# Data downloads ----------------------------------------------------------
# set working directory
dir.create("Rice_GS18", showWarnings = F)
setwd("./Rice_GS18")

# download GeoSurvey data
download("https://www.dropbox.com/s/ugo0q0u9noy7nyt/GIZ_rice_2018.csv.zip?raw=1", "GIZ_rice_2018.csv.zip", mode = "wb")
unzip("GIZ_rice_2018.csv.zip", overwrite = T)
geos <- read.table("GIZ_rice_2018.csv", header = T, sep = ",")

# GS data summary ---------------------------------------------------------
rice <- geos[which(geos$CP == 'Y'), ]
test <- glm(rice~Country-1, family = binomial, rice)  
summary(test)

