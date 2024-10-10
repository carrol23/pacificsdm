#  ---------------------SDM ENVIRONMENTAL COVARIATES (OR PREDICTORS) DATA  -------------------------
# 1 download data from worldclim
# 2 perform correlation and remove highly correlated variables
# 3 write the reduced dataset to a new data folder
# 4 reproject and crop dataset to Pacific Island region

# ----------------- Set Working Directory ----------------
rm(list=ls()) # clean workspace
setwd('C:/Data/sdm/C2_SDM/clim/wc_30')

# ---------------------- INSTALL PACKAGES -----------------
# install.packages(c( 'sp', 'maps', 'tidyverse', 'ggplot2', 'geodata', 'raster', 'devtools', 'corrplot'))

library(dismo)
library(maps)
library(ggplot2)
library(sp)
library(geodata)
library(raster)
library(devtools)
# install mecofun package
devtools::install_git("https://gitup.uni-potsdam.de/macroecology/mecofun.git")
library(mecofun)

# ----------------- DOWNLOAD WORLDCLIM DATA ------------------
# get WORLDCLIM data - download 30 second data
currentEnv = worldclim_global(var = "bio", res = 0.5, path = ".")
elev = elevation_global(res = 0.5, path = "./wc2.1_30s") # same path as above

# convert spatraster to raster
currentEnv.Ras <- stack(as(currentEnv, "Raster"), as(elev, "Raster"))

# rename raster stack with simple names
clim_names = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8",
                           "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15",
                           "bio16", "bio17", "bio18", "bio19", "elevation")
names(currentEnv.Ras) <- clim_names

# check data by plotting the first variable
plot(currentEnv.Ras$bio1)

# write raster stack to file 
writeRaster(stack(currentEnv.Ras), filename = './clim/current_env.tif',
            format = 'GTiff', overwrite = TRUE)

# to load in written raster stack use:
currentEnv <- stack('./clim/current_env.tif')

# -------------------- PERFORM CORRELATION ------------------
library(corrplot)

# load in species presence data
occ = read.csv("C:/Data/sdm/C2_SDM/occurence/at/comb_at.csv.csv", sep = ',')
# extract lat and lon columns only and rename column
occ_df = cbind.data.frame(atulip_occ$decimalLongitude, atulip_occ$decimalLatitude)
names(occ_df)[names(occ_df) == 'occ$decimalLongitude'] <- 'lon'
names(occ_df)[names(occ_df) == 'occ$decimalLatitude'] <- 'lat'
occ_sp = SpatialPoints(occ_df)
crs(occ_sp) <- crs("+init=epsg:4326") # assign crs

# extract values of each global env covariates with each occurrence
current_env <- rast(currentEnv.Ras)
sp_pres <- vect(occ_sp)
crs(sp_pres) <- crs(current_env)
sp_var <- extract(current_env, sp_pres)

# include presence column (and fill with 1)
sp_var$presence <- 1

# get correlation matrix
cor_mat <- cor(sp_var[,c(2:21)], method = "spearman")
corrplot.mixed(cor_mat, tl.pos = 'lt', tl.cex = 1, number.cex = 1, addCoefasPercent = T)

# Run select07() and set threshold (at 0.7)
var_sel <- select07(X = sp_var[,c(2:21)],
                    y = sp_var$presence,
                    threshold = 0.7)

str(var_sel)

# extract names of weakly correlated predictors ordered by univariate variable importance in terms of the AIC
pred_sel <- var_sel$pred_sel
pred_sel
# ---------------- REMOVE CORRELATED VARIABLES ----------
# drop strong correlated variables

clim_stack = dropLayer(currentEnv.Ras, c("bio4","bio5",
                                     "bio6", "bio8", "bio9",
                                     "bio10", "bio11", "bio13",
                                     "bio15","bio16", "bio17",
                                     "bio19"))

# write reduced stack of independent variables to file
# setwd('C:/Data/sdm/C2_SDM/clim/data')
raster::writeRaster(clim_stack, filename = './clim/cor_stack.tif',
            format = 'GTiff', overwrite = TRUE)

# CROP PREDICTOR VARIABLES (CLIMATE) TO PACIFIC
# get pacific region extent first
library(rnaturalearth)
world <- ne_download(type = "countries", scale = "medium")

# select only oceania region
pc <- subset(world, world$REGION_UN == "Oceania")
plot(pc)

# reproject pc first
pdc <- crs("+init=epsg:3832")
pc_rpj <- st_transform(pc, pdc)
plot(pc_rpj, max.plot = 168)

# get extent
pic_ext <- extent(pc_rpj)


# Run this portion once and export to file
# reproject predictor variables to the PDC mercator projection
# NB: faster to process in batch in arcmap or qgis
pic_predictor <- projectRaster(clim_stack, crs = pdc)


# crop reprojected rasters based on extent of pacific region shp
pic_predictor <- crop(pic_predictor, pic_ext)

# write cropped rasters to file (check working directory)
# setwd('C:/Data/sdm/C2_SDM/clim/pic')
writeRaster(pic_predictor, filename = './clim/pic_stack.tif',
            format = 'GTiff', overwrite = TRUE)
