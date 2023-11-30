#  --------------------- SPECIES OCCURENCE DATA  -------------------------
# 1 download data across gbif/ALA
# 2 clean and merge datasets
# 3 perform cleaning based on coordinate cleaner
# 4 create global and native range dataset of occurrence data
# 5 export occurrence data

# ------ Set Working Directory --------
setwd('C:/Data/niche-modelling/SDM_22/merremia/30sec/data')

rm(list=ls()) # clean workspace

# ---------------- INSTALL PACKAGES -------------
install.packages(c('dismo', 'rgbif', 'coordinateCleaner', 'sp', 'maps', 'tidyverse', 'ggplot2'))
library(dismo)
library(maps)
library(ggplot2)
library(sp)
library(rgbif)
library(tidyverse)
library(CoordinateCleaner)


# ---------- DOWNLOAD SPECIES DATA -------
# run occurrence search using alternative name (combined data under Merremia peltata scientific name)
merr_alt = occ_search(scientificName = "Decalobanthus peltatus", limit = 10000, hasCoordinate = TRUE)

# select data
merr_alt <- merr_alt$data # 819 occurrences

# add in ALA data set extracted from local directory (after download)
merr_csv <- read.csv('C:/Data/niche-modelling/SDM_22/merremia/data/ala/records-2022-11-20.csv')

# view data for both imported
View(merr_alt)
View(merr_csv)

# --------------- CLEAN DATA --------------
# scale data down - remove columns not useful (only select the following)
merr_alt <- merr_alt %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)

# remove records without coordinates - records should be the same
merr_alt <- merr_alt %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# identify duplicates based on lat and long coordinates
merr_dup = duplicated(merr_alt[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 256 duplicates (out of 691 observations)
length(merr_dup[merr_dup == TRUE])

# remove duplicates - returns 435 observations
merr_data <- merr_alt[!merr_dup,]

# remove records without coordinates and duplicates for csv loaded
merr_csv <- merr_csv %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

merr_csv_dup <- duplicated(merr_csv[,c("decimalLongitude", "decimalLatitude")])
length(merr_csv_dup[merr_csv_dup == TRUE]) # get number of duplicated = 98
merr_csv <- merr_csv[!merr_csv_dup,]

# scale data down - remove columns not useful (only select the following)
merr_csv <- merr_csv %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         dataResourceUid, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)


# plot remaining occurrences
wrld = borders("world", colour = "gray50", fill = "gray80")

mplot <- ggplot() +  wrld +
  geom_point(data = merr_data, aes(x=decimalLongitude, y=decimalLatitude),
             colour = "yellow", size = 3, alpha = 0.1) +
  geom_point(data = merr_csv, aes(x = decimalLongitude, y = decimalLatitude), colour = "blue", size = 1.5) + 
  coord_sf()
  
mplot

# ---------------------- COORDINATE CLEANER TESTS ---------------------
# use coordinate cleaner to flag problem records
library(countrycode)

# convert country code ISO2 to ISO3
merr_data$countryCode <- countrycode(merr_data$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')
merr_csv$countryCode <- countrycode(merr_csv$countryCode, 
                                       origin = 'iso2c',
                                       destination = 'iso3c')

# flag coordinates not 'clean' based on tests
flags <- clean_coordinates(x = merr_data,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'outliers', 'seas', 'countries'))

# flag coordinates not 'clean' based on tests
flags_csv <- clean_coordinates(x = merr_csv,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'seas')) #countries

summary(flags) # 180 flagged records
summary(flags_csv) # 73 flagged records
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
plot(flags_csv, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude flagged records
merr_occ <- merr_data[flags$.summary,] # 289 records remain
merr_occ_csv <- merr_csv[flags_csv$.summary,] #565 records remain

# merge occurrence datasets
merr_occ = subset(merr_occ, select=-c(gbifID))
merr_occ_csv = subset(merr_occ_csv, select = -c(dataResourceUid))

merr_occ <- rbind(merr_occ, merr_occ_csv)

# identify duplicates based on lat and long coordinates
merr_occ_dup = duplicated(merr_occ[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 3 duplicates (out of 827 observations)
length(merr_occ_dup[merr_occ_dup == TRUE])

# remove duplicates - returns 435 observations
merr_occ <- merr_occ[!merr_occ_dup,]

# ----- METADATA TEST -------
# remove records with low precision
hist(merr_occ$coordinateUncertaintyInMeters / 1000, breaks = 20)

# removed coordinates with uncertainty up to 100km
merr_occ <- merr_occ %>%
  filter(coordinateUncertaintyInMeters/1000 <= 100 | is.na(coordinateUncertaintyInMeters))

# identify basis of records - remove fossil, living & preserved specimen 
table(merr_occ$basisOfRecord)

merr_data_filt <- filter(merr_occ, basisOfRecord == "HUMAN_OBSERVATION" |
                        basisOfRecord == "OBSERVATION" |
                        basisOfRecord == "PRESERVED_SPECIMEN" |
                        basisOfRecord == "OCCURENCE")

table(merr_data_filt$individualCount)

# remove individual counts that equate to 0 and higher than 100
merr_data_filt <-merr_data_filt %>%
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 99 | is.na(individualCount))

# exclude old records- 1925 remove
table(merr_data_filt$year)

merr_occ <- merr_data_filt %>%
  filter(year > 1945)

# check family name
table(merr_occ$family)

# check taxon rank
table(merr_occ$taxonRank)

# ------------------------ CREATE NATIVE RANGE OCCURENCE DATA ----------------------
# Split occurrence data within defined Native Range
nr <- vect("C:/Data/niche-modelling/SDM_22/merremia/30sec/nr_merr.shp")
nr <- as(nr, "Spatial")

# convert to df and add in columns (id and species)
nr_df <- as.data.frame(nr, row.names = T)
nr_df$id <- c("A")
nr_df$species <- c("Decalobanthus peltatus")

# convert native range extent to a spdf (poly)
nr_shp <- SpatialPolygonsDataFrame(nr, nr_df, match.ID = TRUE)
nr <- vect(nr_shp)

range <- ggplot() +
  wrld +
  geom_polygon(data = nr_shp, aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.title = element_blank())
range

# crop spatial points
merr_occ_sdf <- as.data.frame(merr_occ)

# get coordinates and create spatial points dataframe 
merr_occ_coor = as.data.frame(cbind(merr_occ_sdf$decimalLongitude,merr_occ_sdf$decimalLatitude))
merr_occ_sdf <- SpatialPointsDataFrame(merr_occ_coor, merr_occ_sdf)

nr_ext <- extent(nr_shp)
data_fin <- crop(merr_occ_sdf, nr_ext) # crop to get points that fall within the native extent 

data_fin <- as.data.frame(data_fin)
data_fin <- fortify(data_fin)

# plot native and global occurrence points
range +
  geom_point(data = merr_occ, aes(x=decimalLongitude, y = decimalLatitude), color = 'blue') +
  geom_point(data = data_fin, aes(x=decimalLongitude, y=decimalLatitude), color = 'red')

# ------------- Problematic datasets -------------
# ddmm to dd.dd conversion error - 0 records are flagged
out.ddmm <- cd_ddmm(merr_occ, lon = "decimalLongitude",
                    lat = "decimalLatitude", 
                    ds = "species", diagnostic = T, 
                    diff = 1, value = "dataset")

# test for rasterized sampling
# cd_round function identifies datasets with significant proportion of coords
par(mfrow = c(2,2), mar = rep(2, 4))
out.round <- cd_round(data_fin, lon = "decimalLongitude", lat = "decimalLatitude",
                      ds = "species", value = "dataset",
                      T1 = 7, graphs = T)

# export occurrence data as CSV - plot as spatial data
library(sf)
st_write(merr_occ, "merr_global_1.csv", append = TRUE)
st_write(data_fin, "merr_nr_1.csv", append = TRUE)

