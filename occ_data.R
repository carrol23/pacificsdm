# --------- SPECIES OCCURENCE DATA CLEANING -----------
# This script is used to download, import and clean species occurence data
# for the use in species distribution models

# clean out workspace
rm(list=ls())

# set working directory
setwd('C:/Data/sdm/C2_SDM/occurence/at')


# install relevant packages
# uncommented lines of code are present as this script has been rerun

packages = c('devtools', 'rgbif', 'dplyr', 'terra')
install.packages(packages)

library(CoordinateCleaner)
library(rgbif)
library(dplyr)
library(terra)
library(ggplot2)
library(sf)
library(raster)

library(devtools)
install_github("ropensci/CoordinateCleaner")

# -------------- IMPORT SPECIES OCCURRENCE DATA --------------
# Download GBIF Species Occurrence Data for the African Tulip
atulip_data = occ_search(scientificName = "Spathodea campanulata", limit = 10000, hasCoordinate = TRUE)

# select data from list
atulip_data <- atulip_data$data

# scale data down
atulip_occ <- atulip_data %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)

# export occurrence data as CSV - plot as spatial data
st_write(atulip_occ, "gbif_at.csv", append = FALSE)

# ----------------- IMPORT DATA -----------------
# Import CSV files of occurrence data
at_ala = read.csv("ala_records.csv", sep = ',') # ALA

# select columns from the ALA dataset that are relevant for the data cleaning process
at_ala = at_ala %>%
    select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
          dataResourceUid, family, taxonRank, coordinateUncertaintyInMeters, year,
          basisOfRecord, institutionCode, datasetName)

# rename gbif dataset
names(atulip_occ)[names(atulip_occ) == 'gbifID'] <- 'dataResourceUid'

# subset the occurrence dataset with selected columns
atulip_occ = atulip_occ %>%
  select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
          gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
          basisOfRecord, institutionCode, datasetName)

atulip_occ <- rbind(atulip_occ, at_ala)

# -------------- PREPROCESSING STEPS --------------
# remove records without coordinates - records should be the same
atulip_occ <- atulip_occ %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))

# identify duplicates based on lat and long coordinates
atulip_d = duplicated(atulip_occ[,c("decimalLongitude", "decimalLatitude")])

# count number of duplicates - 1003 duplicates
length(atulip_d[atulip_d == TRUE])

# remove duplicates
atulip_occ <- atulip_occ[!atulip_d,]

# convert to dataframe
atulip_occ <- data.frame(atulip_occ)

# plot using ggplot
wrld = borders("world", colour = "gray50", fill = "gray50")
ggplot() +  wrld +
  geom_point(data = atulip_occ, aes(x=decimalLongitude, y=decimalLatitude),
             colour = "red", size = 0.5) +
  coord_sf()

# use coordinate cleaner to flag problem records
install.packages('countrycode')
library(countrycode)
# convert country code ISO2 to ISO3
atulip_occ$countryCode <- countrycode(atulip_occ$countryCode,
                                       origin = 'iso2c',
                                       destination = 'iso3c')

# flag coordinates not 'clean' based on tests
flags <- clean_coordinates(x = atulip_occ,
                           lon = 'decimalLongitude', lat = 'decimalLatitude',
                           countries = 'countryCode', species = 'species',
                           tests = c('capitals', 'centroids', 'duplicates', 'equal', 'gbif',
                                     'institutions', 'zeros', 'outliers', 'seas'))

summary(flags)

plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# exclude flagged records
atulip_occ <- atulip_occ[flags$.summary,]
# flagged records
atulip_flag <- atulip_occ[!flags$.summary,]

# use metadata to remove data
# remove records with low precision
hist(atulip_occ$coordinateUncertaintyInMeters / 1000, breaks = 5)

# removed coordinates with uncertainty up to 1km (using 0.5mins)
atulip_occ <- atulip_occ %>%
  filter(coordinateUncertaintyInMeters/1000 <= 1 | is.na(coordinateUncertaintyInMeters))

# identify basis of records - remove fossil, living & preserved specimen
table(atulip_occ$basisOfRecord)

atulip_occ <- filter(atulip_occ, basisOfRecord == "HUMAN_OBSERVATION" |
                        basisOfRecord == "LIVING_SPECIMEN" |
                        basisOfRecord == "MATERIAL_CITATION" |
                        basisOfRecord == "OBSERVATION" |
                        basisOfRecord == "OCCURENCE")

table(atulip_occ$individualCount)

# remove individual counts that equate to 0 and higher than 100
atulip_occ <- atulip_occ %>%
  filter(individualCount > 0 | is.na(individualCount)) %>%
  filter(individualCount < 99 | is.na(individualCount))

# exclude old records- 1925 remove
table(atulip_occ$year)

# plot records
ggplot(atulip_occ, aes(x = year)) +
  geom_histogram(binwidth=1) +
  ggtitle("Year of African Tulip Occurrences")

atulip_occ <- atulip_occ %>%
  filter(year > 2010)

# check family name
table(atulip_occ$family)

# check taxon rank
table(atulip_occ$taxonRank)

ggplot() +
  geom_point(data = atulip_occ, aes(x=decimalLongitude, y = decimalLatitude), color = 'blue')

# ------------- Problematic datasets -------------
# ddmm to dd.dd conversion error - 0 records
out.ddmm <- cd_ddmm(atulip_occ, lon = "decimalLongitude",
                    lat = "decimalLatitude",
                    ds = "species", diagnostic = T,
                    diff = 1, value = "dataset")

#Flags datasets with periodicity patterns indicative of a rasterized (lattice)
#collection scheme, as often obtain from e.g. atlas data. Using a combination of autocorrelation
#and sliding-window outlier detection to identify periodicity patterns in the data.
#See https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152 for
#further details and a description of the algorithm

par(mfrow = c(2,2), mar = rep(2, 4))
out.round <- cd_round(atulip_occ, lon = "decimalLongitude", lat = "decimalLatitude",
                      ds = "species", value = "dataset",
                      T1 = 7, graphs = T)

# export occurrence data (not thinned) as CSV - plot as spatial data
st_write(atulip_occ, "comb_at.csv", append = FALSE)
