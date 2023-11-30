
# Species Distribution Modelling - Invasive Weeds in the Pacific
Species distribution models (maxent) for three invasive weeds across the Pacific region. The model utilises worldclim variables topographic, FAO global soil organic carbon (GSOC), and terrestrial biome data as environmental predictors and occurence data sourced from GBIF directly accessed in R from R Studio through the RGBIF package.

Occurence or presence data is cleaned using the R package Coordinate Cleaner.

The main code has been published here in three main scripts and is reusable based on the species of study:
1. Download species presence data (example species Merremia), clean and save
2. Download environmental predictors/variables
3. Train and test Maxent model and run evaluation (AUC and CBI)


The scripts have been updated as packages such as maptools, rgdal and enmsdm have been deprecated. Updated on the 30th of November 2023.

This repository is a part of a PhD research project undertaken at the Univeristy of Newcastle.
