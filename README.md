
# Species Distribution Modelling - Invasive Weeds in the Pacific
Species distribution models (maxent) for three invasive weeds across the Pacific region. The model utilises worldclim variables topographic, FAO global soil organic carbon (GSOC) and the 2019 Copernicus landcover dataset as environmental predictors and occurence data sourced from GBIF directly accessed in R from R Studio through the RGBIF package.

Occurence or presence data is cleaned using the R package Coordinate Cleaner.

The main code has been published here in three main scripts and is reusable based on the species of study:
1. <a href = "https://github.com/carrol23/pacificsdm/blob/main/merr_occurence_data_download.R">Download species presence data (example species Merremia), clean and save </a>
2. <a href = "https://github.com/carrol23/pacificsdm/blob/main/merr_predictor_data_download.R">Download environmental predictors/variables </a>
3. <a href = "https://github.com/carrol23/pacificsdm/blob/main/merr_model_prediction_eval.R"> Train and test Maxent model and run evaluation (AUC and CBI) </a>


The scripts have been updated as packages such as maptools, rgdal and enmsdm have been deprecated. Updated on the 30th of November 2023.

This repository is a part of a PhD research project undertaken at the Univeristy of Newcastle with funding support from SPREP MISCCAP project.
