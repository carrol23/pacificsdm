
# Species Distribution Modelling - Invasive Weeds in the Pacific
This repository consists of R code that explored four species distribution models for the African tulip across the Pacific region. The model utilises worldclim climate and elevation as environmental predictors and occurence data sourced from GBIF, RAINBIO, EDDMaps, etc which has been directly downloaded and where possible directly accessed in R from R Studio through the RGBIF package.

Occurence or presence data is cleaned using the R package Coordinate Cleaner.

The main code has been published here in four main scripts and is reusable based on the species of study:
1. download species occurence data (and data cleaning)
2. download environmental covariates (env predictors)
3. sdm data preprocessing
4. sdm model fit, train and evaluation


The scripts have been updated as packages such as maptools, rgdal and enmsdm have been deprecated. Updated on the 1st September 2024.

This repository is a part of a PhD research project undertaken at the Univeristy of Newcastle with funding support from SPREP MISCCAP project.
