# Script: 4-spatial-prediction.R
# Authors: Max T Eyre, Emanuele Giorgi, Peter Diggle 
# Date: 05/05/2021
# Article: 
# Explanation: In this script we make spatial predictions for the fitted parameters at unobserved locations within
# the study area. We do this for two targets: rattiness and the probability of infection. Prediction maps (used
# in the published article) are then produced.
library(plyr)
library(tmap)
library(sf)
library(splancs)
library(raster)
library(pdist)
library(ggmap)
library(mapview)
library(PrevMap)
library(RColorBrewer)
library(tidyverse)
library(paletteer)
library(colorblindr)

## Step 1 - read in data, control parameters and parameter estimates ----

# control is where multiple model arguments are stored
control <- read_csv("Data/control_R1_3Xi.csv")

# read in rat and human data
rat <- read.csv("Data/rat_data.csv")
human <- read.csv("Data/human_data.csv")

# functions
source(file = "Functions/RATT_EPI_alt_geostat_predict_xis_fn.R")
source(file = "Functions/RATT_EPI_alt_model_functions.R")
source(file = "Functions/scaling_fns.R")

# fitted model parameter values
par_hat <- readRDS("Outputs/final_par_R1_3Xi.RDS")
par_hat <- par_hat$`final estimates`$par

## Step 2 - prediction grid set-up ----
cov.pred.human <- control$human[!is.na(control$human)] 
warning("Values for the following model variables must be included in the prediction grid: ", paste0(cov.pred.human,collapse=", "))
# please ensure that all model variables have values for the prediction grid

# predicting for 30 year old male with an household per-capita income of 1 dollars/day who do not work as a travelling salesperson 
# and rarely/never has contact with floodwater.
grid.pred <- read.csv("Data/prediction_grid.csv") %>%
  mutate(age = 30, sex = 1, floodwater_freq_bin = 0, hh_income_pcap_pess30_dailydol = 1, work_salesman = 0) %>%
  mutate(floodwater_freq_bin = factor(floodwater_freq_bin,levels = c("0","1","2")),
         work_salesman = factor(work_salesman,levels = c("0","1")),
         elev_levels = case_when(
           rel_elev <= 6.7 ~ 1,
           rel_elev > 6.7 & rel_elev <= 15.6 ~ 2,
           rel_elev > 15.6 ~ 3,
         ),
         elev_levels = factor(elev_levels, levels = c("1","2","3"))) # create correct levels for prediction categorical variables

# coordinate reference system (CRS)
crs.val <- CRS("+proj=utm +zone=24 +south +datum=WGS84 +units=m +no_defs")

## Step 3 - make predictions ----
## arguments:
# 'n.sim', 'burnin' and 'thin' are MCMC arguments
# 'get.raster' = TRUE - provides outputs in raster form
# 'get.values' = TRUE - provides outputs in dataframe form 
# (beware values include samples from [T|W] and can be ~1-10GB depending on n.pred.samples used)
# 'n.pred.samples' = the number of conditional samples [T|W] used to estimate E[T|W] and Var[T|W]
# 'pred.int' = the Prediction Interval for human predictions c(0.025,0.975) for 95%CI equivalent
# 'pred.human.N' = if human nugget effect in model, then include human nugget effect in prediction target?

# structure of output:
# 1. rasters in predictions.raster:
# r.R.hat - mean predicted rattiness (expectation of all MCMC samples)
# r.R.sd - sd of predicted rattiness
# r.S.rat.hat -  mean predicted S(x) (expectation of all MCMC samples)
# r.human.hat - mean predicted prob. infection (expectation of all MCMC samples)
# r.human.sd - sd of predicted prob. infection
# r.human.hat.li - lower prediction interval for predicted prob. infection
# r.human.hat.ui - upper prediction interval for predicted prob. infection

# 2. dataframe predictions in predictions.val:
# R.pred.cond.samples - samples from predictive distribution of rattiness [T|W] (on linear rattiness scale)
# R.pred.hat - mean predicted rattiness (expectation of all MCMC samples)
# R.pred.sd - sd of predicted rattiness
# S.rat.pred.hat -  mean predicted S(x) (expectation of all MCMC samples)
# human.pred.cond.samples - samples from predictive distribution of [T|W] (transformed from log(p/(1-p)) to p)
# human.pred.hat - mean predicted prob. infection (expectation of all MCMC samples)
# human.pred.sd - sd of predicted prob. infection
# human.pred.hat.int - contains predictions for lower (row 1) and upper (row 2) prediction intervals for predicted prob. infection

pred_out <- rattiness.epi.predict.alt(control, 
                                      rat, human, par_hat,
                                      n.sim = 55000, burnin= 5000, thin = 10,
                                      grid.pred, crs.val,
                                      get.raster=TRUE, get.values=FALSE,
                                      n.pred.samples=5000,
                                      pred.int=c(0.025,0.975),
                                      pred.human.N = TRUE)