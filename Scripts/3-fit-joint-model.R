# Script: 3-fit-joint-model.R
# Authors: Max T Eyre, Emanuele Giorgi, Peter Diggle 
# Date: 05/05/2021
# Article: 
# Explanation: In this script we fit the joint rattiness-infection model to the rat and human data.

library(tidyverse)
library(PrevMap)
set.seed(16000)

# control is where multiple model arguments are stored
control <- read_csv("Data/control_R1_3Xi.csv")

# read in rat and human data
rat <- read.csv("Data/rat_data.csv")
human <- read.csv("Data/human_data.csv")

# functions
source(file = "Functions/RATT_EPI_alt_geostat_model_xis_fn.R")
source(file = "Functions/RATT_EPI_alt_model_functions.R")
source(file = "Functions/scaling_fns.R")
source(file = "Functions/RATT_EPI_alt_geostat_bootstrap_xis_fn.R")

## fitting parameters (iterative method plug-in method)
# 'euclid.norm.method' = TRUE - run until euclid norm of standardised parameter values falls below 'tol'
# 'euclid.norm.method' = FALSE - plug-in and re-run 'n.iter.fit' times
# 'tol' - tolerance for euclid norm of standardised parameter values
# 'n.iter.fit' - use this if euclid.norm.method is set to FALSE, otherwise has no effect.
# 'rel.tol' and 'iter.max' are standard nlminb arguments
# 'n.sim', 'burnin' and 'thin' are MCMC arguments
# 'MCMC.swap' = TRUE - MCMC parameters switch to the second value in 'n.sim' and 'burnin' at iteration = 'n.MCMC.swap'
# 'output.par' = TRUE - output parameter values after every plug-in
# Note: fitting of this complex model is slow and we recommend leaving it to run on a high-end computing network

output <- rattiness.epi.model.alt(control, rat, human,
                                  euclid.norm.method=TRUE, tol=0.1, 
                                  n.iter.fit=10, rel.tol=1e-5, iter.max = 200,
                                  n.sim=c(22000,110000), burnin=c(2000,10000), 
                                  thin=10, MCMC.swap=TRUE, n.MCMC.swap=2, output.par=TRUE)

## Step 2 - example of one bootstrap ----

## fitting parameters (iterative method plug-in method)
# 'euclid.norm.method' = TRUE - run until euclid norm of standardised parameter values falls below 'tol'
# 'euclid.norm.method' = FALSE - plug-in and re-run 'n.iter.fit' times
# 'tol' - tolerance for euclid norm of standardised parameter values
# 'n.iter.fit' - use this if euclid.norm.method is set to FALSE, otherwise has no effect.
# 'rel.tol' and 'iter.max' are standard nlminb arguments
# 'n.sim', 'burnin' and 'thin' are MCMC arguments
# 'MCMC.swap' = TRUE - MCMC parameters switch to the second value in 'n.sim' and 'burnin' at iteration = 'n.MCMC.swap'
output <- rattiness.epi.bootstrap.alt(control, rat, human, par_hat,
                                      euclid.norm.method=TRUE, tol=0.1, 
                                      n.iter.fit=50, rel.tol=1e-5, iter.max = 200,
                                      n.sim=c(11000,110000), burnin=c(1000,10000), 
                                      thin=10, MCMC.swap=TRUE, n.MCMC.swap=4)


