# R-packages 
library(tidyverse)
library(parallel)
library(bnlearn)
# Estimation 
library(nnet) # Multinomial logit
library(PSweight) # AIPW 
# - Non parametric 
library(ranger)
# Figures 
library(gridExtra) 
library(grid) 
library(latex2exp) # For mathematical notation in figures 
# Tables 
library(kableExtra)


# Source files 
source("DGP.R") # Functions for data generation 
source("estimators.R") # Estimators 
source("figures.R") # Creating figures
source("tables.R") # Creating tables
source("run_datalearner.R") # Running data learner and creating all figures and tables

sessionInfo()


#########################
###   PARAMETERS    #####
#########################

### SCENARIO I 
# Marginal probabilities 
marg_prob_dmt <- c(0.19, 0.19, 0.20, 0.25)

### SCENARIO II 
# Marginal probabilities 
marg_prob_smok <- c(0.205, 0.24)

################################
#####   MAIN SIMULATION   #####
###############################

start.time <- Sys.time()

monte_carlo_sim <- run_datalearner(seed=32589189, Rs=1000, n=3470, marg_prob_dmt, marg_prob_smok)

end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


saveRDS(monte_carlo_sim, file="simulation.RData")

save(monte_carlo_sim, file="sim_results.RData")

### SCENARIO I 
# Figure of propensity score distributions under the correct model
grid.draw(monte_carlo_sim$ol_dmt_corr)
# Figure of propensity score distributions under the incorrect model
grid.draw(monte_carlo_sim$ol_dmt_err)
grid.draw(monte_carlo_sim$ol_dmt_np)
# Figure of odds ratio estimates 
grid.draw(monte_carlo_sim$bp_dmt)

# Table of bias and MSE
monte_carlo_sim$res_dmt

# Table of Monte Carlo standard error of bias and MSE
monte_carlo_sim$mcse_dmt


### SCENARIO II
# Figure of propensity score distributions under the correct model
monte_carlo_sim$ol_smok_corr
# Figure of propensity score distributions under the incorrect model
monte_carlo_sim$ol_smok_err
monte_carlo_sim$ol_smok_np
# Figure of odds ratio estimates 
grid.draw(monte_carlo_sim$bp_smok)

# Table of bias and MSE
monte_carlo_sim$res_smok

# Table of Monte Carlo standard error of bias and MSE
monte_carlo_sim$mcse_smok


####################################
###   LARGE SAMPLE PROPERTIES    ###
####################################

ls_tables <- run_datalearner(seed=32589189, Rs=1, n=1000000, marg_prob_dmt, marg_prob_smok, ls=T)

### SCENARIO I 
ls_tables$ls_dmt

### SCENARIO II 
ls_tables$ls_smok

#########################
### SAVE RESULTS      ###
#########################


saveRDS(ls_tables, file="ls_simulation.RData")

#load the results by:

readRDS("simulation.RData")


load("simulation.RData")
