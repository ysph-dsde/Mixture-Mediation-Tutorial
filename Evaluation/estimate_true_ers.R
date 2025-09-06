# This is the code script to run the ERS simulation for Case 1
# ERS case 1 is the case where the more simple form of the ERS score construction is used
# which considers only the exposures and nothing else when running the Elastic Net variable selection

library(tidyverse)
library(MASS)
library(gcdnet)
library(purrr)
library(CMAverse)
library(knitr)
library(kableExtra)
source("functions_ERS.R")

set.seed(1211)
sce <- c("l", "h")

for(level in sce){
  # Read in the data
  list_df <- readRDS(paste0("sim_", level, ".rds"))
  list_medTest <- list()
  
  #### Generate the ERS score ##############################
  
  # The Case 1 of the ERS score calculation
  # includes the exposures (main effects) only
  
  list_df[["ERS_Case1"]] <- ers_Calc(
    data = list_df[["Raw"]],
    exposure = paste0("x", 1:30),
    outcome = "y",
    covar = paste0("c", 1:5),
    include_int = F,
    lambda2_start = exp(seq(log(1e-4), log(1e2), length.out = 100)),
    seed = 1211
  )
  
  post_ERS_data <- list_df[["ERS_Case1"]]$post_ERS_data
  
  #### Test and Estimate the Mediation Effect of the ERS score  #########################
  set.seed(1211)
  
  astar <- quantile(post_ERS_data[, "ERS"], 0.25)
  a <- quantile(post_ERS_data[, "ERS"], 0.75)
  ## Mediation Test on ERS_scores
  list_medTest[["ERS_Case1"]][["CMA Test"]] <- CMAverse::cmest(
    data = post_ERS_data,
    model = "rb",
    full = T,
    EMint = F,
    yreg = "linear",
    mreg =  list("linear"),
    mval = list(1),
    basec = paste0("c", seq(5)),
    outcome = "y",
    exposure = "ERS",
    mediator = "m1",
    a = a,
    astar = astar,
    estimation = "paramfunc",
    inference = "delta"
  )
  
  list_medTest[["ERS_Case1"]][["CMA Summary Table"]] <- cbind(
    list_medTest[["ERS_Case1"]][["CMA Test"]]$effect.pe,
    list_medTest[["ERS_Case1"]][["CMA Test"]]$effect.se,
    list_medTest[["ERS_Case1"]][["CMA Test"]]$effect.ci.low,
    list_medTest[["ERS_Case1"]][["CMA Test"]]$effect.ci.high,
    list_medTest[["ERS_Case1"]][["CMA Test"]]$effect.pval
  )
  colnames(list_medTest[["ERS_Case1"]][["CMA Summary Table"]]) <- 
    c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")
  
  #### Output Results ########
  saveRDS(list_df, paste0("sim_", level, ".rds"))
  saveRDS(list_medTest, paste0("ers1_", level, ".rds"))
}
