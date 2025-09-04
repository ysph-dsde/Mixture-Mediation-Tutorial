# This is the code script to run the ERS simulation for Case 2
# ERS case 2 is the case where the more complex form of the ERS score construction is used
# which considers the interactions and squared terms of each exposures when running the Elastic Net variable selection

library(Rmpi)
library(tidyverse)
library(MASS)
library(gcdnet)
library(purrr)
library(CMAverse)
library(knitr)
library(kableExtra)

source("Functions/Functions_ERS.R")

rank <- mpi.comm.rank(0)
size <- mpi.comm.size(0)

set.seed(1211)
N <- 100
sce <- Sys.getenv("sce")

tasks_per_rank <- split(1:N, cut(1:N, size, labels = FALSE))
my_tasks <- tasks_per_rank[[rank + 1]]

for (j in my_tasks) {
  # Read in the data
  list_df <- readRDS(paste0("sim_", sce, "/sim", j, ".rds"))
  list_medTest <- list()
  
  # The Case 2 of the ERS score calculation
  # includes the exposures (main effects), interactions, and squared terms
  
  # The functions for calculating the ERS score is outlined in the "Functions_ERS_alt.R"
  tmp_df <- list_df[["ERS_Case2"]]
  list_df[["ERS_Case2"]] <- NULL
  list_df[["ERS_Case2"]][["ERS Result"]] <- ers_Calc(
    data = tmp_df,
    exposure = paste0("x", 1:30),
    outcome = "y",
    covar = paste0("c", 1:5),
    include_int = T,
    lambda2_start = exp(seq(log(1e-4), log(1e2), length.out = 100)),
    seed = 1211
  )
  
  list_df[["ERS_Case2"]][["Data"]] <- list_df[["ERS_Case2"]][["ERS Result"]]$post_ERS_data
  
  list_df[["ERS_Case2"]][["ERS Result"]]$ers_fit$lambda # the best lambda value for this
  
  analysis_ids <- which(list_df[["ERS_Case2"]][["Data"]][, "ERS"] != 0)
  
  ############################# Test and Estimate the Mediation Effect of the ERS score  #########################
  set.seed(1211)
  
  astar <- quantile(list_df[["ERS_Case2"]][["Data"]][analysis_ids, "ERS"], 0.25)
  a <- quantile(list_df[["ERS_Case2"]][["Data"]][analysis_ids, "ERS"], 0.75)
  ## Mediation Test on ERS_scores
  list_medTest[["ERS_Case2"]][["CMA Test"]] <- CMAverse::cmest(
    data = list_df[["ERS_Case2"]][["Data"]][analysis_ids, ],
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
  
  list_medTest[["ERS_Case2"]][["CMA Summary Table"]] <- cbind(
    list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.pe,
    list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.se,
    list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.ci.low,
    list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.ci.high,
    list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.pval
  )
  colnames(list_medTest[["ERS_Case2"]][["CMA Summary Table"]]) <- 
    c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")
  
  #### Output Results ########
  saveRDS(list_df, paste0("sim_", sce, "/sim", j, ".rds"))
  saveRDS(list_medTest, paste0(sce, "/ers2_", j, ".rds"))
  print(paste0("Dataset ", j, "/100 Finished"))
}

mpi.barrier(0)
mpi.quit()
