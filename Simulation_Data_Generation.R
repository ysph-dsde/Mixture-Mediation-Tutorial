# This is the code script for Simulation Data Generation
# This would not involve any PCA, ERS, etc. transformations
# The transformations will be carried out in their individual code scripts (PCA_sim.R .... etc.)

library(tidyverse)
library(MASS)
library(CMAverse)
library(cluster)
library(factoextra)
source("Functions/Functions_DataGen.R")

set.seed(1211)

################## Effects Setting for Data Generation ######################

p <- 30
q <- 1
s <- 6

# Generating alpha_c
Alpha_c <- matrix(1, nrow = q, ncol = s)

# Generating beta_c
Beta_c <- matrix(1, nrow = s, ncol = 1)


# Generate effects
## Indirect Effects
Alpha_a <- matrix(0, nrow = 1, ncol = p)
Alpha_a[c(1, 11, 21)] <- 0.3 # weak effect
Alpha_a[c(2, 12, 22)] <- 0.6 # moderate effect
Alpha_a[c(3, 13, 23)] <- 0.9 # strong effect

## Beta_m under only 1 mediation is 1 x 1
Beta_m <- 0.4

## Direct Effect
Beta_a <- rep(c(0.3, 0, 0), times = p / 3) %>% as.matrix()


# Confounder effects on exposures
Theta_c <- matrix(rep(0.1, times = p * (s - 1)), nrow = p)


## Showcase the calculated effects
# Direct Effect
t(Beta_a)

# Indirect Effect
t(t(Alpha_a) %*% Beta_m)

# TE
t(Beta_a + t(Alpha_a) %*% Beta_m)


# Global Proportion Mediated
(t(t(Alpha_a) %*% Beta_m) %>% sum) / (t(Beta_a + t(Alpha_a) %*% Beta_m) %>% sum)

############### The raw data ######################

sce <- c("1000l", "1000h", "2500l", "2500h")
r2_M <- c(0.1, 0.4, 0.1, 0.4)
n_obs <- c(1000, 1000, 2500, 2500)
n_set <- 100

for (j in 1:length(sce)) {
  pb <- txtProgressBar(min = 0, max = n_set, style = 3)
  for (i in 1:n_set) {
    # create a list to store all the data
    
    list_df <- list()
    
    list_df[["Raw"]] <- data_gen(
      n_obs = n_obs[j],
      n_expo = p,
      n_confound = s,
      expo_blockNum =  c(5, 10, 15),
      expo_blockCorr = c(0.4, 0.8, 0.1),
      confound_blockNum = 5,
      confound_blockCorr = 0.2,
      Alpha_a = Alpha_a,
      Alpha_c = Alpha_c,
      Beta_m = Beta_m,
      Beta_a = Beta_a,
      Beta_c = Beta_c,
      Theta_c = Theta_c,
      r2_M = r2_M[j],
      r2_Y = 0.3
    )
    
    #################### PCA ############################
    
    list_df[["PCA"]] <- list_df[["Raw"]]
    
    
    #################### ERS ############################
    
    ## Case 1 of ERS should be the construction of ERS score
    # involves ENET only on the original set of exposures only
    list_df[["ERS_Case1"]] <- list_df[["Raw"]]
    
    
    ## Case 2 of ERS should be the construction of ERS score
    # involves ENET on the expanded matrix of exposures, squared exposures, and interaction terms
    list_df[["ERS_Case2"]] <- list_df[["Raw"]]
    
    
    #################### BKMR ############################
    
    list_df[["BKMR"]] <- list_df[["Raw"]]
    list_df[["BKMR (Hierarchical)"]] <- list_df[["Raw"]]
    
    ################ Output results ########################
    
    write_rds(list_df, paste0("sim_", sce[j], "/sim", i, ".rds"))
    setTxtProgressBar(pb, i)
  }
  close(pb)
}