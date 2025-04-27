# This is the code script for Simulation Data Generation
# This would not involve any PCA, ERS, etc. transformations
# The transformations will be carried out in their individual code scripts (PCA_sim.R .... etc.)

pkgs <- c(
  "tidyverse",
  "MASS",
  "CMAverse",
  "corrplot",
  "cluster",
  "factoextra",
  "ggbiplot",
  "ggfortify",
  "tictoc"
)

invisible(lapply(pkgs, library, character.only = T))

source("Functions/Functions_DataGen.R")

set.seed(1211)

################## Some Preset Effects for Data Generation ######################

n_obs <- 2000
p <- 1
q <- 30
s <- 6

# Generating alpha_c
Alpha_c <- matrix(1, nrow = p, ncol = s)

# Generating beta_c
Beta_c <- matrix(1, nrow = s, ncol = 1)


# Generate effects
## Indirect Effects
Alpha_a <- matrix(0, nrow = 1, ncol = q)
Alpha_a[c(1, 11, 21)] <- 1 # weak effect
Alpha_a[c(2, 12, 22)] <- 4 # moderate effect
Alpha_a[c(3, 13, 23)] <- 8 # strong effect

## Beta_m under only 1 mediation is 1 x 1
Beta_m <- 0.5

## Direct Effect
Beta_a <- rep(c(5, 0, 0), times = q / 3) %>% as.matrix()


# Confounder effects on exposures
Theta_c <- matrix(rep(0.1, times = q * (s - 1)), nrow = 30)


## Showcase the calculated effects
# Direct Effect
t(Beta_a)

# Indirect Effect
t(t(Alpha_a) %*% Beta_m)

# TE
t(Beta_a + t(Alpha_a) %*% Beta_m)


# Global Proportion Mediated
(t(t(Alpha_a) %*% Beta_m) %>% sum) / (t(Beta_a + t(Alpha_a) %*% Beta_m) %>% sum) # around 28%

############### The raw data ######################

# create a list to store all the data
list_df <- list()

list_df[["Raw"]] <- data_gen(
  n_obs = 2000,
  n_expo = 30,
  n_confound = 5,
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
  adjR2_M = 0.3,
  adjR2_Y = 0.3
)

## correlation plot of the raw data

corrplot(
  cor(list_df[["Raw"]][, !(colnames(list_df[["Raw"]]) == "intercept")]),
  order = "original",
  method = "square",
  type = "full",
  number.cex = 0.75,
  diag = T,
  cl.pos = "n",
  # addCoef.col = "black", # this is to add numbers to it
  tl.srt = 90,
  tl.col = "black",
  tl.offset = 0.9,
  tl.cex = 0.8,
  mar = c(0, 0, 3, 0)
)
rect(
  xleft = 3 - 0.5,
  ybottom = 35 + 0.5,
  xright = 7 + 0.5,
  ytop = 31 - 0.5,
  border = "black",
  lwd = 2
)
rect(
  xleft = 17 + 0.5,
  ybottom = 30 + 0.5,
  xright = 8 - 0.5,
  ytop = 21 - 0.5,
  border = "black",
  lwd = 2
)
rect(
  xleft = 32 + 0.5,
  ybottom = 20 + 0.5,
  xright = 18 - 0.5,
  ytop = 6 - 0.5,
  border = "black",
  lwd = 2
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

write_rds(list_df, "RDS/list_df.rds")
