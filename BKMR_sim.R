# This is the code script for running BKMR simulation

library(Rmpi)
library(tidyverse)
library(MASS)
library(bkmr)
library(causalbkmr)
library(purrr)
library(CMAverse)
library(knitr)
library(kableExtra)

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
  
  ################### Single Analysis Screening ################################
  exposures_X <- list_df[["BKMR"]] %>%
    dplyr::select(contains("x")) %>%
    as.matrix()
  
  confounders_C <- list_df[["BKMR"]] %>%
    dplyr::select(starts_with("c")) %>%
    as.matrix()
  
  outcome_Y  <- list_df[["BKMR"]]$y
  mediator_M  <- list_df[["BKMR"]]$m1
  not_sig <- c()
  for (i in 1:ncol(exposures_X)) {
    fit <- lm(outcome_Y ~ exposures_X[, i] + confounders_C)
    p_val <- summary(fit)$coefficients[2, 4]
    if (p_val > 0.1) {
      not_sig <- c(not_sig, i)
    }
  }
  
  exposures_X <- exposures_X[,-not_sig]
  
  ################### BKMR - Component-wise Variable Selection #################
  
  # we need to first define the $Z_M$ and $Z_Y$ matrices for the BKMR mediation and outcome models.
  #Recall that $Z_M$ is the set of exposures and effect modifiers $E_M$
  
  # create Z.M and Z.Y and Z.YM
  #Z.M <- cbind(A, E.M); Z.Y <- cbind(A, E.Y); Zm.Y <- cbind(Z.Y, m)
  
  # we assume no effect modifiers
  E.M <- NULL
  E.Y <- NULL
  
  # create the matrices
  Z.M <- cbind(exposures_X, E.M)
  Z.Y <- cbind(exposures_X, E.Y)
  Zm.Y <- cbind(Z.Y, mediator_M)
  
  # outcome model
  set.seed(1211)
  fit.y <- kmbayes(
    y = outcome_Y,
    Z = Zm.Y,
    X = confounders_C,
    iter = 10000,
    verbose = TRUE,
    varsel = TRUE,
    control.params = list(
      lambda.jump = 0.45,
      r.jump1 = 0.1,
      r.muprop = 0.35,
      r.jump2 = 0.25
    )
  )
  # 9.4 hours
  
  # TE model
  set.seed(1211)
  fit.y.TE <- kmbayes(
    y = outcome_Y,
    Z = Z.Y,
    X = confounders_C,
    iter = 10000,
    verbose = TRUE,
    varsel = TRUE,
    control.params = list(lambda.jump = 0.375, r.jump2 = 0.05)
  )# 11.1 hours
  
  # mediator model
  set.seed(1211)
  fit.m <- kmbayes(
    y = mediator_M,
    Z = Z.M,
    X = confounders_C,
    iter = 10000,
    verbose = TRUE,
    varsel = TRUE,
    control.params = list(lambda.jump = 0.925, r.jump2 = 0.025)
  ) # 8.5 hours
  
  tmp_df <- list_df[["BKMR"]]
  
  list_df[["BKMR"]] <- NULL
  list_df[["BKMR"]][["Data"]] <- tmp_df
  list_df[["BKMR"]][["Screen"]] <- colnames(exposures_X)
  
  list_df[["BKMR"]][["Outcome model"]] <- fit.y
  list_df[["BKMR"]][["TE model"]] <- fit.y.TE
  list_df[["BKMR"]][["Mediator model"]] <- fit.m
  
  ######################### Causal BKMR Code (astar = 25th quantile; a = 75th quantile) #########################
  
  # set the confounders mean level
  X.predict <- matrix(colMeans(confounders_C), nrow = 1)  # take the mean of each confounders
  
  
  # Define the change in exposure levels for mediation effects estimation
  # consider a change in all exposures from their 25th to 75th percentiles
  
  astar <- c(apply(exposures_X, 2, quantile, probs = 0.25)) # the reference level of the exposures at 25th percentile
  a <- c(apply(exposures_X, 2, quantile, probs = 0.75)) # the comparative level at 75th percentile
  
  # The index of the MCMC iterations to be used for inference
  sel <- seq(5001, 10000, by = 10)
  
  #medTest_BKMR
  medTest_BKMR <- mediation.bkmr(
    a = a,
    astar = astar,
    e.m = NULL,
    e.y = NULL,
    fit.m = fit.m,
    # mediation model
    fit.y = fit.y,
    # outcome model
    fit.y.TE = fit.y.TE,
    # TE model
    X.predict.M = X.predict,
    # the mean confounder level for mediation model
    X.predict.Y = X.predict,
    # the mean confounder level for the outcome model
    m.quant = c(0.1, 0.25, 0.5, 0.75),
    # the quantile values of mediator which the CDE(m) is fixed to
    # m.value, # the specific values of the mediator for estimating CDE(m)
    alpha = 0.05,
    # 95% CI
    sel = sel,
    # MCMC iteration indices for making inference
    seed = 1211,
    K = 50
  ) # K is the number of MCMC samples of the mediator
  
  # 1.8 Hours (K = 50)
  
  list_medTest[["BKMR"]][["BKMR_CMA Results"]] <- medTest_BKMR
  list_medTest[["BKMR"]][["Results Table"]] <- medTest_BKMR$est
  
  ######################################## BKMR Hierarchical ##################################################
  
  #######################################################################
  # Hierarchical Clustering of the exposures
  #######################################################################
  
  cor_mat <- cor(exposures_X, method = "pearson")
  
  # hierarchical clustering
  hc <- hclust(as.dist(1 - cor_mat))
  
  # assign groupings
  groups <- cutree(hc, k = 3)
  
  #######################################################################
  # Fitting the BKMR models
  #######################################################################
  
  # The argument `y` is a vector of outcome. `Z` is a $n \times q$ matrix of predictors to be included in the $h(\cdot)$ function.
  # `X` is an $n \times s$ matrix of confounders and should not include the intercept.
  # `iter` is the number of iterations to run the MCMC.
  # `varsel` is a Boolean variable of whether to conduct variable selection on the $Z$ in $h(\cdot)$ function
  # `groups` is an optional vector of group indicators (of length $q$) for fitting the hierarchical variable selection, given `varsel = T`.
  # If `group = NULL` and `varsel = T`, then the component-wise variable selection will be performed
  
  set.seed(1211)
  
  # outcome model
  fit_y_hier <- kmbayes(
    y = outcome_Y,
    Z = Zm.Y,
    X = confounders_C,
    iter = 10000,
    verbose = TRUE,
    varsel = TRUE,
    groups = c(groups, 4),
    control.params = list(lambda.jump = 0.45, r.jump2 = 0.025)
  ) # 9.7 hrs
  
  
  # TE model
  set.seed(1211)
  fit_y_TE_hier <- kmbayes(
    y = outcome_Y,
    Z = Z.Y,
    X = confounders_C,
    iter = 10000,
    verbose = TRUE,
    varsel = TRUE,
    groups = groups,
    control.params = list(lambda.jump = 2, r.jump2 = 0.01)
  )# 8.4 hrs
  
  # mediator model
  set.seed(1211)
  fit_m_hier <- kmbayes(
    y = mediator_M,
    Z = Z.M,
    X = confounders_C,
    iter = 10000,
    verbose = TRUE,
    varsel = TRUE,
    groups = groups,
    control.params = list(
      lambda.jump = 3.25,
      r.jump1 = 0.01,
      r.jump2 = 0.01
    )
  )# 8.4 hrs
  
  list_df[["BKMR (Hierarchical)"]] <- NULL
  
  list_df[["BKMR (Hierarchical)"]][["Outcome model"]] <- fit_y_hier
  list_df[["BKMR (Hierarchical)"]][["TE model"]] <- fit_y_TE_hier
  list_df[["BKMR (Hierarchical)"]][["Mediator model"]] <- fit_m_hier
  
  
  ######################### Causal BKMR Code (astar = 25th quantile; a = 75th quantile) #########################
  
  #######################################################################
  # Causal BKMR Code (astar = 25th quantile; a = 75th quantile)
  ########################################################################
  
  medTest_BKMR <- mediation.bkmr(
    a = a,
    astar = astar,
    e.m = NULL,
    e.y = NULL,
    fit.m = list_df[["BKMR (Hierarchical)"]][["Mediator model"]],
    fit.y = list_df[["BKMR (Hierarchical)"]][["Outcome model"]],
    fit.y.TE = list_df[["BKMR (Hierarchical)"]][["TE model"]],
    X.predict.M = X.predict,
    # the mean confounder level for mediation model
    X.predict.Y = X.predict,
    # the mean confounder level for the outcome model
    m.quant = c(0.1, .25, 0.5, 0.75),
    # the quantile values of mediator which the CDE(m) is fixed to
    # m.value, # the specific values of the mediator for estimating CDE(m)
    alpha = 0.05,
    # 95% CI
    sel = sel,
    # MCMC iteration indices for making inference
    seed = 1211,
    K = 50
  ) # K is the number of MCMC samples of the mediator
  
  # 1.6 HRs (K = 50)
  
  list_medTest[["BKMR (Hierarchical)"]][["BKMR_CMA Results"]] <- medTest_BKMR
  list_medTest[["BKMR (Hierarchical)"]][["Results Table"]] <- medTest_BKMR$est
  
  
  #### Output Results ########
  saveRDS(list_df, paste0("sim_", sce, "/sim", j, ".rds"))
  saveRDS(list_medTest, paste0("~/palmer_scratch/", sce, "/bkmr", j, ".rds"))
  print(paste0("Dataset ", j, "/100 Finished"))
}

mpi.quit()