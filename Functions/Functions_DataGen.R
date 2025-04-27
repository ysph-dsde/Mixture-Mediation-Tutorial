# This is the code script to store the required functions to generate simulation data
# The original code writting process can be found in Data_Gen_Test.R
# The final data generation script is in DataGeneration_Final.R

############## Description ##################

# We will have the following functions in this code script:

## 1. gen_block_corr -- Generate Block Correlation Matrix for the Exposures

## 2. data_gen -- Generate the Simulated Data with arguments of # of obs., # of exposures, exposure block size and numbers,
#                   True values of Alpha_a, Alpha_c, Beta_m, Beta_a, and Beta_c,
#                   and desired adjusted Rsq for the mediation model and the outcome model

######### Function to generate block correlation matrix of exposures #####################

gen_block_corr <- function(exposure_numbers, correlations) {
  if (length(exposure_numbers) != length(correlations)) {
    stop("The lengths of exposure_numbers and correlations must match.")
  }
  
  # Load the Matrix library for bdiag
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "The Matrix package is required. Please install it using install.packages('Matrix')"
    )
  }
  
  # Function to create a single correlation block
  create_correlation_block <- function(size, correlation) {
    block <- matrix(correlation, nrow = size, ncol = size)
    diag(block) <- 1 # Set diagonal elements to 1
    return(block)
  }
  
  # Create each correlation block and store them in a list
  blocks <- lapply(seq_along(exposure_numbers), function(i) {
    create_correlation_block(exposure_numbers[i], correlations[i])
  })
  
  # Combine the blocks into a block diagonal matrix
  correlation_matrix <- Matrix::bdiag(blocks)
  
  # Convert to a regular matrix for compatibility
  correlation_matrix <- as.matrix(correlation_matrix)
  
  return(correlation_matrix)
}

########################## Function to generate data ##############################

data_gen <- function(n_obs,
                     n_expo,
                     n_confound,
                     expo_blockNum,
                     expo_blockCorr,
                     confound_blockNum,
                     confound_blockCorr,
                     Alpha_a,
                     Alpha_c,
                     Beta_m,
                     Beta_a,
                     Beta_c,
                     Theta_c,
                     # theta_c is what confounders contributes to exposures (q times s dim)
                     adjR2_M,
                     adjR2_Y) {
  n_obs <- n_obs
  q <- n_expo
  p <- 1 # we assume 1 mediator for now
  # We assume no confounders
  s <- n_confound + 1 # if s == 1 then it is just the intercept
  
  # Generate intercept and confounders
  if (s == 1) {
    C_i_T <- rep(1, n_obs) %>% as.matrix()
    Sigma_C <- diag(s)
  } else{
    interCept <- rep(1, n_obs) %>% as.matrix()
    
    Sigma_C <- gen_block_corr(exposure_numbers = confound_blockNum, correlations = confound_blockCorr)
    
    conFound <- MASS::mvrnorm(n = n_obs,
                              mu = rep(0, s - 1),
                              Sigma = Sigma_C) # s is the number of confounders
    C_i_T <- cbind(interCept, conFound)
  }
  
  ## Exposures
  Sigma_X <- gen_block_corr(exposure_numbers = expo_blockNum, correlations = expo_blockCorr)
  
  # generate the exposures
  X <- t(Theta_c %*% t(conFound)) + MASS::mvrnorm(n = n_obs,
                                                  mu = rep(0, q),
                                                  Sigma = Sigma_X) #  MASS::mvrnorm(n = n_obs, mu = rep(0, q), Sigma = Sigma_X)
  colnames(X) <- paste0("x", 1:q)
  
  # we can calculate the Sigma_M (a scalar here)
  adjR2_M <- adjR2_M
  r2_M <- 1 - ((n_obs - q - s - 1) / (n_obs - 1)) * (1 - adjR2_M) # set adjusted r-squared to 0.3
  
  
  # some prep
  big_alpha <- cbind(Alpha_a, Alpha_c)
  colnames(big_alpha) <- NULL
  
  V_mat_alpha <- matrix(0, nrow = (q + s), ncol = (q + s))
  V_mat_alpha[1:q, 1:q] <- var(X)
  V_mat_alpha[(q + 1):(q + s), (q + 1):(q + s)] <- var(C_i_T)
  
  
  Sigma_M <- ((1 - r2_M) / (r2_M)) * big_alpha %*% V_mat_alpha %*% t(big_alpha)
  
  # Generate M Mediators
  M <- t(Alpha_a %*% t(X) + Alpha_c %*% t(C_i_T)) + MASS::mvrnorm(n = n_obs, mu = 0, Sigma = Sigma_M) # mu + error
  
  #Generate sigma^2_e (error variance of the outcome model)
  # big_beta
  big_bt <- cbind(Beta_m, t(Beta_a), t(Beta_c))
  colnames(big_bt) <- NULL
  
  # build the V = var covar of M, X
  V_mat <- matrix(0, nrow = (p + q + s), ncol = (p + q + s))
  V_mat[1:p, 1:p] <- var(M)
  V_mat[(p + 1):(p + q), (p + 1):(p + q)] <- var(X)
  V_mat[(p + q + 1):(p + q + s), (p + q + 1):(p + q + s)] <- var(C_i_T)
  
  
  # Calculate the sigma_y
  adjR2_Y <- adjR2_Y
  r2_Y <- 1 - ((n_obs - p - q - s - 1) / (n_obs - 1)) * (1 - adjR2_Y) # set adjusted r-squared to 0.3
  
  # calculate the optimal sigma^2_e for the two cases
  Sigma_Y <- ((1 - r2_Y) / (r2_Y)) * big_bt %*% V_mat %*% t(big_bt)
  
  #Generate Y
  comb_predictors <- cbind(M, X, C_i_T)
  colnames(comb_predictors)[seq(p)] <- paste("m", seq(p))
  colnames(comb_predictors)[p + seq(q)] <- paste("x", seq(q))
  colnames(comb_predictors)[p + q + seq(s)] <- paste(c("intercept", paste("c", seq(s -
                                                                                     1))))
  
  # if there are more than one confounders
  #ifelse(s > 1,   colnames(comb_predictors)[p + q + 1 + seq(s)] <- paste("c", seq(s-1)),
  #      0)
  
  y1 <- t(big_bt %*% t(comb_predictors)) + MASS::mvrnorm(n = n_obs, mu = 0, Sigma = Sigma_Y) # mu + error
  
  #Combine the results
  df_gen <- cbind(y1, M, X, C_i_T) %>% as.data.frame() #%>% dplyr::rename(y = V1, m1 = V2, intercept = V33)
  
  colnames(df_gen) <- c("y",
                        paste0("m", seq(p)),
                        paste0("x", seq(q)),
                        "intercept",
                        paste0("c", seq(s - 1)))
  
  return(df_gen)
}
