############## Function - ers.enet_adapt #####################

# A subfunction inside the wrapper function that
# performs the Elastic net algorithm on the exposure data

ers_enet_adapt = function(x,
                          y,
                          lambda2,
                          nfolds = 5,
                          foldid,
                          pf = rep(1, 35),
                          pf2 = rep(1, 35),
                          method = 'ls',
                          n_confound) {
  # a function to count the number of selected variables
  count_selected_expos <- function(model, n_confound) {
    sum(gcdnet::coef(model) != 0)  - n_confound - 1 # for the intercept
  }
  
  # 5 fold CV Elastic Net over the range of lambda2
  cv_lambda2 <- pbvapply(lambda2, function(lambda) {
    min(
      cv.gcdnet(
        x = x,
        y = y,
        lambda2 = lambda,
        nfolds = nfolds,
        foldid = foldid,
        pf = pf,
        pf2 = pf2,
        method = method
      )$cvm
    )
  }, FUN.VALUE = numeric(1))
  
  # Find the Optimal lambda2 and lambda1
  cv_lambda2_min <- lambda2[which.min(cv_lambda2)]
  cv_lambda1_min <- cv.gcdnet(
    x = x,
    y = y,
    lambda2 = cv_lambda2_min,
    nfolds = nfolds,
    foldid = foldid,
    method = method,
    pf = pf,
    pf2 = pf2
  )$lambda.min
  
  
  best_mod <- gcdnet(
    x = x,
    y = y,
    lambda = cv_lambda1_min,
    lambda2 = cv_lambda2_min,
    pf = pf,
    pf2 = pf2,
    method = method
  )
  
  
  if (count_selected_expos(best_mod, n_confound) < 3) {
    # Case when optimal lambda1 and lambda2 select less than 3 exposures
    
    # Find the Optimal lambda1 that selects at least 3 exposures
    cv_result <- cv.gcdnet(
      x = x,
      y = y,
      lambda2 = cv_lambda2_min,
      nfolds = nfolds,
      foldid = foldid,
      method = method,
      pf = pf,
      pf2 = pf2
    )
    
    lambda1_values <- cv_result$lambda
    optimal_lambda1 <- cv_result$lambda.min
    
    # Initialize progress bar
    pb <- progress_bar$new(
      format = "  Finding optimal lambda 1 values [:bar] :percent eta: :eta",
      total = length(lambda1_values),
      clear = FALSE,
      width = 60
    )
    
    for (lambda1 in lambda1_values) {
      model <- gcdnet(
        x = x,
        y = y,
        lambda = lambda1,
        lambda2 = cv_lambda2_min,
        pf = pf,
        pf2 = pf2,
        method = method
      )
      pb$tick()  # Update progress bar
      
      if (count_selected_expos(model, n_confound) >= 3) {
        optimal_lambda1 <- lambda1
        break
      }
    }
    
    if (is.null(optimal_lambda1)) {
      # If no combination selects at least 3 exposures, use the minimum cross-validated lambda1
      return(best_mod)
    } else{
      mod_atleast_3 <- gcdnet(
        x = x,
        y = y,
        lambda = optimal_lambda1,
        lambda2 = cv_lambda2_min,
        pf = pf,
        pf2 = pf2,
        method = method
      )
      
      # return the best model with at least three variables chosen
      return(mod_atleast_3)
    }
    
  } else{
    # Return the Elastic Net results with optimal lambda settings
    return(best_mod)
  }
}


############## Function - ers.score_adapt #####################

ers_score_adapt = function(data, coef) {
  score <- data %*% coef
  colnames(score) = 'ERS'
  return(score)
}

############## Function - ers_Calc #######################

ers_Calc = function(data,
                    exposure,
                    outcome,
                    covar = NULL,
                    lambda2_start = NULL,
                    include_int = T,
                    method = 'ls',
                    scaled = FALSE,
                    nfolds = 5,
                    seed = NULL,
                    ...) {
  # require the needed pkgs
  pkgs <- c("dplyr", "gcdnet", "magrittr", "progress", "pbapply")
  suppressMessages(lapply(pkgs, require, character.only = T))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  x <- data %>% dplyr::select(all_of(exposure)) %>% as.data.frame()
  y <- data[[outcome]] %>% as.matrix() %>% as.numeric()
  confounders <- covar
  covar <- data %>% dplyr::select(all_of(covar)) %>% as.data.frame()
  n <- length(y)
  
  if (is.null(covar) == F) {
    if (any(!complete.cases(x)) |
        any(!complete.cases(y)) |
        any(!complete.cases(covar))) {
      stop('x, y, or covar contain missing values. This method requires complete data.')
    }
  } else{
    if (any(!complete.cases(x)) | any(!complete.cases(y))) {
      stop('x, y, or covar contain missing values. This method requires complete data.')
    }
  }
  
  foldid <- matrix(data = c(sample(n), rep(1:nfolds, length = n)),
                   nrow = n,
                   ncol = 2)
  foldid <- foldid[order(foldid[, 1]), ]
  foldid <- foldid[, 2]
  
  if (include_int == T) {
    # the complex case where the interactions, square terms are in the ENET var selection
    data_mod <- model.matrix(~ -1 + .^2, data = x)
    x_sq <- x^2
    names(x_sq) <- paste0(names(x), '^2')
    
    if (is.null(covar) == F) {
      pf <- c(rep(1, ncol(data_mod) + ncol(x_sq)), rep(0, ncol(covar)))
      pf2 <- c(rep(1, ncol(data_mod) + ncol(x_sq)), rep(0, ncol(covar)))
      data_mod <- cbind(data_mod, x_sq, covar)
      
      tmp_data <- data_mod
      
      if (!isTRUE(scaled)) {
        data_mod <- as.matrix(apply(data_mod, 2, function(y) {
          scale(y, center = T, scale = T)
        }))
      }
      
    } else {
      pf <- c(rep(1, ncol(data_mod) + ncol(x_sq)))
      pf2 <- c(rep(1, ncol(data_mod) + ncol(x_sq)))
      data_mod <- cbind(data_mod, x_sq)
      
      tmp_data <- data_mod
      
      if (!isTRUE(scaled)) {
        data_mod <- as.matrix(apply(data_mod, 2, function(y) {
          scale(y, center = T, scale = T)
        }))
      }
      
    }
    
  } else {
    # the simple case with just the original exposures
    data_mod <- model.matrix(~ -1 + ., data = x)
    
    if (is.null(covar) == F) {
      pf <- c(rep(1, ncol(data_mod)), rep(0, ncol(covar)))
      pf2 <- c(rep(1, ncol(data_mod)), rep(0, ncol(covar)))
      data_mod <- cbind(data_mod, covar)
      
      tmp_data <- data_mod
      
      if (!isTRUE(scaled)) {
        data_mod <- as.matrix(apply(data_mod, 2, function(y) {
          scale(y, center = T, scale = T)
        }))
      }
      
    } else {
      pf <- c(rep(1, ncol(data_mod)))
      pf2 <- c(rep(1, ncol(data_mod)))
      data_mod <- cbind(data_mod)
      
      tmp_data <- data_mod
      
      if (!isTRUE(scaled)) {
        data_mod <- as.matrix(apply(data_mod, 2, function(y) {
          scale(y, center = T, scale = T)
        }))
      }
      
    }
  }
  
  # ordinary Elastic net
  ers_fit <- ers_enet_adapt(data_mod,
                            y,
                            lambda2_start,
                            nfolds,
                            foldid,
                            pf,
                            pf2,
                            method,
                            ncol(covar))
  
  ers_beta <- as.matrix(coef(ers_fit))
  ers_beta_keep <- ers_beta != 0
  tab <- matrix(0, sum(ers_beta_keep), 1)
  rownames(tab) <- rownames(ers_beta)[ers_beta_keep]
  tab[, 1] <- ers_beta[ers_beta_keep, ]
  
  if (is.null(covar) == F) {
    tab_exposure <- subset(tab, !(row.names(tab) %in% c(
      '(Intercept)', colnames(covar)
    )))
  } else {
    tab_exposure <- subset(tab, !(row.names(tab) %in% c('(Intercept)')))
  }
  
  
  coef_enet <- as.numeric(tab_exposure)
  dat_score <- as.matrix(data_mod[, rownames(tab_exposure)])
  
  # calculate the ERS score for each observation
  ers_scores <- ers_score_adapt(data = dat_score, coef = coef_enet)
  
  if (!is.null(covar)) {
    tmp_data_noExpo <- data %>% dplyr::select(!all_of(exposure)) %>% dplyr::select(!all_of(confounders))
  } else {
    tmp_data_noExpo <- data %>% dplyr::select(!all_of(exposure))
  }
  
  tmp_data <- cbind(tmp_data_noExpo, tmp_data, ers_scores)
  
  ers_obj <- list(
    post_ERS_data = tmp_data,
    ers_scores = ers_scores,
    # constructed ERS score
    ers_fit = ers_fit,
    # ENET result of Y ~ expanded X
    coef = coef_enet,
    # Coefficients of ENET
    dat_score = dat_score # exposures with non-zero effects after ENET
  )
  class(ers_obj) <- 'ers'
  
  return(ers_obj)
}
