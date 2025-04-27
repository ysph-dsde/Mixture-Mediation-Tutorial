# This is the code script to save the useful functions to
# test the mediation effects and conduct causal mediation analysis with the individual exposures.

############################ Function - IndExpo_medTest_NV ################################

## Input
# 1. data -- Raw data with outcome (Y), mediator (M), exposures (Xs), and the confounders (C). The individual exposures are each tested and estimated of mediation effects
#
# 2. nboot -- number of bootstrap simulations to perform for each mediation test (a CMAverse::cmest() argument). The default is 1000 simulations

## The Naive approach means to run the individual testing with one exposure at a time,
## without putting other exposures as the confounders

IndExpo_medTest_NV <- function(data,
                               nboot = 1000,
                               exposures_nm,
                               confounders_nm,
                               mediator_nm,
                               outcome_nm) {
  suppressMessages(invisible(lapply(
    c("CMAverse", "purrr", "rlang"), require, character.only = T
  )))
  
  # extract the list of exposures
  Expo_names <- exposures_nm
  
  if (length(Expo_names) == 0) {
    stop('There are no columns with names in Expo_names')
  }
  
  # a list to store the results
  medTest_res <- list()
  
  # a function for purrr::map()
  run_cmest <- function(Expo_id,
                        data,
                        nboot,
                        confounders_nm = confounders_nm,
                        mediator_nm = mediator_nm,
                        outcome_nm = outcome_nm) {
    if (length(confounders_nm) == 0) {
      y_formula <- as.formula(paste0(outcome_nm, " ~ ", Expo_id, " + ", mediator_nm))
      m_formula <- as.formula(paste0(mediator_nm, " ~ ", Expo_id))
      
      # Evaluate the formulas in the correct environment
      y_model <- eval(bquote(glm(
        .(y_formula), family = gaussian, data = data
      )))
      m_model <- eval(bquote(glm(
        .(m_formula), family = gaussian, data = data
      )))
      
      # the CMA mediation effect testing
      suppressWarnings(
        cma_test <- CMAverse::cmest(
          data = data,
          model = "rb",
          full = T,
          EMint = F,
          yreg = y_model,
          mreg = list(m_model),
          mval = list(0),
          outcome = outcome_nm,
          exposure = Expo_id,
          mediator = mediator_nm,
          inference = "bootstrap",
          nboot = nboot,
          boot.ci.type = "per"
        )
      )
    } else {
      y_formula <- as.formula(paste0(
        outcome_nm,
        " ~ ",
        Expo_id,
        " + ",
        mediator_nm,
        " + ",
        paste0(confounders_nm, collapse = " + ")
      ))
      m_formula <- as.formula(paste0(
        mediator_nm,
        " ~ ",
        Expo_id,
        " + ",
        paste0(confounders_nm, collapse = " + ")
      ))
      
      # Evaluate the formulas in the correct environment
      y_model <- eval(bquote(glm(
        .(y_formula), family = gaussian, data = data
      )))
      m_model <- eval(bquote(glm(
        .(m_formula), family = gaussian, data = data
      )))
      
      # the CMA mediation effect testing
      suppressWarnings(
        cma_test <- CMAverse::cmest(
          data = data,
          model = "rb",
          full = T,
          EMint = F,
          yreg = y_model,
          mreg = list(m_model),
          mval = list(0),
          basec = c(confounders_nm),
          outcome = outcome_nm,
          exposure = Expo_id,
          mediator = mediator_nm,
          inference = "bootstrap",
          nboot = nboot,
          boot.ci.type = "per"
        )
      )
    }
    
    # tidy a summary table
    summary_table <- cbind(
      cma_test$effect.pe,
      cma_test$effect.se,
      cma_test$effect.ci.low,
      cma_test$effect.ci.high,
      cma_test$effect.pval
    )
    
    colnames(summary_table) <- c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")
    
    print(paste("     ", Expo_id, "Bootstrap done"))
    
    list(`CMA Test` = cma_test, `CMA Summary Table` = summary_table)
  }
  
  # run the mediation tests over all the exposures
  medTest_res <- purrr::map(Expo_names,
                            ~ run_cmest(.x, data, nboot, confounders_nm, mediator_nm, outcome_nm))
  names(medTest_res) <- Expo_names
  return(medTest_res)
}


####################### Function - IndExpo_medTest ##########################

## This is the modified version of IndExpo_medTest_NV to consider the coexposures as confounders

IndExpo_medTest <- function(data,
                            nboot = 1000,
                            exposures_nm,
                            confounders_nm,
                            mediator_nm,
                            outcome_nm) {
  suppressMessages(invisible(lapply(
    c("CMAverse", "purrr", "rlang"), require, character.only = T
  )))
  
  # extract the list of exposures
  Expo_names <- exposures_nm
  
  if (length(Expo_names) == 0) {
    stop('There are no columns with names in Expo_names')
  }
  
  # a list to store the results
  medTest_res <- list()
  
  # a function for purrr::map()
  run_cmest <- function(Expo_id,
                        data,
                        nboot,
                        confounders_nm = confounders_nm,
                        mediator_nm = mediator_nm,
                        outcome_nm = outcome_nm) {
    if (length(confounders_nm) == 0) {
      y_formula <- as.formula(paste0(
        outcome_nm,
        " ~ ",
        Expo_id,
        " + ",
        mediator_nm,
        " + ",
        paste0(Expo_names[Expo_names != Expo_id], collapse = " + ")
      ))
      
      m_formula <- as.formula(paste0(
        mediator_nm,
        " ~ ",
        Expo_id,
        " + ",
        paste0(Expo_names[Expo_names != Expo_id], collapse = " + ")
      ))
      
      # Evaluate the formulas in the correct environment
      y_model <- eval(bquote(glm(
        .(y_formula), family = gaussian, data = data
      )))
      m_model <- eval(bquote(glm(
        .(m_formula), family = gaussian, data = data
      )))
      
      # the CMA mediation effect testing
      suppressWarnings(
        cma_test <- CMAverse::cmest(
          data = data,
          model = "rb",
          full = T,
          EMint = F,
          yreg = y_model,
          mreg = list(m_model),
          mval = list(0),
          basec = c(Expo_names[Expo_names != Expo_id]),
          outcome = outcome_nm,
          exposure = Expo_id,
          mediator = mediator_nm,
          inference = "bootstrap",
          nboot = nboot,
          boot.ci.type = "per"
        )
      )
    } else {
      y_formula <- as.formula(paste0(
        outcome_nm,
        " ~ ",
        Expo_id,
        " + ",
        mediator_nm,
        " + ",
        paste0(Expo_names[Expo_names != Expo_id], collapse = " + "),
        " + ",
        paste0(confounders_nm, collapse = " + ")
      ))
      
      m_formula <- as.formula(paste0(
        mediator_nm,
        " ~ ",
        Expo_id,
        " + ",
        paste0(Expo_names[Expo_names != Expo_id], collapse = " + "),
        " + ",
        paste0(confounders_nm, collapse = " + ")
      ))
      
      # Evaluate the formulas in the correct environment
      y_model <- eval(bquote(glm(
        .(y_formula), family = gaussian, data = data
      )))
      m_model <- eval(bquote(glm(
        .(m_formula), family = gaussian, data = data
      )))
      
      # the CMA mediation effect testing
      suppressWarnings(
        cma_test <- CMAverse::cmest(
          data = data,
          model = "rb",
          full = T,
          EMint = F,
          yreg = y_model,
          mreg = list(m_model),
          mval = list(0),
          basec = c(Expo_names[Expo_names != Expo_id], confounders_nm),
          outcome = outcome_nm,
          exposure = Expo_id,
          mediator = mediator_nm,
          inference = "bootstrap",
          nboot = nboot,
          boot.ci.type = "per"
        )
      )
    }
    
    # tidy a summary table
    summary_table <- cbind(
      cma_test$effect.pe,
      cma_test$effect.se,
      cma_test$effect.ci.low,
      cma_test$effect.ci.high,
      cma_test$effect.pval
    )
    
    colnames(summary_table) <- c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")
    
    print(paste("     ", Expo_id, "Bootstrap done"))
    
    list(`CMA Test` = cma_test, `CMA Summary Table` = summary_table)
  }
  
  # run the mediation tests over all the exposures
  medTest_res <- purrr::map(Expo_names,
                            ~ run_cmest(.x, data, nboot, confounders_nm, mediator_nm, outcome_nm))
  names(medTest_res) <- Expo_names
  return(medTest_res)
}
