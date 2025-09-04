# This is the script to conduct PCA mediation analysis using the CMAverse package

require("CMAverse")


########################### Function - PCA_Generate ##########################

## Input
# 1. data -- Raw data with outcome (Y), mediator (M), exposures (Xs). The exposures (variables starting with X) would be PCAed
#
# 2. prop_explained -- the proportion of variance to be explained by the PC components. Default is 80%

PCA_Generate <- function(data, prop_explained = 0.8) {
  if (prop_explained <= 0 | prop_explained > 1) {
    stop("The desired proportion of variance to be explained should be between 0 and 1.")
  }
  
  PC_X_res <- prcomp(data %>% dplyr::select(starts_with("X")))
  # PC_X_res$x # the score matrix of the PCA
  
  # Extract the PC component that cumulatively explains 80% of the variance
  # (var explained/total variance) = proportion explained
  prop_exp_variance <- (PC_X_res$sdev^2) / sum(PC_X_res$sdev^2)
  
  # find the PC component that cumulatively explain 80% variance in Exposures
  id_80prop <- which(cumsum(prop_exp_variance) >= prop_explained)[1]
  
  # Extract and Merge the PC scores to the simulated data
  tmp <- data.frame(PC_X_res$x[, 1:id_80prop])
  colnames(tmp) <- paste0("PCA_", 1:id_80prop)
  
  data <- cbind(data, tmp)
  
  list_output <- list(
    loading_mat = PC_X_res$rotation[, 1:id_80prop],
    # the rotation matrix
    PC_num = id_80prop,
    # the number of PC components to be included in the new data
    post_PCA_data = data
  ) # the new data
  
  return(list_output)
}

###################### Function - PCA_medTest ##############################

PCA_medTest <- function(data,
                        confounders_nm,
                        mediator_nm,
                        outcome_nm) {
  suppressMessages(invisible(lapply(
    c("CMAverse", "purrr", "rlang"), require, character.only = T
  )))
  
  # extract the list of PCA components (names starting with PC)
  PCA_names <- names(data)[grepl("^PC", names(data))]
  
  if (length(PCA_names) == 0) {
    stop('There are no columns with names starting with "PC"')
  }
  
  # set a condition fro the outcome names not in the data
  if (outcome_nm %in% colnames(data)) {
    outcome_nm <- outcome_nm
  } else {
    outcome_nm <- "Y"
  }
  
  # a list to store the results
  medTest_res <- list()
  
  # a function for purrr::map()
  run_cmest <- function(PCA_id,
                        data,
                        confounders_nm = confounders_nm,
                        mediator_nm = mediator_nm,
                        outcome_nm = outcome_nm) {
    if (length(confounders_nm) == 0) {
      # the CMA mediation effect testing
      suppressWarnings(
        cma_test <- CMAverse::cmest(
          data = data,
          model = "rb",
          full = T,
          EMint = F,
          yreg = "linear",
          mreg = list("linear"),
          mval = list(1),
          basec = c(PCA_names[PCA_names != PCA_id]),
          outcome = outcome_nm,
          exposure = PCA_id,
          mediator = mediator_nm,
          estimation = "paramfunc",
          inference = "delta"
        )
      )
    } else {
      # the CMA mediation effect testing
      suppressWarnings(
        cma_test <- CMAverse::cmest(
          data = data,
          model = "rb",
          full = T,
          EMint = F,
          yreg = "linear",
          mreg = list("linear"),
          mval = list(1),
          basec = c(PCA_names[PCA_names != PCA_id], confounders_nm),
          outcome = outcome_nm,
          exposure = PCA_id,
          mediator = mediator_nm,
          estimation = "paramfunc",
          inference = "delta"
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
    
    print(paste("     ", PCA_id, "Calculation done"))
    
    list(`CMA Test` = cma_test, `CMA Summary Table` = summary_table)
  }
  
  # run the mediation tests over all the PC components
  medTest_res <- purrr::map(PCA_names,
                            ~ run_cmest(.x, data, confounders_nm, mediator_nm, outcome_nm))
  names(medTest_res) <- PCA_names
  return(medTest_res)
}
