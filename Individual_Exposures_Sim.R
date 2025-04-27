# This is the R script to analysis the individual exposures directly and do mediation analysis.
# This analysis should serve as a baseline performance on the identification of the (individual) mediation effects

######################## Set up #############################

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

# the simulated data for PCA case 1
list_df <- read_rds("RDS/list_df.rds")

set.seed(1211)

source("Functions/Functions_IndTesting.R")

#################### Run Individual Exposure Mediation Testing WITHOUT co-exposures ###################

# a list to store the mediation testing results
list_medTest <- list()


# set up
expo_nm <- list_df[["Raw"]] %>% dplyr::select(starts_with("x")) %>% colnames()
confounders_nm <- list_df[["Raw"]] %>% dplyr::select(starts_with("c")) %>% colnames()
mediator_nm <- list_df[["Raw"]] %>% dplyr::select(starts_with("m")) %>% colnames()
outcome_nm <- "y"

# run the mediation testing
set.seed(1211)
tictoc::tic("Individual exposure testing (No coexposures)")
list_medTest[["Raw (Naive)"]] <- IndExpo_medTest_NV(
  data = list_df[["Raw"]],
  nboot = 1000,
  exposures_nm = expo_nm,
  mediator_nm = mediator_nm,
  outcome_nm = outcome_nm,
  confounders_nm = confounders_nm
)
tictoc::toc()
# 6 mins

# extract the DE IE and TE of each PC component
# turn them into big dfs and merge to list_medTest
tmp <- paste(c("CDE", "PNDE", "TNDE", "PNIE", "TNIE", "TE", "PM"), "Table")
tidy_res <- list()

for (new_df in tmp) {
  trgt <- str_to_lower(str_replace_all(new_df, " Table", ""))
  
  tidy_res[[new_df]] <- purrr::map_df(names(list_medTest[["Raw (Naive)"]]), function(expo) {
    df <- list_medTest[["Raw (Naive)"]][[expo]][["CMA Summary Table"]] %>% as.data.frame()
    trgt_row <- df[trgt, , drop = FALSE]  # Extract the target effect row
    trgt_row <- cbind(Expo = expo, trgt_row)  # Add a column with the PC name
    trgt_row
  }, .id = NULL)
  # rename each row
  rownames(tidy_res[[new_df]]) <- tidy_res[[new_df]]$Expo
  tidy_res[[new_df]] <- tidy_res[[new_df]] %>% dplyr::select(-Expo)
}

list_medTest[["Raw (Naive)"]][["Tidy Results"]] <- tidy_res

## BH correction for all the effects p-value
correct_trgt <- list_medTest[["Raw (Naive)"]][["Tidy Results"]] %>% names

for (correct_id in correct_trgt) {
  list_medTest[["Raw (Naive)"]][["Tidy Results"]][[correct_id]][["BH.Adj.Pval"]] <- list_medTest[["Raw (Naive)"]][["Tidy Results"]][[correct_id]][["Pval"]] %>%
    p.adjust(method = "BH")
}

#################### Run Individual Exposure Mediation Testing WITH co-exposures ###################
set.seed(1211)
tictoc::tic("Individual exposure testing (W. Coexposure)")
list_medTest[["Raw"]] <- IndExpo_medTest(
  data = list_df[["Raw"]],
  nboot = 1000,
  exposures_nm = expo_nm,
  mediator_nm = mediator_nm,
  outcome_nm = outcome_nm,
  confounders_nm = confounders_nm
)
tictoc::toc()
# 15 mins

# extract the DE IE and TE of each PC component
# turn them into big dfs and merge to list_medTest[["Raw"]][["IE"]]
tmp <- paste(c("CDE", "PNDE", "TNDE", "PNIE", "TNIE", "TE", "PM"), "Table")
tidy_res <- list()

for (new_df in tmp) {
  trgt <- str_to_lower(str_replace_all(new_df, " Table", ""))
  
  tidy_res[[new_df]] <- purrr::map_df(names(list_medTest[["Raw"]]), function(expo) {
    df <- list_medTest[["Raw"]][[expo]][["CMA Summary Table"]] %>% as.data.frame()
    trgt_row <- df[trgt, , drop = FALSE]  # Extract the target effect row
    trgt_row <- cbind(Expo = expo, trgt_row)  # Add a column with the PC name
    trgt_row
  }, .id = NULL)
  # rename each row
  rownames(tidy_res[[new_df]]) <- tidy_res[[new_df]]$Expo
  tidy_res[[new_df]] <- tidy_res[[new_df]] %>% dplyr::select(-Expo)
}

list_medTest[["Raw"]][["Tidy Results"]] <- tidy_res

## BH correction for all the effects p-value
correct_trgt <- list_medTest[["Raw"]][["Tidy Results"]] %>% names

for (correct_id in correct_trgt) {
  list_medTest[["Raw"]][["Tidy Results"]][[correct_id]][["BH.Adj.Pval"]] <- list_medTest[["Raw"]][["Tidy Results"]][[correct_id]][["Pval"]] %>%
    p.adjust(method = "BH")
}

# list_medTest[["Raw"]][["Tidy Results"]]$`PNIE Table` %>% round(2)

############################ Output Results ###########################

write_rds(list_medTest, "RDS/list_medTest.rds")
write_rds(list_df, "RDS/list_df.rds")
