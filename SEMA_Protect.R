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

# Read in the data
list_headZscore <- read_rds("RDS/list_headZscore.rds")

set.seed(1211)

source("Functions/Functions_IndTesting.R")

#################### LTE4 and HEADCIRCUMFERENCEZSCORE (PHTH only) - SEMA WITHOUT co-exposures ###################

# set up
expo_nm <- colnames(list_headZscore[["Exposures"]])
confounders_nm <- colnames(list_headZscore[["Confounders"]])[!grepl("Intercept", colnames(list_headZscore[["Confounders"]]))]
mediator_nm <- "LTE4"
outcome_nm <- colnames(list_headZscore[["Data"]])[1]

tmp_df <- list_headZscore[["Data"]]

# run
set.seed(1211)
tictoc::tic("Individual Exposure Mediation Testing (No coexposures)")
list_headZscore[["medTest"]][["IndExpo (Naive)"]] <- IndExpo_medTest_NV(
  data = tmp_df,
  exposures_nm = expo_nm,
  mediator_nm = mediator_nm,
  outcome_nm = outcome_nm,
  confounders_nm = confounders_nm,
  nboot = 1000
)
tictoc::toc()
# 1.4 mins


# extract the DE IE and TE of each PC component
tmp <- paste(c("CDE", "PNDE", "TNDE", "PNIE", "TNIE", "TE", "PM"), "Table")
tidy_res <- list()

for (new_df in tmp) {
  trgt <- str_to_lower(str_replace_all(new_df, " Table", ""))
  
  tidy_res[[new_df]] <- purrr::map_df(names(list_headZscore[["medTest"]][["IndExpo (Naive)"]]), function(expo) {
    df <- list_headZscore[["medTest"]][["IndExpo (Naive)"]][[expo]][["CMA Summary Table"]] %>% as.data.frame()
    trgt_row <- df[trgt, , drop = FALSE]  # Extract the target effect row
    trgt_row <- cbind(Expo = expo, trgt_row)  # Add a column with the PC name
    trgt_row
  }, .id = NULL)
  # rename each row
  rownames(tidy_res[[new_df]]) <- tidy_res[[new_df]]$Expo
  tidy_res[[new_df]] <- tidy_res[[new_df]] %>% dplyr::select(-Expo)
}

list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]] <- tidy_res

## BH correction for all the effects p-value
correct_trgt <- list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]] %>% names

for (correct_id in correct_trgt) {
  list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]][[correct_id]][["BH.Adj.Pval"]] <- list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]][[correct_id]][["Pval"]] %>%
    p.adjust(method = "BH")
}

list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]][["TNDE Table"]] %>% round(2)
list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]][["PNIE Table"]] %>% round(2)
list_headZscore[["medTest"]][["IndExpo (Naive)"]][["Tidy Results"]][["TE Table"]] %>% round(2)

write_rds(list_headZscore, "RDS/list_headZscore.rds")

#################### LTE4 and HEADCIRCUMFERENCEZSCORE (PHTH only) - Individual Exposures WITH co-exposures ###################

# set up
expo_nm <- colnames(list_headZscore[["Exposures"]])
confounders_nm <- colnames(list_headZscore[["Confounders"]])[!grepl("Intercept", colnames(list_headZscore[["Confounders"]]))]
mediator_nm <- "LTE4"
outcome_nm <- colnames(list_headZscore[["Data"]])[1]

tmp_df <- list_headZscore[["Data"]]

# run
set.seed(1211)
tictoc::tic("Individual Exposure Mediation Testing (W. coexposures)")
list_headZscore[["medTest"]][["IndExpo"]] <- IndExpo_medTest(
  data = tmp_df,
  exposures_nm = expo_nm,
  mediator_nm = mediator_nm,
  outcome_nm = outcome_nm,
  confounders_nm = confounders_nm,
  nboot = 1000
)
tictoc::toc()
# 1.8 mins


# extract the DE IE and TE of each PC component
tmp <- paste(c("CDE", "PNDE", "TNDE", "PNIE", "TNIE", "TE", "PM"), "Table")
tidy_res <- list()

for (new_df in tmp) {
  trgt <- str_to_lower(str_replace_all(new_df, " Table", ""))
  
  tidy_res[[new_df]] <- purrr::map_df(names(list_headZscore[["medTest"]][["IndExpo"]]), function(expo) {
    df <- list_headZscore[["medTest"]][["IndExpo"]][[expo]][["CMA Summary Table"]] %>% as.data.frame()
    trgt_row <- df[trgt, , drop = FALSE]  # Extract the target effect row
    trgt_row <- cbind(Expo = expo, trgt_row)  # Add a column with the PC name
    trgt_row
  }, .id = NULL)
  # rename each row
  rownames(tidy_res[[new_df]]) <- tidy_res[[new_df]]$Expo
  tidy_res[[new_df]] <- tidy_res[[new_df]] %>% dplyr::select(-Expo)
}

list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]] <- tidy_res

## BH correction for all the effects p-value
correct_trgt <- list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]] %>% names

for (correct_id in correct_trgt) {
  list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]][[correct_id]][["BH.Adj.Pval"]] <- list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]][[correct_id]][["Pval"]] %>%
    p.adjust(method = "BH")
}

list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]][["TNDE Table"]] %>% round(2)
list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]][["PNIE Table"]] %>% round(2)
list_headZscore[["medTest"]][["IndExpo"]][["Tidy Results"]][["TE Table"]] %>% round(2)

write_rds(list_headZscore, "RDS/list_headZscore.rds")
