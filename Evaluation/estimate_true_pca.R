# This is the code script to run the PCA simulation

library(tidyverse)
library(MASS)
library(CMAverse)
library(cluster)
library(factoextra)
source("../Functions/Functions_PCA.R")

set.seed(1211)
sce <- c("l", "h")
for(level in sce){
    # Read in the data
    list_df <- readRDS(paste0("sim_", level, ".rds"))
    list_medTest <- list()

    #### Generate PCA components  ######

    ## Perform PCA on the generated X data
    PCA_res <- prcomp(list_df[["Raw"]] %>% dplyr::select(contains("X")),
                    center = TRUE,
                    scale. = TRUE)

    # Check the proportion of variance explained by each PC
    explained_variance <-
    summary(PCA_res)$importance[3,] # the cumulative explained variance
    id_80var <- which.max(explained_variance >= 0.8)

    # combine with the PCs that cumulatively explain 80% of the variance
    list_df[["PCA"]] <- cbind(list_df[["Raw"]], PCA_res$x[, 1:id_80var])

    #### Mediation Effect Testing and estimation  ####

    #setup
    confounders_nm <-
    list_df[["PCA"]] %>% dplyr::select(starts_with("c")) %>% colnames()

    mediator_nm <-
    list_df[["PCA"]] %>% dplyr::select(starts_with("m")) %>% colnames()

    outcome_nm <- "y"

    ## Mediation Test on all selected PCs
    set.seed(1211)
    list_medTest[["PCA"]] <- PCA_medTest(
    data = list_df[["PCA"]],
    outcome_nm = outcome_nm,
    mediator_nm = mediator_nm,
    confounders_nm = confounders_nm
    )

    # extract the DE IE and TE of each PC component
    tmp <-
    paste(c("CDE", "PNDE", "TNDE", "PNIE", "TNIE", "TE", "PM"), "Table")
    tidy_res <- list()

    for (new_df in tmp) {
    trgt <- str_to_lower(str_replace_all(new_df, " Table", ""))
    
    tidy_res[[new_df]] <-
        purrr::map_df(names(list_medTest[["PCA"]]), function(pc) {
        df <-
            list_medTest[["PCA"]][[pc]][["CMA Summary Table"]] %>% as.data.frame()
        trgt_row <-
            df[trgt, , drop = FALSE]  # Extract the target effect row
        trgt_row <-
            cbind(PC = pc, trgt_row)  # Add a column with the PC name
        trgt_row
        }, .id = NULL)
    # rename each row
    rownames(tidy_res[[new_df]]) <- tidy_res[[new_df]]$PC
    tidy_res[[new_df]] <- tidy_res[[new_df]] %>% dplyr::select(-PC)
    }

    list_medTest <- tidy_res

    ## Bonferroni correction for all the effects p-value
    correct_trgt <- list_medTest %>% names

    for (correct_id in correct_trgt) {
    list_medTest[[correct_id]][["BH.Adj.Pval"]] <-
        list_medTest[[correct_id]][["Pval"]] %>%
        p.adjust(method = "BH")
    }

    #### Output Results ########
    saveRDS(list_df, paste0("sim_", level, ".rds"))
    saveRDS(list_medTest, paste0("pca_", level, ".rds"))
}
