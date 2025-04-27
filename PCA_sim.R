# This is the code script to run the PCA simulation

pkgs <- c(
  "tidyverse",
  "MASS",
  "CMAverse",
  "corrplot",
  "cluster",
  "factoextra",
  "ggbiplot",
  "ggfortify"
)

invisible(lapply(pkgs, library, character.only = T))

# Read in the data
list_df <- read_rds("RDS/list_df.rds")
list_medTest <- read_rds("RDS/list_medTest.rds")

set.seed(1211)

source("Functions/Functions_PCA.R")

######################################### Generate PCA components  ############################

# Perform PCA on the generated X data
PCA_res <- prcomp(list_df[["PCA"]] %>% dplyr::select(contains("X")),
                  center = TRUE,
                  scale. = TRUE)

# Check the proportion of variance explained by each PC
explained_variance <- summary(PCA_res)$importance[3, ] # the cumulative explained variance
print(explained_variance) # The 15th PC  cumulatively explain over 80% variance in x

id_80var <- which.max(explained_variance >= 0.8)

# combine with the PCs that cumulatively explain 80% of the variance
list_df[["PCA"]] <- cbind(list_df[["PCA"]], PCA_res$x[, 1:id_80var])

######################## Scree Plots Loadings Plots, Loading Heatmaps  ########################

## Scree plots
plot <- fviz_eig(PCA_res,
                 addlabels = T,
                 choice = "variance",
                 ncp = 16) +
  theme(text = element_text(size = 14)) +
  labs(x = "Principal Components", title = "")
plot
#
# pdf("../pic/pca_sim_scree.pdf", width = 14, height = 10)
# print(plot)
# dev.off()
#
# tiff("../pic/pca_sim_scree.tiff", width = 14, height = 10, units = "in", res = 800, compression = "lzw")
# print(plot)
# dev.off()


## Create a heatmap of the loadings
# Convert loadings matrix to data frame
loadings_PCA <- as.data.frame(PCA_res$rotation[, 1:id_80var]) %>%
  rownames_to_column(var = "Exposures") %>%
  gather(Principal_Component, Loading, -Exposures)
loadings_PCA$Exposures <- gsub("x", "", loadings_PCA$Exposures) %>% as.numeric()
loadings_PCA$Exposures <- factor(loadings_PCA$Exposures, levels = rev(sort(unique(
  loadings_PCA$Exposures
))))

# Draw heat map
plot <- ggplot(loadings_PCA,
               aes(y = Exposures, x = Principal_Component, fill = Loading)) +
  geom_tile(color = "grey") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  scale_y_discrete(
    breaks = seq_along(unique(loadings_PCA$Exposures)),
    labels = unique(loadings_PCA$Exposures),
    expand = c(0, 0)
  ) +
  scale_x_discrete(expand = c(0, 0), limits = paste0("PC", 1:15)) +
  labs(y = "Exposures", x = "Principal Component") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, hjust = 1),
    text = element_text(size = 12),
    #remove plot background
    plot.background = element_blank(),
    #remove plot border
    panel.border = element_blank(),
    # remove the axis ticks
    axis.ticks = element_blank()
  )
plot

# ggsave(
#   "pca_sim_heat.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
#
# tiff("../pic/pca_sim_heat.tiff", width = 7, height = 5, units = "in", res = 800, compression = "lzw")
# print(plot)
# dev.off()



################################ Mediation Effect Testing and estimation  ###############

# In this case, since there are 15 PC components to test
# we will use a single PC component as the exposures and test for mediation effect one at a time
# the other PC components that are not considered as the exposure is considered as confounders

#setup
confounders_nm <- list_df[["PCA"]] %>% dplyr::select(starts_with("c")) %>% colnames()

mediator_nm <- list_df[["PCA"]] %>% dplyr::select(starts_with("m")) %>% colnames()

outcome_nm <- "y"

## Mediation Test on all 15 PCs
set.seed(1211)
tictoc::tic("Mediation Testing with PCA")
list_medTest[["PCA"]] <- PCA_medTest(
  data = list_df[["PCA"]],
  nboot = 1000,
  outcome_nm = outcome_nm,
  mediator_nm = mediator_nm,
  confounders_nm = confounders_nm
)
tictoc::toc()
# 5 mins


# list_medTest[["PCA"]][["PC1"]][["CMA Summary Table"]]

# extract the DE IE and TE of each PC component
tmp <- paste(c("CDE", "PNDE", "TNDE", "PNIE", "TNIE", "TE", "PM"), "Table")
tidy_res <- list()

for (new_df in tmp) {
  trgt <- str_to_lower(str_replace_all(new_df, " Table", ""))
  
  tidy_res[[new_df]] <- purrr::map_df(names(list_medTest[["PCA"]]), function(pc) {
    df <- list_medTest[["PCA"]][[pc]][["CMA Summary Table"]] %>% as.data.frame()
    trgt_row <- df[trgt, , drop = FALSE]  # Extract the target effect row
    trgt_row <- cbind(PC = pc, trgt_row)  # Add a column with the PC name
    trgt_row
  }, .id = NULL)
  # rename each row
  rownames(tidy_res[[new_df]]) <- tidy_res[[new_df]]$PC
  tidy_res[[new_df]] <- tidy_res[[new_df]] %>% dplyr::select(-PC)
}

list_medTest[["PCA"]][["Tidy Results"]] <- tidy_res

## Bonferroni correction for all the effects p-value
correct_trgt <- list_medTest[["PCA"]][["Tidy Results"]] %>% names

for (correct_id in correct_trgt) {
  list_medTest[["PCA"]][["Tidy Results"]][[correct_id]][["BH.Adj.Pval"]] <- list_medTest[["PCA"]][["Tidy Results"]][[correct_id]][["Pval"]] %>%
    p.adjust(method = "BH")
  
}

list_medTest[["PCA"]][["Tidy Results"]]$`PNIE Table` %>% round(2)

########################### Output for future use #########################

write_rds(list_df, "RDS/list_df.rds")
write_rds(list_medTest, "RDS/list_medTest.rds")
