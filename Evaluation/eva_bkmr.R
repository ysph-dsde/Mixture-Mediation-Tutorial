# Load required packages
library(ggplot2)
library(ggh4x)
library(tidyverse)
data <- readRDS("bkmr.rds")

set.seed(1211)

# Global settings
n_dataset <- 100  # Number of simulated datasets per scenario
n_scenario <- 4   # Number of scenarios
sce <- c("1000l", "1000h", "2500l", "2500h")  # Scenario labels
true_signals <-
  c(1, 2, 3, 11, 12, 13, 21, 22, 23)  # True mediating exposures

# Transformation
cutoff_keys <- names(data$bkmr)
method_keys <- c("bkmr", "bkmr_hier")

make_df_from_nested <- function(data, metric = "TPR") {
  do.call(rbind, lapply(cutoff_keys, function(cutoff_key) {
    lapply(method_keys, function(method) {
      mat <- data$bkmr[[cutoff_key]][[method]][[metric]]
      if (is.null(mat)) return(NULL)
      
      df <- do.call(rbind, lapply(1:ncol(mat), function(j) {
        vec <- mat[, j]
        q <- quantile(vec, c(0.025, 0.975), na.rm = TRUE)
        
        # Mapping logic
        sample_size <- ifelse(j <= 2, "1000", "2500")
        effect <- ifelse(j %% 2 == 1, "Weak", "Strong")
        
        data.frame(
          Method = ifelse(method == "bkmr", "BKMR", "BKMR (Hierarchical)"),
          Cutoff = gsub("cutoff_", "", cutoff_key),
          SampleSize = sample_size,
          Effect = effect,
          Mean = mean(vec, na.rm = TRUE),
          Lower = q[1],
          Upper = q[2]
        )
      }))
      df
    })
  }) |> unlist(recursive = FALSE))
}

df_tpr <- make_df_from_nested(data, "TPR")
df_fpr <- make_df_from_nested(data, "FPR")

## TPR
ggplot(df_tpr, aes(x = Cutoff, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(Effect ~ SampleSize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "PIP Threshold", y = "True Positive Rate", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

## FPR
ggplot(df_fpr, aes(x = Cutoff, y = Mean, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper),
                position = position_dodge(0.8), width = 0.2) +
  facet_grid(Effect ~ SampleSize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "PIP Threshold", y = "False Positive Rate", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )
