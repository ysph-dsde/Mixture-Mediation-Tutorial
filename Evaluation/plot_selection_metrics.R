# Load required packages
library(ggplot2)
library(ggh4x)
library(tidyverse)
library(scales)

# Load simulation results
sema <- readRDS("sema.rds")
bkmr <- readRDS("bkmr.rds")

# Format SEMA TPR/FPR
sema_tpr <- sema$TPR
sema_fpr <- sema$FPR
sema_tpr <- sema_tpr |>
  rename(Mean = TPR_mean,
         Lower = TPR_lower,
         Upper = TPR_upper) |>
  mutate(Metric = "TPR", Method = paste0("SE-MA (", Method, ")"))
sema_fpr <- sema_fpr |>
  rename(Mean = FPR_mean,
         Lower = FPR_lower,
         Upper = FPR_upper) |>
  mutate(Metric = "FPR", Method = paste0("SE-MA (", Method, ")"))
sema_combined <- rbind(sema_tpr, sema_fpr)

# Format BKMR TPR/FPR
cutoff_keys <- names(bkmr$bkmr)
method_keys <- c("bkmr", "bkmr_hier")
metric_keys <- c("TPR", "FPR")

make_bkmr_df <- function(metric) {
  do.call(rbind,
          lapply(cutoff_keys, function(cutoff_key) {
            lapply(method_keys, function(method) {
              mat <- bkmr$bkmr[[cutoff_key]][[method]][[metric]]
              if (is.null(mat))
                return(NULL)
              do.call(rbind, lapply(1:ncol(mat), function(j) {
                vec <- mat[, j]
                q <- quantile(vec, c(0.025, 0.975), na.rm = TRUE)
                data.frame(
                  Method = paste0(
                    ifelse(method == "bkmr", "BKMR (C)", "BKMR (H)"),
                    " (",
                    gsub("cutoff_", "", cutoff_key),
                    ")"
                  ),
                  SampleSize = ifelse(j <= 2, "n = 1000", "n = 2500"),
                  Scenario = paste0(ifelse(j %% 2 == 1, "Weak", "Strong"), " indirect effect"),
                  Mean = mean(vec, na.rm = TRUE),
                  Lower = q[1],
                  Upper = q[2],
                  Metric = metric
                )
              }))
            })
          }) |> unlist(recursive = FALSE))
}

bkmr_tpr <- make_bkmr_df("TPR")
bkmr_fpr <- make_bkmr_df("FPR")
bkmr_combined <- rbind(bkmr_tpr, bkmr_fpr)

# Combine all
all_df <- rbind(
  data.frame(
    sema_combined,
    Mean = sema_combined$Mean,
    Lower = sema_combined$Lower,
    Upper = sema_combined$Upper
  )[, c("Method",
        "SampleSize",
        "Scenario",
        "Mean",
        "Lower",
        "Upper",
        "Metric")],
  bkmr_combined
)

# Plot function
plot_metric <- function(metric_label) {
  ggplot(subset(all_df, Metric == metric_label),
         aes(x = Method, y = Mean, fill = Method)) +
    geom_bar(stat = "identity",
             position = position_dodge(0.8),
             width = 0.7) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  position = position_dodge(0.8),
                  width = 0.2) +
    ggh4x::facet_grid2(Scenario ~ SampleSize) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      x = NULL,
      y = ifelse(
        metric_label == "TPR",
        "True Positive Rate",
        "False Positive Rate"
      ),
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      strip.background = element_rect(fill = "white", color = "black"),
      panel.border = element_rect(color = "black", fill = NA)
    )
}

# Draw TPR and FPR plots
plot_metric("TPR")
plot_metric("FPR")

