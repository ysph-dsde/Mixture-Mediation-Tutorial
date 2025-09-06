# Load required packages
library(tidyverse)
library(ggplot2)
library(scales)
library(ggh4x)

set.seed(1211)

# Global settings
n_dataset <- 100  # Number of simulated datasets per scenario
n_scenario <- 4   # Number of scenarios

# Label mapping
labels <- list(
  "1000l" = c("n = 1000", "Weak"),
  "1000h" = c("n = 1000", "Strong"),
  "2500l" = c("n = 2500", "Weak"),
  "2500h" = c("n = 2500", "Strong")
)

# Function to transfer matrices to data frame
build_bias_df <- function(mat, method_label) {
  do.call(rbind, lapply(rownames(mat), function(tag) {
    vec <- mat[tag, ]
    q <- quantile(vec, c(0.025, 0.975), na.rm = TRUE)
    data.frame(
      Method = method_label,
      SampleSize = labels[[tag]][1],
      Scenario = paste0(labels[[tag]][2], " indirect effect"),
      Bias_mean = mean(vec, na.rm = TRUE),
      Bias_lower = q[1],
      Bias_upper = q[2]
    )
  }))
}

# === 1. Load and process SEMA results ===
sema <- readRDS("sema.rds")

true_signals <- c(1, 2, 3, 11, 12, 13, 21, 22, 23)
Alpha_a <- matrix(0, nrow = 1, ncol = 30)
Alpha_a[c(1, 11, 21)] <- 0.3
Alpha_a[c(2, 12, 22)] <- 0.6
Alpha_a[c(3, 13, 23)] <- 0.9
Beta_m <- 0.4
true_ie_sum <- sum(as.vector(t(Alpha_a) %*% Beta_m))

RBias1 <- RBias2 <- matrix(0, nrow = n_scenario, ncol = n_dataset)
rownames(RBias1) <- rownames(RBias2) <- names(labels)

for (i in 1:n_scenario) {
  RBias1[i, ] <- abs((rowSums(sema$IE1[[i]]) - true_ie_sum) / true_ie_sum)
  RBias2[i, ] <- abs((rowSums(sema$IE2[[i]]) - true_ie_sum) / true_ie_sum)
}

sema_bias_df_unadjusted <- build_bias_df(RBias1, "SE-MA (Unadjusted)")
sema_bias_df_adjusted   <- build_bias_df(RBias2, "SE-MA (Adjusted)")

# === 2. Load and process PCA results ===
pca <- readRDS("pca.rds")
res_l <- readRDS("pca_l.rds")
res_h <- readRDS("pca_h.rds")
top_res_h <- readRDS("top_pca_h.rds")
top_res_l <- readRDS("top_pca_l.rds")

true_value_pca <- c(sum(res_l$`PNIE Table`$Estimate),
                    sum(res_h$`PNIE Table`$Estimate))
true_value_first <- c(top_res_l$PCA_first["pnie", "Estimate"], top_res_h$PCA_first["pnie", "Estimate"])
true_value_three <- c(
  sum(
    top_res_l$PCA_three$PC1["pnie", "Estimate"],
    top_res_l$PCA_three$PC2["pnie", "Estimate"],
    top_res_l$PCA_three$PC3["pnie", "Estimate"]
  ),
  sum(
    top_res_h$PCA_three$PC1["pnie", "Estimate"],
    top_res_h$PCA_three$PC2["pnie", "Estimate"],
    top_res_h$PCA_three$PC3["pnie", "Estimate"]
  )
)

RBias_pca <- RBias_first <- RBias_three <- matrix(0, nrow = n_scenario, ncol = n_dataset)
rownames(RBias_pca) <- rownames(RBias_first) <- rownames(RBias_three) <- names(labels)

for (i in 1:n_scenario) {
  true_ie <- ifelse(i %in% c(1, 3), true_value_pca[1], true_value_pca[2])
  true_ie_first <- ifelse(i %in% c(1, 3), true_value_first[1], true_value_first[2])
  true_ie_three <- ifelse(i %in% c(1, 3), true_value_three[1], true_value_three[2])
  for (j in 1:100) {
    RBias_pca[i, j] <- abs((sum(pca$IE[[i]][[j]]) - true_ie) / true_ie)
    RBias_first[i, j] <- abs((pca$IE[[i]][[j]][1] - true_ie_first) / true_ie_first)
    RBias_three[i, j] <- abs((sum(pca$IE[[i]][[j]][1:3]) - true_ie_three) / true_ie_three)
  }
}

pca_bias_df <- build_bias_df(RBias_pca, "PC-MA")
bias_df_first <- build_bias_df(RBias_first, "PC-MA (First PC)")
bias_df_three <- build_bias_df(RBias_three, "PC-MA (Top 3 PCs)")

# === 3. Load and process ERS results ===
ers <- readRDS("ers1.rds")
true_data <- readRDS("true_ers.rds")
true_value_ers <- c(true_data[[1]][4, 1], true_data[[2]][4, 1])

RBias_ers <- matrix(0, nrow = 4, ncol = 100)
rownames(RBias_ers) <- names(labels)

for (i in 1:4) {
  true_ie <- ifelse(i %in% c(1, 3), true_value_ers[1], true_value_ers[2])
  for (j in 1:100) {
    RBias_ers[i, j] <- abs((ers$IE[j, i] - true_ie) / true_ie)
  }
}

ers_bias_df <- build_bias_df(RBias_ers, "ERS")

# === 4. Combine all into a single data frame ===
bias_df_all <- rbind(
  sema_bias_df_unadjusted,
  sema_bias_df_adjusted,
  pca_bias_df,
  bias_df_first,
  bias_df_three,
  ers_bias_df
)
bias_df_all$Method <- factor(
  bias_df_all$Method,
  levels = c(
    "SE-MA (Unadjusted)",
    "SE-MA (Adjusted)",
    "PC-MA",
    "PC-MA (First PC)",
    "PC-MAA (Top 3 PCs)",
    "ERS-MA"
  )
)

# Plot
ggplot(bias_df_all, aes(x = Method, y = Bias_mean, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(
    aes(ymin = Bias_lower, ymax = Bias_upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  ggh4x::facet_grid2(Scenario ~ SampleSize) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = "Percent Relative Bias", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )

# Plot (Log10 scale Y axis)
ggplot(bias_df_all, aes(x = Method, y = Bias_mean, fill = Method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(
    aes(ymin = Bias_lower, ymax = Bias_upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  ggh4x::facet_grid2(Scenario ~ SampleSize) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 0.01),
    labels = scales::label_percent(accuracy = 1),
    breaks = c(0.01, 0.05, 0.1, 0.5, 1, 2, 5)
  ) +
  labs(x = NULL, y = "Percent Relative Bias", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )

