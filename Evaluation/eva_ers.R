# Load required packages
library(tidyverse)

set.seed(1211)

# Global settings
n_dataset <- 100  # Number of simulated datasets per scenario
n_scenario <- 4   # Number of scenarios
sce <- c("1000l", "1000h", "2500l", "2500h")  # Scenario labels
data <- readRDS("../ers1.rds")

## Obtain true value
true_value <- rep(0, 2)
degree <- c("l", "h")
for (i in 1:length(degree)) {
  true_data <- readRDS(paste0("ers1_", degree[i], ".rds"))
  true_value[i] <- true_data$ERS_Case1$`CMA Summary Table`[4, 1]
}
rm(true_data)
gc()

## Relative Bias
RBias <- matrix(0, nrow = n_scenario, ncol = n_dataset)
for (i in 1:n_scenario) {
  true_ie <- ifelse(i <= 2, true_value[1], true_value[2])
  for (j in 1:n_dataset) {
    RBias[i, j] <- abs((data$IE[j, i] - true_ie) /
                         true_ie)
  }
}
rownames(RBias) <- c("1000l", "1000h", "2500l", "2500h")
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
    vec <- mat[tag,]
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

# Construct long-format data frame
bias_df <- build_bias_df(RBias, "ERS")

# Plot
ggplot(bias_df, aes(x = Method, y = Bias_mean, fill = Scenario)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(
    aes(ymin = Bias_lower, ymax = Bias_upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~ SampleSize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Percent Relative Bias", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
