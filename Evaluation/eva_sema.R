# Load required packages
library(ggplot2)
data <- readRDS("sema.rds")

set.seed(1211)

# Global settings
n_dataset <- 100  # Number of simulated datasets per scenario
n_scenario <- 4   # Number of scenarios
sce <- c("1000l", "1000h", "2500l", "2500h")  # Scenario labels
true_signals <-
  c(1, 2, 3, 11, 12, 13, 21, 22, 23)  # True mediating exposures
Alpha_a <- matrix(0, nrow = 1, ncol = 30)
Alpha_a[c(1, 11, 21)] <- 0.3 # weak effect
Alpha_a[c(2, 12, 22)] <- 0.6 # moderate effect
Alpha_a[c(3, 13, 23)] <- 0.9 # strong effect
Beta_m <- 0.4
true_ie <- as.vector(t(Alpha_a) %*% Beta_m)
true_ie_sum <- sum(true_ie)

## Relative Bias
RBias1 <- RBias2 <- matrix(0, nrow = 4, ncol = 100)
for(i in 1:4){
    RBias1[i, ] <- abs((rowSums(data$IE1[[i]]) - true_ie_sum)/true_ie_sum)
    RBias2[i, ] <- abs((rowSums(data$IE2[[i]]) - true_ie_sum)/true_ie_sum)
}
rownames(RBias1) <- rownames(RBias2) <- c("1000l", "1000h", "2500l", "2500h")
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

# Construct long-format data frame
bias_df <- rbind(
  build_bias_df(RBias2, "Adjusted"),
  build_bias_df(RBias1, "Unadjusted")
)

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


## TPR
ggplot(data$TPR, aes(x = Method, y = TPR_mean, fill = Scenario)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(
    aes(ymin = TPR_lower, ymax = TPR_upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~ SampleSize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), ) +
  labs(x = NULL, y = "True Positive Rate", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

## FPR
ggplot(data$FPR, aes(x = Method, y = FPR_mean, fill = Scenario)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7) +
  geom_errorbar(
    aes(ymin = FPR_lower, ymax = FPR_upper),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~ SampleSize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), ) +
  labs(x = NULL, y = "False Positive Rate", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
