# This is the code script for runnning BKMR simulation

pkgs <- c(
  "tidyverse",
  "MASS",
  "corrplot",
  "bkmr",
  "causalbkmr",
  "purrr",
  "CMAverse",
  "knitr",
  "kableExtra"
)
invisible(lapply(pkgs, library, character.only = T))


# Read in the data
list_df <- read_rds("RDS/list_df.rds")
list_medTest <- read_rds("RDS/list_medTest.rds")

set.seed(1211)

################### BKMR - Component-wise Variable Selection #################

# we need to first define the $Z_M$ and $Z_Y$ matrices for the BKMR mediation and outcome models.
#Recall that $Z_M$ is the set of exposures and effect modifiers $E_M$

# create Z.M and Z.Y and Z.YM
#Z.M <- cbind(A, E.M); Z.Y <- cbind(A, E.Y); Zm.Y <- cbind(Z.Y, m)

list_df[["BKMR"]] %>% head

exposures_X <- list_df[["BKMR"]] %>%
  dplyr::select(contains("x")) %>%
  as.matrix()

confounders_C <- list_df[["BKMR"]] %>%
  dplyr::select(starts_with("c")) %>%
  as.matrix()

outcome_Y  <- list_df[["BKMR"]]$y
mediator_M  <- list_df[["BKMR"]]$m1

# we assume no effect modifiers
E.M <- NULL
E.Y <- NULL

# create the matrices
Z.M <- cbind(exposures_X, E.M)
Z.Y <- cbind(exposures_X, E.Y)
Zm.Y <- cbind(Z.Y, mediator_M)


set.seed(1211)

tictoc::tic()
# outcome model
set.seed(1211)
fit.y <- kmbayes(
  y = outcome_Y,
  Z = Zm.Y,
  X = confounders_C,
  iter = 10000,
  verbose = TRUE,
  varsel = TRUE,
  control.params = list(
    lambda.jump = 0.45,
    r.jump1 = 0.1,
    r.muprop = 0.35,
    r.jump2 = 0.25
  )
)
tictoc::toc() # 9.4 hours

# TE model
tictoc::tic()
set.seed(1211)
fit.y.TE <- kmbayes(
  y = outcome_Y,
  Z = Z.Y,
  X = confounders_C,
  iter = 10000,
  verbose = TRUE,
  varsel = TRUE,
  control.params = list(lambda.jump = 0.375, r.jump2 = 0.05)
)
tictoc::toc() # 11.1 hours

# mediator model
set.seed(1211)
tictoc::tic()
fit.m <- kmbayes(
  y = mediator_M,
  Z = Z.M,
  X = confounders_C,
  iter = 10000,
  verbose = TRUE,
  varsel = TRUE,
  control.params = list(lambda.jump = 0.925, r.jump2 = 0.025)
)
tictoc::toc() # 8.5 hours

tmp_df <- list_df[["BKMR"]]

list_df[["BKMR"]] <- NULL
list_df[["BKMR"]][["Data"]] <- tmp_df

list_df[["BKMR"]][["Outcome model"]] <- fit.y
list_df[["BKMR"]][["TE model"]] <- fit.y.TE
list_df[["BKMR"]][["Mediator model"]] <- fit.m

write_rds(list_df, "RDS/list_df.rds")

######################### Examine if the models converge (Component-wise Variable Seelction) #######################

# Beta for Outcome model
par(mfrow = c(3, 1))
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "beta", comp = 1)
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "beta", comp = 2)
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "beta", comp = 3)

# Beta for Outcome model
par(mfrow = c(2, 1))
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "beta", comp = 4)
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "beta", comp = 5)

# sigma_sq epsilon for Outcome model
par(mfrow = c(1, 1))
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "sigsq.eps")

# r_m for Outcome model
par(mfrow = c(3, 1))
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "r", comp = 1)
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "r", comp = 2)
TracePlot(list_df[["BKMR"]][["Outcome model"]], par = "r", comp = 3)

# Univariate functions
pred.resp.univar <- PredictorResponseUnivar(fit = list_df[["BKMR"]][["Outcome model"]])
pred.resp.univar <- pred.resp.univar %>%
  filter(variable %in% c("x4", "x7", "x16", "x18", "x19", "x22", "x24", "x28"))
plot <- ggplot(pred.resp.univar,
               aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable) +
  labs(y = "h(X)", x = "Exposures") +
  theme_bw()
plot

# ggsave(
#   "bkmr_sim_outh.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# 
# tiff(
#   "../pic/bkmr_sim_outh.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

pred.resp.univar <- PredictorResponseUnivar(fit = list_df[["BKMR"]][["TE model"]])
pred.resp.univar <- pred.resp.univar %>%
  filter(variable %in% c("x1", "x4", "x7", "x13", "x16", "x19", "x22", "x25", "x28"))
plot <- ggplot(pred.resp.univar,
               aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable) +
  labs(y = "h(X)", x = "Exposures") +
  theme_bw()
plot

# ggsave(
#   "bkmr_sim_teh.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# 
# tiff(
#   "../pic/bkmr_sim_teh.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

pred.resp.univar <- PredictorResponseUnivar(fit = list_df[["BKMR"]][["Mediator model"]])
pred.resp.univar <- pred.resp.univar %>%
  filter(variable %in% c("x2", "x3", "x13", "x22", "x23"))
plot <- ggplot(pred.resp.univar,
               aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable) +
  labs(y = "h(X)", x = "Exposures") +
  theme_bw()
plot

# ggsave(
#   "bkmr_sim_medh.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# 
# tiff(
#   "../pic/bkmr_sim_medh.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()


######################### Extract and Plot PIPs (Component-wise Variable Selection) #######################

p <- 1
q <- 30
s <- 6
Alpha_c <- matrix(1, nrow = p, ncol = s)
Beta_c <- matrix(1, nrow = s, ncol = 1)
Alpha_a <- matrix(0, nrow = 1, ncol = q)
Alpha_a[c(1, 11, 21)] <- 1 # weak effect
Alpha_a[c(2, 12, 22)] <- 4 # moderate effect
Alpha_a[c(3, 13, 23)] <- 8 # strong effect
Beta_m <- 0.5
Beta_a <- rep(c(5, 0, 0), times = q / 3) %>% as.matrix()
Theta_c <- matrix(rep(0.1, times = q * (s - 1)), nrow = 30)


## TE model
pips <- ExtractPIPs(list_df[["BKMR"]][["TE model"]])

highColor <- which(t(Beta_a + t(Alpha_a) %*% Beta_m) > 0)

pips$exposure <- as.character(pips$variable) %>%
  str_remove(pattern = "x") %>%
  factor(levels = c(1:30))

pips$highlight <- ifelse(as.numeric(as.character(pips$exposure)) %in% highColor,
                         "highlight",
                         "normal")

highColor2 <- ifelse(as.numeric(as.character(pips$exposure)) %in% highColor, "red", "black")

plot <- ggplot(pips, aes(exposure, PIP, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(y = "Posterior Inclusion Probability", x = "Exposures") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  ylim(0, 1) + theme_bw() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14, colour = highColor2),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_sim_tepip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 14,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_tepip.tiff",
#   width = 14,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

## Outcome Model
pips <- ExtractPIPs(list_df[["BKMR"]][["Outcome model"]])

highColor <- c(which(t(Beta_a) > 0), "M")

pips$exposure <- as.character(pips$variable) %>%
  str_remove(pattern = "x") %>% str_remove(pattern = "mediator_") %>%
  factor(levels = c(1:30, "M"))

pips$highlight <- ifelse(as.character(pips$exposure) %in% highColor,
                         "highlight",
                         "normal")

highColor2 <- ifelse(as.character(pips$exposure) %in% highColor, "red", "black")

plot <- ggplot(pips, aes(exposure, PIP, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(y = "Posterior Inclusion Probability", x = "Exposures") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  ylim(0, 1) + theme_bw() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14, colour = highColor2),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_sim_outpip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 14,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_outpip.tiff",
#   width = 14,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

## Mediator Model
pips <- ExtractPIPs(list_df[["BKMR"]][["Mediator model"]])

highColor <- c(which(t(t(Alpha_a) %*% Beta_m) > 0))

pips$exposure <- as.character(pips$variable) %>%
  str_remove(pattern = "x") %>%
  factor(levels = c(1:30))

pips$highlight <- ifelse(as.character(pips$exposure) %in% highColor,
                         "highlight",
                         "normal")

highColor2 <- ifelse(as.character(pips$exposure) %in% highColor, "red", "black")

plot <- ggplot(pips, aes(exposure, PIP, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(y = "Posterior Inclusion Probability", x = "Exposures") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  ylim(0, 1) + theme_bw() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 14, colour = highColor2),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_sim_medpip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 14,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_medpip.tiff",
#   width = 14,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

######################### Causal BKMR Code (astar = 25th quantile; a = 75th quantile) #########################

# set the confounders mean level
X.predict <- matrix(colMeans(confounders_C), nrow = 1)  # take the mean of each confounders


# Define the change in exposure levels for mediation effects estimation
# consider a change in all exposures from their 25th to 75th percentiles

astar <- c(apply(exposures_X, 2, quantile, probs = 0.25)) # the reference level of the exposures at 25th percentile
a <- c(apply(exposures_X, 2, quantile, probs = 0.75)) # the comparative level at 75th percentile

# we can ignore this, we don't have effect modifiers for E_M and E_Y
# Yet if modifiers are considered, you should fix the levels of the modifiers
# e.y10 = quantile(E.Y, probs=0.1)
# e.y90 = quantile(E.Y, probs=0.9)

# The index of the MCMC iterations to be used for inference
sel <- seq(5001, 10000, by = 10)

tictoc::tic()
#medTest_BKMR
medTest_BKMR <- mediation.bkmr(
  a = a,
  astar = astar,
  e.m = NULL,
  e.y = NULL,
  fit.m = fit.m,
  # mediation model
  fit.y = fit.y,
  # outcome model
  fit.y.TE = fit.y.TE,
  # TE model
  X.predict.M = X.predict,
  # the mean confounder level for mediation model
  X.predict.Y = X.predict,
  # the mean confounder level for the outcome model
  m.quant = c(0.1, 0.25, 0.5, 0.75),
  # the quantile values of mediator which the CDE(m) is fixed to
  # m.value, # the specific values of the mediator for estimating CDE(m)
  alpha = 0.05,
  # 95% CI
  sel = sel,
  # MCMC iteration indices for making inference
  seed = 1211,
  K = 50
) # K is the number of MCMC samples of the mediator

tictoc::toc()
# 1.8 Hours (K = 50)

medTest_BKMR$est

list_medTest[["BKMR"]][["BKMR_CMA Results"]] <- medTest_BKMR
list_medTest[["BKMR"]][["Results Table"]] <- medTest_BKMR$est

# forest plot of the estimated effects
plotdf <- as.data.frame(list_medTest[["BKMR"]][["Results Table"]])
plotdf["Effect"] <- rownames(plotdf)

plot <- ggplot(plotdf, aes(Effect, mean, ymin = lower, ymax = upper))  +
  geom_pointrange(position = position_dodge(width = 0.75))  +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    col = "red",
    linewidth = 1
  ) +
  labs(y = "Posterior Mean") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 13))
plot

# ggsave(
#   "bkmr_sim_compeff.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_compeff.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()


############### Output ####################

write_rds(list_medTest, "RDS/list_medTest.rds")
write_rds(list_df, "RDS/list_df.rds")

######################################## BKMR Hierarchical ##################################################


#######################################################################3
# Set up for modeling
########################################################################3

# we need to first define the $Z_M$ and $Z_Y$ matrices for the BKMR mediation and outcome models.
#Recall that $Z_M$ is the set of exposures and effect modifiers $E_M$

# create Z.M and Z.Y and Z.YM
#Z.M <- cbind(A, E.M); Z.Y <- cbind(A, E.Y); Zm.Y <- cbind(Z.Y, m)

list_df[["BKMR (Hierarchical)"]] %>% head

exposures_X <- list_df[["BKMR (Hierarchical)"]] %>%
  dplyr::select(contains("x")) %>%
  as.matrix()

confounders_C <- list_df[["BKMR (Hierarchical)"]] %>%
  dplyr::select(starts_with("c")) %>%
  as.matrix()

outcome_Y  <- list_df[["BKMR (Hierarchical)"]]$y
mediator_M  <- list_df[["BKMR (Hierarchical)"]]$m1

# we assume no effect modifiers
E.M <- NULL
E.Y <- NULL

# create the matrices
Z.M <- cbind(exposures_X, E.M)
Z.Y <- cbind(exposures_X, E.Y)
Zm.Y <- cbind(Z.Y, mediator_M)


set.seed(1211)

#######################################################################3
# Hierarchical Clustering of the exposures
########################################################################3

cor_mat <- cor(exposures_X, method = "pearson")

# hierarchical clustering
hc <- hclust(as.dist(1 - cor_mat))

# Plot the dendrogram
# pdf("../pic/bkmr_sim_hc.pdf",
#     width = 7,
#     height = 5)
plot(
  dendextend::color_branches(hc, k = 3),
  main = "Hierarchical Clustering Dendrogram",
  xlab = "",
  sub = "",
  cex = 0.9
)
# dev.off()

# tiff(
#   "../pic/bkmr_sim_hc.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# plot(
#   dendextend::color_branches(hc, k = 3),
#   main = "Hierarchical Clustering Dendrogram",
#   xlab = "",
#   sub = "",
#   cex = 0.9
# )
# dev.off()

# assign groupings
groups <- cutree(hc, k = 3)

#######################################################################3
# Fitting the BKMR models
########################################################################3

# The argument `y` is a vector of outcome. `Z` is a $n \times q$ matrix of predictors to be included in the $h(\cdot)$ function.
# `X` is an $n \times s$ matrix of confounders and should not include the intercept.
# `iter` is the number of iterations to run the MCMC.
# `varsel` is a Boolean variable of whether to conduct variable selection on the $Z$ in $h(\cdot)$ function
# `groups` is an optional vector of group indicators (of length $q$) for fitting the hierarchical variable selection, given `varsel = T`.
# If `group = NULL` and `varsel = T`, then the component-wise variable selection will be performed

set.seed(1211)

tictoc::tic()
# outcome model
fit_y_hier <- kmbayes(
  y = outcome_Y,
  Z = Zm.Y,
  X = confounders_C,
  iter = 10000,
  verbose = TRUE,
  varsel = TRUE,
  groups = c(groups, 4),
  control.params = list(lambda.jump = 0.45, r.jump2 = 0.025)
)
tictoc::toc() # 9.7 hrs


# TE model
set.seed(1211)
tictoc::tic()
fit_y_TE_hier <- kmbayes(
  y = outcome_Y,
  Z = Z.Y,
  X = confounders_C,
  iter = 10000,
  verbose = TRUE,
  varsel = TRUE,
  groups = groups,
  control.params = list(lambda.jump = 2, r.jump2 = 0.01)
)
tictoc::toc()# 8.4 hrs

# mediator model
set.seed(1211)
tictoc::tic()
fit_m_hier <- kmbayes(
  y = mediator_M,
  Z = Z.M,
  X = confounders_C,
  iter = 10000,
  verbose = TRUE,
  varsel = TRUE,
  groups = groups,
  control.params = list(
    lambda.jump = 3.25,
    r.jump1 = 0.01,
    r.jump2 = 0.01
  )
)
tictoc::toc()# 8.4 hrs

tmp_df <- list_df[["BKMR (Hierarchical)"]]

list_df[["BKMR (Hierarchical)"]] <- NULL
list_df[["BKMR (Hierarchical)"]][["Data"]] <- tmp_df

list_df[["BKMR (Hierarchical)"]][["Outcome model"]] <- fit_y_hier
list_df[["BKMR (Hierarchical)"]][["TE model"]] <- fit_y_TE_hier
list_df[["BKMR (Hierarchical)"]][["Mediator model"]] <- fit_m_hier

write_rds(list_df, "RDS/list_df.rds")

######################### Examine if the models converge (Hierarchical Variable Selection) #######################

#######################################################################
# Trace Plots to examine convergence
########################################################################

# Examine if the models converge
# Beta for Outcome model
par(mfrow = c(3, 1))
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "beta", comp = 1)
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "beta", comp = 2)
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "beta", comp = 3)

# Beta for Outcome model
par(mfrow = c(2, 1))
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "beta", comp = 4)
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "beta", comp = 5)

# sigma_sq epsilon for Outcome model
par(mfrow = c(1, 1))
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "sigsq.eps")

# r_m for Outcome model
par(mfrow = c(3, 1))
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "r", comp = 1)
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "r", comp = 2)
TracePlot(list_df[["BKMR (Hierarchical)"]][["Outcome model"]], par = "r", comp = 3)


######################### Extract and Plot PIPs (Hierarchical Variable Selection) #######################

## TE Model var selection results
## Prep for Group PIP plot
pips <- ExtractPIPs(list_df[["BKMR (Hierarchical)"]][["TE model"]]) %>%
  dplyr::mutate(
    group = case_when(group == 1 ~ "Group 1", group == 2 ~ "Group 2", group == 3 ~ "Group 3"),
    group = factor(group, levels = c("Group 1", "Group 2", "Group 3"))
  )

chem_nms <- colnames(list_df[["BKMR (Hierarchical)"]][["Data"]] %>%
                       dplyr::select(contains("x")))

group_color <- c(rep("#d62d20", times = 5),
                 rep("#008744", times = 10),
                 rep("#0057e7", times = 15))

## prep for conditional PIP plot
pips <- pips %>%
  dplyr::mutate(
    exposure = ifelse(variable == "mediator_M", "M", variable),
    group = as.factor(group),
    group_label = paste(group, "\nGroup PIP:", round(groupPIP, 3))
  )

pips$exposure <- pips$exposure %>% factor(levels = c(chem_nms, "M"))

# Define the order of the levels
levels_order <- c("Group 1 \nGroup PIP: 1",
                  "Group 2 \nGroup PIP: 1",
                  "Group 3 \nGroup PIP: 1")

# Set the order of the new factor
pips$group_label <- factor(pips$group_label, levels = levels_order)

## Conditional PIP showcase (facet by groups)
plot <- ggplot(pips, aes(exposure, condPIP, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "red") +
  facet_wrap( ~ group_label, scales = "free_x") +
  labs(y = "Conditional Posterior Inclusion Probability") +
  scale_color_manual(values = c("#d62d20", "#008744", "#0057e7", "#ffa700"),
                     name = "Groups") +
  ylim(0, 1) + theme_bw() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_sim_hiertepip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 8,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_hiertepip.tiff",
#   width = 8,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

## Outcome Model var selection results
## Prep for Group PIP plot
pips <- ExtractPIPs(list_df[["BKMR (Hierarchical)"]][["Outcome model"]]) %>%
  dplyr::mutate(
    group = case_when(
      group == 1 ~ "Group 1",
      group == 2 ~ "Group 2",
      group == 3 ~ "Group 3",
      group == 4 ~ "Mediator"
    ),
    group = factor(group, levels = c("Group 1", "Group 2", "Group 3", "Mediator"))
  )

chem_nms <- colnames(list_df[["BKMR (Hierarchical)"]][["Data"]] %>%
                       dplyr::select(contains("x")))

group_color <- c(
  rep("#d62d20", times = 5),
  rep("#008744", times = 10),
  rep("#0057e7", times = 15),
  "#ffa700"
)

## prep for conditional PIP plot
pips <- pips %>%
  dplyr::mutate(
    exposure = ifelse(variable == "mediator_M", "M", variable),
    group = as.factor(group),
    group_label = paste(group, "\nGroup PIP:", round(groupPIP, 3))
  )

pips$exposure <- pips$exposure %>% factor(levels = c(chem_nms, "M"))

# Define the order of the levels
levels_order <- c(
  "Group 1 \nGroup PIP: 1",
  "Group 2 \nGroup PIP: 1",
  "Group 3 \nGroup PIP: 1",
  "Mediator \nGroup PIP: 1"
)

# Set the order of the new factor
pips$group_label <- factor(pips$group_label, levels = levels_order)

## Conditional PIP showcase (facet by groups)
plot <- ggplot(pips, aes(exposure, condPIP, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "red") +
  facet_wrap( ~ group_label, scales = "free_x") +
  labs(y = "Conditional Posterior Inclusion Probability") +
  scale_color_manual(values = c("#d62d20", "#008744", "#0057e7", "#ffa700"),
                     name = "Groups") +
  ylim(0, 1) + theme_bw() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_sim_hieroutpip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 8,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_hieroutpip.tiff",
#   width = 8,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

## Mediator Model var selection results
## Prep for Group PIP plot
pips <- ExtractPIPs(list_df[["BKMR (Hierarchical)"]][["Mediator model"]]) %>%
  dplyr::mutate(
    group = case_when(group == 1 ~ "Group 1", group == 2 ~ "Group 2", group == 3 ~ "Group 3"),
    group = factor(group, levels = c("Group 1", "Group 2", "Group 3"))
  )

chem_nms <- colnames(list_df[["BKMR (Hierarchical)"]][["Data"]] %>%
                       dplyr::select(contains("x")))

group_color <- c(rep("#d62d20", times = 5),
                 rep("#008744", times = 10),
                 rep("#0057e7", times = 15))

## prep for conditional PIP plot
pips <- pips %>%
  dplyr::mutate(
    exposure = ifelse(variable == "mediator_M", "M", variable),
    group = as.factor(group),
    group_label = paste(group, "\nGroup PIP:", round(groupPIP, 3))
  )

pips$exposure <- pips$exposure %>% factor(levels = c(chem_nms, "M"))

# Define the order of the levels
levels_order <- c("Group 1 \nGroup PIP: 1",
                  "Group 2 \nGroup PIP: 1",
                  "Group 3 \nGroup PIP: 1")

# Set the order of the new factor
pips$group_label <- factor(pips$group_label, levels = levels_order)

## Conditional PIP showcase (facet by groups)
plot <- ggplot(pips, aes(exposure, condPIP, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "red") +
  facet_wrap( ~ group_label, scales = "free_x") +
  labs(y = "Conditional Posterior Inclusion Probability") +
  scale_color_manual(values = c("#d62d20", "#008744", "#0057e7", "#ffa700"),
                     name = "Groups") +
  ylim(0, 1) + theme_bw() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_sim_hiermedpip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 8,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_hiermedpip.tiff",
#   width = 8,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

######################### Causal BKMR Code (astar = 25th quantile; a = 75th quantile) #########################

#######################################################################
# Causal BKMR Code (astar = 25th quantile; a = 75th quantile)
########################################################################

# set the confounders mean level
X.predict <- matrix(colMeans(confounders_C), nrow = 1)  # take the mean of each confounders


# Define the change in exposure levels for mediation effects estimation
# consider a change in all exposures from their 25th to 75th percentiles

astar <- c(apply(exposures_X, 2, quantile, probs = 0.25)) # the reference level of the exposures at 25th percentile
a <- c(apply(exposures_X, 2, quantile, probs = 0.75)) # the comparative level at 75th percentile

# we can ignore this, we don't have effect modifiers for E_M and E_Y
# Yet if modifiers are considered, you should fix the levels of the modifiers
# e.y10 = quantile(E.Y, probs=0.1)
# e.y90 = quantile(E.Y, probs=0.9)

# The index of the MCMC iterations to be used for inference
sel <- seq(5000, 10000, by = 10)

tictoc::tic()
medTest_BKMR <- mediation.bkmr(
  a = a,
  astar = astar,
  e.m = NULL,
  e.y = NULL,
  fit.m = list_df[["BKMR (Hierarchical)"]][["Mediator model"]],
  fit.y = list_df[["BKMR (Hierarchical)"]][["Outcome model"]],
  fit.y.TE = list_df[["BKMR (Hierarchical)"]][["TE model"]],
  X.predict.M = X.predict,
  # the mean confounder level for mediation model
  X.predict.Y = X.predict,
  # the mean confounder level for the outcome model
  m.quant = c(0.1, .25, 0.5, 0.75),
  # the quantile values of mediator which the CDE(m) is fixed to
  # m.value, # the specific values of the mediator for estimating CDE(m)
  alpha = 0.05,
  # 95% CI
  sel = sel,
  # MCMC iteration indices for making inference
  seed = 1211,
  K = 50
) # K is the number of MCMC samples of the mediator

tictoc::toc()
# 1.6 HRs (K = 50)

medTest_BKMR$est

# write_rds(medTest_BKMR, "RDS/medTest_BKMR.rds")

list_medTest[["BKMR (Hierarchical)"]][["BKMR_CMA Results"]] <- medTest_BKMR
list_medTest[["BKMR (Hierarchical)"]][["Results Table"]] <- medTest_BKMR$est

plotdf <- as.data.frame(list_medTest[["BKMR (Hierarchical)"]][["Results Table"]])
plotdf["Effect"] <- rownames(plotdf)

# forest plot of the estimated effects
plot <- ggplot(plotdf, aes(Effect, mean, ymin = lower, ymax = upper))  +
  geom_pointrange(position = position_dodge(width = 0.75))  +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    col = "red",
    linewidth = 1
  ) +
  labs(y = "Posterior Mean", title = "") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 13))
plot
# ggsave(
#   "bkmr_sim_hier_est.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_sim_hier_est.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()



############### Output 2 ####################

write_rds(list_medTest, "RDS/list_medTest.rds")
write_rds(list_df, "RDS/list_df.rds")
