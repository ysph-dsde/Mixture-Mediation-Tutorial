# This code script is to do BKMR on PROTECT data


####################################### Set up ################################

pkgs <- c(
  "tidyverse",
  "MASS",
  "corrplot",
  "gridExtra",
  "bkmr",
  "causalbkmr",
  "purrr",
  "CMAverse",
  "knitr",
  "kableExtra"
)
invisible(lapply(pkgs, library, character.only = T))

list_headZscore <- read_rds("RDS/list_headZscore.rds")

set.seed(1211)


######################################### LTE4 and Head Zscore (PHTH only) #######################################

################################################3
# Preprocessing
###################################################3

exposures_X <- list_headZscore[["Exposures"]] %>%
  as.matrix()

confounders_C <- list_headZscore[["Confounders"]] %>%
  as.data.frame() %>%
  dplyr::select(-Intercept) %>%
  as.matrix()

outcome_Y  <- list_headZscore[["Data"]]$HEADCIRCUMFERENCEZSCORE
mediator_M  <- list_headZscore[["Mediators"]][, "LTE4"]

# we assume no effect modifiers
E.M <- NULL
E.Y <- NULL

# create the matrices
Z.M <- cbind(exposures_X, E.M)
Z.Y <- cbind(exposures_X, E.Y)
Zm.Y <- cbind(Z.Y, mediator_M)


##############################################3
# Fit BKMR Outcome, Mediator, TE Models
###############################################3

tictoc::tic()
set.seed(1211)
# outcome model
fit_y_headZ <- kmbayes(
  y = outcome_Y,
  Z = Zm.Y,
  X = confounders_C,
  iter = 50000,
  verbose = TRUE,
  varsel = TRUE,
  control.params = list(
    lambda.jump = 0.5,
    r.jump1 = 0.01,
    r.jump2 = 0.1
  )
)
tictoc::toc() # 8.3 mins


# TE model
tictoc::tic()
set.seed(1211)
fit_y_TE_headZ <- kmbayes(
  y = outcome_Y,
  Z = Z.Y,
  X = confounders_C,
  iter = 50000,
  verbose = TRUE,
  varsel = TRUE,
  control.params = list(
    lambda.jump = 0.5,
    r.jump1 = 0.01,
    r.jump2 = 0.1
  )
)
tictoc::toc() # 8.3 mins

# mediator model
tictoc::tic()
set.seed(1211)
fit_m_headZ <- kmbayes(
  y = mediator_M,
  Z = Z.M,
  X = confounders_C,
  iter = 50000,
  verbose = TRUE,
  varsel = TRUE,
  control.params = list(
    lambda.jump = 0.5,
    r.jump1 = 0.01,
    r.jump2 = 0.1
  )
)
tictoc::toc() # 8.2 mins


## store the results
list_headZscore[["BKMR"]][["Outcome model"]] <- fit_y_headZ
list_headZscore[["BKMR"]][["TE model"]] <- fit_y_TE_headZ
list_headZscore[["BKMR"]][["Mediator model"]] <- fit_m_headZ


##################################################3
# BKMR model fitting results
#################################################3

# Examine if the models converge
TracePlot(fit = list_headZscore[["BKMR"]][["Outcome model"]], par = "beta", comp = 1)
TracePlot(fit = list_headZscore[["BKMR"]][["Outcome model"]], par = "sigsq.eps")

TracePlot(fit = list_headZscore[["BKMR"]][["TE model"]], par = "beta", comp = 1)
TracePlot(fit = list_headZscore[["BKMR"]][["TE model"]], par = "sigsq.eps")

TracePlot(fit = list_headZscore[["BKMR"]][["Mediator model"]], par = "beta", comp = 1)
TracePlot(fit = list_headZscore[["BKMR"]][["Mediator model"]], par = "sigsq.eps")

pred.resp.univar <- PredictorResponseUnivar(fit = list_headZscore[["BKMR"]][["Outcome model"]])
pred.resp.univar <- pred.resp.univar %>%
  dplyr::mutate(variable = recode(variable, "mediator_M" = "LTE4"))
plot <- ggplot(pred.resp.univar,
               aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable) +
  labs(y = "h(X)", x = "Exposures", title = "") +
  theme_bw()
plot

# ggsave(
#   "bkmr_data_outh.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_data_outh.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

pred.resp.univar <- PredictorResponseUnivar(fit = list_headZscore[["BKMR"]][["TE model"]])
plot <- ggplot(pred.resp.univar,
               aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable) +
  labs(y = "h(X)", x = "Exposures", title = " ") +
  theme_bw()
plot

# ggsave(
#   "bkmr_data_teh.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_data_teh.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

pred.resp.univar <- PredictorResponseUnivar(fit = list_headZscore[["BKMR"]][["Mediator model"]])
plot <- ggplot(pred.resp.univar,
               aes(z, est, ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  geom_smooth(stat = "identity") +
  facet_wrap( ~ variable) +
  labs(y = "h(X)", x = "Exposures", title = "") +
  theme_bw()
plot

# ggsave(
#   "bkmr_data_medh.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_data_medh.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()


##################################3
## BKMR selected Variables
##################################3

# TE Model var selection results
pips <- ExtractPIPs(list_headZscore[["BKMR"]][["TE model"]])

highColor <- which(pips$PIP > 0.8)

phth_nms <- colnames(list_headZscore[["Exposures"]])

pips$exposure <- pips$variable %>% factor(levels = c(phth_nms))

pips$highlight <- ifelse(pips$exposure %in% pips[highColor, "variable"], "highlight", "normal")

highColor2 <- ifelse(pips$exposure %in% pips[highColor, "variable"], "red", "black")

TE_pip <- ggplot(pips, aes(exposure, PIP, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(y = "", x = "", title = "TE Model") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  ylim(0, 1) + theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = highColor2))

# Outcome Model var selection results
pips <- ExtractPIPs(list_headZscore[["BKMR"]][["Outcome model"]]) %>%
  dplyr::mutate(variable = ifelse(variable == "mediator_M", "LTE4", variable))

highColor <- which(pips$PIP > 0.8)

phth_nms <- colnames(list_headZscore[["BKMR"]][["Outcome model"]]$Z)

pips$exposure <- pips$variable %>% factor(levels = c(phth_nms, "LTE4"))

pips$highlight <- ifelse(pips$exposure %in% pips[highColor, "variable"], "highlight", "normal")

highColor2 <- ifelse(pips$exposure %in% pips[highColor, "variable"], "red", "black")

outcome_pip <- ggplot(pips, aes(exposure, PIP, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(y = "Posterior Inclusion Probability", x = "", title = "Outcome Model") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  ylim(0, 1) + theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = highColor2))

# Mediator Model var selection results
pips <- ExtractPIPs(list_headZscore[["BKMR"]][["Mediator model"]])

highColor <- which(pips$PIP > 0.8)

phth_nms <- colnames(list_headZscore[["BKMR"]][["Mediator model"]]$Z)

pips$exposure <- pips$variable %>% factor(levels = c(phth_nms))

pips$highlight <- ifelse(pips$exposure %in% pips[highColor, "variable"], "highlight", "normal")

highColor2 <- ifelse(pips$exposure %in% pips[highColor, "variable"], "red", "black")

med_pip <- ggplot(pips, aes(exposure, PIP, color = highlight)) +
  geom_point() +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(y = "", x = "", title = "Mediator Model") +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) +
  ylim(0, 1) + theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(colour = highColor2))


plot <- grid.arrange(med_pip, outcome_pip, TE_pip, ncol = 1)
plot

# ggsave(
#   "bkmr_data_pip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 10,
#   height = 8,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_data_pip.tiff",
#   width = 10,
#   height = 8,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()


#####################################################################3
## Causal BKMR Code
##    (astar = 25th quantile; a = 75th quantile)
######################################################################3

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
sel <- seq(25001, 50000, by = 10)

tictoc::tic()
#medTest_BKMR
medTest_BKMR <- mediation.bkmr(
  a = a,
  astar = astar,
  e.m = NULL,
  e.y = NULL,
  fit.m = list_headZscore[["BKMR"]][["Mediator model"]],
  # mediation model
  fit.y = list_headZscore[["BKMR"]][["Outcome model"]],
  # outcome model
  fit.y.TE = list_headZscore[["BKMR"]][["TE model"]],
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
# 53 secs (K = 50)

medTest_BKMR$est


list_headZscore[["medTest"]][["BKMR"]][["BKMR_CMA Result"]] <- medTest_BKMR
list_headZscore[["medTest"]][["BKMR"]][["Results Table"]] <- medTest_BKMR$est

write_rds(list_headZscore, "RDS/list_headZscore.rds")


plotdf <- as.data.frame(list_headZscore[["medTest"]][["BKMR"]][["Results Table"]])
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
#   "bkmr_data_est.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff(
#   "../pic/bkmr_data_est.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# print(plot)
# dev.off()

######################## BKMR Hierarchical Variable Selection (based on correlations) - LTE4 and Head Zscore (PHTH only)  ###############################


####################### Testing Hclust but with PROTECT #######################3

exposures_X <- list_headZscore[["Exposures"]] %>%
  as.matrix()

confounders_C <- list_headZscore[["Confounders"]] %>%
  as.data.frame() %>%
  dplyr::select(-Intercept) %>%
  as.matrix()

outcome_Y  <- list_headZscore[["Data"]]$HEADCIRCUMFERENCEZSCORE
mediator_M  <- list_headZscore[["Mediators"]][, "LTE4"]

# we assume no effect modifiers
E.M <- NULL
E.Y <- NULL

# create the matrices
Z.M <- cbind(exposures_X, E.M)
Z.Y <- cbind(exposures_X, E.Y)
Zm.Y <- cbind(Z.Y, mediator_M)



set.seed(1211)

cor_mat <- cor(exposures_X, method = "spearman")

# hierarchical clustering
hc <- hclust(as.dist(1 - cor_mat))

# Plot the dendrogram
# pdf("../pic/bkmr_data_clust.pdf", width = 7, height = 5)
plot(
  dendextend::color_branches(hc, k = 5),
  main = "",
  xlab = "",
  sub = "",
  cex = 0.9
)
# dev.off()
# tiff(
#   "../pic/bkmr_data_clust.tiff",
#   width = 7,
#   height = 5,
#   units = "in",
#   res = 800,
#   compression = "lzw"
# )
# plot(
#   dendextend::color_branches(hc, k = 5),
#   main = "",
#   xlab = "",
#   sub = "",
#   cex = 0.9
# )
# dev.off()

# Cut the tree to form k clusters
groups <- cutree(hc, k = 5)

corrplot::corrplot(
  cor_mat,
  method = "square",
  type = "full",
  diag = T,
  cl.pos = "n",
  number.cex = 0.8,
  tl.col = "black",
  tl.offset = 0.9,
  tl.cex = 0.8,
  tl.srt = 0.005,
  # arguments for specific for hclust results
  order = "hclust",
  hclust.method = "complete",
  addrect = length(unique(groups)),
  rect.col = "red",
  title = "Spearman Correlation Plot of Exposures",
  mar = c(0, 0, 3, 0)
)

##############################################3
# Fit BKMR Outcome, Mediator, TE Models
###############################################3

tictoc::tic()
# outcome model
set.seed(1211)
fit_y_headZ_hier <- kmbayes(
  y = outcome_Y,
  Z = Zm.Y,
  X = confounders_C,
  iter = 50000,
  verbose = TRUE,
  varsel = TRUE,
  groups = c(groups, 6),
  control.params = list(
    lambda.jump = 0.5,
    r.jump1 = 0.01,
    r.jump2 = 0.1
  )
)
tictoc::toc() # 8.5 mins


# TE model
tictoc::tic()
set.seed(1211)
fit_y_TE_headZ_hier <- kmbayes(
  y = outcome_Y,
  Z = Z.Y,
  X = confounders_C,
  iter = 50000,
  verbose = TRUE,
  varsel = TRUE,
  groups = groups,
  control.params = list(
    lambda.jump = 0.5,
    r.jump1 = 0.01,
    r.jump2 = 0.1
  )
)
tictoc::toc() # 8.9 mins


# mediator model
tictoc::tic()
set.seed(1211)
fit_m_headZ_hier <- kmbayes(
  y = mediator_M,
  Z = Z.M,
  X = confounders_C,
  iter = 50000,
  verbose = TRUE,
  varsel = TRUE,
  groups = groups,
  control.params = list(
    lambda.jump = 0.5,
    r.jump1 = 0.01,
    r.jump2 = 0.1
  )
)
tictoc::toc() # 8.9 mins


## store the results
list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]] <- fit_y_headZ_hier
list_headZscore[["BKMR (Hierarchical)"]][["TE model"]] <- fit_y_TE_headZ_hier
list_headZscore[["BKMR (Hierarchical)"]][["Mediator model"]] <- fit_m_headZ_hier

##################################################3
# BKMR model fitting results
#################################################3

# Examine if the models converge
TracePlot(fit = list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]], par = "beta", comp = 1)
TracePlot(fit = list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]], par = "sigsq.eps")
TracePlot(fit = list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]], par = "r", comp = 1)


##################################3
## BKMR selected Variables
##################################3

# Mediator Model var selection results
ExtractPIPs(list_headZscore[["BKMR (Hierarchical)"]][["Mediator model"]])

# Outcome Model var selection results
ExtractPIPs(list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]])

# Total Effect Model var selection results
ExtractPIPs(list_headZscore[["BKMR (Hierarchical)"]][["TE model"]])



# Outcome Model

## Prep for Group PIP plot
pips <- ExtractPIPs(list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]]) %>%
  dplyr::mutate(
    group = case_when(
      group == 1 ~ "Group 1",
      group == 2 ~ "Group 2",
      group == 3 ~ "Group 3",
      group == 4 ~ "Group 4",
      group == 5 ~ "Group 5",
      group == 6 ~ "LTE4"
    ),
    group = factor(
      group,
      levels = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "LTE4")
    )
  )

chem_nms <- colnames(list_headZscore[["Exposures"]])
group_color <- c(
  "#d62d20",
  "#008744",
  rep("#0057e7", times = 3),
  rep("#ffa700", times = 4),
  "#673ab7",
  "#d62d20",
  "#e91e63"
)

pip_group <- pips %>%
  dplyr::select(c("group", "groupPIP")) %>%
  unique()


## prep for conditional PIP plot
pips <- pips %>%
  dplyr::mutate(
    exposure = ifelse(variable == "mediator_M", "LTE4", variable),
    group = as.factor(group)
  )

pips$exposure <- pips$exposure %>% factor(levels = c(chem_nms, "LTE4"))


## Conditional PIP showcase (facet by groups)
pips <- pips %>%
  mutate(group_label = paste(group, "\nGroup PIP:", round(groupPIP, 3)))

# Define the order of the levels
levels_order <- c(
  "Group 1 \nGroup PIP: 0.493",
  "Group 2 \nGroup PIP: 0.355",
  "Group 3 \nGroup PIP: 0.456",
  "Group 4 \nGroup PIP: 0.429",
  "Group 5 \nGroup PIP: 0.423",
  "LTE4 \nGroup PIP: 0.474"
)

# Set the order of the new factor
pips$group_label <- factor(pips$group_label, levels = levels_order)

plot <- ggplot(pips, aes(exposure, condPIP, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "red") +
  facet_wrap( ~ group_label, scales = "free_x") +
  labs(y = "Conditional Posterior Inclusion Probability", x = "Exposure") +
  scale_color_manual(
    values = c(
      "#d62d20",
      "#008744",
      "#0057e7",
      "#ffa700",
      "#673ab7",
      "#e91e63"
    ),
    name = "Groups"
  ) +
  ylim(0, 1) + theme_bw() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_data_hier_outpip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff("../pic/bkmr_data_hier_outpip.tiff", width = 7, height = 5, units = "in", res = 800, compression = "lzw")
# print(plot)
# dev.off()

# Mediator Model

## Prep for Group PIP plot
pips <- ExtractPIPs(list_headZscore[["BKMR (Hierarchical)"]][["Mediator model"]]) %>%
  dplyr::mutate(
    group = case_when(
      group == 1 ~ "Group 1",
      group == 2 ~ "Group 2",
      group == 3 ~ "Group 3",
      group == 4 ~ "Group 4",
      group == 5 ~ "Group 5"
    ),
    group = factor(
      group,
      levels = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
    )
  )

chem_nms <- colnames(list_headZscore[["Exposures"]])
group_color <- c(
  "#d62d20",
  "#008744",
  rep("#0057e7", times = 3),
  rep("#ffa700", times = 4),
  "#673ab7",
  "#d62d20",
  "#e91e63"
)

pip_group <- pips %>%
  dplyr::select(c("group", "groupPIP")) %>%
  unique()


## prep for conditional PIP plot
pips <- pips %>%
  dplyr::mutate(
    exposure = ifelse(variable == "mediator_M", "LTE4", variable),
    group = as.factor(group)
  )

pips$exposure <- pips$exposure %>% factor(levels = c(chem_nms))


## Conditional PIP showcase (facet by groups)
pips <- pips %>%
  mutate(group_label = paste(group, "\nGroup PIP:", round(groupPIP, 3)))

# Define the order of the levels
levels_order <- c(
  "Group 1 \nGroup PIP: 0.375",
  "Group 2 \nGroup PIP: 0.32",
  "Group 3 \nGroup PIP: 0.579",
  "Group 4 \nGroup PIP: 0.34",
  "Group 5 \nGroup PIP: 0.643"
)

# Set the order of the new factor
pips$group_label <- factor(pips$group_label, levels = levels_order)

plot <- ggplot(pips, aes(exposure, condPIP, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "red") +
  facet_wrap( ~ group_label, scales = "free_x") +
  labs(y = "Conditional Posterior Inclusion Probability", x = "Exposure") +
  scale_color_manual(
    values = c(
      "#d62d20",
      "#008744",
      "#0057e7",
      "#ffa700",
      "#673ab7",
      "#e91e63"
    ),
    name = "Groups"
  ) +
  ylim(0, 1) + theme_bw() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_data_hier_medpip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff("../pic/bkmr_data_hier_medpip.tiff", width = 7, height = 5, units = "in", res = 800, compression = "lzw")
# print(plot)
# dev.off()

# TE Model

## Prep for Group PIP plot
pips <- ExtractPIPs(list_headZscore[["BKMR (Hierarchical)"]][["TE model"]]) %>%
  dplyr::mutate(
    group = case_when(
      group == 1 ~ "Group 1",
      group == 2 ~ "Group 2",
      group == 3 ~ "Group 3",
      group == 4 ~ "Group 4",
      group == 5 ~ "Group 5"
    ),
    group = factor(
      group,
      levels = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
    )
  )

chem_nms <- colnames(list_headZscore[["Exposures"]])
group_color <- c(
  "#d62d20",
  "#008744",
  rep("#0057e7", times = 3),
  rep("#ffa700", times = 4),
  "#673ab7",
  "#d62d20",
  "#e91e63"
)

pip_group <- pips %>%
  dplyr::select(c("group", "groupPIP")) %>%
  unique()


## prep for conditional PIP plot
pips <- pips %>%
  dplyr::mutate(
    exposure = ifelse(variable == "mediator_M", "LTE4", variable),
    group = as.factor(group)
  )

pips$exposure <- pips$exposure %>% factor(levels = c(chem_nms))


## Conditional PIP showcase (facet by groups)
pips <- pips %>%
  mutate(group_label = paste(group, "\nGroup PIP:", round(groupPIP, 3)))

# Define the order of the levels
levels_order <- c(
  "Group 1 \nGroup PIP: 0.514",
  "Group 2 \nGroup PIP: 0.397",
  "Group 3 \nGroup PIP: 0.495",
  "Group 4 \nGroup PIP: 0.514",
  "Group 5 \nGroup PIP: 0.451"
)

# Set the order of the new factor
pips$group_label <- factor(pips$group_label, levels = levels_order)

plot <- ggplot(pips, aes(exposure, condPIP, color = group)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, color = "red") +
  facet_wrap( ~ group_label, scales = "free_x") +
  labs(y = "Conditional Posterior Inclusion Probability", x = "Exposure") +
  scale_color_manual(
    values = c(
      "#d62d20",
      "#008744",
      "#0057e7",
      "#ffa700",
      "#673ab7",
      "#e91e63"
    ),
    name = "Groups"
  ) +
  ylim(0, 1) + theme_bw() +
  theme(
    text = element_text(size = 13),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.position = "none"
  )
plot

# ggsave(
#   "bkmr_data_hier_tepip.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff("../pic/bkmr_data_hier_tepip.tiff", width = 7, height = 5, units = "in", res = 800, compression = "lzw")
# print(plot)
# dev.off()

#####################################################################3
## Causal BKMR Code
##    (astar = 25th quantile; a = 75th quantile)
######################################################################3

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
sel <- seq(25001, 50000, by = 10)

tictoc::tic()
#medTest_BKMR
medTest_BKMR <- mediation.bkmr(
  a = a,
  astar = astar,
  e.m = NULL,
  e.y = NULL,
  fit.m = list_headZscore[["BKMR (Hierarchical)"]][["Mediator model"]],
  # mediation model
  fit.y = list_headZscore[["BKMR (Hierarchical)"]][["Outcome model"]],
  # outcome model
  fit.y.TE = list_headZscore[["BKMR (Hierarchical)"]][["TE model"]],
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
# 57 secs (K = 50)

medTest_BKMR$est


list_headZscore[["medTest"]][["BKMR (Hierarchical)"]][["BKMR_CMA Result"]] <- medTest_BKMR
list_headZscore[["medTest"]][["BKMR (Hierarchical)"]][["Results Table"]] <- medTest_BKMR$est


plotdf <- as.data.frame(list_headZscore[["medTest"]][["BKMR (Hierarchical)"]][["Results Table"]])
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
#   "bkmr_data_hier_est.pdf",
#   path = "../pic/",
#   plot = plot,
#   device = "pdf",
#   width = 7,
#   height = 5,
#   units = "in",
#   dpi = 800
# )
# tiff("../pic/bkmr_data_hier_est.tiff", width = 7, height = 5, units = "in", res = 800, compression = "lzw")
# print(plot)
# dev.off()

### Output
write_rds(list_headZscore, "RDS/list_headZscore.rds")
