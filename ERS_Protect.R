# This is the script for ERS method used on the PROTECT data example

pkgs <- c("tidyverse", "MASS", "CMAverse", "corrplot", "gcdnet")

invisible(lapply(pkgs, library, character.only = T))

list_headZscore <- read_rds("RDS/list_headZscore.rds")

set.seed(1211)

source("Functions/Functions_ERS.R")


######################################### LTE4 and head Zscore (PHTH only) #######################################

########################## Run ERS simple Case ###############################3


expo_nm <- colnames(list_headZscore[["Exposures"]])
confounders_nm <- colnames(list_headZscore[["Confounders"]])[!grepl("Intercept", colnames(list_headZscore[["Confounders"]]))]
mediator_nm <- "LTE4"
outcome_nm <- colnames(list_headZscore[["Data"]])[1]


tmp_df <- list_headZscore[["Data"]]

tictoc::tic()
headZscore_ERS_result <- ers_Calc(
  data = tmp_df,
  exposure = expo_nm,
  outcome = outcome_nm,
  covar = confounders_nm,
  include_int = F,
  lambda2_start = exp(seq(log(1e-4), log(1e2), length.out = 100)),
  seed = 1211
)
tictoc::toc() # 4.36 secs

# optimal lambdas picked up 11 exposures

list_headZscore[["ERS (Simple)"]][["ERS Result"]] <- headZscore_ERS_result
list_headZscore[["ERS (Simple)"]][["ERS_Data"]] <- headZscore_ERS_result$post_ERS_data

list_headZscore[["ERS (Simple)"]][["ERS_Data"]]$ERS

analysis_ids <- which(list_headZscore[["ERS (Simple)"]][["ERS_Data"]][, "ERS"] !=
                        0)

## Correlation Plot of the data
tmp_df <- list_headZscore[["ERS (Simple)"]][["ERS_Data"]] %>% dplyr::select(c(outcome_nm, mediator_nm, "ERS"))

corrplot(
  cor(tmp_df[analysis_ids, ]),
  order = "original",
  method = "square",
  type = "full",
  number.cex = 0.75,
  diag = T,
  cl.pos = "n",
  addCoef.col = "black",
  tl.srt = 0.005,
  tl.col = "black",
  tl.offset = 0.9,
  tl.cex = 1,
  title = "Correlation for PROTECT Data - ERS (Main Effects Only)",
  mar = c(0, 0, 3, 0)
)


## optimal lambda
list_headZscore[["ERS (Simple)"]][["ERS Result"]][["ers_fit"]]$lambda

## Examine the variable selection results
coef_names <- colnames(list_headZscore[["ERS (Simple)"]][["ERS Result"]][["dat_score"]])
names(list_headZscore[["ERS (Simple)"]][["ERS Result"]][["coef"]]) <- coef_names

list_headZscore[["ERS (Simple)"]][["ERS Result"]][["coef"]] %>% round(3)

############################# MedTest on ERS (Simple) #########################3
set.seed(1211)

astar <- quantile(list_headZscore[["ERS (Simple)"]][["ERS_Data"]][analysis_ids, "ERS"], 0.25)
a <- quantile(list_headZscore[["ERS (Simple)"]][["ERS_Data"]][analysis_ids, "ERS"], 0.75)
## Mediation Test on ERS_scores
tictoc::tic()
list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Test"]] <- CMAverse::cmest(
  data = list_headZscore[["ERS (Simple)"]][["ERS_Data"]][analysis_ids, ],
  model = "rb",
  full = T,
  EMint = F,
  yreg = "linear",
  mreg =  list("linear"),
  mval = list(1),
  basec = confounders_nm,
  outcome = outcome_nm,
  exposure = "ERS",
  mediator = mediator_nm,
  a = a,
  astar = astar,
  estimation = "paramfunc",
  inference = "delta"
)
tictoc::toc() # 0.01 sec


list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Summary Table"]] <- cbind(
  list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Test"]]$effect.pe,
  list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Test"]]$effect.se,
  list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Test"]]$effect.ci.low,
  list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Test"]]$effect.ci.high,
  list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Test"]]$effect.pval
)
colnames(list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Summary Table"]]) <- c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")

list_headZscore[["medTest"]][["ERS (Simple)"]][["CMA Summary Table"]]


########################## Run ERS Complex Case ###############################3

expo_nm <- colnames(list_headZscore[["Exposures"]])
confounders_nm <- colnames(list_headZscore[["Confounders"]])[!grepl("Intercept", colnames(list_headZscore[["Confounders"]]))]
mediator_nm <- "LTE4"
outcome_nm <- colnames(list_headZscore[["Data"]])[1]


tmp_df <- list_headZscore[["Data"]]

tictoc::tic()
headZscore_ERS_result_comp <- ers_Calc(
  data = tmp_df,
  exposure = expo_nm,
  outcome = outcome_nm,
  covar = confounders_nm,
  include_int = T,
  lambda2_start = exp(seq(log(1e-4), log(1e2), length.out = 100)),
  seed = 1211
)
tictoc::toc() # 1m 6s

list_headZscore[["ERS (Complex)"]][["ERS Result"]] <- headZscore_ERS_result_comp
list_headZscore[["ERS (Complex)"]][["ERS_Data"]] <- headZscore_ERS_result_comp$post_ERS_data

list_headZscore[["ERS (Complex)"]][["ERS_Data"]]$ERS

analysis_ids <- which(list_headZscore[["ERS (Complex)"]][["ERS_Data"]][, "ERS"] !=
                        0)

## Correlation Plot of the data
tmp_df <- list_headZscore[["ERS (Complex)"]][["ERS_Data"]] %>% dplyr::select(c(outcome_nm, mediator_nm, "ERS"))

corrplot(
  cor(tmp_df[analysis_ids, ]),
  order = "original",
  method = "square",
  type = "full",
  number.cex = 0.75,
  diag = T,
  cl.pos = "n",
  addCoef.col = "black",
  # this is to add numbers to it
  tl.srt = 0.005,
  tl.col = "black",
  tl.offset = 0.9,
  tl.cex = 1,
  title = "Correlation for PROTECT Data - ERS (Main Effects and Interactions)",
  mar = c(0, 0, 3, 0)
)


## optimal lambda
list_headZscore[["ERS (Complex)"]][["ERS Result"]][["ers_fit"]]$lambda

## Examine the variable selection results
coef_names <- colnames(list_headZscore[["ERS (Complex)"]][["ERS Result"]][["dat_score"]])
names(list_headZscore[["ERS (Complex)"]][["ERS Result"]][["coef"]]) <- coef_names

list_headZscore[["ERS (Complex)"]][["ERS Result"]][["coef"]]

list_headZscore[["ERS (Complex)"]][["ERS Result"]][["ers_fit"]]$beta

############################# MedTest on ERS (Complex) #########################3
set.seed(1211)

astar <- quantile(list_headZscore[["ERS (Complex)"]][["ERS_Data"]][analysis_ids, "ERS"], 0.25)
a <- quantile(list_headZscore[["ERS (Complex)"]][["ERS_Data"]][analysis_ids, "ERS"], 0.75)
## Mediation Test on ERS_scores
tictoc::tic()
list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Test"]] <- CMAverse::cmest(
  data = list_headZscore[["ERS (Complex)"]][["ERS_Data"]][analysis_ids, ],
  model = "rb",
  full = T,
  EMint = F,
  yreg = "linear",
  mreg =  list("linear"),
  mval = list(1),
  basec = confounders_nm,
  outcome = outcome_nm,
  exposure = "ERS",
  mediator = mediator_nm,
  a = a,
  astar = astar,
  estimation = "paramfunc",
  inference = "delta"
)
tictoc::toc() # 0.03 secs

list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Summary Table"]] <- cbind(
  list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Test"]]$effect.pe,
  list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Test"]]$effect.se,
  list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Test"]]$effect.ci.low,
  list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Test"]]$effect.ci.high,
  list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Test"]]$effect.pval
)
colnames(list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Summary Table"]]) <- c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")

list_headZscore[["medTest"]][["ERS (Complex)"]][["CMA Summary Table"]]

### Output ###3
write_rds(list_headZscore, "RDS/list_headZscore.rds")
