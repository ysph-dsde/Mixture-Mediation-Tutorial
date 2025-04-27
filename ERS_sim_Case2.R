# This is the code script to run the ERS simulation for Case 2
# ERS case 2 is the case where the more complex form of the ERS score construction is used
# which considers the interactions and squared terms of each exposures when running the Elastic Net variable selection

pkgs <- c("tidyverse",
          "MASS",
          "corrplot",
          "gcdnet",
          "purrr",
          "CMAverse",
          "knitr",
          "kableExtra")
invisible(lapply(pkgs, library, character.only = T))


# Read in the data
list_df <- read_rds("RDS/list_df.rds")
list_medTest <- read_rds("RDS/list_medTest.rds")

set.seed(1211)

source("Functions/Functions_ERS.R")

############################### Generate the ERS score  ##############################

# The Case 2 of the ERS score calculation
# includes the exposures (main effects), interactions, and squared terms

# The functions for calculating the ERS score is outlined in the "Functions_ERS_alt.R"
tmp_df <- list_df[["ERS_Case2"]]
list_df[["ERS_Case2"]] <- NULL

tictoc::tic()
list_df[["ERS_Case2"]][["ERS Result"]] <- ers_Calc(
  data = tmp_df,
  exposure = paste0("x", 1:30),
  outcome = "y",
  covar = paste0("c", 1:5),
  include_int = T,
  lambda2_start = exp(seq(log(1e-4), log(1e2), length.out = 100)),
  seed = 1211
)
tictoc::toc() # 3 hrs

list_df[["ERS_Case2"]][["Data"]] <- list_df[["ERS_Case2"]][["ERS Result"]]$post_ERS_data

list_df[["ERS_Case2"]][["ERS Result"]]$ers_fit$lambda # the best lambda value for this

analysis_ids <- which(list_df[["ERS_Case2"]][["Data"]][, "ERS"] != 0)

## Correlation Plot of the data
tmp_df <- list_df[["ERS_Case2"]][["Data"]] %>% dplyr::select(!c(contains("x"), "intercept"))

# ERS Case 2
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
  title = "Correlation for Generated Data - ERS (Main Effects and Interaction)",
  mar = c(0, 0, 3, 0)
)


############################# Test and Estimate the Mediation Effect of the ERS score  #########################
set.seed(1211)

astar <- quantile(list_df[["ERS_Case2"]][["Data"]][analysis_ids, "ERS"], 0.25)
a <- quantile(list_df[["ERS_Case2"]][["Data"]][analysis_ids, "ERS"], 0.75)
## Mediation Test on ERS_scores
tictoc::tic()
list_medTest[["ERS_Case2"]][["CMA Test"]] <- CMAverse::cmest(
  data = list_df[["ERS_Case2"]][["Data"]][analysis_ids, ],
  model = "rb",
  full = T,
  EMint = F,
  yreg = "linear",
  mreg =  list("linear"),
  mval = list(1),
  basec = paste0("c", seq(5)),
  outcome = "y",
  exposure = "ERS",
  mediator = "m1",
  a = a,
  astar = astar,
  inference = "bootstrap",
  nboot = 1000,
  boot.ci.type = "per"
)
tictoc::toc() # 13 seconds

summary(list_medTest[["ERS_Case2"]][["CMA Test"]])

list_medTest[["ERS_Case2"]][["CMA Summary Table"]] <- cbind(
  list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.pe,
  list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.se,
  list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.ci.low,
  list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.ci.high,
  list_medTest[["ERS_Case2"]][["CMA Test"]]$effect.pval
)
colnames(list_medTest[["ERS_Case2"]][["CMA Summary Table"]]) <- c("Estimate", "SE", "CI_Low", "CI_Upper", "Pval")

list_medTest[["ERS_Case2"]][["CMA Summary Table"]]


#################### Output for Future use ######################

write_rds(list_df, "RDS/list_df.rds")
write_rds(list_medTest, "RDS/list_medTest.rds")
