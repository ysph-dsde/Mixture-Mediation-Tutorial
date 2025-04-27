# This is the R code script to first do the data cleaning for the PROTECT data
# This code is adopted from the PROTECT data cleaning code written by Dr. Jonathan Boss

pkgs <- c("tidyverse", "Rcpp", "RcppArmadillo",
          "corrplot", "infinitefactor", "coda",
          "doParallel", "hdi")

suppressMessages(invisible(lapply(pkgs, library, character.only = T)))

# read in the data
exps <- read.csv("Data/dat.exp.csv")
out <- read.csv("Data/dat.out.csv")
eic <- read.csv("Data/eic.data.csv")
eic.groupings <- read.csv("Data/eic.names.csv")
extra.ids <- read.csv("Data/IDs_missing_info.csv")

################# General Set up ######################

# set up start
all.outcomes <- c("FINALGA_BEST", "WEIGHTZSCORE", "HEADCIRCUMFERENCEZSCORE")

all.covar <- list(c("MAGE_CAT", "EDU_CAT"),
                  c("MAGE_CAT", "EDU_CAT", "PREBMI_CAT"),
                  c("MAGE_CAT", "EDU_CAT", "PREBMI_CAT"))

all.C_frm <- list(as.formula("~ factor(MAGE_CAT) + factor(EDU_CAT)"),
                  as.formula("~ factor(MAGE_CAT) + factor(EDU_CAT) + factor(PREBMI_CAT)"),
                  as.formula("~ factor(MAGE_CAT) + factor(EDU_CAT) + factor(PREBMI_CAT)"))



############################ Data Cleaning (HEADCIRCUMFERENCEZSCORE) - PHTH only ##################################

# assign arrayid
arrayid <- 3

outcome <- all.outcomes[arrayid]
covariate_vec <- all.covar[[arrayid]]
C_frm <- all.C_frm[[arrayid]]


#Exposure Data Cleaning
# Replace specific values with 0 in  (blod_OH) columns
col_2_replace_97 <- exps %>% dplyr::select(contains("blod_OH")) %>% names

exps[col_2_replace_97] <- lapply(exps[col_2_replace_97], function(x) {
  replace(x, x == 97, 0)
})

# Replace specific values with NA in (`metal`_urine) columns
col_2_replace_min999 <- exps %>% dplyr::select(ends_with("urine")) %>% names

exps[col_2_replace_min999] <- lapply(exps[col_2_replace_min999], function(x) {
  replace(x, x == -999, NA)
})

#Percent Below LOD
prop.exp.blod <- apply(exps %>% dplyr::select(all_of(names(exps)[grepl("BLOD", toupper(names(exps)))])), 2,
                       function(x){sum(x, na.rm = TRUE)})/apply(exps %>% dplyr::select(all_of(names(exps)[grepl("BLOD", toupper(names(exps)))])), 2,
                                                                function(x){sum(!is.na(x))})

exp.above.lod <- names(prop.exp.blod[prop.exp.blod < 0.5])
exp.above.lod.upper <- gsub("(^BLOD_|\\.BLOD)", "", toupper(exp.above.lod))
exp.above.lod.upper <- exp.above.lod.upper[!(exp.above.lod.upper %in% c("AMPA","GLY"))]

names(exps) <- toupper(names(exps))
exps <- exps %>% dplyr::select(ID, VISIT, SG, all_of(exp.above.lod.upper))

exp.subset <- exps %>% dplyr::select(-c(PFHXS, PFNA, PFOS, SE.URINE, SN.BLOOD, SB.BLOOD, MO.BLOOD, MECPTP, MONP, MEHHTP,
                                        ZN.BLOOD, PB.BLOOD, AS.BLOOD, CD.BLOOD, CO.BLOOD, CS.BLOOD, CU.BLOOD, HG.BLOOD,
                                        MN.BLOOD, NI.BLOOD)) %>%
  dplyr::select(-c(MHBP, MHIBP, BPS)) %>% na.omit()

#Specific Gravity Adjustment
exp.nms.for.sg <- names(exp.subset)[!(names(exp.subset) %in% c("ID","VISIT","SG"))]
for(i in 1:length(exp.nms.for.sg)){
  exp.subset[[exp.nms.for.sg[i]]] <- exp.subset[[exp.nms.for.sg[i]]]*((median(exp.subset$SG)-1)/(exp.subset$SG-1))
}

#Eicosanoid Cleaning
eic <- eic %>% filter(visitid == 3)
eic.subset <- eic %>% dplyr::select(-all_of(eic.groupings$Variable.Name[eic.groupings$Pathway.Group == "Cyclooxygenase Pathway"])) %>%
  dplyr::select(-c(X, visitid, LPXB4, LTB3, LTC4.ME, X19HETE, X9HETE)) %>%
  dplyr::select(-all_of(eic.groupings$Variable.Name[eic.groupings$Pathway.Group == "Parent Compound"])) %>% 
  na.omit()

names(eic.subset) <- gsub("^STUDYID$", "ID", toupper(names(eic.subset)))

#A few values of LTE4 are equal to zero, need to set to a small, postive value to avoid issues with log-transformation
eic.subset$LTE4[eic.subset$LTE4 == 0] <- min(eic.subset$LTE4[eic.subset$LTE4 != 0])/2

#Merge Eicosanoid Data with Exposure Data
exp.eic.df <- merge(exp.subset, eic.subset, by = "ID")

#Merge Outcome and Covariate Data
names(out) <- toupper(names(out))
names(extra.ids) <- toupper(names(extra.ids))

out <- out %>% dplyr::select(ID, HEADCIRCUMFERENCEZSCORE, WEIGHTZSCORE, FINALGA_BEST, 
                             MAGE_CAT, EDU_CAT, PREBMI_CAT) 

extra.ids <- extra.ids %>% dplyr::select(ID, HEADCIRCUMFERENCEZSCORE, WEIGHTZSCORE, FINALGA_BEST,
                                         MAGE_CAT, EDU_CAT, PREBMI_CAT) 
out <- rbind(out, extra.ids)

clean.df <- merge(out, exp.eic.df, by = "ID") %>% 
  filter(!is.na(HEADCIRCUMFERENCEZSCORE) & !is.na(FINALGA_BEST) & !is.na(PREBMI_CAT))

##########################################################################3
## Correlation Heat Map
##########################################################################3

A.nms <- names(exp.subset)
A.nms <- A.nms[!(A.nms %in% c("ID", "VISIT", "SG"))]
A.nms <- A.nms[1:11] # include phthalates only

# correlation plot for exposures
exp.cors <- cor(clean.df %>% dplyr::select(all_of(A.nms)), method = "spearman")

corrplot(exp.cors, method = "color", tl.col = "black")


M.nms <- names(eic.subset)[names(eic.subset) != "ID"]

eic.groupings$Abbreviation <- gsub(" ", "", eic.groupings$Abbreviation)
eic.groupings$Variable.Name <- toupper(eic.groupings$Variable.Name)

eic.groupings <- eic.groupings %>% filter(Variable.Name %in% M.nms)

clean.df <- cbind(clean.df %>% dplyr::select(-eic.groupings$Variable.Name),
                  clean.df %>% dplyr::select(eic.groupings$Variable.Name))
M.nms <- eic.groupings$Variable.Name

# correlation plot for mediators
exp.eics <- cor(clean.df %>% dplyr::select(all_of(M.nms)), method = "spearman")

corrplot(exp.eics, method = "color", tl.col = "black")


#Remove missing values for relevant covariates and outcomes
clean.df.spec.out <- clean.df %>% filter(!is.na(.data[[outcome]])) # outcome here is HEADCIRCUMFERENCEZSCORE
for(i in 1:length(covariate_vec)){
  clean.df.spec.out <- clean.df.spec.out %>% filter(!is.na(.data[[covariate_vec[i]]])) # extract the complete cases
}

## These are the things you need!
C_confound <- model.matrix(C_frm, data = clean.df.spec.out %>% filter(VISIT == 1)) # 175 obs 6 confounders (including intercept)
A <- scale(as.matrix(log(clean.df.spec.out %>% filter(VISIT == 1) %>%  # 175 obs 11 exposures (40 if no picking)
                           dplyr::select(all_of(A.nms)))))
A.iqr <- apply(A, 2, IQR)
M <- scale(as.matrix(log(clean.df.spec.out %>% filter(VISIT == 1) %>% 
                           dplyr::select(all_of(M.nms))))) # 175 obs and 15 mediators (if no picking, then 32 mediators)

Y <- (clean.df.spec.out %>% filter(VISIT == 1))[[outcome]] # 175 obs 1 outcome (HEADCIRCUMFERENCEZSCORE)


## Tidy the cleaned data
colnames(C_confound) <- c("Intercept", "MAGE_CAT_fac2", "MAGE_CAT_fac3", "MAGE_CAT_fac4",
                          "EDU_CAT_fac2", "EDU_CAT_fac3", "PRE_BMI_fac2", "PRE_BMI_fac3")

df_protect <- cbind(Y, M, A, C_confound) %>% as.data.frame() %>% 
  dplyr::rename(HEADCIRCUMFERENCEZSCORE = Y)


list_headZscore <- list(Data = df_protect,
                     Outcome = Y,
                     Mediators = M, # assign arrayid
                     Exposures = A,
                     Confounders = C_confound)

write_rds(list_headZscore, "RDS/list_headZscore.rds")