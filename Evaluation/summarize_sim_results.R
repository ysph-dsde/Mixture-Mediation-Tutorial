# Load required packages
library(zeallot)
library(ggplot2)
library(bkmr)

# Global settings
n_dataset <- 100  # Number of simulated datasets per scenario
n_scenario <- 4   # Number of scenarios
sce <- c("1000l", "1000h", "2500l", "2500h")  # Scenario labels
true_signals <-
  c(1, 2, 3, 11, 12, 13, 21, 22, 23)  # True mediating exposures

#### SEMA --------

compute_tpr_fpr <- function(fdr_vec, true_signals, cutoff = 0.05) {
  p <- length(fdr_vec)
  selected <- which(fdr_vec <= cutoff)
  TP <- sum(selected %in% true_signals)
  FP <- sum(!(selected %in% true_signals))
  FN <- sum(!(true_signals %in% selected))
  TN <- p - length(true_signals) - FP
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  return(c(TPR = TPR, FPR = FPR))
}

TPR1 <- TPR2 <- FPR1 <- FPR2 <- matrix(0, n_dataset, n_scenario)
colnames(TPR1) <-
  colnames(TPR2) <- colnames(FPR1) <- colnames(FPR2) <- sce
IE1 <- IE2 <- list()

for (j in 1:n_scenario) {
  pb <- txtProgressBar(min = 0, max = n_dataset, style = 3)
  
  for (i in 1:n_dataset) {
    data <- readRDS(paste0(sce[j], "/sema", i, ".rds"))
    nie1 <- data[["Raw (Naive)"]][["Tidy Results"]]$`PNIE Table`
    nie2 <- data[["Raw"]][["Tidy Results"]]$`PNIE Table`
    
    # Calculate TRP and FPR
    c(TPR1[i, j], FPR1[i, j]) %<-% compute_tpr_fpr(nie1$BH.Adj.Pval, true_signals)
    c(TPR2[i, j], FPR2[i, j]) %<-% compute_tpr_fpr(nie2$BH.Adj.Pval, true_signals)
    
    # Save NIE for each exposure
    if (i == 1) {
      IE1[[j]] <- nie1[, 1]
      IE2[[j]] <- nie2[, 2]
    } else{
      IE1[[j]] <- rbind(IE1[[j]], nie1[, 1])
      IE2[[j]] <- rbind(IE2[[j]], nie2[, 1])
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
}

tpr_list <- list(
  "Unadjusted|n = 1000|Strong" = TPR1[, 2],
  "Unadjusted|n = 1000|Weak"  = TPR1[, 1],
  "Unadjusted|n = 2500|Strong" = TPR1[, 4],
  "Unadjusted|n = 2500|Weak"  = TPR1[, 3],
  "Adjusted|n = 1000|Strong"  = TPR2[, 2],
  "Adjusted|n = 1000|Weak"   = TPR2[, 1],
  "Adjusted|n = 2500|Strong"  = TPR2[, 4],
  "Adjusted|n = 2500|Weak"   = TPR2[, 4]
)

tpr_df <- do.call(rbind, lapply(names(tpr_list), function(name) {
  parts <- strsplit(name, "\\|")[[1]]
  vec <- tpr_list[[name]]
  q <- quantile(vec, c(0.025, 0.975))
  data.frame(
    Method = parts[1],
    SampleSize = parts[2],
    Scenario = paste0(parts[3], " indirect effect"),
    TPR_mean = mean(vec),
    TPR_lower = q[1],
    TPR_upper = q[2]
  )
}))

fpr_list <- list(
  "Unadjusted|n = 1000|Strong" = FPR1[, 2],
  "Unadjusted|n = 1000|Weak"  = FPR1[, 1],
  "Unadjusted|n = 2500|Strong" = FPR1[, 4],
  "Unadjusted|n = 2500|Weak"  = FPR1[, 3],
  "Adjusted|n = 1000|Strong"  = FPR2[, 2],
  "Adjusted|n = 1000|Weak"   = FPR2[, 1],
  "Adjusted|n = 2500|Strong"  = FPR2[, 4],
  "Adjusted|n = 2500|Weak"   = FPR2[, 4]
)

fpr_df <- do.call(rbind, lapply(names(tpr_list), function(name) {
  parts <- strsplit(name, "\\|")[[1]]
  vec <- fpr_list[[name]]
  q <- quantile(vec, c(0.025, 0.975))
  data.frame(
    Method = parts[1],
    SampleSize = parts[2],
    Scenario = paste0(parts[3], " indirect effect"),
    FPR_mean = mean(vec),
    FPR_lower = q[1],
    FPR_upper = q[2]
  )
}))

saveRDS(list(
  TPR = tpr_df,
  FPR = fpr_df,
  IE1 = IE1,
  IE2 = IE2
), file = "sema.rds")


#### PCA ----

PCA_IE <- vector("list", n_scenario)
for (j in 1:n_scenario) {
  PCA_IE[[j]] <- vector("list", n_dataset)
  pb <- txtProgressBar(min = 0, max = n_dataset, style = 3)
  
  for (i in 1:n_dataset) {
    data <- readRDS(paste0(sce[j], "/pca", i, ".rds"))
    PCA_IE[[j]][[i]] <- data[["PCA"]][["Tidy Results"]]$`PNIE Table`$Estimate
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

saveRDS(list(IE = PCA_IE), file = "pca.rds")


#### ERS1 ----

ERS_IE <- matrix(0, ncol = n_scenario, nrow = n_dataset)
for (j in 1:n_scenario) {
  pb <- txtProgressBar(min = 0, max = n_dataset, style = 3)
  
  for (i in 1:n_dataset) {
    data <- readRDS(paste0(sce[j], "/ers1_", i, ".rds"))
    ERS_IE[i, j] <- data[["ERS_Case1"]][["CMA Summary Table"]][4, 1]
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

saveRDS(list(IE = ERS_IE), file = "ers1.rds")

#### BKMR ----

# Define PIP thresholds to evaluate
pip_thresholds <- c(0.5, 0.3, 0.1)
n_cutoff <- length(pip_thresholds)
true_exposure <- paste0("x", true_signals)

# List to hold exposure-level performance results across cutoffs
TPR_bkmr_cw_list <-
  FPR_bkmr_cw_list <- vector("list", n_cutoff)
TPR_bkmr_hier_list <-
  FPR_bkmr_hier_list <- vector("list", n_cutoff)

# Initialize storage matrices for each cutoff
for (k in 1:n_cutoff) {
  TPR_bkmr_cw_list[[k]] <- matrix(NA, n_dataset, n_scenario)
  FPR_bkmr_cw_list[[k]] <- matrix(NA, n_dataset, n_scenario)
  
  TPR_bkmr_hier_list[[k]] <- matrix(NA, n_dataset, n_scenario)
  FPR_bkmr_hier_list[[k]] <- matrix(NA, n_dataset, n_scenario)
}

# Function to compute exposure-level TPR/FPR from PIPs
compute_pip_tpr_fpr <-
  function(pips_outcome,
           true_exposure,
           threshold = 0.5,
           selection = "c") {
    p <- 30
    true_len <- length(true_exposure)
    pips <- subset(pips_outcome, variable != "mediator_M")
    selected <- if (selection == "c") {
      pips$variable[pips$PIP > threshold]
    } else {
      pips$variable[(pips$groupPIP * pips$condPIP) > threshold]
    }
    
    TP <- sum(selected %in% true_exposure)
    FP <- sum(!(selected %in% true_exposure))
    FN <- true_len - TP
    TN <- p - true_len - FP
    TPR <- if ((TP + FN) == 0)
      NA
    else
      TP / true_len
    FPR <- if ((FP + TN) == 0)
      NA
    else
      FP / (FP + TN)
    return(c(
      TPR = TPR,
      FPR = FPR
    ))
  }

# Loop over scenarios and datasets
for (j in 1:n_scenario) {
  pb <- txtProgressBar(min = 0, max = n_dataset, style = 3)
  
  for (i in 1:n_dataset) {
    data <- readRDS(paste0(sce[j], "/bkmr", i, ".rds"))
    sim <-
      readRDS(paste0("~/project/tutorial/sim_", sce[j], "/sim", i , ".rds"))
    
    ### Component-wise BKMR
    pips_outcome_cw <-
      ExtractPIPs(sim[["BKMR"]][["Outcome model"]])
    
    ### Hierarchical BKMR
    pips_outcome_hier <-
      ExtractPIPs(sim[["BKMR (Hierarchical)"]][["Outcome model"]])
    
    ### Loop over cutoff values
    for (k in 1:n_cutoff) {
      cutoff <- pip_thresholds[k]
      
      res_cw <-
        compute_pip_tpr_fpr(pips_outcome_cw,
                            true_exposure,
                            threshold = cutoff)
      TPR_bkmr_cw_list[[k]][i, j] <- res_cw["TPR"]
      FPR_bkmr_cw_list[[k]][i, j] <- res_cw["FPR"]
      
      res_hier <-
        compute_pip_tpr_fpr(pips_outcome_hier,
                            true_exposure,
                            threshold = cutoff)
      TPR_bkmr_hier_list[[k]][i, j] <- res_hier["TPR"]
      FPR_bkmr_hier_list[[k]][i, j] <- res_hier["FPR"]
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
}

# Assign names for easier lookup
names(TPR_bkmr_cw_list) <-
  names(FPR_bkmr_cw_list) <-
  names(TPR_bkmr_hier_list) <-
  names(FPR_bkmr_hier_list) <-
  paste0("cutoff_", pip_thresholds)

# Create result list across all cutoffs
bkmr_cutoff_sensitivity <- list()
for (k in 1:n_cutoff) {
  key <- paste0("cutoff_", pip_thresholds[k])
  bkmr_cutoff_sensitivity[[key]] <- list(
    bkmr = list(
      TPR = TPR_bkmr_cw_list[[k]],
      FPR = FPR_bkmr_cw_list[[k]]
    ),
    bkmr_hier = list(
      TPR = TPR_bkmr_hier_list[[k]],
      FPR = FPR_bkmr_hier_list[[k]]
    )
  )
}
for (cutoff_key in names(bkmr_cutoff_sensitivity)) {
  for (type in c("bkmr", "bkmr_hier")) {
    for (metric in c("TPR", "FPR")) {
      colnames(bkmr_cutoff_sensitivity[[cutoff_key]][[type]][[metric]]) <-
        sce
    }
  }
}

#### Save Results----

list_res <- list(
  bkmr = bkmr_cutoff_sensitivity
)

saveRDS(list_res, file = "bkmr.rds")
