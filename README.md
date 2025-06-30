# A Tutorial on Conducting Mediation Analysis with Exposure Mixtures
This repository provides code and functions for conducting mediation analysis in the context of exposure mixtures, focusing on both simulated data and real-world applications from the PROTECT cohort study. The methods implemented include single exposure mediation analysis (SEMA), principal component analysis (PCA)-based mediation, environmental risk score (ERS)-based mediation, and Bayesian kernel machine regression causal mediation analysis (BKMR-CMA).

## Overview

In environmental health research, individuals are often exposed to complex mixtures of chemicals rather than isolated exposures. This project illustrates practical strategies for mediation analysis under such settings, addressing challenges such as multicollinearity, sparsity, and nonlinearity among exposures.

The codebase supports:
- Simulation of exposure mixtures and outcomes
- Application of mediation methods on simulated datasets
- Application of mediation methods on real-world PROTECT cohort data
- Comparison of different mediation analysis strategies

## Data Availability

Due to the size of the datasets, they are hosted separately.  
You can download the datasets here:  
🔗 [Download Data from Google Drive](https://drive.google.com/drive/folders/1ok6uO5dyDF9X2qel9ZBfBH9omQtXsEw_?usp=sharing)

After downloading, place the data files in the appropriate directories as indicated in the scripts.

## Repository Structure

- **Functions/**
  - **Functions_DataGen.R**: Functions for simulation data generation
  - **Functions_ERS.R**: Functions for ERS-based mediation analysis
  - **Functions_IndTesting.R**: Functions for single exposure mediation analysis (SEMA)
  - **Functions_PCA.R**: Functions for PCA-based mediation analysis
- **Evaluation/**
  - **data_sim.R**: Master script for simulating 1 dataset with 100,000 observation under each degree of effect
  - **estimate_true_ers.R**: Estimates "true" ERS effects using a large simulated dataset (n = 100,000)
  - **estimate_true_pca.R**: Estimates "true" PCA effects using a large simulated dataset (n = 100,000)
  - **eva_bkmr.R**: Computes variable selection accuracy for BKMR-CMA
  - **eva_ers.R**: Computes relative bias for ERS
  - **eva_pca.R**: Computes relative bias for PCA
  - **eva_sema.R**: Computes relative bias and variable selection accuracy for SEMA
  - **functions_ERS.R**: Functions for ERS-based mediation analysis
  - **plot_relative_bias.R**: Generates comparative plots of percent relative bias across methods
  - **plot_selection_metrics.R**: Generates comparative plots of true/false positive rates across methods
  - **summarize_sim_results.R**: Extracts and stores relevant quantities from each simulated dataset
- **BKMR_Protect.R**: BKMR-CMA analysis on PROTECT dataset
- **BKMR_sim.R**: BKMR-CMA analysis on simulated dataset
- **Data_Clean_PROTECT.R**: Data cleaning and preprocessing for PROTECT dataset
- **ERS_Protect.R**: ERS-based mediation analysis on PROTECT dataset
- **ERS_sim_Case1.R**: ERS-based mediation analysis on simulated data (case 1)
- **ERS_sim_Case2.R**: ERS-based mediation analysis on simulated data (case 2)
- **SEMA_Protect.R**: Single exposure mediation analysis (SEMA) on PROTECT dataset
- **SEMA_sim.R**: Single exposure mediation analysis (SEMA) on simulated dataset
- **PCA_Protect.R**: PCA-based mediation analysis on PROTECT dataset
- **PCA_sim.R**: PCA-based mediation analysis on simulated dataset
- **Simulation_Data_Generation.R**: Script for generating simulated datasets
- **LICENSE**: License information
- **README.md**: Project documentation (this file)

## Notes
- Supporting functions are located in the Functions/ directory and sourced within scripts.
- For computationally intensive methods such as BKMR, runtime may be substantial depending on your computing environment.

## License
This project is licensed under the terms specified in the LICENSE file.
