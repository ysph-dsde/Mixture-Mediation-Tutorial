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

## Notes
- Supporting functions are located in the Functions/ directory and sourced within scripts.
- For computationally intensive methods such as BKMR, runtime may be substantial depending on your computing environment.

## License
This project is licensed under the terms specified in the LICENSE file.
