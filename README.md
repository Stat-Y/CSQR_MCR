# CSQR_MCR
R Code and Datasets for "Convolution smoothed quantile regression for multiple correlated responses".

## Main Functions

[LSSQMR_Pilot_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Pilot_Example.R) : Pilot example for fitting the proposed method using simulated dataset.

[LSSQMR_Simulation_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Simulation_Example.R) : Compares the proposed method with other methods through a simulation study, following the setup in the main paper's "Simulation Study" section.

[LSSQMR_Real_Application_CCLE_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Real_Application_CCLE_Example.R) : Applies the proposed method using the CCLE dataset, following the guidelines outlined in the "Real data applications" section of the main paper and the supplementary materials' "Details of real datasets" section.

[LSSQMR_Real_Application_GDSC_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Real_Application_GDSC_Example.R) : Applies the proposed method using the GDSC dataset, following the guidelines outlined in the "Real data applications" section of the main paper and the supplementary materials' "Details of real datasets" section.

[LSSQMR_Real_Application_Yeast_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Real_Application_Yeast_Example.R) : Applies the proposed method using the Yeast cell type dataset, following the guidelines outlined in the "Real data applications" section of the main paper and the supplementary materials' "Details of real datasets" section.

## Directory

### [Functions](https://github.com/Stat-Y/CSQR_MCR/tree/main/Functions)

- [LSSQMR_Basic_Functions.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/Functions/LSSQMR_Basic_Functions.R) : Contains fundamental functions utilized in various steps of the ADMM algorithm.

- [LSSQMR_ADMM_Functions.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/Functions/LSSQMR_ADMM_Functions.R) : Implements the ADMM algorithm for fitting the proposed method. 

- [LSSQMR_Simulation_Data_Functions.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/Functions/LSSQMR_Simulation_Data_Functions.R) : Provides functions for generating simulation datasets, which are employed in either [LSSQMR_Pilot_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Pilot_Example.R) or [LSSQMR_Simulation_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Simulation_Example.R).

- [LSSQMR_Simulation_Functions.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/Functions/LSSQMR_Simulation_Functions.R) : Functions for simulation that compares the proposed method with other methods, used in [LSSQMR_Simulation_Example.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Simulation_Example.R).

- [LSSQMR_Real_Application_Functions.R](https://github.com/Stat-Y/CSQR_MCR/blob/main/Functions/LSSQMR_Real_Application_Functions.R) : Functions used in real data applications for [CCLE](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Real_Application_CCLE_Example.R), [GDSC](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Real_Application_GDSC_Example.R), and [Yeast cell type](https://github.com/Stat-Y/CSQR_MCR/blob/main/LSSQMR_Real_Application_CCLE_Example.R) datasets.

### [Simulation_Parameters](https://github.com/Stat-Y/CSQR_MCR/tree/main/Simulation_Parameters)

- [beta0.RData](https://github.com/Stat-Y/CSQR_MCR/blob/main/Simulation_Parameters/beta0.RData) : 

- [Beta_p800_q20_s20_r4.RData](https://github.com/Stat-Y/CSQR_MCR/blob/main/Simulation_Parameters/Beta_p800_q20_s20_r4.RData) : 

- [Beta_p800_q20_s40_r4.RData](https://github.com/Stat-Y/CSQR_MCR/blob/main/Simulation_Parameters/Beta_p800_q20_s40_r4.RData) :  

### [Data](https://github.com/Stat-Y/CSQR_MCR/tree/main/Data)

- [CCLE.RData](https://github.com/Stat-Y/CSQR_MCR/blob/main/Data/CCLE.RData) : https://sites.broadinstitute.org/ccle/ gene expressions and IC50 see the main paper for details for the data

- [GDSC.RData](https://github.com/Stat-Y/CSQR_MCR/blob/main/Data/GDSC.RData) : https://www.cancerrxgene.org gene expressions and IC50 see the main paper for details for the data

- [Yeast.RData](https://github.com/Stat-Y/CSQR_MCR/blob/main/Data/Yeast.RData) : R package spls details in the main paper

## Contact
youngjin@vt.edu
