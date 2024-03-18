# Seed
set.seed(777)

# Library
library(progress)
library(rootSolve)
library(MASS)
library(Matrix)
library(gplots)
library(ald)
library(parallel)
library(doParallel)
library(foreach)
library(parallelly)
library(plyr)

# Functions
setwd("~/LSSQMR_Git/Functions")
source("LSSQMR_Basic_Functions.R")
source("LSSQMR_ADMM_Functions.R")
source("LSSQMR_Simulation_Data_Functions.R")
source("LSSQMR_Simulation_Functions.R")

# Load B and beta0
setwd("~/LSSQMR_Git/Simulation_Parameters")
load("~/LSSQMR_Git/Simulation_Parameters/beta0.RData")
load("~/LSSQMR_Git/Simulation_Parameters/Beta_p800_q20_s40_r4.RData")

# Simulation Example
Mult_Run_p800_q20_s40_r4_Gaussian_075=Simulation_Multiple_Run(n=400, Beta, beta0, tau=0.75, X_AR_par=0.7, 
                                                              E_noise_level=0.2, E_error_type="Gaussian", 
                                                              rho=1, init_option="Lasso",
                                                              M1_lambda_1_range=c(2,14), M1_lambda_2_range=c(4,14), 
                                                              M1_lambda_1_dense_level=11, M1_lambda_2_dense_level=11, 
                                                              M2_lambda_1_range=c(2,14), M2_lambda_2_range=c(4,14),
                                                              M2_lambda_1_dense_level=11, M2_lambda_2_dense_level=11, 
                                                              M3_lambda_1_range=c(8,22), M3_lambda_2_range=c(4,14),
                                                              M3_lambda_1_dense_level=11, M3_lambda_2_dense_level=11, 
                                                              Replicate=100)
Mult_Run_p800_q20_s40_r4_Gaussian_075$"Performance Mean"