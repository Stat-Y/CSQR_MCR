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
setwd("~/CSQR_MCR/Functions")
source("LSSQMR_Basic_Functions.R")
source("LSSQMR_ADMM_Functions.R")
source("LSSQMR_Real_Application_Functions.R")

# Data
load("~/CSQR_MCR/Data/GDSC.RData")

# GIC
## Skipping this step is acceptable since optimal tuning parameters from the same GIC selection were already loaded by above line.
## If the subsequent procedure updates the optimal tuning parameters, the outcome would be identical to the ones already loaded.
GIC_025=LSSQMR_Real_GIC(Sparsity_option="row-wise", Penalty_option="SCAD", rho=1, Y=Y, X=X, tau=0.25,
                        lambda_1_range=c(0,12), lambda_2_range=c(0,24), lambda_1_dense_level=36, lambda_2_dense_level=72, 
                        rweight=0.007, sweight=0.007)

GIC_050=LSSQMR_Real_GIC(Sparsity_option="row-wise", Penalty_option="SCAD", rho=1, Y=Y, X=X, tau=0.5,
                        lambda_1_range=c(0,12), lambda_2_range=c(0,24), lambda_1_dense_level=36, lambda_2_dense_level=72, 
                        rweight=0.007, sweight=0.007)

GIC_075=LSSQMR_Real_GIC(Sparsity_option="row-wise", Penalty_option="SCAD", rho=1, Y=Y, X=X, tau=0.75,
                        lambda_1_range=c(0,12), lambda_2_range=c(0,24), lambda_1_dense_level=36, lambda_2_dense_level=72, 
                        rweight=0.007, sweight=0.007)

## Optimal Regularization Parameters by GIC
optimal_025=GIC_025[order(GIC_025$GIC),][1,1:2]
names(optimal_025)=c("Lambda1", "Lambda2")

optimal_050=GIC_050[order(GIC_050$GIC),][1,1:2]
names(optimal_050)=c("Lambda1", "Lambda2")

optimal_075=GIC_075[order(GIC_075$GIC),][1,1:2]
names(optimal_075)=c("Lambda1", "Lambda2")

# Cross Validation
CV_025=LSSQMR_Real_CVFit(Sparsity_option="row-wise", Penalty_option="SCAD", Y=Y, X=X, tau=0.25,
                         lambda_1=optimal_025$Lambda1, lambda_2=optimal_025$Lambda2, CVgrp=CV10order)

CV_050=LSSQMR_Real_CVFit(Sparsity_option="row-wise", Penalty_option="SCAD", Y=Y, X=X, tau=0.5,
                         lambda_1=optimal_050$Lambda1, lambda_2=optimal_050$Lambda2, CVgrp=CV10order)

CV_075=LSSQMR_Real_CVFit(Sparsity_option="row-wise", Penalty_option="SCAD", Y=Y, X=X, tau=0.75,
                         lambda_1=optimal_075$Lambda1, lambda_2=optimal_075$Lambda2, CVgrp=CV10order)

# Final Fit
## Initial values for tau = 0.25, 0.5, 0.75
Base_fit_025=LSSQMR_ADMM_real(Sparsity_option="row-wise", Penalty_option="Lasso", Y=Y, X=X,
                              int_init=matrix(0,dim(Y)[2]), init=matrix(0,dim(X)[2],dim(Y)[2]), 
                              res_init=matrix(0,dim(Y)[1],dim(Y)[2]), 
                              tau=0.25, lambda_1=0, lambda_2=0.0001)
Init_par_025=list(Base_fit_025$D, Base_fit_025$`beta 0`,
                  Y-matrix(1,dim(Y)[1])%*%t(Base_fit_025$`beta 0`)-X%*%Base_fit_025$D)
names(Init_par_025)=c("Beta:Init", "beta0:Init", "R:Init")

Base_fit_050=LSSQMR_ADMM_real(Sparsity_option="row-wise", Penalty_option="Lasso", Y=Y, X=X,
                              int_init=matrix(0,dim(Y)[2]), init=matrix(0,dim(X)[2],dim(Y)[2]), 
                              res_init=matrix(0,dim(Y)[1],dim(Y)[2]), 
                              tau=0.5, lambda_1=0, lambda_2=0.0001)
Init_par_050=list(Base_fit_050$D, Base_fit_050$`beta 0`,
                  Y-matrix(1,dim(Y)[1])%*%t(Base_fit_050$`beta 0`)-X%*%Base_fit_050$D)
names(Init_par_050)=c("Beta:Init", "beta0:Init", "R:Init")

Base_fit_075=LSSQMR_ADMM_real(Sparsity_option="row-wise", Penalty_option="Lasso", Y=Y, X=X,
                              int_init=matrix(0,dim(Y)[2]), init=matrix(0,dim(X)[2],dim(Y)[2]), 
                              res_init=matrix(0,dim(Y)[1],dim(Y)[2]), 
                              tau=0.75, lambda_1=0, lambda_2=0.0001)
Init_par_075=list(Base_fit_075$D, Base_fit_075$`beta 0`,
                  Y-matrix(1,dim(Y)[1])%*%t(Base_fit_075$`beta 0`)-X%*%Base_fit_075$D)
names(Init_par_075)=c("Beta:Init", "beta0:Init", "R:Init")

## Fit the proposed methods for tau = 0.25, 0.5, 0.75
Fit_SCAD_025=LSSQMR_ADMM_real(Sparsity_option="row-wise", Penalty_option="SCAD", Y=Y, X=X, 
                              int_init=Init_par_025$`beta0:Init`, init=Init_par_025$`Beta:Init`, 
                              res_init=Init_par_025$`R:Init`, 
                              tau=0.25, lambda_1=optimal_025$Lambda1, lambda_2=optimal_025$Lambda2)

Fit_SCAD_050=LSSQMR_ADMM_real(Sparsity_option="row-wise", Penalty_option="SCAD", Y=Y, X=X, 
                              int_init=Init_par_050$`beta0:Init`, init=Init_par_050$`Beta:Init`, 
                              res_init=Init_par_050$`R:Init`, 
                              tau=0.5, lambda_1=optimal_050$Lambda1, lambda_2=optimal_050$Lambda2)

Fit_SCAD_075=LSSQMR_ADMM_real(Sparsity_option="row-wise", Penalty_option="SCAD", Y=Y, X=X, 
                              int_init=Init_par_075$`beta0:Init`, init=Init_par_075$`Beta:Init`, 
                              res_init=Init_par_075$`R:Init`, 
                              tau=0.75, lambda_1=optimal_075$Lambda1, lambda_2=optimal_075$Lambda2)
