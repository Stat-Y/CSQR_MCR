# Seed
set.seed(777)

# Library
library(rootSolve)
library(MASS)
library(Matrix)

# Functions
setwd("~/CSQR_MCR/Functions")
source("LSSQMR_Basic_Functions.R")
source("LSSQMR_ADMM_Functions.R")
source("LSSQMR_Simulation_Data_Functions.R")

# Load B and beta0
setwd("~/CSQR_MCR/Simulation_Parameters")
load("~/CSQR_MCR/Simulation_Parameters/beta0.RData")
load("~/CSQR_MCR/Simulation_Parameters/Beta_p800_q20_s40_r4.RData")

# Generate Simulation Data
Simul_Data=Simulation_Data(n=400, Beta, beta0, tau=0.75, 
                           X_AR_par=0.7, E_noise_level=0.5, 
                           E_error_type="Gaussian")

# Generate Initial Estimators
Base_fit=LSSQMR_ADMM(Sparsity_option="row-wise", Penalty_option="Lasso", rho=1, Y=Simul_Data$response, X=Simul_Data$input, 
                     int_init=matrix(0,dim(beta0)[1]), init=matrix(0,dim(Beta)[1],dim(Beta)[2]), 
                     res_init=matrix(0,dim(Simul_Data$response)[1],dim(Beta)[2]), 
                     tau=Simul_Data$tau, lambda_1=0, lambda_2=0.001, toll=0.01, maxiter=25, 
                     Beta_maxit=500, Beta_tol=0.2, Beta_stepsize=0.0003)
Init_par=list(Base_fit$D, Base_fit$`beta 0`, 
             Simul_Data$response-matrix(1,dim(Simul_Data$response)[1])%*%t(Base_fit$`beta 0`)-Simul_Data$input%*%Base_fit$D)
names(Init_par)=c("Beta:Init", "beta0:Init", "R:Init")

# Fit the proposed method
## Note that the regularization parameters are arbitrarily chosen in this pilot example.
## But for all simulation and real application, they are carefully chosen by GIC.
Fit_SCAD=LSSQMR_ADMM(Sparsity_option="row-wise", Penalty_option="SCAD", rho=1, Y=Simul_Data$response, X=Simul_Data$input, 
                     int_init=Init_par$`beta0:Init`, init=Init_par$`Beta:Init`, res_init=Init_par$`R:Init`, 
                     tau=Simul_Data$tau, lambda_1=10, lambda_2=8, toll=0.01, maxiter=25, 
                     Beta_maxit=500, Beta_tol=0.2, Beta_stepsize=0.0003)

# Compare true B and estimated B
heatmap.2(Beta[apply(Beta^2,1,sum)!=0,],Rowv=FALSE,Colv=FALSE,dendrogram="none",trace="none",col=bluered(1000), main="True B")
heatmap.2(Fit_SCAD$D[apply(Fit_SCAD$D^2,1,sum)!=0,],Rowv=FALSE,Colv=FALSE,dendrogram="none",trace="none",col=bluered(1000), main="Fitted B")
