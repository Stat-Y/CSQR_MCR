# Sparsity_option : "element-wise" : estimated B is element-wise sparse matrix
# or "row-wise" : estimated B is row-wise sparse matrix.
# Penalty_option : SCAD or Lasso (Use SCAD for the proposed method).
# Y : n by q response matrix.
# X : n by p covariate matrix.
# int_init : q by 1 matrix, which is initial value for intercept beta0.
# init : p by q matrix, which is initial value for B. 
# res_init : n by q matrix, which is initial value for residual matrix.
# tau : quantile, which should be in (0,1).
# lambda_1 : regularization parameter for low-rank structure.
# lambda_2 : regularization paramter for sparsity structure.

## ADMM Function
LSSQMR_ADMM=function(Sparsity_option, Penalty_option, rho, Y, X, int_init, init, res_init, 
                     tau, lambda_1, lambda_2, toll, maxiter, Beta_maxit, Beta_tol, Beta_stepsize, 
                     R_normalize=TRUE, hweight=1,
                     Beta_toll_weight=1, beta_0_toll_weight=1, 
                     C_toll_weight=1, D_toll_weight=1, R_toll_weight=1,
                     M1_toll_weight=1, M2_toll_weight=1, M3_toll_weight=1){
  ### Restriction Bound for R update
  Rbd=5000*max(abs(res_init))+50000
  ### Iteration Record
  iter=1
  ### Parameters
  Beta_new=init
  beta_0_new=int_init
  C_new=init
  D_new=init
  R_new=res_init
  M1_new=init
  M2_new=init
  M3_new=res_init
  ### Bandwidth h
  h=hweight*max(0.05,sqrt(tau*(1-tau))*(log(dim(X)[2])/dim(X)[1])^0.25)
  ### Error Record
  Beta_error=c(Inf,numeric(maxiter))
  beta_0_error=c(Inf,numeric(maxiter))
  C_error=c(Inf,numeric(maxiter))
  D_error=c(Inf,numeric(maxiter))
  R_error=c(Inf,numeric(maxiter))
  M1_error=c(Inf,numeric(maxiter))
  M2_error=c(Inf,numeric(maxiter))
  M3_error=c(Inf,numeric(maxiter))
  ### Beta Step Record
  Beta_iter=c(Inf,numeric(maxiter))
  Beta_grad=c(Inf,numeric(maxiter))
  ### Time Record
  Beta_time=c(Inf,numeric(maxiter))
  beta_0_time=c(Inf,numeric(maxiter))
  C_time=c(Inf,numeric(maxiter))
  D_time=c(Inf,numeric(maxiter))
  R_time=c(Inf,numeric(maxiter))
  M1_time=c(Inf,numeric(maxiter))
  M2_time=c(Inf,numeric(maxiter))
  M3_time=c(Inf,numeric(maxiter))
  ### Stopping Criteria
  while ( (iter <= maxiter) && (iter >= 1) && ( (Beta_error[iter]>toll*Beta_toll_weight)||
                                                (beta_0_error[iter]>toll*beta_0_toll_weight)||
                                                (C_error[iter]>toll*C_toll_weight)||
                                                (D_error[iter]>toll*D_toll_weight)||
                                                (R_error[iter]>toll*R_toll_weight)||
                                                (M1_error[iter]>toll*M1_toll_weight)||
                                                (M2_error[iter]>toll*M2_toll_weight)||
                                                (M3_error[iter]>toll*M3_toll_weight) ) ){
    #### Define Old
    Beta_old=Beta_new
    beta_0_old=beta_0_new
    C_old=C_new
    D_old=D_new
    R_old=R_new
    M1_old=M1_new
    M2_old=M2_new
    M3_old=M3_new
    #### Update New
    Beta_time_start=Sys.time()
    if (Penalty_option=="SCAD"){
      updated_Beta_Beta=Beta_update_GD(sparsity_option=Sparsity_option, Beta_old, rho, Y, X, lambda_1, lambda_2, 
                                       M1_old, M2_old, M3_old, C_old, D_old, R_old, beta_0_old, 
                                       Beta_maxit, Beta_tol, Beta_stepsize)
      Beta_new=updated_Beta_Beta$Beta_updated
      Beta_iter[iter+1]=updated_Beta_Beta$Beta_iterated
      Beta_grad[iter+1]=updated_Beta_Beta$Beta_gradgradgrad
    } else{
      A_old=(C_old+D_old+M1_old/rho+M2_old/rho)/2
      Beta_new=solve(t(X)%*%X+2*diag(rep(1,dim(X)[2])))%*%(2*A_old-t(X)%*%(M3_old/rho+R_old-Y+matrix(1,dim(Y)[1])%*%t(beta_0_old)))
    }
    beta_0_time_start=Sys.time()
    beta_0_new=t(Y-R_old-M3_old/rho-X%*%Beta_new)%*%matrix(1,dim(Y)[1])/dim(Y)[1]
    C_time_start=Sys.time()
    C_new=SVT_update(rho, lambda_1, Beta_new, M1_old)
    D_time_start=Sys.time()
    if (Sparsity_option=="element-wise"){
      D_new=ST_update(rho, lambda_2, Beta_new, M2_old)
    } else{
      D_new=BST_update(rho, lambda_2, Beta_new, M2_old)
    }
    R_time_start=Sys.time()
    if (R_normalize==TRUE){
      R_new=R_update_univ(Rbd, R_old ,tau, rho, h, Y, X, Beta_new, beta_0_new, M3_old)
    } else{
      R_new=R_update_univ_unnormalized(Rbd, R_old ,tau, rho, h, Y, X, Beta_new, beta_0_new, M3_old)
    }
    M1_time_start=Sys.time()
    M1_new=M1_old+rho*(C_new-Beta_new)
    M2_time_start=Sys.time()
    M2_new=M2_old+rho*(D_new-Beta_new)
    M3_time_start=Sys.time()
    M3_new=M3_old+rho*(R_new-Y+matrix(1,dim(Y)[1])%*%t(beta_0_new)+X%*%Beta_new)
    M3_time_end=Sys.time()
    #### Update Error
    Beta_error[iter+1]=sqrt(mean((Beta_new-Beta_old)^2))
    beta_0_error[iter+1]=sqrt(mean((beta_0_new-beta_0_old)^2))
    C_error[iter+1]=sqrt(mean((C_new-C_old)^2))
    D_error[iter+1]=sqrt(mean((D_new-D_old)^2))
    R_error[iter+1]=sqrt(mean((R_new-R_old)^2))
    M1_error[iter+1]=sqrt(mean((M1_new-M1_old)^2))
    M2_error[iter+1]=sqrt(mean((M2_new-M2_old)^2))
    M3_error[iter+1]=sqrt(mean((M3_new-M3_old)^2))
    ### Update Time
    Beta_time[iter+1]=-as.numeric(difftime(time1 = Beta_time_start, time2 = beta_0_time_start, units = "secs"))
    beta_0_time[iter+1]=-as.numeric(difftime(time1 = beta_0_time_start, time2 = C_time_start, units = "secs"))
    C_time[iter+1]=-as.numeric(difftime(time1 = C_time_start, time2 = D_time_start, units = "secs"))
    D_time[iter+1]=-as.numeric(difftime(time1 = D_time_start, time2 = R_time_start, units = "secs"))
    R_time[iter+1]=-as.numeric(difftime(time1 = R_time_start, time2 = M1_time_start, units = "secs"))
    M1_time[iter+1]=-as.numeric(difftime(time1 = M1_time_start, time2 = M2_time_start, units = "secs"))
    M2_time[iter+1]=-as.numeric(difftime(time1 = M2_time_start, time2 = M3_time_start, units = "secs"))
    M3_time[iter+1]=-as.numeric(difftime(time1 = M3_time_start, time2 = M3_time_end, units = "secs"))
    #### Update Iteration
    iter=iter+1
  }
  ### Output
  outout=list(iter-1, Beta_new, beta_0_new, C_new, D_new, R_new, M1_new, M2_new, M3_new,  
              Beta_error[-1], beta_0_error[-1], C_error[-1], D_error[-1], R_error[-1], 
              M1_error[-1], M2_error[-1], M3_error[-1],
              Beta_time[-1], beta_0_time[-1], C_time[-1], D_time[-1], R_time[-1], 
              M1_time[-1], M2_time[-1], M3_time[-1], Beta_iter[-1], Beta_grad[-1])
  names(outout)=c("Iteration", "Beta", "beta 0", "C", "D", "R", "M1", "M2", "M3", 
                  "Beta Error", "beta 0 Error", "C Error", "D Error", "R Error", 
                  "M1 Error", "M2 Error", "M3 Error",
                  "Beta Time", "beta 0 Time", "C Time", "D Time", "R Time", 
                  "M1 Time", "M2 Time", "M3 Time", "Beta Iterations", "Beta Final Gradients")
  ### Return Output
  return(outout)
}