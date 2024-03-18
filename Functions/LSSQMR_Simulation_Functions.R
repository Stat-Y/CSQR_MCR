# GIC function : calculates GIC value for the fitted model.
GIC=function(Fitted_Model, Y, X, tauu, lambda_1, lambda_2, rweight, sweight, hhweight=1){
  hh=hhweight*max(0.05,sqrt(tauu*(1-tauu))*(log(dim(X)[2])/dim(X)[1])^0.25)
  if ((lambda_1!=0) && (lambda_2!=0)){
    # Used for the proposed method
    RR=Y-matrix(1,dim(Y)[1])%*%t(Fitted_Model$`beta 0`)-X%*%Fitted_Model$D
    Q_loss=sum(exp(-(RR^2)/(2*hh^2))*hh/sqrt(2*pi)-RR*pnorm(-RR/hh, mean=0, sd=1)+RR*tauu)/dim(Y)[1]
    shat=sum(sqrt(apply((Fitted_Model$D)^2,1,mean))!=0)
    rhat=as.numeric(rankMatrix(Fitted_Model$C))
    GIC_value=log(Q_loss)+(rweight*rhat*(shat+dim(Y)[2])+sweight*shat*dim(Y)[2])*(log(log(dim(Y)[1]))*log(dim(X)[2]))/dim(Y)[1]
  } else if ((lambda_1!=0) && (lambda_2==0)){
    RR=Y-matrix(1,dim(Y)[1])%*%t(Fitted_Model$`beta 0`)-X%*%Fitted_Model$C
    Q_loss=sum(exp(-(RR^2)/(2*hh^2))*hh/sqrt(2*pi)-RR*pnorm(-RR/hh, mean=0, sd=1)+RR*tauu)/dim(Y)[1]
    shat=dim(Fitted_Model$D)[1]
    rhat=as.numeric(rankMatrix(Fitted_Model$C))
    GIC_value=log(Q_loss)+(rhat)*(as.numeric(rankMatrix(X))+dim(Y)[2])*rweight*(log(log(dim(Y)[1]))*log(dim(X)[2]))/dim(Y)[1]
  } else if ((lambda_1==0) && (lambda_2!=0)){
    RR=Y-matrix(1,dim(Y)[1])%*%t(Fitted_Model$`beta 0`)-X%*%Fitted_Model$D
    Q_loss=sum(exp(-(RR^2)/(2*hh^2))*hh/sqrt(2*pi)-RR*pnorm(-RR/hh, mean=0, sd=1)+RR*tauu)/dim(Y)[1]
    shat=sum(sqrt(apply((Fitted_Model$D)^2,1,mean))!=0)
    rhat=min(dim(Fitted_Model$C))
    GIC_value=log(Q_loss)+(shat*dim(Y)[2])*sweight*(log(log(dim(Y)[1]))*log(dim(X)[2]))/dim(Y)[1]
  } else{
    RR=Y-matrix(1,dim(Y)[1])%*%t(Fitted_Model$`beta 0`)-X%*%Fitted_Model$C
    Q_loss=sum(exp(-(RR^2)/(2*hh^2))*hh/sqrt(2*pi)-RR*pnorm(-RR/hh, mean=0, sd=1)+RR*tauu)/dim(Y)[1]
    shat=dim(Fitted_Model$D)[1]
    rhat=min(dim(Fitted_Model$C))
    GIC_value=log(Q_loss)
  } 
  return(GIC_value)
}

# For each regularization parameters in the given range, calculate GIC, shat, rhat, relative MSE, prediction error, U, C, O.
# Select regularization parameters with smallest GIC
LSSQMR_ADMM_GIC=function(n, Sparsity_option, Penalty_option, rho, Y, X, beta0, Beta, 
                         Y_test, X_test, int_init, init, res_init, tau,
                         lambda_1_range, lambda_2_range, lambda_1_dense_level, lambda_2_dense_level, 
                         toll, maxiter, selection_threshold, Beta_maxit, Beta_tol, Beta_stepsize, rweight, sweight, 
                         hweightt=1){
  GIC_record=expand.grid(seq(lambda_1_range[1], lambda_1_range[2], length=lambda_1_dense_level), 
                         seq(lambda_2_range[1], lambda_2_range[2], length=lambda_2_dense_level))
  colnames(GIC_record)=c("Lambda 1", "Lambda 2")
  GIC_record$GIC=0
  GIC_record$Penalty_option=rep(Penalty_option, dim(GIC_record)[1])
  GIC_record$Sparsity_option=rep(Sparsity_option,dim(GIC_record)[1])
  GIC_record$shat=0
  GIC_record$rhat=0
  GIC_record$UF=0
  GIC_record$CF=0
  GIC_record$OF=0
  GIC_record$Rel_MSE=0
  GIC_record$Pred_Q_Error=0
  
  
  cores=parallelly::availableCores()
  cl <- makeCluster(cores-1) 
  clusterExport(cl, c("scad", "check_loss", "scad_tilde_grad", "Beta_gradient", "Beta_update_GD",
                      "SVT_update", "ST_update", "BST_update", 
                      "R_equation_univ", "R_update_univ", "GIC", 
                      "LSSQMR_ADMM", "progress_bar", "multiroot", "rankMatrix"))
  registerDoParallel(cl)
  chunk=foreach(i = 1:dim(GIC_record)[1], .combine= 'cbind') %dopar% {
    Temp_Fit=LSSQMR_ADMM(Sparsity_option=GIC_record[i,5], Penalty_option=GIC_record[i,4], 
                         rho, Y, X, int_init, init, res_init, tau, 
                         lambda_1=GIC_record[i,1], lambda_2=GIC_record[i,2], toll, maxiter, Beta_maxit, Beta_tol, 
                         Beta_stepsize, hweight=hweightt)
    Temp_out=numeric(8)
    Temp_out[1]=GIC(Fitted_Model=Temp_Fit, Y, X, tauu=tau, 
                    lambda_1=GIC_record[i,1], lambda_2=GIC_record[i,2], 
                    rweight=rweight, sweight=sweight, hhweight=hweightt)
    Temp_out[2]=sum(sqrt(apply((Temp_Fit$D)^2,1,sum))!=0)
    Temp_out[3]=sum(svd(Temp_Fit$D)$d>selection_threshold)
    Temp_out[4]=as.numeric(sum((sqrt(apply(Temp_Fit$D^2,1,sum))!=0)[1:sum(sqrt(apply(Beta^2,1,sum))!=0)])
                           !=sum(sqrt(apply(Beta^2,1,sum))!=0))
    Temp_out[5]=as.numeric(sum((sqrt(apply(Beta^2,1,sum))!=0)!=(sqrt(apply(Temp_Fit$D^2,1,sum))!=0))==0)
    Temp_out[6]=as.numeric(sum((sqrt(apply(Temp_Fit$D^2,1,sum))!=0)[1:sum(sqrt(apply(Beta^2,1,sum))!=0)])
                           ==sum(sqrt(apply(Beta^2,1,sum))!=0))*(1-Temp_out[5])
    Temp_out[7]=sqrt(sum((Beta-Temp_Fit$D)^2)/(sum((Beta)^2)))
    RRRR=Y_test-matrix(1,n)%*%t(Temp_Fit$`beta 0`)-X_test%*%Temp_Fit$D
    Temp_out[8]=mean(check_loss(tau,RRRR))
    Temp_out
  }
  stopCluster(cl)
  GIC_record[,-c(1,2,4,5)]=t(chunk)
  return(GIC_record)
}

# Simulation single run function
## M1 : SCAD for low rank, Group SCAD for Group sparse 
## M2 : Lasso for low rank and Group Lasso for Group sparse
## M3 : Lasso for low rank and Lasso for sparse
Simulation_Single_Run=function(n, Beta, beta0, tau, X_AR_par, E_noise_level, E_error_type, rho, init_option, 
                               M1_lambda_1_range, M1_lambda_2_range, M1_lambda_1_dense_level, M1_lambda_2_dense_level, 
                               M2_lambda_1_range, M2_lambda_2_range, M2_lambda_1_dense_level, M2_lambda_2_dense_level, 
                               M3_lambda_1_range, M3_lambda_2_range, M3_lambda_1_dense_level, M3_lambda_2_dense_level, 
                               toll, maxiter, selection_threshold, Beta_maxit, Beta_tol, Beta_stepsize, rweight, sweight){
  Simul_Data=Simulation_Data(n, Beta, beta0, tau, X_AR_par, E_noise_level, E_error_type)
  Test_Data=Simulation_Data(n, Beta, beta0, tau, X_AR_par, E_noise_level, E_error_type)
  
  if (init_option=="LSE"){
    LSE_par=LSE_generator(Simul_Data)
  } else {
    Base_fit=LSSQMR_ADMM(Sparsity_option="row-wise", Penalty_option="Lasso", rho, Y=Simul_Data$response, X=Simul_Data$input, 
                         int_init=matrix(0,dim(beta0)[1]), init=matrix(0,dim(Beta)[1],dim(Beta)[2]), res_init=matrix(0,n,dim(Beta)[2]), 
                         tau, lambda_1=0, lambda_2=0.001, toll, maxiter, Beta_maxit, Beta_tol, Beta_stepsize)
    LSE_par=list(Base_fit$D, Base_fit$`beta 0`, 
                 Simul_Data$response-matrix(1,n)%*%t(Base_fit$`beta 0`)-Simul_Data$input%*%Base_fit$D)
    names(LSE_par)=c("Beta:LSE", "beta0:LSE", "R:LSE")
  }
  
  GIC_recordd_SCAD=LSSQMR_ADMM_GIC(n, Sparsity_option="row-wise", Penalty_option="SCAD", rho, 
                                   Y=Simul_Data$response, X=Simul_Data$input, beta0, Beta, 
                                   Y_test=Test_Data$response, X_test=Test_Data$input, 
                                   int_init=LSE_par$`beta0:LSE`, init=LSE_par$`Beta:LSE`, res_init=LSE_par$`R:LSE`, tau,
                                   M1_lambda_1_range, M1_lambda_2_range, M1_lambda_1_dense_level, M1_lambda_2_dense_level, 
                                   toll, maxiter, selection_threshold, Beta_maxit, Beta_tol, Beta_stepsize, rweight, sweight)
  
  GIC_recordd_Lasso=LSSQMR_ADMM_GIC(n, Sparsity_option="row-wise", Penalty_option="Lasso", rho, 
                                    Y=Simul_Data$response, X=Simul_Data$input, beta0, Beta, 
                                    Y_test=Test_Data$response, X_test=Test_Data$input, 
                                    int_init=LSE_par$`beta0:LSE`, init=LSE_par$`Beta:LSE`, res_init=LSE_par$`R:LSE`, tau,
                                    M2_lambda_1_range, M2_lambda_2_range, M2_lambda_1_dense_level, M2_lambda_2_dense_level, 
                                    toll, maxiter, selection_threshold, Beta_maxit, Beta_tol, Beta_stepsize, rweight, sweight)
  
  GIC_recordd_element_lasso=LSSQMR_ADMM_GIC(n, Sparsity_option="element-wise", Penalty_option="Lasso", rho, 
                                            Y=Simul_Data$response, X=Simul_Data$input, beta0, Beta, 
                                            Y_test=Test_Data$response, X_test=Test_Data$input, 
                                            int_init=LSE_par$`beta0:LSE`, init=LSE_par$`Beta:LSE`, res_init=LSE_par$`R:LSE`, tau,
                                            M3_lambda_1_range, M3_lambda_2_range, M3_lambda_1_dense_level, M3_lambda_2_dense_level, 
                                            toll, maxiter, selection_threshold, Beta_maxit, Beta_tol, Beta_stepsize, rweight, sweight)
  
  GIC_recordd=rbind(GIC_recordd_SCAD, GIC_recordd_Lasso, GIC_recordd_element_lasso)
  Perfornamce_table=rbind(arrange(subset(GIC_recordd, Penalty_option=="SCAD" & Sparsity_option=="row-wise"), GIC)[1,],
                          arrange(subset(GIC_recordd, Penalty_option=="Lasso" & Sparsity_option=="row-wise"), GIC)[1,],
                          arrange(subset(GIC_recordd, Penalty_option=="Lasso" & Sparsity_option=="element-wise"), GIC)[1,])
  outoutoutoutout=list(GIC_recordd, Perfornamce_table)
  names(outoutoutoutout)=c("GIC", "Performance")
  return(outoutoutoutout)
}

# Simulation multiple run function
Simulation_Multiple_Run=function(n, Beta, beta0, tau, X_AR_par, E_noise_level, E_error_type, rho, init_option, 
                                 M1_lambda_1_range, M1_lambda_2_range, M1_lambda_1_dense_level=11, M1_lambda_2_dense_level=11, 
                                 M2_lambda_1_range, M2_lambda_2_range, M2_lambda_1_dense_level=11, M2_lambda_2_dense_level=11, 
                                 M3_lambda_1_range, M3_lambda_2_range, M3_lambda_1_dense_level=11, M3_lambda_2_dense_level=11, 
                                 toll=0.01, maxiter=25, selection_threshold=1, Beta_maxit=500, Beta_tol=0.2, Beta_stepsize=0.0003,
                                 rweight=1, sweight=1, Replicate){
  GIC_record_array=array(0,c(M1_lambda_1_dense_level*M1_lambda_2_dense_level
                             +M2_lambda_1_dense_level*M2_lambda_2_dense_level
                             +M3_lambda_1_dense_level*M3_lambda_2_dense_level,10,Replicate))
  Performance_array=array(0,c(3,10,Replicate))
  Timerec=numeric(Replicate)
  pbb=progress_bar$new(total = Replicate)
  
  for(r in 1:Replicate){
    Current.Replicate.Rate<<- r/Replicate
    tempstarttime=Sys.time()
    tempsinglerunrun=Simulation_Single_Run(n, Beta, beta0, tau, X_AR_par, E_noise_level, E_error_type, rho, init_option, 
                                           M1_lambda_1_range, M1_lambda_2_range, M1_lambda_1_dense_level, M1_lambda_2_dense_level, 
                                           M2_lambda_1_range, M2_lambda_2_range, M2_lambda_1_dense_level, M2_lambda_2_dense_level, 
                                           M3_lambda_1_range, M3_lambda_2_range, M3_lambda_1_dense_level, M3_lambda_2_dense_level, 
                                           toll, maxiter, selection_threshold, Beta_maxit, Beta_tol, Beta_stepsize, rweight, sweight)
    tempendtime=Sys.time()
    print(paste("Result for Replicate ", r, sep=""))
    print(tempsinglerunrun[[2]])
    GIC_record_array[,,r]=as.matrix(tempsinglerunrun[[1]][,-c(4,5)])
    Performance_array[,,r]=as.matrix(tempsinglerunrun[[2]][,-c(4,5)])
    Timerec[r]=-as.numeric(difftime(time1 = tempstarttime, time2 = tempendtime, units = "secs"))
    Previous.Replicate.Time<<- Timerec[r]
    pbb$tick()
    r=r+1
  }
  
  Performance_mean=apply(Performance_array,c(1,2),mean)
  row.names(Performance_mean)=paste("M",1:3,sep="")
  colnames(Performance_mean)=colnames(tempsinglerunrun[[1]])[-c(4,5)]
  
  outoutoutout=list(GIC_record_array, Timerec, Performance_array, Performance_mean)
  names(outoutoutout)=c("GIC Record", "Time Record", "Performance", "Performance Mean")
  return(outoutoutout)
}