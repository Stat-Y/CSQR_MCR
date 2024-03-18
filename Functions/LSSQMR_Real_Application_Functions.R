# Real Data Analysis : ADMM
LSSQMR_ADMM_real=function(Sparsity_option, Penalty_option, Y, X, 
                          int_init, init, res_init, 
                          tau, lambda_1, lambda_2){
  LSSQMR_ADMM(Sparsity_option, Penalty_option, rho=1, Y, X, int_init, init, res_init, 
              tau, lambda_1, lambda_2, toll=0.01, maxiter=25, Beta_maxit=2000, Beta_tol=0.2, Beta_stepsize = 0.001, 
              R_normalize=FALSE, hweight=1,
              Beta_toll_weight=1, beta_0_toll_weight=1, 
              C_toll_weight=1, D_toll_weight=1, R_toll_weight=1,
              M1_toll_weight=1, M2_toll_weight=1, M3_toll_weight=1)
}

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

# Real Data Analysis : GIC calculation
LSSQMR_Real_GIC=function(Sparsity_option, Penalty_option, rho, Y, X, tau,
                         lambda_1_range, lambda_2_range, lambda_1_dense_level, lambda_2_dense_level, 
                         toll=0.01, maxiter=25, Beta_maxit=2000, Beta_tol=0.2, Beta_stepsize=0.001, rweight, sweight){
  GIC_record=expand.grid(seq(lambda_1_range[1], lambda_1_range[2], length=lambda_1_dense_level), 
                         seq(lambda_2_range[1], lambda_2_range[2], length=lambda_2_dense_level))
  colnames(GIC_record)=c("Lambda 1", "Lambda 2")
  GIC_record$GIC=0
  GIC_record$Penalty_option=rep(Penalty_option, dim(GIC_record)[1])
  GIC_record$Sparsity_option=rep(Sparsity_option,dim(GIC_record)[1])
  GIC_record$shat=0
  GIC_record$rhat=0
  
  hhhhh=max(0.05,sqrt(tau*(1-tau))*(log(dim(X)[2])/dim(X)[1])^0.25)
  Base_fit=LSSQMR_ADMM(Sparsity_option="row-wise", Penalty_option="Lasso", rho=1, Y=Y, X=X, 
                       int_init=matrix(0,dim(Y)[2]), init=matrix(0,dim(X)[2],dim(Y)[2]), 
                       res_init=matrix(0,dim(Y)[1],dim(Y)[2]), 
                       tau=tau, lambda_1=0, lambda_2=0.0001, toll=toll, maxiter=maxiter, 
                       Beta_maxit=Beta_maxit, Beta_tol=Beta_tol, Beta_stepsize=Beta_stepsize, R_normalize=FALSE)
  LSE_par=list(Base_fit$D, Base_fit$`beta 0`, 
               Y-matrix(1,dim(Y)[1])%*%t(Base_fit$`beta 0`)-X%*%Base_fit$D)
  names(LSE_par)=c("Beta:Init", "beta0:Init", "R:Init")
  
  cores=parallelly::availableCores()
  cl <- makeCluster(cores-1) #not to overload your computer
  clusterExport(cl, c("scad", "check_loss", "scad_tilde_grad", "Beta_gradient", "Beta_update_GD",
                      "SVT_update", "ST_update", "BST_update", 
                      "R_equation_univ", "R_update_univ",
                      "R_equation_univ_unnormalized", "R_update_univ_unnormalized", "GIC", 
                      "LSSQMR_ADMM", "progress_bar", "multiroot", "rankMatrix"))
  registerDoParallel(cl)
  chunk=foreach(i = 1:dim(GIC_record)[1], .combine= 'cbind') %dopar% {
    Temp_Fit=LSSQMR_ADMM(Sparsity_option=GIC_record[i,5], Penalty_option=GIC_record[i,4], rho, 
                         Y, X, int_init=LSE_par$`beta0:Init`, init=LSE_par$`Beta:Init`, res_init=LSE_par$`R:Init`, tau, 
                         lambda_1=GIC_record[i,1], lambda_2=GIC_record[i,2], 
                         toll, maxiter, Beta_maxit, Beta_tol, Beta_stepsize, R_normalize=FALSE)
    Temp_out=numeric(3)
    Temp_out[1]=GIC(Fitted_Model=Temp_Fit, Y, X, tauu=tau, 
                    lambda_1=GIC_record[i,1]/dim(Y)[1], lambda_2=GIC_record[i,2]/dim(Y)[1], rweight=rweight, sweight=sweight)
    Temp_out[2]=sum(sqrt(apply((Temp_Fit$D)^2,1,sum))!=0)
    Temp_out[3]=as.numeric(rankMatrix(Temp_Fit$C))
    Temp_out
  }
  stopCluster(cl)
  GIC_record[,-c(1,2,4,5)]=t(chunk)
  return(GIC_record)
}

# Real Data Analysis : Cross Validation Calculation
LSSQMR_Real_CVFit=function(Sparsity_option, Penalty_option, rho=1, Y, X, tau,
                           lambda_1, lambda_2, 
                           toll=0.01, maxiter=25, Beta_maxit=2000, Beta_tol=0.2, Beta_stepsize=0.001, CVgrp){
  ordering_vec=arrange(data.frame(1:dim(Y)[1],CVgrp),CVgrp)[,1]
  
  cores=length(unique(CVgrp))+1
  cl <- makeCluster(cores-1)
  clusterExport(cl, c("scad", "check_loss", "scad_tilde_grad", "Beta_gradient", "Beta_update_GD",
                      "SVT_update", "ST_update", "BST_update", 
                      "R_equation_univ", "R_update_univ",
                      "R_equation_univ_unnormalized", "R_update_univ_unnormalized", "GIC", 
                      "LSSQMR_ADMM", "progress_bar", "multiroot", "rankMatrix"))
  registerDoParallel(cl)
  
  chunk=foreach(i = 1:(cores-1), .combine= 'cbind') %dopar% {
    Ytr=Y[CVgrp!=i,]
    Xtr=X[CVgrp!=i,]
    Ytst=Y[CVgrp==i,]
    Xtst=X[CVgrp==i,]
    
    Base_fit_tr=LSSQMR_ADMM(Sparsity_option="row-wise", Penalty_option="Lasso", rho=rho, Y=Ytr, X=Xtr, 
                            int_init=matrix(0,dim(Ytr)[2]), init=matrix(0,dim(Xtr)[2],dim(Ytr)[2]), 
                            res_init=matrix(0,dim(Ytr)[1],dim(Ytr)[2]), 
                            tau=tau, lambda_1=0, lambda_2=0.0001, toll=toll, maxiter=maxiter, 
                            Beta_maxit=Beta_maxit, Beta_tol=Beta_tol, Beta_stepsize=Beta_stepsize, R_normalize=FALSE)
    LSE_par_tr=list(Base_fit_tr$D, Base_fit_tr$`beta 0`, 
                    Ytr-matrix(1,dim(Ytr)[1])%*%t(Base_fit_tr$`beta 0`)-Xtr%*%Base_fit_tr$D)
    names(LSE_par_tr)=c("Beta:Init", "beta0:Init", "R:Init")
    
    Tmpfit=LSSQMR_ADMM(Sparsity_option=Sparsity_option, Penalty_option=Penalty_option, 
                       rho=rho, Y=Ytr, X=Xtr, 
                       int_init=LSE_par_tr$`beta0:Init`, init=LSE_par_tr$`Beta:Init`, res_init=LSE_par_tr$`R:Init`,
                       tau=tau, lambda_1=lambda_1, lambda_2=lambda_2, toll=toll, maxiter=maxiter, 
                       Beta_maxit=Beta_maxit, Beta_tol=Beta_tol, Beta_stepsize=Beta_stepsize, R_normalize=FALSE)
    hhhhh=max(0.05,sqrt(tau*(1-tau))*(log(dim(Xtr)[2])/dim(Xtr)[1])^0.25)
    RRtst=Ytst-matrix(1,dim(Ytst)[1])%*%t(Tmpfit$`beta 0`)-Xtst%*%Tmpfit$D
    SQSQ=exp(-(RRtst^2)/(2*hhhhh^2))*hhhhh/sqrt(2*pi)-RRtst*pnorm(-RRtst/hhhhh, mean=0, sd=1)+RRtst*tau
    QQ=check_loss(tau, RRtst)
    rbind(t(QQ),t(SQSQ),as.numeric(rankMatrix(Tmpfit$C)),sum(apply((Tmpfit$D)^2,1,sum)!=0))
  }
  stopCluster(cl)
  
  Q_CV_error=t(chunk[1:dim(Y)[2],])
  SQ_CV_error=t(chunk[(dim(Y)[2]+1):(2*dim(Y)[2]),])
  
  Q_CV_error=as.matrix(arrange(data.frame(ordering_vec,Q_CV_error),ordering_vec)[,-1])
  SQ_CV_error=as.matrix(arrange(data.frame(ordering_vec,SQ_CV_error),ordering_vec)[,-1])
  
  colnames(Q_CV_error)=colnames(Y)
  colnames(SQ_CV_error)=colnames(Y)
  rownames(Q_CV_error)=rownames(Y)
  rownames(SQ_CV_error)=rownames(Y)
  
  meanr=mean(chunk[(dim(chunk)[1]-1),])
  means=mean(chunk[dim(chunk)[1],])
  
  outsouts=list(Q_CV_error, SQ_CV_error, meanr, means)
  names(outsouts)=c("Check Loss", "Smoothed Quantile Loss", "Mean Rank", "Mean Sparsity")
  return(outsouts)
}