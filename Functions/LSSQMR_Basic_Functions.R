# SCAD function
scad=function(x, lambda){
  if (abs(x) <= lambda){
    lambda*abs(x)
  } else if (abs(x) > lambda & abs(x) <= 3.7*lambda){
    (2*3.7*lambda*abs(x)-x^2-lambda^2)/(2*(3.7-1))
  } else {
    ((3.7+1)*lambda^2)/2
  }
}

# Gradient of SCAD function
scad_tilde_grad=function(x, lambda){
  if (x <= lambda){
    0
  } else if (x > lambda & x <= 3.7*lambda){
    lambda-lambda*(3.7-x/lambda)/(3.7-1)
  } else {
    lambda
  } 
}

# Check Loss function
check_loss=function(tau, r){
  r*(tau-as.numeric(r<0))
}

# Smoothed Quantile Check Loss function
SQ_check_loss=function(tau,Resp,Covar,intercept_beta,slope_beta){
  h_h=max(0.05,sqrt(tau*(1-tau))*(log(dim(Covar)[2])/dim(Covar)[1])^0.25)
  Res=Resp-matrix(1,dim(Resp)[1])%*%t(intercept_beta)-Covar%*%slope_beta
  return(exp(-(Res^2)/(2*h_h^2))*h_h/sqrt(2*pi)-Res*pnorm(-Res/h_h, mean=0, sd=1)+Res*tau)
}

# Steps in ADMM
## B Step in ADMM
Beta_gradient=function(sparsity_option, BetaBeta, rho, Y, X, lambda_1, lambda_2, 
                       M1_old, M2_old, M3_old, C_old, D_old, R_old, beta_0_old){
  svdBeta=svd(BetaBeta)
  svdBeta_scad_tilde_grad=sapply(svdBeta$d, scad_tilde_grad, lambda=lambda_1)
  A_old=(C_old+D_old+M1_old/rho+M2_old/rho)/2
  if (sparsity_option=="element-wise"){
    return(rho*t(X)%*%(M3_old/rho+R_old-Y+matrix(1,dim(Y)[1])%*%t(beta_0_old)+X%*%BetaBeta)+2*rho*(BetaBeta-A_old)
           -(svdBeta$u)%*%diag(svdBeta_scad_tilde_grad)%*%t(svdBeta$v)
           -sapply(abs(BetaBeta), scad_tilde_grad, lambda=lambda2)*BetaBeta/abs(BetaBeta)
    )
  } else{
    BetaEucnorm=sqrt(apply(BetaBeta^2,1,sum))
    BetaEucnorm_scad_tilde_grad=sapply(BetaEucnorm, scad_tilde_grad, lambda=lambda_2)
    return(rho*t(X)%*%(M3_old/rho+R_old-Y+matrix(1,dim(Y)[1])%*%t(beta_0_old)+X%*%BetaBeta)+2*rho*(BetaBeta-A_old)
           -(svdBeta$u)%*%diag(svdBeta_scad_tilde_grad)%*%t(svdBeta$v)
           -diag(BetaEucnorm_scad_tilde_grad/BetaEucnorm)%*%BetaBeta
    )
  }
}
Beta_update_GD=function(sparsity_option, Beta_old, rho, Y, X, lambda_1, lambda_2, 
                        M1_old, M2_old, M3_old, C_old, D_old, R_old, beta_0_old, Beta_maxit, Beta_tol, Beta_stepsize){
  B_iter=1
  BetaBeta_n=Beta_old
  BetaBeta_gradgrad=c(Inf,numeric(Beta_maxit))
  while ( (B_iter <= Beta_maxit) && (B_iter >= 1) && (BetaBeta_gradgrad[B_iter]>Beta_tol) ){
    BetaBeta_o=BetaBeta_n
    BetaBeta_direction=Beta_gradient(sparsity_option, BetaBeta=BetaBeta_o, rho, Y, X, lambda_1, lambda_2, 
                                     M1_old, M2_old, M3_old, C_old, D_old, R_old, beta_0_old)
    BetaBeta_n=BetaBeta_o-Beta_stepsize*BetaBeta_direction
    BetaBeta_gradgrad[B_iter+1]=sqrt(mean(BetaBeta_direction^2))
    B_iter=B_iter+1
  }
  BetaBeta_returned=list(BetaBeta_n, B_iter-1, BetaBeta_gradgrad[B_iter])
  names(BetaBeta_returned)=c("Beta_updated", "Beta_iterated", "Beta_gradgradgrad")
  return(BetaBeta_returned)
}

## Singular Value Thresholding (SVT) in ADMM
SVT_update=function(rho, lambda_1, Beta_new, M1_old){
  AA=svd(Beta_new-M1_old/rho)
  aa=AA$d-lambda_1/rho
  aa[aa<0]=0
  return(AA$u%*%diag(aa)%*%t(AA$v))
}

## Soft Thresholding (ST) in ADMM
ST_update=function(rho, lambda_2, Beta_new, M2_old){
  BB=Beta_new-M2_old/rho
  CC=abs(BB)-lambda_2/rho
  CC[CC<0]=0
  return(CC*sign(BB))
}

## Block Soft Thresholding (BST) in ADMM
BST_update=function(rho, lambda_2, Beta_new, M2_old){
  BBB=Beta_new-M2_old/rho
  normBBB=sqrt(apply(BBB^2,1,sum))
  CCC=normBBB-lambda_2/rho
  CCC[CCC<0]=0
  return(matrix(CCC/normBBB,dim(Beta_new)[1],dim(Beta_new)[2])*BBB)
}

## R Step in ADMM
R_equation_univ=function(rr, vecindex, tau, rho, h, Y, X, Beta_new, beta_0_new, M3_old){
  matindex=arrayInd(vecindex,dim(Y))
  r_i=matindex[1]
  r_k=matindex[2]
  return(tau/dim(Y)[1]-pnorm(-rr/h)/dim(Y)[1]+M3_old[r_i,r_k]+rho*(rr-Y[r_i,r_k]+beta_0_new[r_k]+X[r_i,]%*%Beta_new[,r_k]))
}
R_update_univ=function(Rbd, R_old ,tau, rho, h, Y, X, Beta_new, beta_0_new, M3_old){
  unisolverR=function(vecindex){
    uniroot(f=R_equation_univ,vecindex=vecindex, tau=tau, rho=rho, h=h, 
            Y=Y, X=X, Beta_new=Beta_new, beta_0_new=beta_0_new, M3_old=M3_old, interval = c(-Rbd,Rbd))$root
  }
  solver=sapply(1:prod(dim(Y)),unisolverR)
  return(matrix(solver, dim(Y)[1], dim(Y)[2]))
}
R_equation_univ_unnormalized=function(rr, vecindex, tau, rho, h, Y, X, Beta_new, beta_0_new, M3_old){
  matindex=arrayInd(vecindex,dim(Y))
  r_i=matindex[1]
  r_k=matindex[2]
  return(tau-pnorm(-rr/h)+M3_old[r_i,r_k]+rho*(rr-Y[r_i,r_k]+beta_0_new[r_k]+X[r_i,]%*%Beta_new[,r_k]))
}
R_update_univ_unnormalized=function(Rbd, R_old ,tau, rho, h, Y, X, Beta_new, beta_0_new, M3_old){
  unisolverR=function(vecindex){
    uniroot(f=R_equation_univ_unnormalized,vecindex=vecindex, tau=tau, rho=rho, h=h, 
            Y=Y, X=X, Beta_new=Beta_new, beta_0_new=beta_0_new, M3_old=M3_old, interval = c(-Rbd,Rbd))$root
  }
  solver=sapply(1:prod(dim(Y)),unisolverR)
  return(matrix(solver, dim(Y)[1], dim(Y)[2]))
}

# LSE Generator : for initial values in ADMM, only used when p <= n
LSE_generator=function(Simul_Data){
  X=Simul_Data$input
  Y=Simul_Data$response
  n=dim(Y)[1]
  XX=cbind(rep(1,n),X)
  LSE=solve(t(XX)%*%XX)%*%t(XX)%*%Y
  Beta_LSE=LSE[-1,]
  beta0_LSE=t(t(LSE[1,]))
  R_LSE=Y-XX%*%LSE
  outoutout=list(Beta_LSE, beta0_LSE, R_LSE)
  names(outoutout)=c("Beta:LSE", "beta0:LSE", "R:LSE")
  return(outoutout)
}