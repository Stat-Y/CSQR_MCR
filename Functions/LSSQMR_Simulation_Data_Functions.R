# Generating B
Simulation_Beta=function(p, q, r, s, Sparsity_option, min_signal, max_signal){
  if (Sparsity_option=="element-wise"){
    Betaa=matrix(0,p,q)
    B1=matrix(0,p,r)
    for(i in 1:r){
      B1[r*(0:(p/r-1))+i,i]=((-1)^rbinom(p/r,1,0.5))*runif(p/r,min_signal,max_signal)
    }
    for(j in 1:r){
      Betaa[,r*(((1:(q/r))-1))+j]=replicate(q/r,B1[,j])*t(matrix((-1)^rbinom(q/r,1,0.5),q/r,p))
    }
  } else{
    B1=matrix((-1)^rbinom(s*q,1,0.5)*runif(s*q,min_signal,max_signal),s,q)
    B1=svd(B1)$u%*%diag(c(svd(B1)$d[1:r],rep(0,min(s,q)-r)))%*%t(svd(B1)$v)
    Betaa=rbind(B1,matrix(0,p-s,q))
  }
  return(Betaa)
}

# Generating beta0
Simulation_beta0=function(q, min_signal, max_signal){
  beta0=matrix((-1)^rbinom(q,1,0.5)*runif(q,min_signal,max_signal),q)
  return(beta0)
}

# Simulation Data Generation
Simulation_Data=function(n, Beta, beta0, tau, X_AR_par, E_noise_level, E_error_type){
  p=dim(Beta)[1]
  q=dim(Beta)[2]
  X=mvrnorm(n,rep(0,p),X_AR_par^(as.matrix(dist(1:p))))
  if(E_error_type=="ALD"){
    E=matrix(rALD(n*q,0,E_noise_level,tau),n,q)
  } else if(E_error_type=="Gaussian"){
    E=matrix(rnorm(n*q,0,E_noise_level),n,q)
    E=E-matrix(qnorm(tau,0,E_noise_level),n,q)
  } else if(E_error_type=="Gaussian_Hetero"){
    noisevecc=runif(n*q,0.5*E_noise_level,2*E_noise_level)
    E=matrix(rnorm(n*q,0,noisevecc),n,q)
    E=E-qnorm(tau,0,matrix(noisevecc,n,q))
  } else if(E_error_type=="Gaussian_Inc"){
    Xsize=rank(apply(X^2,1,mean))
    E=t(mvrnorm(q,rep(0,n),2*E_noise_level*diag(Xsize)/n))
  } else {
    # Gaussian_Correlated
    E=mvrnorm(n,rep(0,q),(E_noise_level^2)*(0.7^(as.matrix(dist(1:q)))))
    E=E-matrix(qnorm(tau,0,E_noise_level),n,q)
  }
  Y=matrix(1,n)%*%t(beta0)+X%*%Beta+E
  # return : X, Y, tau
  outout=list(X, Y, tau)
  names(outout)=c("input", "response", "tau")
  return(outout)
}