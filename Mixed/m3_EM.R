
library(matrixStats)
library(EnvStats)

### Checking of validated parameters in binary comparison (i.e. at least one agreement/disagreement happen)
check1 <- function(m,N){
  min.par = 1/N
  if (length(which(m > (1-min.par))) > 0) {
    m[which(m > (1-min.par))] <- 1-min.par
  }
  if (length(which(m < min.par)) > 0) {
    m[which(m < min.par)] <- min.par
  }
  
  return(m)
  
}

##### Sub-function for computing starting points
start_hurdle_gammaK <- function(X,K,p0M,nB){
  lambdaM = nB/nrow(X)
  lambdaU = 1 - lambdaM
  
  p0 = matrix(0, nrow = 2, ncol =K)
  alpha = matrix(0, nrow =2, ncol = K)
  beta = matrix(0, nrow =2, ncol = K)
  for (k in 1:K){
    x = X[,k]
    
    x0 = x[x==0]
    x1 = sort(x[x>0])
    xM = x1[1:round((1-p0M)*nB)]
    xU = x1[-(1:round((1-p0M)*nB))]
    
    p0U = (length(x0)-nB*p0M)/(length(x)-nB)
    
    if (min(xM)<max(xM)){
      fitM = egamma(xM)
      paramM = fitM$parameters
    }else{
      paramM = c(1,1)
    }
    
    fitU = egamma(xU)
    paramU = fitU$parameters
    
    p0[,k] =  c(p0M, p0U)
    alpha[,k] = c(paramM[1], paramU[1])
    beta[,k] = c(paramM[2], paramU[2])
  }
  
  return(list(lambda = c(lambdaM, lambdaU), p0 = p0, alpha = alpha, beta = beta))
}

########## Sub-functions for E step and M step
# Comparison vector: (K1 binary value, K2 3-categorical, K3  continuous)

# For K1 binary comparison values
E1 <- function(m11, m10, u11, u10, N, c11, c10){
  m11.mat = matrix(rep(m11,N), nrow =N, byrow = TRUE)
  m10.mat = matrix(rep(m10,N), nrow =N, byrow = TRUE)
  
  u11.mat = matrix(rep(u11,N), nrow =N, byrow = TRUE)
  u10.mat = matrix(rep(u10,N), nrow =N, byrow = TRUE)
  
  pM = rowProds(m11.mat^(c11)*m10.mat^(c10))
  pU = rowProds(u11.mat^(c11)*u10.mat^(c10))
  
  return(cbind(pM, pU))
}

M1 <- function(g, K1, c11,c10){
  N = length(g)
 
  
  m11 = g%*%c11/sum(g)
  m11 = check1(m11,N) 
  
  m10 = 1-m11
  
  u11 = (1-g)%*%c11/sum(1-g)
  u11 = check1(u11,N) 
  
  u10 = 1 - u11
  
  return(list(m11=m11, m10 = m10, u11 = u11, u10 = u10))
}

# For K2 3-categorical comparison values


E2 <- function(m22, m21, m20, u22,u21, u20, N, c22, c21, c20){
  m22.mat = matrix(rep(m22,N), nrow =N, byrow = TRUE)
  m21.mat = matrix(rep(m21,N), nrow =N, byrow = TRUE)
  m20.mat = matrix(rep(m20,N), nrow =N, byrow = TRUE)
  
  u22.mat = matrix(rep(u22,N), nrow =N, byrow = TRUE)
  u21.mat = matrix(rep(u21,N), nrow =N, byrow = TRUE)
  u20.mat = matrix(rep(u20,N), nrow =N, byrow = TRUE)
  
  
  p2M = rowProds(m22.mat^(c22)*m21.mat^(c21)*m20.mat^(c20))
  p2U = rowProds(u22.mat^(c22)*u21.mat^(c21)*u20.mat^(c20))
  
  return(cbind(p2M, p2U))
}

M2 <- function(g,K2,c22,c21,c20){
  N = length(g)
  g2.mat = matrix(rep(g,K2),ncol = K2)
  
  m22 = colSums(g2.mat*c22)/sum(g)
  m22 = check1(m22,N)
  m21 = colSums(g2.mat*c21)/sum(g)
  m21 = check1(m21,N)
  m20 = 1-m22-m21 #check 0,1
  
  if (length(which(m20 < 1/N)) > 0){
    m22[which(m20 < 1/N)] = m22[which(m20 < 1/N)] - 1/(2*N)
    m21[which(m20 < 1/N)] = m21[which(m20 < 1/N)] - 1/(2*N)
    m20[which(m20 < 1/N)] = m20[which(m20 < 1/N)] + 1/N 
  }
  
  u22 = colSums((1-g2.mat)*c22)/sum(1-g)
  u22 = check1(u22,N)
  u21 = colSums((1-g2.mat)*c21)/sum(1-g)
  u21 = check1(u21,N)
  u20 = 1-u22-u21 #check 0,1
  
  if (length(which(u20 < 1/N)) > 0){
    u22[which(u20 < 1/N)] = u22[which(u20 < 1/N)] - 1/(2*N)
    u21[which(u20 < 1/N)] = u21[which(u20 < 1/N)] - 1/(2*N)
    u20[which(u20 < 1/N)] = u20[which(u20 < 1/N)] + 1/N 
  }
  return(list(m22 = m22,m21 = m21,m20=m20,u22=u22,u21=u21,u20=u20))
}

# For K3 continuous comparison values
E3 <- function(p0, alpha, beta, X){
  #p0, alpha, beta: matrix( 2 rows, K3 columns) contains parameters for M (first row) and U (second row)
  K3 = ncol(X)
  
  dhgamma= function(x,p0,shape,scale){
    
    ind0 = which(x==0)
    ind1 = which(x>0)
    
    if (length(ind0)>0 & p0 == 0){
      cat("WARNING! p0 should > 0!", "\n")}
    
    y = rep(0, length(x))
    
    y[ind0] = p0
    y[ind1] = (1-p0)* dgamma(x[ind1],shape = shape,scale = scale)
    
    return(y)
  }
  
  #lambda = c(pM, pU)
  dens <- function(X, K, p0, alpha, beta) {
    ncomp <- 2
    
    L = 1
    for (k in 1:K){
      temp = sapply(1:ncomp, function(j) dhgamma(X[,k], p0 = p0[j,k], shape = alpha[j,k], 
                                                 scale = beta[j,k]))
      L = L*temp
    }
    return(L)
  }
  
  dens3 <- dens(X, K3,  p0 , alpha , beta )
  
  return(dens3)
}



M3 <- function(g, X, p0, alpha, beta){
  
  
  fn.alpha <- function(alpha, beta, z, x) -log(beta) + sum(z * log(x))/sum(z) - digamma(alpha)
  fn.alpha.2 <- function(alpha, beta, z, x) (-log(beta) + sum(z * log(x))/sum(z) - digamma(alpha))^2
  fn.beta <- function(z, x, alpha) sum(z * x)/(sum(z) * alpha)
  
  K3 = ncol(X)
  N = nrow(X)
  z = cbind(g, 1-g)
  
  scale.mle = beta
  shape.mle = alpha
  old.scale.mle = scale.mle
  old.shape.mle = shape.mle
  p0.mle = p0
  
  ncomp = 2
  
  ##### M step
  maximize <- function(k){
    ind0 = which(X[,k]==0)
    ind1 = which(X[,k]>0)
    
    z0 <- z[ind0,]
    p0.mle[,k] <- colSums(z0)/colSums(z)
    p0.mle[p0.mle < 1/N] = 1/N
    
    temp <- try(sapply(1:ncomp, function(i) uniroot(fn.alpha, interval = c(1e-06, 10000), beta = old.scale.mle[i,k], 
                                                    z = z[ind1, i], x = X[ind1,k])$root), silent = TRUE)
    if (class(temp) == "try-error"){
      temp <- sapply(1:ncomp, function(i) nlminb(old.shape.mle[i], 
                                                 fn.alpha.2, lower = 0, beta = scale.mle[i,k], 
                                                 z = z[ind1, i], x = X[ind1,k])$par)}
    shape.mle[,k] = temp
    
    scale.mle[,k] <- sapply(1:ncomp, function(i) fn.beta(z = z[ind1,i], x = X[ind1,k], alpha = shape.mle[i,k]))
    
    return(c(shape.mle[,k] , scale.mle[,k] , p0.mle[,k]))
  }
  
  max_temp <- sapply(1:K3, FUN = maximize)
  shape.mle <- max_temp[1:2,]
  scale.mle <- max_temp[3:4,]
  p0.mle <- max_temp[5:6,]

  return(list(p0 = p0.mle, alpha = shape.mle, beta = scale.mle))
}

################ Implementation of ECM algorithm for the proposed model
EM_mixed = function(comp_mat, datA, datB, K1, K2, K3, tol = 1e-6, maxits = 500){

  # initializations
  nA = nrow(datA)
  nB = nrow(datB)
  
  N <- nrow(comp_mat)
  
  c11 <- matrix(as.numeric(comp_mat[,1:K1]==1), ncol = K1)
  c10 <- matrix(as.numeric(comp_mat[,1:K1]==0), ncol = K1)
  
  c22 <- matrix(as.numeric(comp_mat[,(K1+1):(K1+K2)]==2), ncol = K2)
  c21 <- matrix(as.numeric(comp_mat[,(K1+1):(K1+K2)]==1), ncol = K2)
  c20 <- matrix(as.numeric(comp_mat[,(K1+1):(K1+K2)]==0), ncol = K2)
  
  c3 <- matrix(comp_mat[,(K1+K2+1):(K1+K2+K3)], ncol = K3)
  
  # Starting points
  p = 1/nA
  
  m11 = rep(0.9,K1)  #length K1
  m10 = 1-m11
  
  u11 = rep(0.1,K1) #length K1
  u10 = 1-u11
  
  m22 = rep(0.05,K2) #length K2
  m21 = rep(0.85,K2) #length K2
  m20 = 1 - m22 - m21
  
  u22 = rep(0.1,K2) #length K2
  u21 = rep(0.5,K2) #length K2
  u20 = 1 - u22 - u21
  
  
  start = start_hurdle_gammaK(c3, K3, p0M = 0.8, nB = nB)
  p0    = start$p0
  alpha = start$alpha
  beta = start$beta
  
  
  loglikelihood <- function(p,N,m11, m10, u11, u10, c11, c10,
                            m22, m21, m20, u22,u21, u20, c22, c21, c20,
                            p0, alpha, beta, c3) {
    p1 = E1(m11, m10, u11, u10, N, c11, c10)
    p2 = E2(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    p3 = E3(p0, alpha, beta, c3)
    
    pM = p1[,1]*p2[,1]*p3[,1]
    pU = p1[,2]*p2[,2]*p3[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #log-likelihood not complete log-likelihood
    return(ll)
  }
  
  

  iter <- 0
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m11, m10, u11, u10, c11, c10,
                         m22, m21, m20, u22,u21, u20, c22, c21, c20,
                         p0, alpha, beta, c3)
  ll= old.ll
  converge = TRUE
  while (diff > tol && iter < maxits){ 
    
    p.old = p
    
    m11.old = m11 #length K1
    m10.old = m10
    
    u11.old = u11 #length K1
    u10.old = u10
    
    m22.old = m22 #length K2
    m21.old = m21 #length K2
    m20.old = m20
    
    u22.old = u22 #length K2
    u21.old = u21 #length K2
    u20.old = u20
    
    p0.old = p0 #length K3
    alpha.old = alpha 
    beta.old = beta
    
    
    ############################ E step
    # Compute partly probabilities
    p1 = E1(m11, m10, u11, u10, N, c11, c10)
    p2 = E2(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    p3 = E3(p0, alpha, beta, c3)
    
    pM = p1[,1]*p2[,1]*p3[,1]
    pU = p1[,2]*p2[,2]*p3[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################ M step
    p = sum(g)/N
    
    ###  Binary comparison values
    max1 = M1(g,K1,c11,c10)
    m11 = max1$m11
    m10 = max1$m10
    u11 = max1$u11
    u10 = max1$u10
    
    ### 3-categorical comparison values
    max2 = M2(g,K2,c22,c21,c20)
    m22 = max2$m22
    m21 = max2$m21
    m20 = max2$m20
    u22 = max2$u22
    u21 = max2$u21
    u20 = max2$u20
    
    ### Continuous comparison values
    max3 = M3(g, c3, p0, alpha, beta)
    p0 = max3$p0
    alpha = max3$alpha
    beta = max3$beta

    
    ### Stopping 
    new.ll = loglikelihood(p,N,m11, m10, u11, u10, c11, c10,
                           m22, m21, m20, u22,u21, u20, c22, c21, c20,
                           p0, alpha, beta,  c3)
    diff <- abs(new.ll - old.ll)
    
    old.ll <- new.ll
    ll <- c(ll, old.ll)
    
    iter = iter + 1
    
    if (iter == maxits) {
      converge = FALSE
      cat("WARNING! NOT CONVERGENT!", "\n")}
     
  }
  
  return(list(g=g, loglikelihood = ll, converge = converge))
}


######################### Implementation of EM algorithm for standard model
EM_binary <- function(comp_mat, datA, datB, K, tol = 1e-5, maxits = 500){
  nA = nrow(datA)
  nB = nrow(datB)
  # Starting point
  
  N <- nrow(comp_mat)
  p = nB/N
  m = rep(0.9, K)
  u = rep(0.1, K)
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  
  g = rep(0,N) # probability of being in Match  for each pair l
  it = 0
  converge = FALSE
  
  
  while ((!converge) & (it < maxits)){ 
    
    pOld = p
    mOld = m
    uOld = u
    ### E
    # Compute expectation
    
    m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
    u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
    
    
    probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
    probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
    
    
    g = p*probM/(p*probM+(1-p)*probU)
    
    ### Maximization
    g_mat = matrix(rep(g,K),ncol = K)
    
    p = sum(g)/N # p fix
    m = colSums(g_mat*comp_mat)/sum(g)
    u = colSums((1-g_mat)*comp_mat)/sum(1-g)
    
    if (length(which(m > 0.99999)) > 0) {
      m[which(m > 0.99999)] <- 0.99999
    }
    if (length(which(m < 1e-05)) > 0) {
      m[which(m < 1e-05)] <- 1e-05
    }
    if (length(which(u > 0.99999)) > 0) {
      u[which(u > 0.99999)] <- 0.99999
    }
    if (length(which(u < 1e-05)) > 0) {
      u[which(u < 1e-05)] <- 1e-05
    }
    
    
    parmlistold = c(mOld, uOld)
    parmlistcurrent = c(m,u)
    
    it = it + 1
    
    converge = (max(abs(parmlistold - parmlistcurrent)/abs(parmlistold)) <= tol)
  }
  m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
  u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
  
  
  probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
  probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
  
  
  g = p*probM/(p*probM+(1-p)*probU)
  
  
  return(list(g=g, p=p, m=m, u = u, it = it, converge = converge))
}



