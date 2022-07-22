rhgamma <- function(n,p0,shape,scale){
  
  X = rbinom(n, 1, 1-p0)
  gammaH = X* rgamma(n=n, shape = shape, scale = scale)
  
  return(gammaH)
}


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

library(EnvStats)

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

#### EM algorithm for mixture of hurdle gamma distribution when K=1
EM_hurdle_gamma = function(x, start, tol = 1e-6, maxits = 500){
  
  ind0 = which(x==0)
  ind1 = which(x>0)
  
  lambda = start$lambda
  p0    = start$p0
  alpha = start$alpha
  beta = start$beta
  
  fn.alpha <- function(alpha, beta, z, x) -log(beta) + sum(z * log(x))/sum(z) - digamma(alpha)
  fn.alpha.2 <- function(alpha, beta, z, x) (-log(beta) + sum(z * log(x))/sum(z) - digamma(alpha))^2
  fn.beta <- function(z, x, alpha) sum(z * x)/(sum(z) * alpha)
  
  
  lambda.mle <- lambda
  p0.mle <- p0
  shape.mle <- alpha
  scale.mle <- beta
  
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
  
  dens <- function(x, lambda, p0, alpha, beta) {
    ncomp <- length(lambda)
    temp = sapply(1:ncomp, function(j) dhgamma(x, p0 = p0[j], shape = alpha[j], 
                                               scale = beta[j]))
    temp = t(lambda * t(temp))
    temp
  }
  
  iter <- 0
  diff <- tol + 1
  dens1 <- dens(x = x,  lambda = lambda.mle, p0 = p0.mle, alpha = shape.mle, beta = scale.mle)
  
  old.obs.ll <- sum(log(apply(dens1, 1, sum)))
  ll <- old.obs.ll
  ncomp <- length(lambda)
  
  while (diff > tol && iter < maxits){
    old.p0.mle <- p0.mle
    old.shape.mle <- shape.mle
    old.scale.mle <- scale.mle
    old.lambda.mle <- lambda.mle
    
    z <- dens1/apply(dens1, 1, sum)
    
    z0 <- z[ind0,]
    p0.mle <- colSums(z0)/colSums(z)
    
    shape.mle <- try(sapply(1:ncomp, function(i) uniroot(fn.alpha, interval = c(1e-06, 10000), beta = old.scale.mle[i], 
                                                         z = z[ind1, i], x = x[ind1])$root), silent = TRUE)
    if (class(shape.mle) == "try-error"){
      shape.mle <- sapply(1:ncomp, function(i) nlminb(old.shape.mle[i], 
                                                      fn.alpha.2, lower = 0, beta = scale.mle[i], 
                                                      z = z[ind1, i], x = x[ind1])$par)}
    
    scale.mle <- sapply(1:ncomp, function(i) fn.beta(z = z[ind1,i], x = x[ind1], alpha = shape.mle[i]))
    
    lambda.mle <- apply(z, 2, mean)
    
    dens1 <- dens(x = x, p0 = p0.mle, lambda = lambda.mle, alpha = shape.mle, beta = scale.mle)
    new.obs.ll <- sum(log(apply(dens1, 1, sum)))
    diff <- new.obs.ll - old.obs.ll
    
    old.obs.ll <- new.obs.ll
    ll <- c(ll, old.obs.ll)
    iter = iter + 1
    
    if (iter == maxits){
      cat("WARNING! NOT CONVERGENT!", "\n")}
  }
  return(list(g=z[,1], lambda = lambda.mle, p0 = p0.mle, shape = shape.mle, scale = scale.mle, iter = iter))
}

#### EM algorithm for mixture of hurdle gamma distribution when K>1

EM_hurdle_gammaK = function(X, K, nB, tol = 1e-6, maxits = 500){
  
  N =  nrow(X)
  start = start_hurdle_gammaK(X, K, p0M = 0.8, nB = nB)
  lambda = start$lambda
  p0    = start$p0
  alpha = start$alpha
  beta = start$beta
  
  fn.alpha <- function(alpha, beta, z, x) -log(beta) + sum(z * log(x))/sum(z) - digamma(alpha)
  fn.alpha.2 <- function(alpha, beta, z, x) (-log(beta) + sum(z * log(x))/sum(z) - digamma(alpha))^2
  fn.beta <- function(z, x, alpha) sum(z * x)/(sum(z) * alpha)
  
  
  lambda.mle <- lambda
  p0.mle <- p0
  p0.mle[p0.mle < 1/N] = 1/N
  shape.mle <- alpha
  scale.mle <- beta
  
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
  
  dens <- function(X, K, lambda, p0, alpha, beta) {
    ncomp <- length(lambda)
    
    L = 1
    for (k in 1:K){
      temp = sapply(1:ncomp, function(j) dhgamma(X[,k], p0 = p0[j,k], shape = alpha[j,k], 
                                                 scale = beta[j,k]))
      L = L*temp
    }
    
    density = t(lambda * t(L))
    density
  }
  
  iter <- 0
  diff <- tol + 1
  dens1 <- dens(X=X, K=K,  lambda = lambda.mle, p0 = p0.mle, alpha = shape.mle, beta = scale.mle)
  
  old.obs.ll <- sum(log(apply(dens1, 1, sum)))
  ll <- old.obs.ll
  converge = FALSE
  ncomp <- length(lambda)
  
  while (!converge && iter < maxits){
    old.p0.mle <- p0.mle
    old.shape.mle <- shape.mle
    old.scale.mle <- scale.mle
    old.lambda.mle <- lambda.mle
    
    z <- dens1/apply(dens1, 1, sum)
    
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
    
    max_temp <- sapply(1:K, FUN = maximize)
    shape.mle <- max_temp[1:2,]
    scale.mle <- max_temp[3:4,]
    p0.mle <- max_temp[5:6,]
    
    lambda.mle <- apply(z, 2, mean)
    
    dens1 <- dens(X=X, K=K,  lambda = lambda.mle, p0 = p0.mle, alpha = shape.mle, beta = scale.mle)
    
    new.obs.ll <- sum(log(apply(dens1, 1, sum)))
    diff <- abs((new.obs.ll - old.obs.ll)/old.obs.ll)
    
    
    converge <- (abs(old.lambda.mle - lambda.mle)/old.lambda.mle< tol) && 
      all(abs(old.shape.mle - shape.mle)/old.shape.mle < tol) && 
      all(abs(old.scale.mle - scale.mle)/old.scale.mle < tol) &&
      all(abs(old.p0.mle - p0.mle)/old.p0.mle < tol) 
    
    
    old.obs.ll <- new.obs.ll
    ll <- c(ll, old.obs.ll)
    iter = iter + 1
    
    
    if (iter == maxits){
      cat("WARNING! NOT CONVERGENT!", "\n")}
  }
  return(list(g=z[,1], lambda = lambda.mle, p0 = p0.mle, shape = shape.mle, scale = scale.mle, iter = iter, converge = converge))
}


library(EnvStats)
### This function assumes the true matching status are observed
EM_hurdle_gammaK_obs = function(comp_mat, K, nB){
  
  M = comp_mat[comp_mat[,K+1]==1,1:K]
  
  U = comp_mat[comp_mat[,K+1]==0,1:K]
  
  N =  nrow(comp_mat)
  
  p = nrow(M)/N
  
  lambda.mle <- c(p,1-p)
  
  p0.mle <- rbind(colMeans(M==0), colMeans(U==0))
  p0.mle[p0.mle < 1/N] = 1/N
  
  shape.mle = matrix(0, ncol =K, nrow = 2)
  scale.mle = matrix(0, ncol =K, nrow = 2)
  for (k in 1:K){
    Mk = M[M[,k]>0,k]
    fitM = egamma(Mk)
    shape.mle[1,k] = fitM$parameters[1]
    scale.mle[1,k] = fitM$parameters[2]
    
    Uk = U[U[,k]>0,k]
    fitU = egamma(Uk)
    shape.mle[2,k] = fitU$parameters[1]
    scale.mle[2,k] = fitU$parameters[2]
  }
  
  
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
  
  dens <- function(X, K, lambda, p0, alpha, beta) {
    ncomp <- length(lambda)
    
    L = 1
    for (k in 1:K){
      temp = sapply(1:ncomp, function(j) dhgamma(X[,k], p0 = p0[j,k], shape = alpha[j,k], 
                                                 scale = beta[j,k]))
      L = L*temp
    }
    
    density = t(lambda * t(L))
    density
  }
  
  X = comp_mat[,1:K]
  dens1 <- dens(X=X, K=K,  lambda = lambda.mle, p0 = p0.mle, alpha = shape.mle, beta = scale.mle)
  
  z <- dens1/apply(dens1, 1, sum)
    
   

  return(list(g=z[,1]))
}

#### EM algorithm for binary comparison
library(matrixStats)
EM_binary <- function(comp_mat, nB, K,  tol = 1e-6, maxits = 500){
  
  # Starting point
  N <- nrow(comp_mat)
  p = nB/N
  
  
  u = rep(0.1,K)
  m = rep(0.8,K)
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  
  g = rep(0,N) # probability of being in Match  for each pair l
  it = 0
  converge = FALSE
  
  
  while ((!converge) & (it < maxits)){ 
    
    p.old = p
    m.old = m
    u.old = u
    ### E
    # Compute expectation
    
    m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
    u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
    
    
    probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
    probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
    
    
    g = p*probM/(p*probM+(1-p)*probU)
    
    ### Maximization
    g_mat = matrix(rep(g,K),ncol = K)
    
    p = sum(g)/N
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
    
    
    it = it + 1
    
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) && 
      all(abs(u.old - u)/u.old < tol)
    
    if (it == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }
  }
  
  
  
  return(list(g=g, p=p, m=m, u = u, it = it, converge = converge))
}

## This version assumes the true matching statuses are observed
EM_binary_obs <- function(comp_mat,nB, K){
  
  M = comp_mat[comp_mat[,K+1]==1,1:K]
  U = comp_mat[comp_mat[,K+1]==0,1:K]
  
  m = colMeans(M)
  u = colMeans(U)
  
  N = nrow(comp_mat)
  p = nrow(M)/N
  
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  
  m_mat = matrix(rep(m,N), nrow =N, byrow = TRUE)
  u_mat = matrix(rep(u,N), nrow =N, byrow = TRUE)
  
  
  probM = rowProds(m_mat^(comp_mat)*(1-m_mat)^(1-comp_mat))
  probU = rowProds(u_mat^(comp_mat)*(1-u_mat)^(1-comp_mat))
  
  
  g = p*probM/(p*probM+(1-p)*probU)
  
  
  return(list(g=g))
}



################ 3 categorical
library(clue)
library(matrixStats)
library(klaR)


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




E3 <- function(m2, m1, m0, u2,u1, u0, N, c2, c1, c0){
  m2.mat = matrix(rep(m2,N), nrow =N, byrow = TRUE)
  m1.mat = matrix(rep(m1,N), nrow =N, byrow = TRUE)
  m0.mat = matrix(rep(m0,N), nrow =N, byrow = TRUE)
  
  u2.mat = matrix(rep(u2,N), nrow =N, byrow = TRUE)
  u1.mat = matrix(rep(u1,N), nrow =N, byrow = TRUE)
  u0.mat = matrix(rep(u0,N), nrow =N, byrow = TRUE)
  
  
  p2M = rowProds(m2.mat^(c2)*m1.mat^(c1)*m0.mat^(c0))
  p2U = rowProds(u2.mat^(c2)*u1.mat^(c1)*u0.mat^(c0))
  
  return(cbind(p2M, p2U))
}



M3 <- function(g,K2,c2,c1,c0){
  N = length(g)
  g2.mat = matrix(rep(g,K2),ncol = K2)

  m2 = colSums(g2.mat*c2)/sum(g)
  m2 = check1(m2,N)
  m1 = colSums(g2.mat*c1)/sum(g)
  m1 = check1(m1,N)
  m0 = 1-m2-m1 #check 0,1

  if (length(which(m0 < 1/N)) > 0){
    m2[which(m0 < 1/N)] = m2[which(m0 < 1/N)] - 1/(2*N)
    m1[which(m0 < 1/N)] = m1[which(m0 < 1/N)] - 1/(2*N)
    m0[which(m0 < 1/N)] = m0[which(m0 < 1/N)] + 1/N
  }

  u2 = colSums((1-g2.mat)*c2)/sum(1-g)
  u2 = check1(u2,N)
  u1 = colSums((1-g2.mat)*c1)/sum(1-g)
  u1 = check1(u1,N)
  u0 = 1-u2-u1 #check 0,1

  if (length(which(u0 < 1/N)) > 0){
    u2[which(u0 < 1/N)] = u2[which(u0 < 1/N)] - 1/(2*N)
    u1[which(u0 < 1/N)] = u1[which(u0 < 1/N)] - 1/(2*N)
    u0[which(u0 < 1/N)] = u0[which(u0 < 1/N)] + 1/N
  }
  return(list(m2 = m2,m1 = m1,m0=m0,u2=u2,u1=u1,u0=u0))
}


################### EM algorithm for three categorical comparison

EM3 <- function(comp_mat, nB, K, tol = 1e-6, maxits = 500){
  
  # Starting point
  N <- nrow(comp_mat)
  p = nB/N
  
  
  m2 = rep(0.1,K)
  m1 = rep(0.1,K)
  m0 = rep(0.8,K)
  
  
  u2 = rep(0.7,K)
  u1 = rep(0.2,K)
  u0 = rep(0.1,K)
  
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c0 <- as.numeric(comp_mat==0)
  c1 <- as.numeric(comp_mat==1)
  c2 <- as.numeric(comp_mat==2)
  
  loglikelihood <- function(p,N,m2, m1, m0, u2,u1, u0, c2, c1, c0){
    
    p2 = E3(m2, m1, m0, u2,u1, u0, N, c2, c1, c0)
    pM = p2[,1]
    pU = p2[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #likelihood not complete likelihood
    return(ll)
  }
  
  iter <- 0
  converge = FALSE
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m2, m1, m0, u2,u1, u0, c2, c1, c0)
  ll = old.ll
  
  
  while (!converge && iter < maxits){ 
    
    p.old = p
    
    m2.old = m2 #length K2
    m1.old = m1 #length K2
    m0.old = m0
    
    u2.old = u2 #length K2
    u1.old = u1 #length K2
    u0.old = u0
    
    m.old = c(m0.old,m1.old,m2.old)
    u.old = c(u0.old,u1.old,u2.old)
    
    ############################ E step
    p2 = E3(m2, m1, m0, u2,u1, u0, N, c2, c1, c0)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################### M step
    p = sum(g)/N
    
    max2 = M3(g,K,c2,c1,c0)
    m2 = max2$m2
    m1 = max2$m1
    m0 = max2$m0
    u2 = max2$u2
    u1 = max2$u1
    u0 = max2$u0
    
    m = c(m0,m1,m2)
    u = c(u0,u1,u2)
    
    ### Stopping 
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) &&
      all(abs(u.old - u)/u.old < tol)
    
    new.ll = loglikelihood(p,N,m2, m1, m0, u2,u1, u0, c2, c1, c0)
    
    diff <- new.ll - old.ll
    
    old.ll <- new.ll
    ll <- c(ll, old.ll)
    
    iter = iter + 1
    
    
    if (iter == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = FALSE
    }
    
  }
  
  
  return(list(g=g, loglikelihood = ll, iter = iter, converge = converge))
}


## This version assumes the true matching statuses are observed
EM3_obs <- function(comp_mat, nB, K){
  
  M = comp_mat[comp_mat[,K+1]==1,1:K]
  U = comp_mat[comp_mat[,K+1]==0,1:K]
  
  M20 <- (M==0)+0
  M21 <- (M==1)+0
  M22 <- (M==2)+0
  
  
  U20 <- (U==0)+0
  U21 <- (U==1)+0
  U22 <- (U==2)+0
  
  
  # initialization
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- (comp_mat==0)+0
  c21 <- (comp_mat==1)+0
  c22 <- (comp_mat==2)+0
  
  
  p = nrow(M)/N
  
  m22 = colMeans(M22)
  m21 = colMeans(M21)
  m20 = colMeans(M20)
  
  u22 = colMeans(U22)
  u21 = colMeans(U21)
  u20 = colMeans(U20)
  
  
  ############################ E step
  p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
  
  pM = p2[,1]
  pU = p2[,2]
  
  
  # Expectations
  g = p*pM/(p*pM+(1-p)*pU)
  
  
  
  return(list(g=g))
}

