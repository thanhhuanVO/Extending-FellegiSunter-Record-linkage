library(clue)
library(matrixStats)
library(klaR)

############ Implementation of EM algorithm
check <- function(param, N){
  if (length(which(param < 1/N)) > 0) {
    param[which(param < 1/N)] <- 1/N
  }
  
  if (length(which(param > (1-1/N))) > 0) {
    param[which(param > (1-1/N))] <- 1-1/N
  }
  return(param)
}

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


start_3 <- function(datA, datB,K){
  e = 0.01 #is q in your documents
  nA = nrow(datA)
  nB = nrow(datB)
  p = 1/nA
  
  # This is the empirical prevalence that I use for estimate starting point
  # To compute the true, you replace it with true parameters
  prev =  colMeans(datA[,1:K])
  
  p0M = (1-e)*(1-prev) #(0,0)
  p1M = rep(e,K) #(0,1) + (1,0)
  p2M = 1-p1M-p0M # (1,1)
  
  p0U = (1-prev)*((1-e)*(1-prev)+e*prev) #(0,0)
  p1U =   (1-prev)*((1-e)*prev + e*(1-prev)) + prev*((1-e)*(1-prev) + e*prev) #(0,1) + (1,0)
  p2U = 1- p1U - p0U # (1,1)
  
  return(c(p, p2M, p1M, p0M, p2U, p1U, p0U ))
  
}


start_4 <- function(datA, datB,K){
  e = 0.01 #is q in your documents
  nA = nrow(datA)
  nB = nrow(datB)
  p = 1/nA
  
  # This is the empirical prevalence that I use for estimate starting point
  # To compute the true, you replace it with true parameters
  prev =  colMeans(datA[,1:K])
  
  # To compute the true, you replace it with true parameters
  p0M = (1-e)*(1-prev) #(0,0)
  p1M = e*(1-prev) #(0,1) 
  p2M = e*prev #(1,0)
  p3M = 1-p2M-p1M-p0M # (1,1)
  
  p0U = (1-prev)*((1-e)*(1-prev)+e*prev) #(0,0)
  p1U = (1-prev)*((1-e)*prev + e*(1-prev))  #(0,1)
  p2U = prev*((1-e)*(1-prev) + e*prev)  #(1,0)
  p3U = 1- p2U- p1U - p0U # (1,1)
  
  
  return(c(p, p3M, p2M, p1M, p0M, p3U, p2U, p1U, p0U ))
  
}

E3 <- function(m22, m21, m20, u22,u21, u20, N, c22, c21, c20){
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

E4 <- function(m23, m22, m21, m20, u23, u22,u21, u20, N, c23, c22, c21, c20){
  m23.mat = matrix(rep(m23,N), nrow =N, byrow = TRUE)
  m22.mat = matrix(rep(m22,N), nrow =N, byrow = TRUE)
  m21.mat = matrix(rep(m21,N), nrow =N, byrow = TRUE)
  m20.mat = matrix(rep(m20,N), nrow =N, byrow = TRUE)
  
  u23.mat = matrix(rep(u23,N), nrow =N, byrow = TRUE)
  u22.mat = matrix(rep(u22,N), nrow =N, byrow = TRUE)
  u21.mat = matrix(rep(u21,N), nrow =N, byrow = TRUE)
  u20.mat = matrix(rep(u20,N), nrow =N, byrow = TRUE)
  
  
  p2M = rowProds(m23.mat^(c23)*m22.mat^(c22)*m21.mat^(c21)*m20.mat^(c20))
  p2U = rowProds(u23.mat^(c23)*u22.mat^(c22)*u21.mat^(c21)*u20.mat^(c20))
  
  return(cbind(p2M, p2U))
}

M3 <- function(g,K2,c22,c21,c20){
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

M4 <- function(g,K2,c23,c22,c21,c20){
  N = length(g)
  g2.mat = matrix(rep(g,K2),ncol = K2)
  
  m23 = colSums(g2.mat*c23)/sum(g)
  m23 = check1(m23,N)
  m22 = colSums(g2.mat*c22)/sum(g)
  m22 = check1(m22,N)
  m21 = colSums(g2.mat*c21)/sum(g)
  m21 = check1(m21,N)
  m20 = 1-m23-m22-m21 #check 0,1
  
  if (length(which(m20 < 1/N)) > 0){
    m23[which(m20 < 1/N)] = m23[which(m20 < 1/N)] - 1/(3*N)
    m22[which(m20 < 1/N)] = m22[which(m20 < 1/N)] - 1/(3*N)
    m21[which(m20 < 1/N)] = m21[which(m20 < 1/N)] - 1/(3*N)
    m20[which(m20 < 1/N)] = m20[which(m20 < 1/N)] + 1/N 
  }
  
  u23 = colSums((1-g2.mat)*c23)/sum(1-g)
  u23 = check1(u23,N)
  u22 = colSums((1-g2.mat)*c22)/sum(1-g)
  u22 = check1(u22,N)
  u21 = colSums((1-g2.mat)*c21)/sum(1-g)
  u21 = check1(u21,N)
  u20 = 1-u22-u21 #check 0,1
  
  if (length(which(u20 < 1/N)) > 0){
    u23[which(u20 < 1/N)] = u23[which(u20 < 1/N)] - 1/(3*N)
    u22[which(u20 < 1/N)] = u22[which(u20 < 1/N)] - 1/(3*N)
    u21[which(u20 < 1/N)] = u21[which(u20 < 1/N)] - 1/(3*N)
    u20[which(u20 < 1/N)] = u20[which(u20 < 1/N)] + 1/N 
  }
  return(list(m23 = m23,m22 = m22,m21 = m21,m20=m20, u23 = u23, u22=u22,u21=u21,u20=u20))
}


EM3_ll <- function(comp_mat, datA, datB, K, tol = 1e-6, maxits = 300){
  
  # Starting point
  #start <- find_start(comp_mat, datA, datB, K)
  start <- start_3(datA, datB, K)
  
  p = start[1]
  
  m22 = start[2:(K+1)]
  m21 = start[(K+2):(2*K+1)]
  m20 = start[(2*K+2):(3*K+1)]
  
  
  u22 = start[(3*K+2):(4*K+1)]
  u21 = start[(4*K+2):(5*K+1)]
  u20 = start[(5*K+2):(6*K+1)]
  

  # initializations
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- as.numeric(comp_mat==0)
  c21 <- as.numeric(comp_mat==1)
  c22 <- as.numeric(comp_mat==2)
  
  loglikelihood <- function(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20){
    
    p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    pM = p2[,1]
    pU = p2[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #likelihood not complete likelihood
    return(ll)
  }
  
  iter <- 0
  converge = 1
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20)
  ll = old.ll


  while (diff > tol && iter < maxits){ 
    
    p.old = p
  
    m22.old = m22 #length K2
    m21.old = m21 #length K2
    m20.old = m20
    
    u22.old = u22 #length K2
    u21.old = u21 #length K2
    u20.old = u20

    
    ############################ E step
    p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################### M step
    p = sum(g)/N
    
    max2 = M3(g,K,c22,c21,c20)
    m22 = max2$m22
    m21 = max2$m21
    m20 = max2$m20
    u22 = max2$u22
    u21 = max2$u21
    u20 = max2$u20
    
    m = rbind(m20,m21,m22)
    u = rbind(u20,u21,u22)
    
    ### Stopping 
    new.ll = loglikelihood(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20)
    
    diff <- abs(new.ll - old.ll)/abs(old.ll)
    
    old.ll <- new.ll
    ll <- c(ll, old.ll)
    
    iter = iter + 1
    
    
    if (iter == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
      converge = 0
      }
      
  }

  
  return(list(g=g, loglikelihood = ll, iter = iter, converge = converge))
}

EM3 <- function(comp_mat, datA, datB, K, tol = 1e-6, maxits = 500){
  
  # Starting point
  #start <- find_start(comp_mat, datA, datB, K)
  start <- start_3(datA, datB, K)
  
  p = start[1]
  
  m22 = start[2:(K+1)]
  m21 = start[(K+2):(2*K+1)]
  m20 = start[(2*K+2):(3*K+1)]
  
  
  u22 = start[(3*K+2):(4*K+1)]
  u21 = start[(4*K+2):(5*K+1)]
  u20 = start[(5*K+2):(6*K+1)]
  
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- as.numeric(comp_mat==0)
  c21 <- as.numeric(comp_mat==1)
  c22 <- as.numeric(comp_mat==2)
  
  loglikelihood <- function(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20){
    
    p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    pM = p2[,1]
    pU = p2[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #likelihood not complete likelihood
    return(ll)
  }
  
  iter <- 0
  converge = FALSE
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20)
  ll = old.ll
  
  
  while (!converge && iter < maxits){ 
    
    p.old = p
    
    m22.old = m22 #length K2
    m21.old = m21 #length K2
    m20.old = m20
    
    u22.old = u22 #length K2
    u21.old = u21 #length K2
    u20.old = u20
    
    m.old = c(m20.old,m21.old,m22.old)
    u.old = c(u20.old,u21.old,u22.old)
    
    ############################ E step
    p2 = E3(m22, m21, m20, u22,u21, u20, N, c22, c21, c20)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################### M step
    p = sum(g)/N
    
    max2 = M3(g,K,c22,c21,c20)
    m22 = max2$m22
    m21 = max2$m21
    m20 = max2$m20
    u22 = max2$u22
    u21 = max2$u21
    u20 = max2$u20
    
    m = c(m20,m21,m22)
    u = c(u20,u21,u22)
    
    ### Stopping 
    converge <- (abs(p.old - p)/p.old < tol) && 
               all(abs(m.old - m)/m.old < tol) &&
               all(abs(u.old - u)/u.old < tol)
    
     new.ll = loglikelihood(p,N,m22, m21, m20, u22,u21, u20, c22, c21, c20)
    
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

EM4 <- function(comp_mat, datA, datB, K, tol = 1e-6, maxits = 300){
  
  # Starting point
  #start <- find_start(comp_mat, datA, datB, K)
  start <- start_4(datA, datB, K)
  
  p = start[1]
  
  m23 = start[2:(K+1)]
  m22 = start[(K+2):(2*K+1)]
  m21 = start[(2*K+2):(3*K+1)]
  m20 = start[(3*K+2):(4*K+1)]
  
  u23 = start[(4*K+2):(5*K+1)]
  u22 = start[(5*K+2):(6*K+1)]
  u21 = start[(6*K+2):(7*K+1)]
  u20 = start[(7*K+2):(8*K+1)]
  
  # initializations
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- as.numeric(comp_mat==0)
  c21 <- as.numeric(comp_mat==1)
  c22 <- as.numeric(comp_mat==2)
  c23 <- as.numeric(comp_mat==3)
  
  loglikelihood <- function(p,N,m23,m22, m21, m20,u23, u22,u21, u20, c23, c22, c21, c20){
    
    p2 = E4(m23, m22, m21, m20, u23, u22,u21, u20, N, c23, c22, c21, c20)
    pM = p2[,1]
    pU = p2[,2]
    
    ll = sum(log(p*pM+(1-p)*pU)) #likelihood not complete likelihood
    return(ll)
  }
  
  iter <- 0
  converge = FALSE
  diff <- tol + 1
  old.ll = loglikelihood(p,N,m23,m22, m21, m20, u23, u22,u21, u20, c23, c22, c21, c20)
  ll = old.ll
  
  
  while (!converge && iter < maxits){ 
    
    p.old = p
    
    m23.old = m23
    m22.old = m22 #length K2
    m21.old = m21 #length K2
    m20.old = m20
    
    u23.old = u23
    u22.old = u22 #length K2
    u21.old = u21 #length K2
    u20.old = u20
    
    m.old = c(m20.old,m21.old,m22.old, m23.old )
    u.old = c(u20.old,u21.old,u22.old, u23.old)
    ############################ E step
    p2 = E4(m23,m22, m21, m20, u23, u22,u21, u20, N, c23, c22, c21, c20)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
    ############################### M step
    p = sum(g)/N
    
    max2 = M4(g,K,c23,c22,c21,c20)
    m23 = max2$m23
    m22 = max2$m22
    m21 = max2$m21
    m20 = max2$m20
    
    u23 = max2$u23
    u22 = max2$u22
    u21 = max2$u21
    u20 = max2$u20
    

    m = c(m20,m21,m22, m23)
    u = c(u20,u21,u22, u23)
    ### Stopping 
    converge <- (abs(p.old - p)/p.old < tol) && 
      all(abs(m.old - m)/m.old < tol) &&
      all(abs(u.old - u)/u.old < tol)
    
    new.ll = loglikelihood(p,N,m23,m22, m21, m20, u23, u22,u21, u20, c23, c22, c21, c20)
    
    diff <- new.ll - old.ll
    
    old.ll <- new.ll
    ll <- c(ll, old.ll)
    
    iter = iter + 1
    
    
    if (iter == maxits) {
      cat("WARNING! NOT CONVERGENT!", "\n")
    converge = FALSE
    }
  }
  
  
  return(list(g=g, loglikelihood = ll, iter = iter, converge=converge))
}

EM_binary <- function(comp_mat, datA, datB, K, e = 0.01, tol = 1e-5, maxits = 200){
  
  # Starting point
  nB =  nrow(datB)
  N <- nrow(comp_mat)
  p = nB/N
  
  prev = (colMeans(datA[,1:K]))
  
  u = ((1-e)*(1-prev) + e*prev)*(1-prev) + prev*((1-e)*prev+e*(1-prev))
  m = rep(1 - e,K)
  
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


EM4_obs <- function(comp_mat, datA, datB, K){
  
  M = comp_mat[comp_mat[,K+1]==1,1:K]
  U = comp_mat[comp_mat[,K+1]==0,1:K]
  
  M20 <- (M==0)+0
  M21 <- (M==1)+0
  M22 <- (M==2)+0
  M23 <- (M==3)+0
  
  U20 <- (U==0)+0
  U21 <- (U==1)+0
  U22 <- (U==2)+0
  U23 <- (U==3)+0
  
  # initialization
  comp_mat  <- comp_mat[,1:K]
  N <- nrow(comp_mat)
  c20 <- (comp_mat==0)+0
  c21 <- (comp_mat==1)+0
  c22 <- (comp_mat==2)+0
  c23 <- (comp_mat==3)+0
  
    p = nrow(M)/N
    
    m23 = colMeans(M23)
    m22 = colMeans(M22)
    m21 = colMeans(M21)
    m20 = colMeans(M20)
    
    u23 = colMeans(U23)
    u22 = colMeans(U22)
    u21 = colMeans(U21)
    u20 = colMeans(U20)
    

    ############################ E step
    p2 = E4(m23,m22, m21, m20, u23, u22,u21, u20, N, c23, c22, c21, c20)
    
    pM = p2[,1]
    pU = p2[,2]
    
    
    # Expectations
    g = p*pM/(p*pM+(1-p)*pU)
    
  
  
  return(list(g=g))
}


EM3_obs <- function(comp_mat, datA, datB, K){
  
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

EM_binary_obs <- function(comp_mat, datA, datB, K){
  
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
