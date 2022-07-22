### Generate two database A (nA units) and database B (nB units) 
# with K continuous matching variables
# We consider different combination distribution for matching data and errors


#Exponential data and exponential error
generate_data_exp_exp = function(nA, nB, K, lambdaK, error, lambdaE, round = FALSE){
  if (round == FALSE){
    XA = matrix(rexp(nA*K,rate = lambdaK),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(rexp(nB*K, rate = lambdaE),ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
  
    
  }else if (round == TRUE){
    XA = matrix(ceiling(rexp(nA*K,rate = lambdaK)),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(ceiling(rexp(nB*K, rate = lambdaE)), ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
  }
  
  #####
  return(list(dataA=datA, dataB = datB))
}


#Uniform data and exponential error
generate_data_unif_exp = function(nA, nB, K, maxK, error, lambdaE, round = FALSE){
  if (round == FALSE){
    XA = matrix(runif(nA*K, min = 0 , max = maxK),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(rexp(nB*K, rate = lambdaE),ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
    
    
  }else if (round == TRUE){
    XA = matrix(ceiling(runif(nA*K, min = 0 , max = maxK)),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(ceiling(rexp(nB*K, rate = lambdaE)), ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
  }
  
  #####
  return(list(dataA=datA, dataB = datB))
}

#Exponential data and normal error
generate_data_exp_norm = function(nA, nB, K, lambdaK, error, meanE, sdE, round = FALSE){
  if (round == FALSE){
    XA = matrix(rexp(nA*K,rate = lambdaK),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(rnorm(nB*K, mean = meanE, sd = sdE),ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
    
    
  }else if (round == TRUE){
    XA = matrix(ceiling(rexp(nA*K,rate = lambdaK)),ncol = K)
    #S= sample(1:nA,nB)
    XA1 = XA[1:nB,]
    XA2  = XA[-(1:nB),]
    
    X = matrix(rbinom(nB*K, 1, error), ncol =K)
    XB = (1-X)* XA1 + X*(XA1 + matrix(ceiling(rnorm(nB*K, mean = meanE, sd = sdE)), ncol =K))
    
    datA = matrix(0,nrow = nA, ncol = K+1)
    datB = matrix(0,nrow = nB, ncol = K+1)
    datA[,1:K] = XA
    datA[,K+1] = 1:nA #id
    datB[,1:K] = XB
    datB[,K+1] = 1:nB #id
  }
  
  #####
  return(list(dataA=datA, dataB = datB))
}