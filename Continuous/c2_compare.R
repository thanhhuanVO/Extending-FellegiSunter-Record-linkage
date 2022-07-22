compare_abs <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = abs(temp[,1]-temp[,2])
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

compare_sqr <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = (temp[,1]-temp[,2])^2
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

compare_binary <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)

    gamma.k = as.numeric(temp[,1]==temp[,2])

    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

compare3 <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = abs(temp[,1]-temp[,2])
      gamma.k[(gamma.k <= 3) & (gamma.k > 0)] = 1
      gamma.k[(gamma.k > 3) ] = 2
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}
