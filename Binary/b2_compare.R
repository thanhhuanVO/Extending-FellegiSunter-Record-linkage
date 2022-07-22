library(doParallel)

# The proposed 3 categorical comparison for binary matching variables
compare3 <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = temp[,1]+temp[,2]
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}

# The proposed 4 categorical comparison for binary matching variables
compare4 <- function(datA, datB, K){
  compare1 <- function(k, datA, datB, K){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    if (k == (K+1)){
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else{
      gamma.k = 2*temp[,1]+temp[,2]
    }
    
    return(gamma.k)
  }
  
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K =K)
  return(comp_mat)
}


# # The standard binary comparison for binary matching variables
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
