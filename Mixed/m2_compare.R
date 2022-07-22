library(doParallel)

# 1. This function compute comparison vector for the proposed model 
# Assume that there are K1 categorical, K2 low prevalence binary and K3 continuous matching variables
compare <- function(datA, datB, K1, K2, K3){
  # Compare the k-th matching variables 
  compare1 <- function(k, datA, datB, K1, K2, K3){
    XA.k = datA[,k]
    XB.k = datB[,k]
    temp = expand.grid(XA.k, XB.k)
    
    if (k %in% 1:K1){ # For categorical matching variables
      gamma.k = as.numeric(temp[,1]==temp[,2])
    }else if (k %in% ((K1+1):(K1+K2))){ #For low prevalence binary matching variables
      gamma.k = temp[,1]+temp[,2]
    }else if (k %in% ((K1+K2+1):(K1+K2+K3))){ # For continuous matching variables
      gamma.k = abs(temp[,1]-temp[,2])
    }else if (k==K1+K2+K3+1){
      gamma.k = as.numeric(temp[,1]==temp[,2]) #Matching indicator
    }
    
      return(gamma.k)
  }
  K = K1+K2+K3
  comp_mat = sapply(1:(K+1), FUN = compare1, datA = datA, datB = datB, K1=K1, K2=K2, K3=K3)
  
  return(comp_mat)
}

# 2. This function computes binary comparison vectors for standard Fellegi-Sunter model
# Binary comparison for all variables
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
 