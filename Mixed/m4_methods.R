library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)
library(RecordLinkage)


######## Record linkage implementation of proposed model
FS_mixed <- function(datA, datB, K1, K2, K3){
  nA = nrow(datA)
  nB =  nrow(datB)
  K = K1 + K2 + K3
  
  comp_mat = compare(datA, datB, K1, K2, K3)
  fit = EM_mixed(comp_mat, datA, datB, K1, K2, K3)
  g = fit$g
  converge = fit$converge
  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 1
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  return(c(TPR_FS5, PPV_FS5, converge))

}

### Record linkage implementation of standard model
FS_binary <- function(datA, datB, K1, K2, K3){
  
  nA = nrow(datA)
  nB =  nrow(datB)
  K = K1 + K2 +K3
  
  comp_mat = compare_binary(datA, datB, K)
  fit = EM_binary(comp_mat, datA, datB, K)
  
  g = fit$g
  converge = fit$converge
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 1
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  return(c(TPR_FS5, PPV_FS5, converge))
}

### Record linkage implementation of standard model by using package RecordLinkage
FS_RecordLinkage <- function(datA, datB, K1, K2, K3){
  
  K <- K1+K2+K3
  rpairs <- compare.linkage(datA[,1:K], datB[,1:K], identity1 = datA[, K+1], identity2 = datB[,K+1])
  weights <- emWeights(rpairs)
  
  p = 1/nrow(datA)
  q = 0.5
  # Relation between matching probability threshold and FS score
  threshold <- log2((1-p)/(1-q)*q/p)
  
  results <- emClassify(weights, threshold)
  measures <- getErrorMeasures(results)
  TPR <- measures$sensitivity
  PPV <- measures$ppv
  
  return(c(TPR, PPV))
}

