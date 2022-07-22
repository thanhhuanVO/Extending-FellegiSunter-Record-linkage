library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)


FS <- function(datA, datB, K,tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB =  nrow(datB)
  
  comp_mat <- compare_binary(datA, datB, K)
  
  fit <- EM_binary(comp_mat, datA, datB, K,tol = tol, maxits = maxits)
  g = fit$g
  converge = fit$converge
  
  
  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_FS11= nTruePredict/nB
  # PPV_FS11 = nTruePredict/nPredict
  
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



# Fellegi-Sunter with 3 categorical comparison
FS3 <- function(datA, datB, K, tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB = nrow(datB)
  comp_mat <- compare3(datA, datB, K=K)
  
  
  ## Using EM with the above estimated starting point
  fit = EM3(comp_mat, datA, datB, K, tol=tol, maxits = maxits)
  g = fit$g
  converge = fit$converge
  

  

  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_FS11= nTruePredict/nB
  # PPV_FS11 = nTruePredict/nPredict

  
  
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


bayesian <- function(datA, datB, K){
  nB = nrow(datB)
  nA = nrow(datA)
  bayes = recordLink(datB[,1:K], datA[,1:K], eps_plus =0.01, eps_minus = 0.01,use_diff = FALSE)
  g = as.vector(t(bayes))

  # Gmat = bayes
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_bayes11= nTruePredict/nB
  # PPV_bayes11 = nTruePredict/nPredict
  
  #############################
  threshold = 0.5
  upper = which(bayes >= threshold)
  indexB = upper%%nB
  indexB[which(indexB==0)] = nB
  indexA = ceiling(upper/nB)
  idB = datB[indexB, K+1]
  idA = datA[indexA, K+1]
  
  TPR_bayes5 = sum(idA==idB)/nB
  if (length(idB) == 0){
    PPV_bayes5 = 1
  }else{
    PPV_bayes5 = sum(idA==idB)/length(idB)
  }
  
  
 
 
  return(c(TPR_bayes5, PPV_bayes5, 1))
}

### Fellegi-Sunter with 4 categorical comparison
FS4 <- function(datA, datB, K, comp_mat4, tol=1e-6, maxits = 500){
  nA = nrow(datA)
  nB = nrow(datB)
  comp_mat <- compare4(datA, datB, K=K)
  
  
  ## Using EM
  fit = EM4(comp_mat, datA, datB, K, tol=tol, maxits = maxits)
  
  g = fit$g
  converge = fit$converge

  
  
  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_FS11= nTruePredict/nB
  # PPV_FS11 = nTruePredict/nPredict
  
  threshold = 0.5
  index = which(g>= threshold)
  TPR_FS5 = sum(comp_mat[index,K+1])/nB
  if (length(index)==0){
    PPV_FS5 = 0
  }else{
    PPV_FS5 = sum(comp_mat[index,K+1])/length(index)
  }
  
  
  return(c(TPR_FS5, PPV_FS5, converge))
}

