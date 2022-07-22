library(clue)

## Different record linkage methods corresponding to proposed comparison approach
FS_gamma <- function(datA, datB, K,  tol = 1e-6, maxits = 500){
  
  comp_mat = compare_abs(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  
  fit = EM_hurdle_gammaK(X = comp_mat[,1:K],K = K,nB = nB, tol = tol, maxits = maxits)
  
  g = fit$g
  converge = fit$converge
  
  ## Threshold
  
  
  indM = which(comp_mat[,K+1]==1)
  gm = g[indM]
  
  n_true = sum(gm>0.5)
  n_predict = sum(g>0.5)
  n_matches = nB
  
  TPR_5 = n_true/n_matches
  if (n_predict >0){
    PPV_5 = n_true/n_predict
  }else{
    PPV_5 = 1
  }
  
  
  return(c(TPR_5, PPV_5, converge))
}


FS_gamma_sqr <- function(datA, datB, K,  tol = 1e-6, maxits = 500){

  comp_mat = compare_sqr(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  fit = EM_hurdle_gammaK(X = comp_mat[,1:K],K = K,nB = nB, tol = tol, maxits = maxits)
  
  g = fit$g
  converge = fit$converge

  ## Threshold
  
  
  indM = which(comp_mat[,K+1]==1)
  gm = g[indM]
  
  n_true = sum(gm>0.5)
  n_predict = sum(g>0.5)
  n_matches = nB
  
  TPR_5 = n_true/n_matches
  if (n_predict >0){
    PPV_5 = n_true/n_predict
  }else{
    PPV_5 = 1
  }
  
  
  return(c(TPR_5, PPV_5, converge))
}

FS_binary <- function(datA, datB, K,  tol = 1e-6, maxits = 500){
  
  comp_mat = compare_binary(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  fit = EM_binary(comp_mat = comp_mat[,1:K], nB =  nB, K = K, tol = tol, maxits = maxits)
  
  g = fit$g
  converge = fit$converge

  ## Threshold
  
  
  indM = which(comp_mat[,K+1]==1)
  gm = g[indM]
  
  n_true = sum(gm>0.5)
  n_predict = sum(g>0.5)
  n_matches = nB
  
  TPR_5 = n_true/n_matches
  if (n_predict >0){
    PPV_5 = n_true/n_predict
  }else{
    PPV_5 = 0
  }
  
  
  
  return(c(TPR_5, PPV_5, converge))
}

FS3 <- function(datA, datB, K,  tol = 1e-6, maxits = 500){
  
  comp_mat = compare3(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  fit = EM3(comp_mat = comp_mat[,1:K], nB =  nB, K = K, tol = tol, maxits = maxits)
  
  g = fit$g
  converge = fit$converge
  

  ## Threshold
  
  
  indM = which(comp_mat[,K+1]==1)
  gm = g[indM]
  
  n_true = sum(gm>0.5)
  n_predict = sum(g>0.5)
  n_matches = nB
  
  TPR_5 = n_true/n_matches
  if (n_predict >0){
    PPV_5 = n_true/n_predict
  }else{
    PPV_5 = 0
  }
  
  
  return(c(TPR_5, PPV_5, converge))
}
