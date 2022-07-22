library(clue)

FS_gamma <- function(datA, datB, K,  tol = 1e-6, maxits = 500){
  #start = list(lambda = c(0.002,0.998), p0 = matrix(rep(c(0.3,0.02),K), ncol =K), 
  #             alpha = matrix(rep(c(1,1),K), ncol =K), 
  #             beta = matrix(rep(c(1,50),K), ncol =K))
  
  comp_mat = compare_abs(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  
  fit = EM_hurdle_gammaK(X = comp_mat[,1:K],K = K,nB = nB, tol = tol, maxits = maxits)
  
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
  
  # TPR_11= nTruePredict/nB
  # if (nPredict > 0){
  #   PPV_11 = nTruePredict/nPredict
  # }else{
  #   PPV_11 = 0
  # }
  
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

FS_gamma_start <- function(datA, datB, K,  tol = 1e-6, maxits = 500){
  #start = list(lambda = c(0.002,0.998), p0 = matrix(rep(c(0.3,0.02),K), ncol =K), 
  #             alpha = matrix(rep(c(1,1),K), ncol =K), 
  #             beta = matrix(rep(c(1,50),K), ncol =K))
  
  comp_mat = compare_abs(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  fit = EM_hurdle_gammaK_start(X = comp_mat[,1:K],K = K,nB = nB, tol = tol, maxits = maxits)
  
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
  
  # TPR_11= nTruePredict/nB
  # if (nPredict > 0){
  #   PPV_11 = nTruePredict/nPredict
  # }else{
  #   PPV_11 = 0
  # }
  
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
  #start = list(lambda = c(0.002,0.998), p0 = matrix(rep(c(0.3,0.02),K), ncol =K), 
  #             alpha = matrix(rep(c(1,1),K), ncol =K), 
  #             beta = matrix(rep(c(1,50),K), ncol =K))
  
  comp_mat = compare_sqr(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  fit = EM_hurdle_gammaK(X = comp_mat[,1:K],K = K,nB = nB, tol = tol, maxits = maxits)
  
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
  # TPR_11= nTruePredict/nB
  # if (nPredict > 0){
  #   PPV_11 = nTruePredict/nPredict
  # }else{
  #   PPV_11 = 0
  # }
  
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
  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_11= nTruePredict/nB
  # if (nPredict > 0){
  #   PPV_11 = nTruePredict/nPredict
  # }else{
  #   PPV_11 = 0
  # }
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
  
  # Gmat = matrix(g, nrow = nB, byrow = TRUE)
  # Gmat[Gmat<0.5] = 0
  # opti_cate <- solve_LSAP(Gmat, maximum = TRUE)
  # predict_cate = cbind(seq_along(opti_cate), opti_cate)
  # temp = (predict_cate[,1]-1)*nA+ predict_cate[,2]
  # # Percentage of correct link
  # nPredict = sum((g[temp]>=0.5))
  # nTruePredict = sum((predict_cate[,2] == datB[,K+1])&(g[temp]>=0.5))
  # 
  # TPR_11= nTruePredict/nB
  # if (nPredict > 0){
  #   PPV_11 = nTruePredict/nPredict
  # }else{
  #   PPV_11 = 0
  # }
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
