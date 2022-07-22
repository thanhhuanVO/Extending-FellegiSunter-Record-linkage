library(clue)

tpr_ppv <- function(g, indM, n_matches){
  #thresholds = seq(0,1,length.out = 11)
  
  thresholds = c(0,
                 1e-100,1e-20,1e-5,1e-2,0.05,
                 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
                 0.95,0.975,0.995,0.999,
                 1)
  
  one_threshold <- function(s,g, indM, n_matches){
    gm = g[indM]
    n_true = sum(gm>s)
    n_predict = sum(g>s)
    
    TPR = n_true/n_matches
    if (n_predict >0){
      PPV = n_true/n_predict
    }else{
      PPV = 1
    }
    return(c(TPR, PPV))
  }
  
  results = t(sapply(thresholds,FUN = one_threshold, g = g, indM = indM, n_matches = n_matches))
  
}

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
  ## Thresholds
  indM = which(comp_mat[,K+1]==1)
  res = tpr_ppv(g=g, indM = indM, n_matches = nB)
  #################### Observed
  fit = EM_hurdle_gammaK_obs(comp_mat = comp_mat,K = K,nB = nB)
  
  g = fit$g
  
  ## Thresholds
  indM = which(comp_mat[,K+1]==1)
  res_obs = tpr_ppv(g=g, indM = indM, n_matches = nB)
  

  
  return(c(converge,as.vector(res), as.vector(res_obs)))
}



FS_gamma_sqr <- function(datA, datB, K,  tol = 1e-6, maxits = 500){
  #start = list(lambda = c(0.002,0.998), p0 = matrix(rep(c(0.3,0.02),K), ncol =K), 
  #             alpha = matrix(rep(c(1,1),K), ncol =K), 
  #             beta = matrix(rep(c(1,50),K), ncol =K))
  
  comp_mat = compare_abs(datA = datA,datB = datB, K =K)
  nB =  nrow(datB)
  nA =  nrow(datA)
  
  fit = EM_hurdle_gammaK(X = comp_mat[,1:K],K = K,nB = nB, tol = tol, maxits = maxits)
  
  g = fit$g
  
  ## Thresholds
  indM = which(comp_mat[,K+1]==1)
  tpr_ppv = tpr_ppv(g=g, indM = indM, n_matches = nB)
  #################### Observed
  fit = EM_hurdle_gammaK_obs(comp_mat = comp_mat,K = K,nB = nB)
  
  g = fit$g
  
  
  ## Thresholds
  indM = which(comp_mat[,K+1]==1)
  tpr_ppv_obs = tpr_ppv(g=g, indM = indM, n_matches = nB)
  
  return(c( as.vector(tpr_ppv),  as.vector(tpr_ppv_obs)))
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
  res = tpr_ppv(g=g, indM = indM, n_matches = nB)
  
  ###################### Observed
  fit = EM_binary_obs(comp_mat = comp_mat, nB =  nB, K = K)
  
  g = fit$g
  
  ## Threshold
  indM = which(comp_mat[,K+1]==1)
  res_obs = tpr_ppv(g=g, indM = indM, n_matches = nB)
  
  return(c(converge, as.vector(res),  as.vector(res_obs)))
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
  res = tpr_ppv(g=g, indM = indM, n_matches = nB)
  ################ Observed
  fit = EM3_obs(comp_mat = comp_mat, nB =  nB, K = K)
  
  g = fit$g
  
  ## Threshold
  indM = which(comp_mat[,K+1]==1)
  res_obs = tpr_ppv(g=g, indM = indM, n_matches = nB)
  
  return(c(converge, as.vector(res),  as.vector(res_obs)))
}
