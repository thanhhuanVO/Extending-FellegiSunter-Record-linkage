library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(EnvStats)


doOne_exp_exp <- function(nA, nB, K, lambdaK, error, lambdaE, round ){
  library(tictoc)
  source("c1_generate_data.R")
  source("c2_compare.R")
  source("c3_EM.R")
  source("c4_methods.R")
  
  
  data = generate_data_exp_exp(nA = nA, nB=nB,K=K, lambdaK = lambdaK, error = error, lambdaE = lambdaE, round = round)
  datA = data$dataA
  datB = data$dataB
  

  tic()
  res_FS_gamma = FS_gamma(datA, datB, K, tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS_gamma <- as.numeric(temp$toc-temp$tic)
  
  tic.clear()
  tic()
  res_FS3 = FS3(datA, datB, K, tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS3 <- as.numeric(temp$toc-temp$tic)
  
  tic.clear()
  tic()
  res_FS_binary = FS_binary(datA, datB, K, tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS_binary <- as.numeric(temp$toc-temp$tic)
  
  return(c(res_FS_gamma, res_FS3,  res_FS_binary, time_FS_gamma, time_FS3, time_FS_binary))
}

doOne_exp_norm <- function(nA, nB, K, lambdaK, error, meanE, sdE, round ){
  source("1c_generate_data.R")
  source("2c_compare.R")
  source("3c_EM.R")
  source("4c_methods.R")
  
  data = generate_data_exp_norm(nA = nA, nB=nB,K=K, lambdaK = lambdaK, error = error, meanE = meanE, sdE = sdE, round = round)
  datA = data$dataA
  datB = data$dataB
  
  res_FS_gamma_sqr = FS_gamma_sqr(datA, datB, K)
  res_FS_gamma = FS_gamma(datA, datB, K)
  res_FS3 = FS3(datA, datB, K)
  res_FS_binary = FS_binary(datA, datB, K)
  
  
  return(c(res_FS_gamma_sqr, res_FS_gamma, res_FS3, res_FS_binary))
}

doOne_unif_exp <- function(nA, nB, K, maxK, error, lambdaE, round ){
  source("1c_generate_data.R")
  source("2c_compare.R")
  source("3c_EM.R")
  source("4c_methods.R")
  
  data = generate_data_unif_exp(nA = nA, nB=nB,K=K, maxK = maxK, error = error, lambdaE = lambdaE, round = round)
  datA = data$dataA
  datB = data$dataB
  
  res_FS_gamma = FS_gamma(datA, datB, K)
  res_FS3 = FS3(datA, datB, K)
  res_FS_binary = FS_binary(datA, datB, K)
  
  
  return(c(res_FS_gamma, res_FS3, res_FS_binary))
}

vlis_exp_exp1 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    lambdaK = list(type="frozen", value = 0.02),
    error = list(type="grid", value = c(0.1,0.2,0.3)),
    lambdaE = list(type="grid", value = c(1/2,1/3,1/4)), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_exp_exp1_test <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    lambdaK = list(type="frozen", value = 0.02),
    error = list(type="grid", value = c(0.2,0.3)),
    lambdaE = list(type="grid", value = c(1/3,1/4)), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_exp_exp2 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    lambdaK = list(type="grid", value = c(0.01,0.02,0.03)),
    error = list(type="frozen", value = 0.2),
    lambdaE = list(type="frozen", value = 1/3), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_exp_exp3 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="grid", value = c(400,800,1200)),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    lambdaK = list(type="frozen", value = 0.02),
    error = list(type="frozen", value = 0.2),
    lambdaE = list(type="frozen", value = 1/3), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_exp_exp4 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K = list(type="grid", value = c(2,3,4)),
    lambdaK = list(type="grid", value = c(0.005, 0.002)),
    error = list(type="frozen", value = 0.2),
    lambdaE = list(type="frozen", value = 1/3), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_unif_exp <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    maxK = list(type="frozen", value = 100),
    error = list(type="grid", value = c(0.1,0.2,0.3)),
    lambdaE = list(type="grid", value = c(1/2,1/3,1/4)), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_exp_norm <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    lambdaK = list(type="frozen", value = 0.02),
    error = list(type="grid", value = c(0.1,0.2,0.3)),
    meanE = list(type="grid", value = c(1,2,3)), 
    sdE = list(type="frozen", value = 1), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

runSims <- function(vList=vlis1(),doOne = doOne1, seedList=NULL){
  res <- simsalapar::doForeach(vList,  doOne = doOne,  cluster=makeCluster(8, type="PSOCK"), seed = seedList)
  return(res)
}

#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")
#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Continuous\\Programs")

nsim = 100

res_exp_exp1 <- runSims(vlis_exp_exp1(nsim), doOne = doOne_exp_exp)
save(res_exp_exp1, file="0623_100s_res_exp_exp1.RData")

#res_exp_norm <- runSims(vlis_exp_norm(nsim), doOne = doOne_exp_norm)
#save(res_exp_norm, file="0329_1000s_res_exp_norm.RData")
 
#res_unif_exp <- runSims(vlis_unif_exp(nsim), doOne = doOne_unif_exp)
#save(res_unif_exp, file="0410_1000s_res_unif_exp.RData")

#res_exp_exp2 <- runSims(vlis_exp_exp2(nsim), doOne = doOne_exp_exp)
#save(res_exp_exp2, file="0421_1000s_res_exp_exp2.RData")
 
#res_exp_exp3 <- runSims(vlis_exp_exp3(nsim), doOne = doOne_exp_exp)
#save(res_exp_exp3, file="0422_1000s_res_exp_exp3.RData")

#res_exp_exp4 <- runSims(vlis_exp_exp4(nsim), doOne = doOne_exp_exp)
#save(res_exp_exp4, file="0420_1000s_res_exp_exp4.RData")


