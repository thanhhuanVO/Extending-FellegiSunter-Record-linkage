library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(EnvStats)

# Simulation for computing precison-recall with diferent threshold
doOne_exp_exp <- function(nA, nB, K, lambdaK, error, lambdaE, round ){
  source("c1_generate_data.R")
  source("c2_compare.R")
  source("c3_EM.R")
  source("cprc4_methods.R")
  
  data = generate_data_exp_exp(nA = nA, nB=nB,K=K, lambdaK = lambdaK, error = error, lambdaE = lambdaE, round = round)
  datA = data$dataA
  datB = data$dataB
  
  FS_gamma = FS_gamma(datA, datB, K)
  FS3 = FS3(datA, datB, K)
  FS_binary = FS_binary(datA, datB, K)
  
  return(c(FS_gamma, FS3, FS_binary))
}

doOne_exp_norm <- function(nA, nB, K, lambdaK, error, meanE, sdE, round ){
  source("1c_generate_data.R")
  source("2c_compare.R")
  source("3c_EM.R")
  source("4prc_methods.R")
  
  data = generate_data_exp_norm(nA = nA, nB=nB,K=K, lambdaK = lambdaK, error = error, meanE = meanE, sdE = sdE, round = round)
  datA = data$dataA
  datB = data$dataB
  
  FS_gamma_sqr = FS_gamma_sqr(datA, datB, K)
  FS_gamma = FS_gamma(datA, datB, K)
  FS_binary = FS_binary(datA, datB, K)
  FS3 = FS3(datA, datB, K)
  
  return(c(FS_gamma_sqr, FS_gamma, FS3, FS_binary))
}



vlis_exp_exp1 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 3),
    lambdaK = list(type="frozen", value = 0.02),
    error = list(type="frozen", value = 0.2),
    lambdaE = list(type="grid", value = c(1/2,1/3)), 
    round = list(type="frozen", value = TRUE)
  )
  return(vList)
}

vlis_exp_exp1_testK2 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K = list(type="grid", value = c(1,2)),
    lambdaK = list(type="frozen", value = 0.02),
    error = list(type="frozen", value = 0.2),
    lambdaE = list(type="frozen", value = 1/3), 
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

res1_exp_exp_prc_testK2 <- runSims(vlis_exp_exp1_testK2(nsim), doOne = doOne_exp_exp)
save(res1_exp_exp_prc_testK2, file="0415_res1_exp_exp_prc_testK2.RData")

