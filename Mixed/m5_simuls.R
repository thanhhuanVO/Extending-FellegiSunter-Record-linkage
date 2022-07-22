library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(ludic)
library(maxLik)
library(Rfast)
library(mvtnorm)
library(mixtools)
library(EnvStats)
library(fastLink)
library(mixR)

# Function for parallel simulation runnings
doOne_mixed <- function(nA, nB, K1, K2, K3, error1, error2, error3, prev2, lambda3, mu3){
  library(tictoc)
  source("m1_generate_data.R")
  source("m2_compare.R")
  source("m3_EM.R")
  source("m4_methods.R")
  
  
  data = generate_data(nA = nA, nB = nB, K1 = K1, K2 = K2, K3 = K3,
                       error1 = error1, error2 = error2, error3 = error3, 
                       prev2 = prev2, lambda3 = lambda3, mu3 = mu3)
  
  datA = data$dataA
  datB = data$dataB
  
  tic()
  FS_mixed <- FS_mixed(datA, datB, K1, K2, K3)
  temp = toc(quiet = TRUE)
  time_FS_mixed <- temp$toc-temp$tic
  
  tic()
  FS_binary <- FS_binary(datA, datB, K1, K2, K3)
  temp = toc(quiet = TRUE)
  time_FS_binary = temp$toc - temp$tic
  
  tic()
  FS_RecordLinkage <- FS_RecordLinkage(datA, datB, K1, K2, K3)
  temp = toc(quiet = TRUE)
  time_FS_RecordLinkage = temp$toc - temp$tic
  
  results <- c(FS_mixed, FS_binary, FS_RecordLinkage, time_FS_mixed, time_FS_binary, time_FS_RecordLinkage)
  
  return(results)
}

# Simulation parameters
vlis1 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500),
    nB = list(type="frozen", value = 200),
    K1 = list(type="frozen", value = 2),
    K2 = list(type="frozen", value = 2),
    K3 = list(type="frozen", value = 2),
    error1 = list(type="frozen", value = 0.1),
    error2 = list(type="frozen", value = 0.04),
    error3 = list(type="grid", value = c(0.1,0.2,0.3)),
    prev2 = list(type="frozen", value = 0.01),
    lambda3 = list(type="frozen", value = 0.01),
    mu3 = list(type="frozen", value = 2)
  )
  return(vList)
}


runSims <- function(vList=vlis1(),doOne = doOne1, seedList=NULL){
  res <- simsalapar::doForeach(vList,  doOne = doOne, cluster = makeCluster(8, type="PSOCK"), seed = seedList)
  return(res)
}


nsim = 1000

res <- runSims(vlis1(nsim), doOne = doOne_mixed, seedList = 1:nsim)

