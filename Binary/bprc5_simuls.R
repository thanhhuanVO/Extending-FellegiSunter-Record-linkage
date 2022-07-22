
library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)


doOne1 <- function(nA, nB, K, prev, error){
  source("b1_generate_data.R")
  source("b2_compare.R")
  source("b3_EM.R")
  source("bprc4_methods.R")
  
  nA = 500
  nB = 200
  K = 40
  error = 0.04
  prev = 0.2
  

  prevalence = rep(prev, K)
  data = generate_data(nA = nA, nB = nB, K = K, prevalence = prevalence, error = error)
  
  datA = data$dataA
  datB = data$dataB
  prev = data$prev
  
  g <- true_prob(datA, datB, K, prev, error)
  g_true = g$g_true
  indM = g$indM
  indU = g$indU
  comp_mat4 = g$comp_mat4
 
  
  FS <- FS(datA, datB, K, g_true, indM, indU, tol=1e-6, maxits = 500)
  FS3 <- FS3(datA, datB, K,  g_true, indM, indU, tol=1e-6, maxits = 500)
  FS4 <- FS4(datA, datB, K, comp_mat4, g_true, indM, indU, tol=1e-6, maxits = 500)
  Bayesian <- bayesian(datA, datB, K,  g_true, indM, indU)
  
  results <- c(FS, FS3,FS4, Bayesian)
  return(results)
}


vlis1 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="grid", value = c(40,50)),
    error = list(type="frozen", value = 0.04),
    prev = list(type="frozen", value = 0.2)
  )
  return(vList)
}


#####################################
runSims <- function(vList=vlis1(),doOne = doOne1, seedList=NULL){
  res <- simsalapar::doForeach(vList,  doOne = doOne,  cluster=makeCluster(8, type="PSOCK"), seed = seedList)
  return(res)
}

#setwd("C:\\Users\\thanhvo\\Dropbox\\R_program\\First_paper\\Binary\\Programs")
setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Programs")
nsim = 1000

res1 <- runSims(vlis1(nsim), doOne = doOne1)
save(res1, file="0402_1000s_res1_binary_prc.RData")






