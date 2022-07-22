
library(simsalapar)
library(doParallel)
library(clue)
library(matrixStats)
library(klaR)
library(ludic)


doOne1 <- function(nA, nB, K, prev, error){
  library(tictoc)
  source("b1_generate_data.R")
  source("b2_compare.R")
  source("b3_EM.R")
  source("b4_methods.R")
  

  prevalence = rep(prev, K)
  data = generate_data(nA = nA, nB = nB, K = K, prevalence = prevalence, error = error)
  
  datA = data$dataA
  datB = data$dataB

  tic()
  FS <- FS(datA, datB, K,  tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS <- temp$toc-temp$tic
  
  tic()
  FS3 <- FS3(datA, datB, K,  tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS3 <- temp$toc-temp$tic
  
  tic()
  FS4 <- FS4(datA, datB, K,  tol=1e-6, maxits = 500)
  temp = toc(quiet = TRUE)
  time_FS4 <- temp$toc-temp$tic
  
  tic()
  Bayesian <- bayesian(datA, datB, K)
  temp = toc(quiet = TRUE)
  time_Bayesian <- temp$toc-temp$tic
  
  results <- c(FS, FS3,FS4, Bayesian, time_FS, time_FS3, time_FS4, time_Bayesian)
  return(results)
}


vlis1 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="grid", value = c(30,40,50)),
    error = list(type="grid", value = c(0.02,0.04,0.06)),
    prev = list(type="frozen", value = 0.2)
  )
  return(vList)
}

vlis2 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="frozen", value = 500 ),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 40),
    error = list(type="frozen", value = 0.04),
    prev = list(type="grid", value = c(0.1,0.2,0.3))
  )
  return(vList)
}

vlis3 <- function(nsim=100) {
  vList <- simsalapar::varlist(
    n.sim = list(value = nsim, expr = quote(N[sim])), # , type = "N"
    nA = list(type="grid", value = c(400,800,1200) ),
    nB = list(type="frozen", value = 200),
    K = list(type="frozen", value = 40),
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
#setwd("C:\\Users\\thhvo.BCOM\\Dropbox\\R_program\\First_paper\\Binary\\Programs")


nsim = 100

res1 <- runSims(vlis1(nsim), doOne = doOne1, seedList = 1:nsim)
save(res1, file="0722_100s_res1_binary.RData")

#res2 <- runSims(vlis2(nsim), doOne = doOne1)
#save(res2, file="0324_1000s_res2_binary.RData")

#res3 <- runSims(vlis3(nsim), doOne = doOne1)
#save(res3, file="0324_1000s_res3_binary.RData")





