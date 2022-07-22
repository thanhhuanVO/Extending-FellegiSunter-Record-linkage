
# Make error for one binary vector x
makeError <- function(x,error){
  #x is a column of B
  #error is a proportion of error
  nE  = round(length(x)*error)
  index = sample(1:length(x), nE)
  x[index] = 1-x[index]
  return(x)
}

generate_data <- function(nA, nB, K, prevalence, error, min_prev = 0.01){
  #prevalence = rep(c(0.05,0.1,0.2,0.3,0.4), K/5)

  # First database A
  datA = matrix(0, nrow = nA, ncol = K+1)
  datA[,K+1] = 1:nA #id
  
  conditionA = TRUE
  while (conditionA){
    datA[,1:K] = sapply(prevalence, function(x){rbinom(n=nA, size = 1,prob = x)})
    conditionA = (sum(colSums(datA[,1:K]/nA) >= min_prev) < K)
  }
  
  datA = data.frame(datA)
  colnames(datA)=c(paste("R", 1:K, sep = ""),"id") 
  
  conditionB = TRUE
  while (conditionB) {
    # Second database B
    idAB <- sample(1:nA,nB) #ident in A appearing in B
    datB <- datA[idAB,]
    
    # Make error for B
    datB[,1:K]= apply(datB[,1:K], MARGIN = 2, FUN = makeError, error = error)
    
    conditionB = (sum(colSums(datB[,1:K]) >= 1) < K)
  }
  
  
  return(list(dataA=datA, dataB = datB, prev = prevalence))
}

