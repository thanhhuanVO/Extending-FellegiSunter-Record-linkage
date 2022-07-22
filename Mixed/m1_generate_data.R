# These functions aim to generate two database for linkage with mix-typed matching variables (categorical, low prevalence binary and continuous)

###################### 1. Categorical
# This function can be used to generate a database of nA units, K categorical matching variables with 30/50 categories
generate_cate <- function(nA, K){
  # First database A
  datA = matrix(0, nrow = nA, ncol = K)
  for (k in (1:K)){
    n_cate <- sample(c(30,50), size = 1)
    datA[,k] = sample(1:n_cate, size = nA, replace = TRUE)
  }
  
  return(datA)
}

###################### 2. Low prevalence binary
# This function can be used to generate a dataset of nA units with K binary matching variables
# prevalence: proportion of low prevalence value
# min_prev: minimum bound for the frequency of low prevalence value (e.g. to avoid a column of all 0s)

generate_binary <- function(nA,  K, prevalence, min_prev){
  # First database A
  datA = matrix(0, nrow = nA, ncol = K)
  
  conditionA = TRUE
  while (conditionA){
    datA[,1:K] = sapply(prevalence, function(x){rbinom(n=nA, size = 1,prob = x)})
    if (K==1){
      conditionA = sum(datA[,1]) < 1
    }else{
      conditionA = (sum(colSums(datA[,1:K]/nA) >= min_prev) < K)
    }
  }
  
  return(datA)
}

######################### 3. Continuous (dates)
# This function is used to generate a dataset of nA units with K continuous matching variables 
# following exponential distribution with parameter lambdaK

generate_cont_exp = function(nA, K, lambdaK = 0.02){
  
  datA = matrix(ceiling(rexp(nA*K,rate = lambdaK)),ncol = K)

  return(datA)
}
######################### 4. Mixed data
# This function aggregates all to generate two database for linkage
# nA: number of units in database A
# nB: number of units in database B. We assume B belongs to A
# K1: number of categorical matching variables
# K2: number of low prevalence binary matching variables
# K3: number of continuous matching variables
# error1, error2, error3: proportion of errors for the K1, K2, K3 matching variables, respectively
# prev2: prevalence for K2 binary matching variables
# lambda3: K3 continuous matching variables ~ Exponential(lambda3)
# mu3: mean of time lag for continuous matching variable

generate_data = function(nA, nB, K1, K2, K3,  error1, error2,  error3, prev2, lambda3 , mu3){
  min_prev = 0.01
  #First database
  datA = matrix(0, nrow = nA, ncol = K1+K2+K3+1)
  
  datA[,K1+K2+K3+1] = 1:nA #id
  
  # K1 categorical variable
  datA[,1:K1] = generate_cate(nA, K1)
  # K2 low prevalence binary
  datA[,(K1+1):(K1+K2)] = generate_binary(nA, K2, prevalence = rep(prev2, K2), min_prev)
  # K3 dates
  datA[,(K1+K2+1):(K1+K2+K3)] = generate_cont_exp(nA, K3, lambda3 )
  
  datA = data.frame(datA)
  colnames(datA)=c(paste("R", 1:(K1+K2+K3), sep = ""),"id") 
  
  ## Second database B
  idAB <- sample(1:nA,nB) #ident in A appearing in B
  datB <- datA[idAB,]
  ## Making error for B
  makeError_cate <- function(x,error){
    #x is a column of B
    #error is a proportion of error
    nE  = round(length(x)*error)
    index = sample(1:length(x), nE)
    x[index] = sapply(index, function(i) sample(setdiff(1:length(unique(x)),x[i]), size = 1, replace = TRUE))
    return(x)
  }
  makeError_binary <- function(x,error){
    #x is a column of B
    #error is a proportion of error
    nE  = round(length(x)*error)
    index = sample(1:length(x), nE)
    x[index] = 1-x[index]
    return(x)
  }
  
  makeError_continuous <- function(x,error, mu){
    #x is a column of B
    #error is a proportion of error
    nE  = round(length(x)*error)
    index = sample(1:length(x), nE)
    x[index] = x[index] + ceiling(rexp(nE, rate = 1/mu))
    return(x)
  }
  
  ## Error for 1st part
  if (K1==1){
    datB[,1:K1] = makeError_cate(datB[,1], error1)
  }else{
    datB[,1:K1]= apply(datB[,1:K1], MARGIN = 2, FUN = makeError_cate, error = error1)
  }
  
  #Error for 2nd part
  conditionB = TRUE
  while (conditionB) {
    # Make error for B
    if (K2==1){
      datB[,(K1+1):(K1+K2)] = makeError_binary(datB[,(K1+1):(K1+K2)],error2)
      conditionB = sum(datB) < 1
    }else{
    datB[,(K1+1):(K1+K2)]= apply(datB[,(K1+1):(K1+K2)], MARGIN = 2, FUN = makeError_binary, error = error2)
    conditionB = (sum(colSums(datB[,(K1+1):(K1+K2)]) >= 1) < K2)
    }
  }
  
  #Error for 3rd part
  if (K3==1){
    datB[,K1+K2+1] = makeError_continuous(datB[,K1+K2+1], error = error3, mu = mu3)
  }else{
    datB[,(K1+K2+1):(K1+K2+K3)]= apply(datB[,(K1+K2+1):(K1+K2+K3)], MARGIN = 2, FUN = makeError_continuous, error = error3, mu = mu3)
  }
  
  return(list(dataA=datA, dataB = datB))
}



