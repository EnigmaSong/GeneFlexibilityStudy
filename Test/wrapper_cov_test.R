## Wrapper for covariance tests

library(equalCovs) #R package by Jun Li et al.
Rcpp::sourceCpp("Test/CLX_cpp/CLX_cpp.cpp")
library(heplots)
library(covTestR)

#wrapper for equalCovs()
LC_wrapper<- function(X,Y){
  n_X = dim(X)[1]
  n_Y = dim(Y)[1]
  res= equalCovs(X,Y, n_X, n_Y)
  
  return(list(statistic = res[1],
              p.value = res[2]))
}
#wrapper for CLX_cpp()
CLX_wrapper<- function(X,Y){
  res= CLX_cpp(X,Y)
  
  return(list(statistic = res$TSvalue,
              p.value = res$pvalue))
}

#wrapper for boxM()
BoxM_wrapper<- function(X,Y){
  n_X = dim(X)[1]
  n_Y = dim(Y)[1]
  group = c(rep(1, n_X), rep(2,n_Y))
  
  res = boxM(rbind(X,Y), group)
  return(list(statistic = res$statistic,
              p.value = res$p.value))
}

#wrapper for Sristava()
Srivastava_wrapper <- function(X,Y){
  res=Srivastava2007(list(X,Y))
  
  return(list(statistic = res$statistic,
              p.value = res$p.value))
}

#wrapper for Schott2001()
Schott2001_wrapper <- function(X,Y){
  res=Schott2001(list(X,Y))
  
  return(list(statistic = res$statistic,
              p.value = res$p.value))
}

#wrapper for Schott2001()
Schott2007_wrapper <- function(X,Y){
  res=Schott2007(list(X,Y))
  
  return(list(statistic = res$statistic,
              p.value = res$p.value))
}