Rcpp::sourceCpp("Test/CL_cpp/CL_stat_cpp.cpp")
library(compiler)
enableJIT(3)
#x,y: data matrix
#FDR_level: pre-determined false discovery rate level
#norm.approx: TRUE : use normal approximation, (9) in CL.
#             FALSE: Use bootstrap, (12) in CL.

cor_CL<-function(x,y = NA, FDR_level = 0.05,norm.approx= TRUE, B = 2000, verbose=FALSE){
  n_x = dim(x)[1]
  n_y = dim(y)[1]
  p = dim(x)[2]
  
  onesample_flag <- is.na(y)
  ############################
  # Test statistic computation
  ####
  #Kurtosis like estimates 
  kappa_x = 1/(3*p) * sum(n_x*apply(x,2, function(x) sum((x-mean(x))^4))/apply(x,2, function(x) sum((x-mean(x))^2)^2))
  if(!onesample_flag){
    kappa_y = 1/(3*p) * sum(n_y*apply(y,2, function(x) sum((x-mean(x))^4))/apply(y,2, function(x) sum((x-mean(x))^2)^2))
  }
  # if(verbose) print(c(kappa_x,kappa_y))
  ##################
  # Need large amount of memory to conduct p by p matrix computation!
  
  # #Sample correlation
  # cor_x = cor(x)
  # cor_y = cor(y)
  # 
  # #Thresholded version of the sample correlation
  # tcor_x = cor_x*(abs(cor_x)/sqrt(kappa_x/n_x*(1-cor_x^2)) > threshold[1])
  # tcor_y = cor_x*(abs(cor_y)/sqrt(kappa_y/n_y*(1-cor_x^2)) > threshold[1])
  
  #p by p test statistics. Individual entries are asymptotically normal.
  # T_stat= (cor_x-cor_y)/sqrt(kappa_x/n_x*(1-tcor_x^2)^2 + kappa_y/n_y*(1-tcor_y^2)^2)
  if(onesample_flag){
    T_stat = CL_stat_onesample(x, kappa_x)
  }else{
    T_stat = CL_stat_twosample(x,y, kappa_x,kappa_y)
  }
  rownames(T_stat) = colnames(x)
  colnames(T_stat) = colnames(x)
  ########################
  # Find threshold
  if(norm.approx){
    #Using normal approximation (see equation 9 in Cai and Liu, 2016)
    t_hat = threshold(T_stat[upper.tri(T_stat)], alpha=FDR_level)
  }else{
    if(p>500) warning(paste("p is larger than 500, 
                            make sure there is enough memory for saving bootstrap results (need",
                      round((216+p*(p-1)/2*B)/1024^3,2)," GB)."))
    #Using bootstrap (see equation 12 in Cai and Liu, 2016)
    t_hat = threshold_B(T_stat[upper.tri(T_stat)], x, y, 
                        kappa_x, kappa_y, alpha=FDR_level, B = B)
  }
  adj_mat = (abs(T_stat)>t_hat)
  diag(adj_mat) <- FALSE #Set diagonal = FALSE in default
  
  return(list(stat= T_stat, threshold = t_hat, 
              fdr_level = FDR_level, decision = adj_mat, one.sample = onesample_flag))
}

# #This is the pure-R version of CL test based on the normal approximation
# cor_CL_R<-function(x,y, FDR_level = 0.05, verbose=FALSE){
#   #x,y : n by p matrix
#   #number of sample size
#   n_x = dim(x)[1]
#   n_y = dim(y)[1]
#   #number of features
#   p = dim(x)[2]
#   
#   threshold = 2*sqrt(log(p)/c(n_x,n_y))
#   if(verbose){
#     print(paste("n_x=",n_x, ", n_y=",n_y, ", p=",p, sep=""))
#     print(paste("threshold=", threshold[1], " ", threshold[2], sep=""))
#   }
#   ############################
#   # Test statistic computation
#   ####
#   #Kurtosis like estimates 
#   #See page 7 line 8.
#   kappa_x = 1/(3*p) * sum(n_x*apply(x,2, function(x) sum((x-mean(x))^4))/apply(x,2, function(x) sum((x-mean(x))^2)^2))
#   kappa_y = 1/(3*p) * sum(n_y*apply(y,2, function(x) sum((x-mean(x))^4))/apply(y,2, function(x) sum((x-mean(x))^2)^2))
#   if(verbose) print(c(kappa_x,kappa_y))
#   ##################
#   # Need large amount of memory to conduct p by p matrix computation!
#   
#   #Sample correlation
#   cor_x = cor(x)
#   cor_y = cor(y)
#   if(verbose) print(cor_x[1:5,1:5])
#   
#   #Thresholded version of the sample correlation
#   #See page 7 line 11
#   tcor_x = cor_x*(abs(cor_x)/sqrt(kappa_x/n_x*(1-cor_x^2)^2) >= threshold[1])
#   tcor_y = cor_y*(abs(cor_y)/sqrt(kappa_y/n_y*(1-cor_y^2)^2) >= threshold[2])
#   if(verbose){
#     print(pmax(tcor_x[1:5,1:5]^2,tcor_y[1:5,1:5]^2))
#     print(c(sum(tcor_x!=0),sum(tcor_y!=0)))
#   } 
#   #p by p test statistics. Individual entries are asymptotically normal.
#   # T_stat= (cor_x-cor_y)/sqrt(kappa_x/n_x*(1-tcor_x^2)^2 + kappa_y/n_y*(1-tcor_y^2)^2)
#   T_stat= (cor_x-cor_y)/sqrt((kappa_x/n_x + kappa_y/n_y)*(1-pmax(tcor_x^2,tcor_y^2))^2)
#   # T_stat = CL_stat_cpp(x,y, kappa_x,kappa_y)
#   ########################
#   # Find threshold
#   t_hat = threshold(T_stat[upper.tri(T_stat)], alpha=FDR_level)
#   
#   adj_mat = (abs(T_stat)>t_hat)
#   
#   return(list(stat= T_stat, threshold = t_hat, fdr_level = FDR_level, decision = adj_mat))
# }

#Threshold finder by normal approximation
threshold<-function(stat, alpha=0.05){
  #length of stat should be p(p-1)/2 = l, p = (1/2) * (1+sqrt(1+8l))
  p = 0.5*(1+sqrt(1+8*length(stat))) 
  t_ub = sqrt(4*log(p)-2*log(log(p)))
  #If there exist
  t = tryCatch(uniroot(function(t,p,alpha) (2*pnorm(t,lower.tail = FALSE)*p*(p-1)/2)/
                max(sum(abs(stat)>=t),1) - alpha, p = p, alpha = alpha,
              interval= c(0, t_ub))$root,
              error = function(e) return(NA))
  if(is.na(t)){
    return(sqrt(4*log(p)))
  }else{
    return(t)
  }
}

#Threshold finder by bootstrap
threshold_B<-cmpfun(function(stat, x, y, kappa_x, kappa_y, 
                      alpha=0.05, B = 2000){
  two_sample_flag = !is.na(y)
  #length of stat should be p(p-1)/2 = l, p = (1/2) * (1+sqrt(1+8l))
  p = 0.5*(1+sqrt(1+8*length(stat))) 
  t_ub = sqrt(4*log(p)-2*log(log(p)))
  #Compute bootstrap samples
  b_stat = matrix(NA, length(stat), B)
  
  n_x = dim(x)[1]
  if(two_sample_flag){
    n_y = dim(y)[1]
  }
  upper_diag_elem = upper.tri(matrix(0,p,p))
  if(two_sample_flag){
    for(b in 1:B){
      ind_x = sample(1:n_x, replace = TRUE)
      ind_y = sample(1:n_y, replace = TRUE)
      
      b_stat_temp = CL_stat_twosample(x[ind_x,],y[ind_y,], kappa_x,kappa_y)
      b_stat[,b]= b_stat_temp[upper_diag_elem]
    }
  }else{
    for(b in 1:B){
      ind_x = sample(1:n_x, replace = TRUE)
      
      b_stat_temp = CL_stat_onesample(x[ind_x,], kappa_x)
      b_stat[,b]= b_stat_temp[upper_diag_elem]
    }
  }
  
  #If there exist threshold
  t = tryCatch(uniroot(function(t,B,p,alpha, abs_stat,abs_b_stat){
                          ((2/(B*p*(p-1))* sum(abs_b_stat>t,na.rm=TRUE))*p*(p-1)/2)/max(sum(abs_stat>=t),1) - alpha
    }, B = B,p = p, alpha = alpha,
    abs_stat = abs(stat), abs_b_stat = abs(b_stat),
                       interval= c(0, t_ub))$root,
               error = function(e) return(NA))
  if(is.na(t)){
    return(sqrt(4*log(p)))
  }else{
    return(t)
  }
})
