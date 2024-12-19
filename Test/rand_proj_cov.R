require(SetTest)
require(heplots)
require(equalCovs)

rand_proj_cov<- function(x,y, proj.dim = 3, method = 'LC', n_proj = 2000L){
  n1 = dim(x)[1]
  n2 = dim(y)[1]
  p = dim(x)[2]
  
  if(proj.dim==1){
    proj_vec = matrix(rnorm(p*n_proj*proj.dim), nrow = p, ncol = n_proj*proj.dim)
    proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
    
    x = scale(x, center=TRUE, scale=FALSE)
    y = scale(y, center=TRUE, scale=FALSE)
    XtR = x%*%proj_vec
    YtR = y%*%proj_vec
    
    F_stat = (colSums(XtR*XtR)/n1) / (colSums(YtR*YtR)/n2)
    p_val = pf(pmax(F_stat,1/F_stat), n1,n2, lower.tail = FALSE)*2
  
    T_stat = max(log(F_stat)/sqrt(2/n1+2/n2))
    var_F = (2*n2^2*(n1+n2-2))/(n1*(n2-2)^2*(n2-4))
    return(list(p.value = p_val,F.stat = F_stat, T.stat = T_stat, p.val_T = 1-pnorm(T_stat, sd=sqrt(var_F))^p))
  }else{
    x = scale(x, center=TRUE, scale=FALSE)
    y = scale(y, center=TRUE, scale=FALSE)

    temp_test_stats = rep(NA, n_proj)
    if(method %in% c("LC","LC_max")){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[1]
      }
      test_stat =  max(temp_test_stats)
    }else if(method == "boxM"){
      group= c(rep(1,n1),rep(2,n2))
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test_stats[i] = heplots::boxM(rbind(XtR, YtR),
                                           group)$statistics
      }
      test_stat =  max(temp_test_stats)
    }else if(method=='FC'){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        LC_pval = equalCovs::equalCovs(XtR, YtR, n1, n2)[2]
        
        #Use independent projection for CLX (2022-01-01)
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        CLX_pval = CLX_cpp(XtR,YtR)$pvalue
        temp_test_stats[i] = -2*log(LC_pval) - 2*log(CLX_pval)
      }
      test_stat =  max(temp_test_stats)
    }else if(method%in%c('LC_sq_max',"LC_max_sq")){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[1]
      }
      test_stat =  max(temp_test_stats^2)
    }else if(method%in%c('LC_mean')){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[1]
      }
      #Even we call it 'mean', it is not actual mean. sum(T_k)/sqrt(K)
      test_stat =  sum(temp_test_stats)/sqrt(n_proj)
    }else if(method%in%c('LC_abs_mean',"LC_mean_abs")){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[1]
      }
      #Even we call it 'mean', it is not actual mean. sum(T_k)/sqrt(K)
      test_stat =  sum(abs(temp_test_stats))/sqrt(n_proj)
    }else if(method %in% c('LC_sq_mean','LC_mean_sq')){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[1]
      }
      test_stat =  sum(temp_test_stats^2)/sqrt(n_proj)
    }else if(method %in% c('LC_avg_pvals')){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[2]
      }
      test_stat =  mean(temp_test_stats)
    }else if(method %in% c('LC_HC')){
      for(i in 1:n_proj){
        proj_vec = matrix(rnorm(p*proj.dim), nrow = p, ncol = proj.dim)
        proj_vec = apply(proj_vec, 2, function(x) x/sqrt(sum(x^2)))
        XtR = x%*%proj_vec
        YtR = y%*%proj_vec
        
        temp_test = equalCovs::equalCovs(XtR, YtR, n1, n2)
        temp_test_stats[i] = temp_test[2]
      }
      test_res =  SetTest::test.hc(temp_test_stats, M = diag(1,n_proj,n_proj), k0=1, k1=n_proj/2, LS = T, onesided = TRUE)
      return(list(p.value = test_res$pvalue, test.stat = test_res$hcstat))
    }
  }
  return(list(p.value = NA, test.stat = test_stat))
}


