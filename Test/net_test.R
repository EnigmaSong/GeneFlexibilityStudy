#This is an implementation of Shao et al. (2022).
#The code is from version 1.2.10 of the personal R package ysong.
#Requires Rcpp 
net_test <- function(net1, net2, method = "edgeworth", alternative = "two.sided",
                     motif = "triangle", n_boot, subsample_sizes){
  if(!(motif %in% c("triangle", "vshape", "threestar"))) stop("motif must be one of triangle, vshape, or threestar")
  alternative <- alternative[1]
  if(!(alternative %in% c("two.sided", "less", "greater"))) stop("alternative must be one of two.sided, less, or greater")
  
  result<- switch(method,
                  edgeworth = net_test_edgeworth(net1,net2,motif, alternative),
                  subsample = net_test_subsample(net1,net2,motif, alternative, n_boot, subsample_sizes),
                  stop("Error: method must be one of edgeworth and subsample"))
  return(result)
}
net_test_edgeworth <- function(net1, net2, motif, alternative){
  if(!(motif %in% c("triangle", "vshape", "threestar"))) stop("motif must be one of triangle, vshape, or threestar")
  net1_args <- network_hashing(net1, motif)
  net2_args <- network_hashing(net2, motif)
  
  test_stat <- get_net_teststat(net1_args, net2_args)
  p_value <- get_pvalue_nettest(test_stat, net1_args, net2_args, alternative)
  
  result <- list(test_stat = test_stat,
                 p_value = p_value,
                 method = "edgeworth",
                 alternative = alternative)
  result <- structure(result, class="net_test")
  return(result)
}
net_test_subsample <- function(net1, net2, motif,  
                               alternative, n_boot, subsample_sizes){
  if(!(motif %in% c("triangle", "vshape", "threestar"))) stop("motif must be one of triangle, vshape, threestar")
  if(any(missing(subsample_sizes))){
    warning(paste0("subsample_sizes is not specified."))
    subsample_sizes <- floor(sqrt(c(n1,n2)))
  }
  if(length(subsample_sizes)==1){ 
    warning(paste("length(subsample_sizes) == 1: Use", subsample_sizes, "for subsample_sizes of two networks."))
    subsample_sizes <- rep(subsample_sizes, 2)
  }
  
  n1 <- dim(net1)[1]
  n2 <- dim(net2)[2]
  
  
  net1_args <- network_hashing(net1, motif)
  net2_args <- network_hashing(net2, motif)
  test_stat <- get_net_teststat_subsample(net1_args, net2_args)
  
  boot_stat <- rep(NA, n_boot)
  for(i in 1:n_boot){
    net1_subind <- sample.int(n1, subsample_sizes[1])
    net2_subind <- sample.int(n2, subsample_sizes[2])
    
    net1_args <- network_hashing(net1[net1_subind,net1_subind], motif)
    net2_args <- network_hashing(net2[net2_subind,net2_subind], motif)
    
    boot_stat[i] <- get_net_teststat_subsample(net1_args, net2_args)
  }
  result <- list(test_stat = test_stat,
                 p_value = get_boot_pvalue_nettest(test_stat,boot_stat, alternative),
                 method = "subsample",
                 alternative = alternative,
                 n_boot = n_boot,
                 subsample_sizes = subsample_sizes,
                 boot_stat = boot_stat)
  attr(result, "hidden") <- c("n_boot", "subsample_sizes","boot_stat")
  result <- structure(result, class="net_test")
  return(result)
}

get_alpha_fun_args <- function(net, motif){
  num_node <- dim(net)[1]
  rho_hat <- get_edge_density(net)
  motif_num_nodes <- switch(motif,
                            vshape = 3,
                            triangle = 3,
                            threestar = 4)
  motif_num_edges <- switch(motif,
                            vshape = 2,
                            triangle = 3,
                            threestar = 4)
  net_motif_moment <- count_motif(net, motif)/max_count_motif(num_node, motif)
  
  #(S11)
  g_A1_hat <- get_g_A1_hat(net, motif, net_motif_moment, num_node, motif_num_nodes)
  #(S12)
  g_rhoA1_hat <- get_g_rhoA1_hat(net, rho_hat, num_node)
  #(S13)
  g_A2_hat <- get_g_A2_hat(net, motif, net_motif_moment, 
                           num_node, motif_num_nodes, 
                           g_A1_hat)
  #(S14)
  g_rhoA2_hat <- get_g_rhoA2_hat(net, rho_hat, g_rhoA1_hat)
  #(S15)
  xi_A_square <- mean(g_A1_hat^2)
  xi_A_rhoA <- mean(g_A1_hat*g_rhoA1_hat)
  xi_rhoA_square <- mean(g_rhoA1_hat^2)
  
  #merge all args
  alpha_funs_args <- list(num_node = num_node,
                          rho_hat_inv = ginv(rho_hat),
                          r = motif_num_nodes,
                          s = motif_num_edges,
                          net_motif_moment = net_motif_moment,
                          g_A1_hat = g_A1_hat,
                          g_rhoA1_hat = g_rhoA1_hat,
                          g_A2_hat = g_A2_hat,
                          g_rhoA2_hat = g_rhoA2_hat,
                          xi_A_square = xi_A_square,
                          xi_A_rhoA = xi_A_rhoA,
                          xi_rhoA_square = xi_rhoA_square)
  
  return(alpha_funs_args)
}

network_hashing<-function(net, motif){
  alpha_funs_args <- get_alpha_fun_args(net, motif)
  
  #get alphas
  alpha0 <- get_alpha0(alpha_funs_args$r,
                       alpha_funs_args$s,
                       alpha_funs_args$rho_hat_inv, 
                       alpha_funs_args$net_motif_moment,
                       alpha_funs_args$xi_rhoA_square, 
                       alpha_funs_args$xi_A_rhoA)
  alpha1 <- get_alpha1(alpha_funs_args$r,
                       alpha_funs_args$s,
                       alpha_funs_args$rho_hat_inv, 
                       alpha_funs_args$net_motif_moment,
                       alpha_funs_args$g_A1_hat, 
                       alpha_funs_args$g_rhoA1_hat)
  alpha2 <- get_alpha2(alpha_funs_args$r,
                       alpha_funs_args$s,
                       alpha_funs_args$rho_hat_inv, 
                       alpha_funs_args$net_motif_moment,
                       alpha_funs_args$g_A1_hat, 
                       alpha_funs_args$g_A2_hat,
                       alpha_funs_args$g_rhoA1_hat, 
                       alpha_funs_args$g_rhoA2_hat)
  alpha3 <- get_alpha3(alpha_funs_args$r,
                       alpha_funs_args$s,
                       alpha_funs_args$rho_hat_inv, 
                       alpha_funs_args$net_motif_moment,
                       alpha_funs_args$g_A1_hat, 
                       alpha_funs_args$g_rhoA1_hat, 
                       alpha_funs_args$xi_A_square,
                       alpha_funs_args$xi_A_rhoA,
                       alpha_funs_args$xi_rhoA_square)
  alpha4 <- get_alpha4(alpha_funs_args$r,
                       alpha_funs_args$s,
                       alpha_funs_args$rho_hat_inv, 
                       alpha_funs_args$net_motif_moment,
                       alpha_funs_args$g_A1_hat, 
                       alpha_funs_args$g_A2_hat,
                       alpha_funs_args$g_rhoA1_hat, 
                       alpha_funs_args$g_rhoA2_hat)
  
  #moments of alphas
  xi_alpha_square <- mean(alpha1^2)
  e_a1_cubic <- mean(alpha1^3)
  e_a1a3 <- mean(alpha1*alpha3)
  e_a4a1 <- .E_a4a1(alpha1,alpha4)
  e_a1a1a2 <- .E_uutA_uppertri(alpha1, alpha2)
  
  return(structure(list(num_node = alpha_funs_args$num_node,
                        rho_hat_inv = alpha_funs_args$rho_hat_inv,
                        r = alpha_funs_args$r,
                        s = alpha_funs_args$s,
                        net_motif_moment = alpha_funs_args$net_motif_moment,
                        alpha0 = alpha0,
                        xi_alpha_square = xi_alpha_square,
                        e_a1_cubic = e_a1_cubic,
                        e_a1a3 = e_a1a3,
                        e_a4a1 = e_a4a1,
                        e_a1a1a2 = e_a1a1a2
  ),
  class = "net_test_expterms")
  )
}


get_edge_density<-function(adj_mat){
  n <- dim(adj_mat)[1]
  return(sum(adj_mat)/(n*(n-1)))
}
get_alpha0<- function(r,s,rho_hat_inv, net_motif_moment,
                      xi_rhoA_square, xi_A_rhoA){
  #(S1)
  alpha0_hat <- 2*s*(s+1)*(rho_hat_inv^(s+2))*net_motif_moment*xi_rhoA_square
  alpha0_hat <- alpha0_hat - 2*r*s*(rho_hat_inv^(s+1))*xi_A_rhoA
  
  return(alpha0_hat)
}
get_alpha1<-function(r,s,rho_hat_inv, net_motif_moment,
                     g_A1_hat, g_rhoA1_hat){
  alpha1_hat <- (rho_hat_inv^s)*(r*g_A1_hat - 2*(s*rho_hat_inv)*net_motif_moment*g_rhoA1_hat)
  return(alpha1_hat)
}
get_alpha2<-function(r,s,rho_hat_inv, net_motif_moment,
                     g_A1_hat, g_A2_hat,
                     g_rhoA1_hat, g_rhoA2_hat){
  ##(S3)
  alpha2_hat <- (r*(r-1)/2)*rho_hat_inv^(s)*(g_A2_hat)
  alpha2_hat <- alpha2_hat - s*rho_hat_inv^(s+1)*net_motif_moment*g_rhoA2_hat
  alpha2_hat <- alpha2_hat + 2*rho_hat_inv^(s+2)*s*(s+1)*net_motif_moment*outer(g_rhoA1_hat,g_rhoA1_hat,"*")
  alpha2_hat <- alpha2_hat - 2*rho_hat_inv^(s+1)*s*r*outer(g_rhoA1_hat,g_A1_hat,"*")
  
  return(alpha2_hat)
}
get_alpha3<-function(r,s,rho_hat_inv, net_motif_moment,
                     g_A1_hat, g_rhoA1_hat, 
                     xi_A_square,xi_A_rhoA,xi_rhoA_square){
  #(S4)
  alpha3_hat <- -4*r^2*s*rho_hat_inv^(2*s+1)*xi_A_square*g_rhoA1_hat
  alpha3_hat <- alpha3_hat + r^2*(rho_hat_inv^(2*s))*(g_A1_hat^2-xi_A_square)
  alpha3_hat <- alpha3_hat - 16*(s^2)*(s+1)*(rho_hat_inv^(2*s+3))*(net_motif_moment^2)*xi_rhoA_square*g_rhoA1_hat
  alpha3_hat <- alpha3_hat + 8*r*(s^2)*(rho_hat_inv^(2*s+2))*net_motif_moment*g_A1_hat*xi_rhoA_square
  alpha3_hat <- alpha3_hat + 4*(s^2)*(rho_hat_inv^(2*s+2))*(net_motif_moment^2)*((g_rhoA1_hat^2)-xi_rhoA_square)
  alpha3_hat <- alpha3_hat - 4*r*s*(-(4*s+2)*(rho_hat_inv^(2*s+2))*g_rhoA1_hat*net_motif_moment*xi_A_rhoA
                                    + r*(rho_hat_inv^(2*s+1))*g_A1_hat*xi_A_rhoA
                                    + (rho_hat_inv^(2*s+1))*net_motif_moment*(g_A1_hat*g_rhoA1_hat - xi_A_rhoA)
  )
  return(alpha3_hat)
}
get_alpha4 <- function(r,s,rho_hat_inv, net_motif_moment,
                       g_A1_hat, g_A2_hat,
                       g_rhoA1_hat, g_rhoA2_hat){
  #(S5)
  alpha4_hat <- 2*(r^2)*(r-1)*(rho_hat_inv^(2*s))*g_A1_hat*g_A2_hat
  alpha4_hat <- alpha4_hat + 8*(s^2)*(rho_hat_inv^(2*s+2))*(net_motif_moment^2)*g_rhoA1_hat*g_rhoA2_hat
  alpha4_hat <- alpha4_hat - 4*r*(r-1)*s*(rho_hat_inv^(2*s+1))*net_motif_moment*g_rhoA1_hat*g_A2_hat
  alpha4_hat <- alpha4_hat - 4*r*s*(rho_hat_inv^(2*s+1))*net_motif_moment*g_A1_hat*g_rhoA2_hat
  
  return(alpha4_hat)
}
get_g_A1_hat <- function(net, motif, net_motif_moment, num_nodes, r){
  result <- count_rooted_motif(net, motif)/max_root_count_motif(num_nodes, motif) - net_motif_moment
  return(result)
}
get_g_A2_hat <- function(net, motif, net_motif_moment, num_nodes, r, g_A1_hat){
  result <- count_rooted_2V_motif(net, motif)/max_root_2V_count_motif(num_nodes, motif)
  result <- result - outer(g_A1_hat, 
                           g_A1_hat, 
                           "+")
  result <- result - net_motif_moment
  diag(result) <- 0
  return(result)
}

get_g_rhoA1_hat <- function(net, rho_hat, num_nodes){
  return(colSums(net)/(num_nodes-1) - rho_hat)
}
get_g_rhoA2_hat <- function(net, rho_hat, g_rhoA1_hat){
  result <- net - outer(g_rhoA1_hat,g_rhoA1_hat,"+") - rho_hat
  diag(result) <- 0
  return(result)
}
get_net_teststat <- function(net1_args, net2_args, c_delta = 0.01){
  motif_num_edges <- net1_args$s
  net1_motif_moment <- net1_args$net_motif_moment
  net2_motif_moment <- net2_args$net_motif_moment
  rho_hat_inv1 <- net1_args$rho_hat_inv
  rho_hat_inv2 <- net2_args$rho_hat_inv
  net1_num_node <- net1_args$num_node
  net2_num_node <- net2_args$num_node
  xi_alpha_square <- net1_args$xi_alpha_square
  xi_beta_square <- net2_args$xi_alpha_square
  
  D_mn_hat <- net1_motif_moment*(rho_hat_inv1^motif_num_edges) - net2_motif_moment*(rho_hat_inv2^motif_num_edges)
  S_mn_hat <- sqrt(xi_alpha_square/net1_num_node + xi_beta_square/net2_num_node)
  test_stat <- ginv(S_mn_hat)*(D_mn_hat) + 
    rnorm(1, 0, sqrt(c_delta*(log(net1_num_node)/net1_num_node+log(net2_num_node)/net2_num_node)))
  
  return(test_stat)
}
get_net_teststat_subsample <- function(net1_args, net2_args){
  motif_num_edges <- net1_args$s
  net1_motif_moment <- net1_args$net_motif_moment
  net2_motif_moment <- net2_args$net_motif_moment
  rho_hat_inv1 <- net1_args$rho_hat_inv
  rho_hat_inv2 <- net2_args$rho_hat_inv
  net1_num_node <- net1_args$num_node
  net2_num_node <- net2_args$num_node
  xi_alpha_square <- net1_args$xi_alpha_square
  xi_beta_square <- net2_args$xi_alpha_square
  
  D_mn_hat <- net1_motif_moment*(rho_hat_inv1^motif_num_edges) - net2_motif_moment*(rho_hat_inv2^motif_num_edges)
  S_mn_hat <- sqrt(xi_alpha_square/net1_num_node + xi_beta_square/net2_num_node)
  test_stat <- ginv(S_mn_hat)*(D_mn_hat)
  
  return(test_stat)
}
get_boot_pvalue_nettest <- function(test_stat,boot_stat, alternative){
  if(any(is.nan(boot_stat)|is.na(boot_stat))){
    warning(paste0(sum(is.nan(boot_stat))," NaN or NA produced from bootstrap"))
  }
  p_value <- mean(test_stat>boot_stat, na.rm=TRUE)
  p_value <- switch(match(alternative, c("two.sided","less","greater")),
                    2*min(p_value, 1-p_value),
                    p_value,
                    1-p_value
  )
  return(p_value)
}

get_pvalue_nettest <- function(test_stat, net1_args, net2_args, alternative){
  p_value <- p_Gmn(test_stat, net1_args, net2_args)
  p_value <- max(0, min(1, p_value))
  p_value <- switch(match(alternative, c("two.sided","less","greater")),
                    2*min(p_value, 1-p_value),
                    p_value,
                    1-p_value
  )
  
  return(p_value)
}

ginv <- function(x){
  non_zero_ind <- (x!=0)
  x[non_zero_ind] <- 1/x[non_zero_ind]
  return(x)
}

p_Gmn<- function(x, net1_args, net2_args){
  net1_motif_moment <- net1_args$net_motif_moment
  net2_motif_moment <- net2_args$net_motif_moment
  rho_hat_inv1 <- net1_args$rho_hat_inv
  rho_hat_inv2 <- net2_args$rho_hat_inv
  net1_num_node <- net1_args$num_node
  net2_num_node <- net2_args$num_node
  alpha0 <- net1_args$alpha0
  beta0 <- net2_args$alpha0
  
  xi_alpha_square <- net1_args$xi_alpha_square
  xi_beta_square <- net2_args$xi_alpha_square
  net1_e_a1_cubic = net1_args$e_a1_cubic
  net2_e_a1_cubic = net2_args$e_a1_cubic
  net1_e_a1a3 = net1_args$e_a1a3
  net2_e_a1a3 = net2_args$e_a1a3
  net1_e_a4a1 = net1_args$e_a4a1
  net2_e_a4a1 = net2_args$e_a4a1
  net1_e_a1a1a2 = net1_args$e_a1a1a2
  net2_e_a1a1a2 = net2_args$e_a1a1a2
  
  S_mn_hat <- sqrt(xi_alpha_square/net1_num_node 
                   + xi_beta_square/net2_num_node)
  S_mn_hat_inv <- ginv(S_mn_hat)
  
  I_zero <- S_mn_hat_inv*(alpha0/net1_num_node - beta0/net2_num_node)
  Q_mn_first <- (- net1_num_node^(-2)*(net1_e_a4a1+net1_e_a1a3)
                 + net2_num_node^(-2)*(net2_e_a4a1+net2_e_a1a3))
  Q_mn_first <- 0.5*(S_mn_hat_inv^3)*Q_mn_first
  
  Q_mn_second_first <- (net1_num_node^(-2))*(net1_e_a1_cubic/6 + net1_e_a1a1a2)
  Q_mn_second_first <- Q_mn_second_first - (net2_num_node^(-2))*(net2_e_a1_cubic/6 + net2_e_a1a1a2) 
  Q_mn_second_first <- (S_mn_hat_inv^3)*Q_mn_second_first
  
  Q_mn_second_second <- -(net1_num_node^(-3)*xi_alpha_square + 
                            net1_num_node^(-2)*net2_num_node^(-1)*xi_beta_square)*(net1_e_a4a1 + net1_e_a1a3)
  Q_mn_second_second <- Q_mn_second_second + (net2_num_node^(-3)*xi_beta_square + 
                                                net2_num_node^(-2)*net1_num_node^(-1)*xi_alpha_square)*(net2_e_a4a1 + net2_e_a1a3)
  Q_mn_second_second <- (0.5*S_mn_hat_inv^5)*Q_mn_second_second
  Q_mn_second <- Q_mn_second_first + Q_mn_second_second
  
  res <- pnorm(x) - dnorm(x) * (Q_mn_first + Q_mn_second*(x^2-1) + I_zero)
  
  return (res)
}

