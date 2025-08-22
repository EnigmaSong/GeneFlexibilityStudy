####
# Load test source
source("Test/cor_CL.R")

# library(WGCNA)
library(openxlsx)
library(dplyr)
library(foreach)
library(doSNOW)
library(parallel)

# Setting ----
file.save = FALSE #Save results in a file
Use_hashed_value <- TRUE #Network hasing and use hashed value for test
n_cores <- 4 # Number of cores to use for parallel processing
####
# Load sub-modules

#Load behaviral data
source("../DE_github/Data_preprocessing/pre_processing_genome.R")

####
# Step 0: Pre-processing----
# Follow WGCNA pre-processing to reduce number of genes used in covariance test
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-03-Preprocessing.pdf

###
# Load the count and behaviroal data
#Data pre-processing for Guppy data
#Network inference
## Data loading
## Load normalized counts for genes used in DE analysis
norm_counts_Aripo_genome = read.csv("../outputs/norm_counts_aripo_genome.csv",row.names = 1)
norm_counts_Quare_genome = read.csv("../outputs/norm_counts_quare_genome.csv",row.names = 1)

options(stringsAsFactors = FALSE)

# Step 0-1: Identification of outlying sample
# The tutorial keeps the samples containing less than 500 missing entries, 
# i.e. in the transcriptome list, no missing data. Hence, we can skip this step.
# meanExpressionByArray=apply(counts_genome_Aripo[rownames(norm_counts_Aripo_genome),],2,mean, na.rm=T)
# NumberMissingByArray=apply( is.na(data.frame(counts_genome_Aripo[rownames(norm_counts_Aripo_genome),])),2, sum)

# Step 0-2: Handling missing data and zero variance in probe profiles
# No missing in the count matrix
# Filter the low variance.
norm_counts_Aripo = norm_counts_Aripo_genome
norm_counts_Quare = norm_counts_Quare_genome

var_Aripo=as.vector(apply(norm_counts_Aripo,1,var, na.rm=T)) 
var_Quare=as.vector(apply(norm_counts_Quare,1,var, na.rm=T)) 

KeepGenes_Aripo = (var_Aripo>100)#13447
KeepGenes_Quare = (var_Quare>100)#14379

names(KeepGenes_Aripo)<- rownames(norm_counts_Aripo)
names(KeepGenes_Quare)<- rownames(norm_counts_Quare)

norm_counts_Aripo_Keep = norm_counts_Aripo[KeepGenes_Aripo,]
norm_counts_Quare_Keep = norm_counts_Quare[KeepGenes_Quare,]

####
#Load DE analysis results from GLMM
Aripo_DE_g_042920 = read.xlsx("../outputs/DE_GLMM_genome_042920.xlsx",sheet = 2)
Quare_DE_g_042920 = read.xlsx("../outputs/DE_GLMM_genome_042920.xlsx",sheet = 3)

print(head(Aripo_DE_g_042920))

####
# Load Id, pop, rear, interaction DE flags from DE analysis results file
#Aripo
Aripo_filtered_summary=inner_join(Aripo_DE_g_042920[,c(1,3,7,11)], 
                                  data.frame("Id"=names(KeepGenes_Aripo), "Kept"=KeepGenes_Aripo,
                                             stringsAsFactors = FALSE), 
                                  by = "Id")
###
# Num of genes (used vs unused, DE vs NDE)
# population main effect
table("pop_DE"=Aripo_filtered_summary$pop_DE_flag,
      "Kept"=Aripo_filtered_summary$Kept)
# rearing main effect
table("rear_DE"=Aripo_filtered_summary$rear_DE_flag,
      "Kept"=Aripo_filtered_summary$Kept)
# Interaction effect
table("inter_DE"=Aripo_filtered_summary$inter_DE_flag,
      "Kept"=Aripo_filtered_summary$Kept)

#Quare
Quare_filtered_summary=inner_join(Quare_DE_g_042920[,c(1,3,7,11)], 
                                  data.frame("Id"=names(KeepGenes_Quare), 
                                             "Kept"=KeepGenes_Quare,
                                             stringsAsFactors = FALSE), 
                                  by = "Id")
# population main effect
table("pop_DE"=Quare_filtered_summary$pop_DE_flag,
      "Kept"=Quare_filtered_summary$Kept)
# rearing main effect
table("rear_DE"=Quare_filtered_summary$rear_DE_flag,
      "Kept"=Quare_filtered_summary$Kept)
# Interaction effect
table("inter_DE"=Quare_filtered_summary$inter_DE_flag,
      "Kept"=Quare_filtered_summary$Kept)

####
# Step 1: Load residuals from csv file------
# Originally computed by 
#Load residuals from lmer results
resid_Aripo = read.table(file= "../outputs/network/residual_aripo.csv",sep=",", row.names=1, header=T) %>% as.matrix()
resid_Quare = read.table(file= "../outputs/network/residual_quare.csv",sep=",", row.names=1, header=T) %>% as.matrix()

print(head(resid_Aripo))

####
# Save the size of each group for future usage
length_LPP_Aripo = length(LPP_Aripo)
length_LPNP_Aripo = length(LPNP_Aripo)
length_HPP_Aripo = length(HPP_Aripo)
length_HPNP_Aripo = length(HPNP_Aripo)
length_LPP_Quare = length(LPP_Quare)
length_LPNP_Quare = length(LPNP_Quare)
length_HPP_Quare = length(HPP_Quare)
length_HPNP_Quare = length(HPNP_Quare)

####
##LC----
#Using all genes
na_ind = apply(resid_Aripo, 1, function(x) any(is.na(x)))#Only Aripo has the NA issue

# Step 5: network comparison based on Correlation "networks"  ------
#load extra libraries
library(ggplot2)
library(nethist)
library(ysong)

set.seed(2021)
#Population DE lists
Aripo_popkept_DE_id = setdiff(unlist(subset(Aripo_filtered_summary, select = "Id", subset = (pop_DE_flag == 1)&(Kept == 1))), names(which(na_ind)))
Quare_popkept_DE_id = unlist(subset(Quare_filtered_summary, select = "Id", subset = (pop_DE_flag == 1)&(Kept == 1)))

## Setup ----
# cor_fdr_levels = c(0.01,0.05,0.1)
fdr_level <- 0.05

### pairs to be compared----
#### Within each drainage, group-wise comparison of DE networks, those for NDE networks
sim_settings <- data.frame(
  drainage1 = c(rep("Aripo",8),rep("Quare",8)),
  group1 = rep(c("HPP", "HPNP", "HPP","LPP"),4),
  DE1 = rep(c(rep("DE",4),rep("NDE",4)),2),
  drainage2 = c(rep("Aripo",8),rep("Quare",8)),
  group2 = rep(c("HPNP", "LPNP", "LPP", "LPNP"),2),
  DE2 = rep(c(rep("DE",4),rep("NDE",4)),2)
)

#### Between drainage, in the same group and DE/NDE flags (e.g. Aripo LPP DE vs Quare LPP DE)
sim_settings <- rbind(sim_settings,
                      data.frame(drainage1 = c(rep("Aripo",8)),
                                 group1 = rep(c("HPP","HPNP","LPNP","LPP"),2),
                                 DE1 = c(rep("DE",4), rep("NDE",4)),
                                 drainage2 = c(rep("Quare",8)),
                                 group2 = rep(c("HPP","HPNP","LPNP","LPP"),2),
                                 DE2 = c(rep("DE",4), rep("NDE",4))
                      )
)
sim_settings$drain_group1 <- paste(sim_settings$drainage1,
                                   sim_settings$group1,sep=" ")
sim_settings$drain_group2 <- paste(sim_settings$drainage2,
                                   sim_settings$group2,sep=" ")
## Correlation network construction----
cor_nets <- vector("list",length = 8L)
message("Network construction by Cai and Liu (2016)")
begin_time <- Sys.time()
### Aripo HPP----
message("Aripo HPP")
### correlation network, group info, DE gene index info
cor_nets[[1]]$net <- cor_CL(x=t(resid_Aripo[!na_ind,HPP_Aripo]), FDR_level = fdr_level)$decision
cor_nets[[1]]$group <- "Aripo HPP"
cor_nets[[1]]$DE_ind <- match(Aripo_popkept_DE_id, rownames(cor_nets[[1]]$net))

###Aripo HPNP----
message("Aripo HPNP")
cor_nets[[2]]$net <- cor_CL(x=t(resid_Aripo[!na_ind,HPNP_Aripo]),
                            FDR_level = fdr_level)$decision
cor_nets[[2]]$group <- "Aripo HPNP"
cor_nets[[2]]$DE_ind <- match(Aripo_popkept_DE_id, rownames(cor_nets[[2]]$net))

###Aripo LPP----
message("Aripo LPP")
#two genes has 0 variance. excluding the as well
zero_var_ind <- (apply(resid_Aripo[,LPP_Aripo],1,var)==0)
cor_nets[[3]]$net <- cor_CL(x=t(resid_Aripo[!(na_ind|zero_var_ind),LPP_Aripo]), FDR_level = fdr_level)$decision
cor_nets[[3]]$group <- "Aripo LPP"
cor_nets[[3]]$DE_ind <- match(Aripo_popkept_DE_id, rownames(cor_nets[[3]]$net))

###Aripo LPNP----
message("Aripo LPNP")
cor_nets[[4]]$net <- cor_CL(x=t(resid_Aripo[!na_ind,LPNP_Aripo]), FDR_level = fdr_level)$decision
cor_nets[[4]]$group <- "Aripo LPNP"
cor_nets[[4]]$DE_ind <- match(Aripo_popkept_DE_id, rownames(cor_nets[[4]]$net))

###Quare HPP----
message("Quare HPP")
cor_nets[[5]]$net <- cor_CL(x=t(resid_Quare[,HPP_Quare]), FDR_level = fdr_level)$decision
cor_nets[[5]]$group <- "Quare HPP"
cor_nets[[5]]$DE_ind <- match(Quare_popkept_DE_id, rownames(cor_nets[[5]]$net))

###Quare HPNP----
message("Quare HPNP")
cor_nets[[6]]$net <- cor_CL(x=t(resid_Quare[,HPNP_Quare]), FDR_level = fdr_level)$decision
cor_nets[[6]]$group <- "Quare HPNP"
cor_nets[[6]]$DE_ind <- match(Quare_popkept_DE_id, rownames(cor_nets[[6]]$net))

###Quare LPP----
message("Quare LPP")
cor_nets[[7]]$nets <- cor_CL(x=t(resid_Quare[,LPP_Quare]), FDR_level = fdr_level)$decision
cor_nets[[7]]$group <- "Quare LPP"
cor_nets[[7]]$DE_ind <- match(Quare_popkept_DE_id, rownames(cor_nets[[7]]$net))

###Quare LPNP----
message("Quare LPNP")
cor_nets[[8]]$net <- cor_CL(x=t(resid_Quare[,LPNP_Quare]), FDR_level = fdr_level)$decision
cor_nets[[8]]$group <- "Quare LPNP"
cor_nets[[8]]$DE_ind <- match(Quare_popkept_DE_id, rownames(cor_nets[[8]]$net))
#End of network construction
print("Network construction running time: ")
print(round(difftime(Sys.time(),begin_time,units = "mins"),2))
print(paste0("length of DE_ind=", lapply(cor_nets, function(x) length(x$DE_ind))))
## Network Hashing----
motifs <- c("vshape", "triangle", "threestar")
if(Use_hashed_value){
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(cor_nets)*2, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  print("Network hashing start")
  begin_time <- Sys.time()
  ### Hashing the networks
  hash_nets <- foreach(i = 1:length(cor_nets),.options.snow = opts) %:%
    foreach(DE_flag = c("DE","NDE")) %dopar% {
      message(paste("Hashing network:", i, DE_flag))
      if(DE_flag=="DE"){
        gene_ind <- cor_nets[[i]]$DE_ind
      }else if(DE_flag=="NDE"){
        gene_ind <- -cor_nets[[i]]$DE_ind
      }
      # gene_ind <- 1:100 #for test
      net_hash <- list()  
      net_hash$group <- cor_nets[[i]]$group
      net_hash$DE_flag <- DE_flag
      net_hash$hashed$vshape <- ysong:::network_hashing(cor_nets[[i]]$net[gene_ind,gene_ind], "vshape")
      net_hash$hashed$triangle <- ysong:::network_hashing(cor_nets[[i]]$net[gene_ind,gene_ind], "triangle")
      net_hash$hashed$threestar <- ysong:::network_hashing(cor_nets[[i]]$net[gene_ind,gene_ind], "threestar")
      
      return(net_hash)
  }
  #End of network hashing
  print("Network hashing running time: ")
  print(round(difftime(Sys.time(),begin_time,units = "mins"),2))
  stopCluster(cl)
  print(paste("Group matches?- ",all(lapply(hash_nets, function(x) print(x[[1]]$group==x[[2]]$group)))))
  rm(cor_nets)
}

##Test----
test_result <- data.frame(
  drain_group1 = sim_settings$drain_group1,
  DE_flag_group1 = sim_settings$DE1,
  drain_group2 = sim_settings$drain_group2,
  DE_flag_group2 = sim_settings$DE2,
  #vshape
  test_stat_vshape =  rep(NA,nrow(sim_settings)),
  p_value_vshape_twosided = rep(NA,nrow(sim_settings)),
  p_value_vshape_less = rep(NA,nrow(sim_settings)),
  p_value_vshape_greater = rep(NA,nrow(sim_settings)),
  #triangle
  test_stat_triangle =  rep(NA,nrow(sim_settings)),
  p_value_triangle_twosided = rep(NA,nrow(sim_settings)),
  p_value_triangle_less = rep(NA,nrow(sim_settings)),
  p_value_triangle_greater = rep(NA,nrow(sim_settings)),
  #threestar
  test_stat_threestar =  rep(NA,nrow(sim_settings)),
  p_value_threestar_twosided = rep(NA,nrow(sim_settings)),
  p_value_threestar_less = rep(NA,nrow(sim_settings)),
  p_value_threestar_greater = rep(NA,nrow(sim_settings))
)
begin_time <- Sys.time()
print("Test start")
if(Use_hashed_value){
  ### If we use hashed values----
  ## Takes shorter time
  for(i in 1:nrow(sim_settings)){
    ind1 <- match(sim_settings$drain_group1[i], 
                  sapply(hash_nets, function(x) x[[1]]$group))
    ind2 <- match(sim_settings$drain_group2[i], 
                  sapply(hash_nets, function(x) x[[1]]$group))
    if(sim_settings$DE1[i]=="DE"){
      hashed1 <- hash_nets[[ind1]][[1]]$hashed
    }else if(sim_settings$DE1[i]=="NDE"){
      hashed1 <- hash_nets[[ind1]][[2]]$hashed
    }
    if(sim_settings$DE2[i]=="DE"){
      hashed2 <- hash_nets[[ind2]][[1]]$hashed
    }else if(sim_settings$DE2[i]=="NDE"){
      hashed2 <- hash_nets[[ind2]][[2]]$hashed
    }
    
    ## vshape
    test_result$test_stat_vshape[i] <- ysong:::get_net_teststat(hashed1$vshape, hashed2$vshape)
    test_result$p_value_vshape_twosided[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_vshape[i], hashed1$vshape, hashed2$vshape, 
                                                     "two.sided")
    test_result$p_value_vshape_less[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_vshape[i], hashed1$vshape, hashed2$vshape, 
                                                                "less")
    test_result$p_value_vshape_greater[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_vshape[i], hashed1$vshape, hashed2$vshape, 
                                                                 "greater")
    
    ## triangle
    test_result$test_stat_triangle[i] <- ysong:::get_net_teststat(hashed1$triangle, hashed2$triangle)
    test_result$p_value_triangle_twosided[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_triangle[i], hashed1$triangle, hashed2$triangle, 
                                                     "two.sided")
    test_result$p_value_triangle_less[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_triangle[i], hashed1$triangle, hashed2$triangle, 
                                                                "less")
    test_result$p_value_triangle_greater[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_triangle[i], hashed1$triangle, hashed2$triangle, 
                                                                 "greater")
    ## threestar
    test_result$test_stat_threestar[i] <- ysong:::get_net_teststat(hashed1$threestar, hashed2$threestar)
    test_result$p_value_threestar_twosided[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_threestar[i], hashed1$threestar, hashed2$threestar, 
                                                     "two.sided")
    test_result$p_value_threestar_less[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_threestar[i], hashed1$threestar, hashed2$threestar, 
                                                                "less")
    test_result$p_value_threestar_greater[i] <- ysong:::get_pvalue_nettest(test_result$test_stat_threestar[i], hashed1$threestar, hashed2$threestar, 
                                                                 "greater")
  }
}else{
  ### If we use networks themselves ----
  ## Takes longer
  for(i in 1:nrow(sim_settings)){
    ind1 <- match(sim_settings$drain_group1[i],
                  sapply(cor_nets, function(x) x$group))
    ind2 <- match(sim_settings$drain_group2[i],
                  sapply(cor_nets, function(x) x$group))
    if(sim_settings$DE1[i]=="DE"){
      gene_ind1 <- cor_nets[[ind1]]$DE_ind
    }else if(sim_settings$DE1[i]=="NDE"){
      gene_ind1 <- -cor_nets[[ind1]]$DE_ind
    }
    if(sim_settings$DE2[i]=="DE"){
      gene_ind2 <- cor_nets[[ind2]]$DE_ind
    }else if(sim_settings$DE2[i]=="NDE"){
      gene_ind2 <- -cor_nets[[ind2]]$DE_ind
    }
    #two-sided test
    ## vshape
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "vshape", alternative = "two.sided")
    test_result$p_value_vshape_twosided[i] <- temp_res$p_value
    test_result$test_stat_vshape_twosided[i] <- temp_res$test_stat
    ## triangle
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                                cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "triangle", alternative = "two.sided")
    test_result$p_value_triangle_twosided[i] <- temp_res$p_value
    test_result$test_stat_triangle_twosided[i] <- temp_res$test_stat
    ## threestar
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                                cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "threestar", alternative = "two.sided")
    test_result$p_value_threestar_twosided[i] <- temp_res$p_value
    test_result$test_stat_threestar_twosided[i] <- temp_res$test_stat
  
    #less
    ## vshape
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "vshape", alternative = "less")
    test_result$p_value_vshape_less[i] <- temp_res$p_value
    test_result$test_stat_vshape_less[i] <- temp_res$test_stat
    ## triangle
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "triangle", alternative = "less")
    test_result$p_value_triangle_less[i] <- temp_res$p_value
    test_result$test_stat_triangle_less[i] <- temp_res$test_stat
    ## threestar
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "threestar", alternative = "less")
    test_result$p_value_threestar_less[i] <- temp_res$p_value
    test_result$test_stat_threestar_less[i] <- temp_res$test_stat
    #greater
    ## vshape
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "vshape", alternative = "greater")
    test_result$p_value_vshape_greater[i] <- temp_res$p_value
    test_result$test_stat_vshape_greater[i] <- temp_res$test_stat
    ## triangle
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "triangle", alternative = "greater")
    test_result$p_value_triangle_greater[i] <- temp_res$p_value
    test_result$test_stat_triangle_greater[i] <- temp_res$test_stat
    ## threestar
    temp_res <- net_test(cor_nets[[ind1]]$net[gene_ind1,gene_ind1],
                         cor_nets[[ind2]]$net[gene_ind2,gene_ind2],
                         motif = "threestar", alternative = "greater")
    test_result$p_value_threestar_greater[i] <- temp_res$p_value
    test_result$test_stat_threestar_greater[i] <- temp_res$test_stat
  }
}
print("Network Comparison test running time: ")
print(round(difftime(Sys.time(),begin_time,units = "mins"),2))

# #Save the results----
if(file.save){
  write.csv(file= "../outputs/network/For revision/net_test_extra_results.csv", x = test_result, row.names = F)
  save.image("../outputs/network/For revision/net_test_extra_results.RData")
}
print(test_result)
message("Finished!")
