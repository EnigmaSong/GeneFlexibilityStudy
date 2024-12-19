#Code for Fischer et al. (2024)
getwd()
####
# Step 0: Load libraries, pre-coded functions, settings ----

## Load Libraries ----
library(openxlsx) #Load/save Excel files
library(dplyr) 
library(ggplot2)
# install.packages("devtools")
# devtools::install_github("EnigmaSong/nethist")
library(nethist) #For violin summary plot (Maugis et al., 2017)
library(Rcpp)
library(RcppArmadillo)

## Load source files for hypothesis tests----
source("Test/wrapper_cov_test.R")
source("Test/rand_proj_cov.R")
source("Test/cor_CL.R")
source("Test/net_test.R")
source("Test/count_motif_utils.R")
sourceCpp("src/net_test_util.cpp")

## Setting ----
file.save = F #Save results in a file if T
###
### parameter setups for random projection test ----
rand_project_test_stat = c("LC_max")
q = 20 #dimension of projection
K = 500 #Number of projections

###
## Load the count and behavioral data ----
# 
### Data pre-processing from Fischer et al. (2021) ----
source("../DE_github/Data_preprocessing/pre_processing_genome.R")

## Data loading
### Load normalized counts for genes used in DE analysis ----
norm_counts_Aripo_genome = read.csv("../outputs/norm_counts_aripo_genome.csv",row.names = 1)
norm_counts_Quare_genome = read.csv("../outputs/norm_counts_quare_genome.csv",row.names = 1)

options(stringsAsFactors = FALSE)

####
### Load DE analysis results from GLMM (Fischer et al., 2021) ----
Aripo_DE_g_042920 = read.xlsx("../outputs/DE_GLMM_genome_042920.xlsx",sheet = 2)
Quare_DE_g_042920 = read.xlsx("../outputs/DE_GLMM_genome_042920.xlsx",sheet = 3)

print(head(Aripo_DE_g_042920))

# End of the code for loading results from Fischer et al. (2021)

# Step 1: Handling missing data and zero variance in probe profiles ----
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
# Load Id, pop, rear, interaction DE flags from DE analysis results file
# to extract summary tables (Table 3 in supporting information)
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

#Quare
Quare_filtered_summary=inner_join(Quare_DE_g_042920[,c(1,3,7,11)], 
                                  data.frame("Id"=names(KeepGenes_Quare), 
                                             "Kept"=KeepGenes_Quare,
                                             stringsAsFactors = FALSE), 
                                  by = "Id")
# population main effect
table("pop_DE"=Quare_filtered_summary$pop_DE_flag,
      "Kept"=Quare_filtered_summary$Kept)

# Obtain residuals from lmer using log-transformed normalized counts.
# (Do not run the next commented line unless you want to overwrite residual files)
# source("get_residuals_from_lmer.R":)
resid_Aripo = read.table(file= "../outputs/network/residual_aripo.csv",sep=",", row.names=1, header=T) %>% as.matrix()
resid_Quare = read.table(file= "../outputs/network/residual_quare.csv",sep=",", row.names=1, header=T) %>% as.matrix()

print(head(resid_Aripo))

####
#Step 2 : Covaraince test -----
# LC (2012), CLX (2014),
#Critical values from the 95th percentile of the simulated data.
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
na_ind = apply(resid_Aripo, 1, function(x) any(is.na(x)))#Only Aripo has the NA issue, removed, so using 13446 genes

### Aripo----
  res_LC = matrix(nrow=4,ncol=2, 
                  dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","LPNP_v_HPNP","HPP_v_HPNP"),
                                  c("Aripo","Quare")))
  test_stat_LC = matrix(nrow=4,ncol=2, 
                        dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","LPNP_v_HPNP","HPP_v_HPNP"),
                                        c("Aripo","Quare")))
  asym_pval_LC = matrix(nrow=4,ncol=2, 
                        dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","LPNP_v_HPNP","HPP_v_HPNP"),
                                        c("Aripo","Quare")))
  
  temp_res = equalCovs(t(resid_Aripo[!na_ind,LPP_Aripo]),t(resid_Aripo[!na_ind,LPNP_Aripo]),
                       size1 = length_LPP_Aripo, size2 = length_LPNP_Aripo)
  asym_pval_LC["LPP_v_LPNP","Aripo"] = temp_res[2]
  test_stat_LC["LPP_v_LPNP","Aripo"] = temp_res[1]
  
  temp_res = equalCovs(t(resid_Aripo[!na_ind,LPP_Aripo]),t(resid_Aripo[!na_ind,HPP_Aripo]),
                       size1 = length_LPP_Aripo, size2 = length_HPP_Aripo)
  asym_pval_LC["LPP_v_HPP","Aripo"] = temp_res[2]
  test_stat_LC["LPP_v_HPP","Aripo"] = temp_res[1]
  
  temp_res = equalCovs(t(resid_Aripo[!na_ind,LPNP_Aripo]),t(resid_Aripo[!na_ind,HPNP_Aripo]),
                       size1 = length_LPNP_Aripo, size2 = length_HPNP_Aripo)
  asym_pval_LC["LPNP_v_HPNP","Aripo"] = temp_res[2]
  test_stat_LC["LPNP_v_HPNP","Aripo"] = temp_res[1]
  
  temp_res = equalCovs(t(resid_Aripo[!na_ind,HPP_Aripo]),t(resid_Aripo[!na_ind,HPNP_Aripo]),
                       size1 = length_HPP_Aripo, size2 = length_HPNP_Aripo)
  asym_pval_LC["HPP_v_HPNP","Aripo"] = temp_res[2]
  test_stat_LC["HPP_v_HPNP","Aripo"] = temp_res[1]
  
  ### Quare----
  
  temp_res = equalCovs(t(resid_Quare[,LPP_Quare]),t(resid_Quare[,LPNP_Quare]),
                       size1 = length_LPP_Quare, size2 = length_LPNP_Quare)
  asym_pval_LC["LPP_v_LPNP","Quare"] = temp_res[2]
  test_stat_LC["LPP_v_LPNP","Quare"] = temp_res[1]
  
  temp_res = equalCovs(t(resid_Quare[,LPP_Quare]),t(resid_Quare[,HPP_Quare]),
                       size1 = length_LPP_Quare, size2 = length_HPP_Quare)
  asym_pval_LC["LPP_v_HPP","Quare"] = temp_res[2]
  test_stat_LC["LPP_v_HPP","Quare"] = temp_res[1]
  
  temp_res = equalCovs(t(resid_Quare[,LPNP_Quare]),t(resid_Quare[,HPNP_Quare]),
                       size1 = length_LPNP_Quare, size2 = length_HPNP_Quare)
  asym_pval_LC["LPNP_v_HPNP","Quare"] = temp_res[2]
  test_stat_LC["LPNP_v_HPNP","Quare"] = temp_res[1]
  
  temp_res = equalCovs(t(resid_Quare[,HPP_Quare]),t(resid_Quare[,HPNP_Quare]),
                       size1 = length_HPP_Quare, size2 = length_HPNP_Quare)
  asym_pval_LC["HPP_v_HPNP","Quare"] = temp_res[2]
  test_stat_LC["HPP_v_HPNP","Quare"] = temp_res[1]
  
  
  #CLX
  res_CLX = matrix(nrow=4,ncol=2, 
                   dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","LPNP_v_HPNP","HPP_v_HPNP"),
                                   c("Aripo","Quare")))
  test_stat_CLX = matrix(nrow=4,ncol=2, 
                         dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","LPNP_v_HPNP","HPP_v_HPNP"),
                                         c("Aripo","Quare")))
  asym_pval_CLX = matrix(nrow=4,ncol=2, 
                         dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","LPNP_v_HPNP","HPP_v_HPNP"),
                                         c("Aripo","Quare")))
  
  ## CLX (2014)----
  ### Aripo----
  temp_res = CLX_cpp(t(resid_Aripo[!na_ind,LPP_Aripo]),t(resid_Aripo[!na_ind,LPNP_Aripo]),b_size = 400)
  asym_pval_CLX["LPP_v_LPNP","Aripo"] = temp_res$pvalue
  test_stat_CLX["LPP_v_LPNP","Aripo"] = temp_res$TSvalue
  
  temp_res = CLX_cpp(t(resid_Aripo[!na_ind,LPP_Aripo]),t(resid_Aripo[!na_ind,HPP_Aripo]),b_size = 400)
  asym_pval_CLX["LPP_v_HPP","Aripo"] = temp_res$pvalue
  test_stat_CLX["LPP_v_HPP","Aripo"] = temp_res$TSvalue
  
  temp_res = CLX_cpp(t(resid_Aripo[!na_ind,LPNP_Aripo]),t(resid_Aripo[!na_ind,HPNP_Aripo]),b_size = 400)
  asym_pval_CLX["LPNP_v_HPNP","Aripo"] = temp_res$pvalue
  test_stat_CLX["LPNP_v_HPNP","Aripo"] = temp_res$TSvalue
  
  temp_res = CLX_cpp(t(resid_Aripo[!na_ind,HPP_Aripo]),t(resid_Aripo[!na_ind,HPNP_Aripo]),b_size = 400)
  asym_pval_CLX["HPP_v_HPNP","Aripo"] = temp_res$pvalue
  test_stat_CLX["HPP_v_HPNP","Aripo"] = temp_res$TSvalue
  ### Quare----
  temp_res = CLX_cpp(t(resid_Quare[!na_ind,LPP_Quare]),t(resid_Quare[,LPNP_Quare]),b_size = 400)
  asym_pval_CLX["LPP_v_LPNP","Quare"] = temp_res$pvalue
  test_stat_CLX["LPP_v_LPNP","Quare"] = temp_res$TSvalue
  
  temp_res = CLX_cpp(t(resid_Quare[!na_ind,LPP_Quare]),t(resid_Quare[,HPP_Quare]),b_size = 400)
  asym_pval_CLX["LPP_v_HPP","Quare"] = temp_res$pvalue
  test_stat_CLX["LPP_v_HPP","Quare"] = temp_res$TSvalue
  
  temp_res = CLX_cpp(t(resid_Quare[!na_ind,LPNP_Quare]),t(resid_Quare[,HPNP_Quare]),b_size = 400)
  asym_pval_CLX["LPNP_v_HPNP","Quare"] = temp_res$pvalue
  test_stat_CLX["LPNP_v_HPNP","Quare"] = temp_res$TSvalue
  
  temp_res = CLX_cpp(t(resid_Quare[!na_ind,HPP_Quare]),t(resid_Quare[,HPNP_Quare]),b_size = 400)
  asym_pval_CLX["HPP_v_HPNP","Quare"] = temp_res$pvalue
  test_stat_CLX["HPP_v_HPNP","Quare"] = temp_res$TSvalue
  
  
  #Test statistics
  print("Test statistics:")
  print("LC")
  print(test_stat_LC)
  print("CLX")
  print(test_stat_CLX)
  
  print("p-values from asymptotic distribution:")
  print("LC")
  print(asym_pval_LC)
  print("CLX")
  print(asym_pval_CLX)
  
  
####
# Step 3: run the random projection tests ## ----
  
resid_Aripo_NA_cleaned = resid_Aripo[!na_ind,]

p_Aripo = dim(resid_Aripo_NA_cleaned)[1]
p_Quare = dim(resid_Quare)[1]

set.seed(2021)

temp_stat_rand_proj =matrix(NA, nrow = 4, ncol = 2, 
                            dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","HPP_v_HPNP","LPNP_v_HPNP"), c("Aripo","Quare")))

begin.time = Sys.time()
## use all genes----
####
### LPP vs LPNP----
#### Aripo----
X_1 = t(resid_Aripo_NA_cleaned[, LPP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[, LPNP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
print(temp_result)
temp_stat_rand_proj["LPP_v_LPNP","Aripo"] = temp_result$test.stat

#### Quare----
Y_1 = t(resid_Quare[, LPP_Quare])
Y_2 = t(resid_Quare[, LPNP_Quare])

temp_result = rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_LPNP","Quare"] = temp_result$test.stat

####
### LPP vs HPP----
#### Aripo----
X_1 = t(resid_Aripo_NA_cleaned[, LPP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[, HPP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_HPP","Aripo"] = temp_result$test.stat

#### Quare----
Y_1 = t(resid_Quare[, LPP_Quare])
Y_2 = t(resid_Quare[, HPP_Quare])

###LC
temp_result = rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_HPP","Quare"] = temp_result$test.stat

####
### LPNP vs HPNP----
#### Aripo----
X_1 = t(resid_Aripo_NA_cleaned[, LPNP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[, HPNP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPNP_v_HPNP","Aripo"] = temp_result$test.stat

#### Quare----
Y_1 = t(resid_Quare[, LPNP_Quare])
Y_2 = t(resid_Quare[, HPNP_Quare])

temp_result = rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPNP_v_HPNP","Quare"] = temp_result$test.stat

####
### HPP vs HPNP----
#### Aripo----
X_1 = t(resid_Aripo_NA_cleaned[, HPP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[, HPNP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["HPP_v_HPNP","Aripo"] = temp_result$test.stat

####Quare----
Y_1 = t(resid_Quare[, HPP_Quare])
Y_2 = t(resid_Quare[, HPNP_Quare])

###LC
temp_result =  rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["HPP_v_HPNP","Quare"] = temp_result$test.stat
# }
end.time = Sys.time()
print(paste("Projection test takes", round(difftime(end.time, begin.time, units='mins'),4), "minutes"))

res_rand_proj = temp_stat_rand_proj

print(paste("Test type:", rand_project_test_stat , ", K:", K, ", q: ",q))
print("Test stats:")
print(round(temp_stat_rand_proj,4))

## DE Genes only----
temp_stat_rand_proj <- matrix(NA, nrow = 4, ncol = 2, 
                              dimnames = list(c("LPP_v_LPNP","LPP_v_HPP","HPP_v_HPNP","LPNP_v_HPNP"), c("Aripo","Quare")))

#DE gene names
Aripo_DE_gene_name <- intersect(rownames(resid_Aripo_NA_cleaned), subset(Aripo_filtered_summary, select = Id, subset = (pop_DE_flag == TRUE))%>%sapply(as.character))
Quare_DE_gene_name <- intersect(rownames(resid_Quare), subset(Quare_filtered_summary, select = Id, subset = (pop_DE_flag == TRUE))%>%sapply(as.character))

###LPP vs LPNP----
####Aripo----
X_1 <- t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, LPP_Aripo])
X_2 <- t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, LPNP_Aripo])

temp_result <- rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_LPNP","Aripo"] <- temp_result$test.stat

####Quare----
Y_1 = t(resid_Quare[Quare_DE_gene_name, LPP_Quare])
Y_2 = t(resid_Quare[Quare_DE_gene_name, LPNP_Quare])

temp_result = rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_LPNP","Quare"] = temp_result$test.stat

####
###LPP vs HPP----
####Aripo----
X_1 = t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, LPP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, HPP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_HPP","Aripo"] = temp_result$test.stat

####Quare----
Y_1 = t(resid_Quare[Quare_DE_gene_name, LPP_Quare])
Y_2 = t(resid_Quare[Quare_DE_gene_name, HPP_Quare])

###LC
temp_result = rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPP_v_HPP","Quare"] = temp_result$test.stat

####
###LPNP vs HPNP----
####Aripo----
X_1 = t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, LPNP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, HPNP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPNP_v_HPNP","Aripo"] = temp_result$test.stat

####Quare----
Y_1 = t(resid_Quare[Quare_DE_gene_name, LPNP_Quare])
Y_2 = t(resid_Quare[Quare_DE_gene_name, HPNP_Quare])

temp_result = rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["LPNP_v_HPNP","Quare"] = temp_result$test.stat

####
###HPP vs HPNP----
####Aripo----
X_1 = t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, HPP_Aripo])
X_2 = t(resid_Aripo_NA_cleaned[Aripo_DE_gene_name, HPNP_Aripo])

temp_result = rand_proj_cov(X_1,X_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["HPP_v_HPNP","Aripo"] = temp_result$test.stat

####Quare----
Y_1 = t(resid_Quare[Quare_DE_gene_name, HPP_Quare])
Y_2 = t(resid_Quare[Quare_DE_gene_name, HPNP_Quare])

###LC
temp_result =  rand_proj_cov(Y_1,Y_2, proj.dim = q, n_proj = K, method = rand_project_test_stat)
temp_stat_rand_proj["HPP_v_HPNP","Quare"] = temp_result$test.stat

print(paste("Test type:", rand_project_test_stat , ", K:", K, ", q: ",q))
print("Test stats:")
print(round(temp_stat_rand_proj,4))


# Step 4: network comparison based on Correlation "networks"  ------
#load extra libraries
  
set.seed(2021)
#Population DE lists
Aripo_popkept_DE_id = setdiff(unlist(subset(Aripo_filtered_summary, select = "Id", subset = (pop_DE_flag == 1)&(Kept == 1))), names(which(na_ind)))
Quare_popkept_DE_id = unlist(subset(Quare_filtered_summary, select = "Id", subset = (pop_DE_flag == 1)&(Kept == 1)))

violinplot_save = FALSE
# cor_fdr_levels = c(0.01,0.05,0.1) #To draw violin summary figures 
cor_fdr_levels = c(0.05) #To have test results.

violin_summary_stats <- data.frame()
cor_net_summary_stats <- data.frame()

for(fdr_level in cor_fdr_levels){
  ##Aripo HPP----
  #correlation "networks" from 
  message("Aripo HPP")
  cor_nets <- cor_CL(x=t(resid_Aripo[!na_ind,HPP_Aripo]), FDR_level = fdr_level)$decision
  #To compare DE vs NDE networks extract index
  Aripo_popkept_DE_ind <- match(Aripo_popkept_DE_id, rownames(cor_nets))
  nettest_res <- list()
  size_subgraph <- 200 #Used for violin summary plot (Maugis et al., 2017)
  
  ###Triangle----
  nettest_res$triangle$Aripo$HPP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                             cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]
  )
  nettest_res$triangle$Aripo$HPP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                       cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Aripo$HPP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                  cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                  alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Aripo$HPP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                           cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                           motif = "vshape"
  )
  nettest_res$vshape$Aripo$HPP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                  cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                  motif = "vshape",
                                                  alternative = "less"
  )
  nettest_res$vshape$Aripo$HPP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     motif = "vshape",
                                                     alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Aripo$HPP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                              cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                              motif = "threestar"
  )
  nettest_res$threestar$Aripo$HPP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                motif = "threestar",
                                                alternative = "less"
  )
  nettest_res$threestar$Aripo$HPP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "greater"
  )
  ### violin summary plots----
  #### DE----
  temp_violin_summary_stats <- violin_netsummary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                 subsample_sizes = size_subgraph, R=2000, max_cycle_order=7, 
                                                 y.max = 0.3, 
                                                 save.plot=violinplot_save,
                                                 filename = paste0("../outputs/network/netsummary/Aripo_HPP_DE_", fdr_level*100,".pdf"))
  
  violin_summary_stats <- cbind(dataset = "Aripo", group = "HPP" , DE_NDE = "DE", temp_violin_summary_stats)
  print(head(temp_violin_summary_stats))
  #### NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.3, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Aripo_HPP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "HPP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  print(head(temp_violin_summary_stats))
  print(head(violin_summary_stats))
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Aripo_HPP_DE = network_summary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind]),
                                 Aripo_HPP_NDE = network_summary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]))
  print(cor_net_summary_stats)
  ##Aripo HPNP----
  message("Aripo HPNP")
  cor_nets <- cor_CL(x=t(resid_Aripo[!na_ind,HPNP_Aripo]), 
                     FDR_level = fdr_level)$decision
  
  ###Triangle----
  nettest_res$triangle$Aripo$HPNP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                       cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]
  )
  nettest_res$triangle$Aripo$HPNP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                  cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Aripo$HPNP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Aripo$HPNP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Aripo$HPNP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Aripo$HPNP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Aripo$HPNP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                        cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Aripo$HPNP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Aripo$HPNP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                      cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  ### violin summary----
  #### DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.3, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Aripo_HPNP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "HPNP" , DE_NDE = "DE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  #### NDE----
  temp_violin_summary_stats<- violin_netsummary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                                y.max = 0.3, 
                                                save.plot=violinplot_save,
                                                filename = paste0("../outputs/network/netsummary/Aripo_HPNP_NDE_",fdr_level*100,".pdf"))
  
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "HPNP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Aripo_HPNP_DE = network_summary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind]),
                                 Aripo_HPNP_NDE = network_summary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]))
  
  
  ##Aripo LPP----
  message("Aripo LPP")
  #two genes has 0 variance. excluding the as well
  zero_var_ind = (apply(resid_Aripo[,LPP_Aripo],1,var)==0)
  cor_nets = cor_CL(x=t(resid_Aripo[!(na_ind|zero_var_ind),LPP_Aripo]), FDR_level = fdr_level)$decision
  
  ###Triangle----
  nettest_res$triangle$Aripo$LPP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                       cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]
  )
  nettest_res$triangle$Aripo$LPP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                  cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Aripo$LPP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Aripo$LPP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Aripo$LPP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Aripo$LPP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Aripo$LPP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                        cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Aripo$LPP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Aripo$LPP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                      cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  ###violin summary----
  #### DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.3, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Aripo_LPP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "LPP" , DE_NDE = "DE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Aripo_HPP_DE = network_summary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind]),
                                 Aripo_HPP_NDE = network_summary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]))
  #### NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.3, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Aripo_LPP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "LPP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Aripo_HPP_DE = network_summary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind]),
                                 Aripo_HPP_NDE = network_summary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]))
  
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Aripo_LPP_DE = network_summary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind]),
                                 Aripo_LPP_NDE = network_summary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]))
  
  
  ##Aripo LPNP----
  message("Aripo LPNP")
  cor_nets = cor_CL(x=t(resid_Aripo[!na_ind,LPNP_Aripo]), FDR_level = fdr_level)$decision
  
  ###Triangle----
  nettest_res$triangle$Aripo$LPNP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                       cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]
  )
  nettest_res$triangle$Aripo$LPNP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                  cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Aripo$LPNP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Aripo$LPNP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                     cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Aripo$LPNP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Aripo$LPNP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Aripo$LPNP$two.sided <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                        cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Aripo$LPNP$less <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                   cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Aripo$LPNP$greater <- net_test(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                                      cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  ###violin summary----
  #### DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.3, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Aripo_LPNP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "LPNP" , DE_NDE = "DE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  ####NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.3, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Aripo_LPNP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Aripo", group = "LPNP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  
  
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Aripo_LPNP_DE = network_summary(cor_nets[Aripo_popkept_DE_ind,Aripo_popkept_DE_ind]),
                                 Aripo_LPNP_NDE = network_summary(cor_nets[-Aripo_popkept_DE_ind,-Aripo_popkept_DE_ind]))
  
  
  
  ##Quare HPP----
  #correlation "networks" from 
  message("Quare HPP")
  cor_nets <- cor_CL(x=t(resid_Quare[,HPP_Quare]), FDR_level = fdr_level)$decision
  
  #Compare DE vs NDE networks
  Quare_popkept_DE_ind <- match(Quare_popkept_DE_id, rownames(cor_nets))
  size_subgraph <- 200
  ###Triangle----
  nettest_res$triangle$Quare$HPP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                       cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
  )
  nettest_res$triangle$Quare$HPP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                  cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Quare$HPP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Quare$HPP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Quare$HPP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Quare$HPP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Quare$HPP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                        cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Quare$HPP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Quare$HPP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                      cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  ###violin summary----
  ####DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_HPP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "HPP" , DE_NDE = "DE", temp_violin_summary_stats)
  print(head(temp_violin_summary_stats))
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  ####NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_HPP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "HPP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  
  
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Quare_HPP_DE = network_summary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind]),
                                 Quare_HPP_NDE = network_summary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind]))
  
  
  ##Quare HPNP----
  message("Quare HPNP")
  cor_nets <- cor_CL(x=t(resid_Quare[,HPNP_Quare]), FDR_level = fdr_level)$decision
  ###Triangle----
  nettest_res$triangle$Quare$HPNP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                       cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
  )
  nettest_res$triangle$Quare$HPNP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                  cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Quare$HPNP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Quare$HPNP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Quare$HPNP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Quare$HPNP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Quare$HPNP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                        cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Quare$HPNP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Quare$HPNP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                      cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  ### violin summary
  #### DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_HPNP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "HPNP" , DE_NDE = "DE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  ####NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_HPNP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "HPNP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Quare_HPNP_DE = network_summary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind]),
                                 Quare_HPNP_NDE = network_summary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind]))
  
  
  ##Quare LPP----
  message("Quare LPP")
  cor_nets <- cor_CL(x=t(resid_Quare[,LPP_Quare]), FDR_level = fdr_level)$decision
  ###Triangle----
  nettest_res$triangle$Quare$LPP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                       cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
  )
  nettest_res$triangle$Quare$LPP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                  cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Quare$LPP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Quare$LPP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Quare$LPP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Quare$LPP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Quare$LPP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                        cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Quare$LPP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Quare$LPP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                      cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  
  ### violin summary----
  #### DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_LPP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "LPP" , DE_NDE = "DE", temp_violin_summary_stats)
  
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  ####NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_LPP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "LPP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Quare_LPP_DE = network_summary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind]),
                                 Quare_LPP_NDE = network_summary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind]))
  
  ##Quare LPNP----
  message("Quare LPNP")
  cor_nets = cor_CL(x=t(resid_Quare[,LPNP_Quare]), FDR_level = fdr_level)$decision
  ###Triangle----
  nettest_res$triangle$Quare$LPNP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                       cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
  )
  nettest_res$triangle$Quare$LPNP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                  cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                  alternative = "less"
  )
  nettest_res$triangle$Quare$LPNP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     alternative = "greater"
  )
  ###v-shape----
  nettest_res$vshape$Quare$LPNP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                     cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                     motif = "vshape"
  )
  nettest_res$vshape$Quare$LPNP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                motif = "vshape",
                                                alternative = "less"
  )
  nettest_res$vshape$Quare$LPNP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "vshape",
                                                   alternative = "greater"
  )
  ###threestar----
  nettest_res$threestar$Quare$LPNP$two.sided <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                        cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                        motif = "threestar"
  )
  nettest_res$threestar$Quare$LPNP$less <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                   cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                   motif = "threestar",
                                                   alternative = "less"
  )
  nettest_res$threestar$Quare$LPNP$greater <- net_test(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                                      cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                                      motif = "threestar",
                                                      alternative = "greater"
  )
  ### vioin summary ----
  #### DE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_LPNP_DE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "LPNP" , DE_NDE = "DE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  ####NDE----
  temp_violin_summary_stats<-violin_netsummary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind],
                                               subsample_sizes = size_subgraph, R=2000, max_cycle_order=7,
                                               y.max = 0.6, 
                                               save.plot=violinplot_save,
                                               filename = paste0("../outputs/network/netsummary/Quare_LPNP_NDE_",fdr_level*100,".pdf"))
  temp_violin_summary_stats <- cbind(dataset = "Quare", group = "LPNP" , DE_NDE = "NDE", temp_violin_summary_stats)
  violin_summary_stats <- rbind(violin_summary_stats, temp_violin_summary_stats)
  cor_net_summary_stats <- rbind(cor_net_summary_stats,
                                 Quare_LPNP_DE = network_summary(cor_nets[Quare_popkept_DE_ind,Quare_popkept_DE_ind]),
                                 Quare_LPNP_NDE = network_summary(cor_nets[-Quare_popkept_DE_ind,-Quare_popkept_DE_ind]))
}
rm(cor_nets)

#### Print test results ----
message("print test (triangle)")
print(nettest_res$triangle$Aripo)
print(nettest_res$triangle$Quare)
lapply(nettest_res$triangle$Aripo, print)
lapply(nettest_res$triangle$Quare, print)
message("Aripo net comparison test p-values (triangle)")
lapply(nettest_res$triangle$Aripo, function(x) print(x$p_value))
message("Quare net comparison test p-values (triangle)")
lapply(nettest_res$triangle$Quare, function(x) print(x$p_value))

message("print test (v-shape)")
print(nettest_res$vshape$Aripo)
print(nettest_res$vshape$Quare)
lapply(nettest_res$vshape$Aripo, print)
lapply(nettest_res$vshape$Quare, print)
message("Aripo net comparison test p-values (v-shape)")
lapply(nettest_res$vshape$Aripo, function(x) print(x$p_value))
message("Quare net comparison test p-values (v-shape)")
lapply(nettest_res$vshape$Quare, function(x) print(x$p_value))

message("print test (threestar)")
print(nettest_res$threestar$Aripo)
print(nettest_res$threestar$Quare)
lapply(nettest_res$threestar$Aripo, print)
lapply(nettest_res$threestar$Quare, print)
message("Aripo net comparison test p-values (threestar)")
lapply(nettest_res$threestar$Aripo, function(x) print(x$p_value))
message("Quare net comparison test p-values (threestar)")
lapply(nettest_res$threestar$Quare, function(x) print(x$p_value))

if(file.save){
  write.csv(violin_summary_stats, file = paste0("../outputs/network/netsummary/violin_summary_",Sys.time(),".csv"))
  save.image(file = paste0("coexpression_",Sys.time(),".RData"))
}
message("Finished!")
