# Step 1: Use residual to test
library(lme4)
#Aripo
p = dim(norm_counts_Aripo_Keep)[1]
n = dim(norm_counts_Aripo_Keep)[2]
# behave_Aripo$pop =factor(behave_Aripo$pop,levels=c("LP","HP"))
# behave_Aripo$rear=factor(behave_Aripo$rear,levels=c("P","NP"))
# behave_Aripo$week=factor(behave_Aripo$week)
# behave_Aripo$family=factor(behave_Aripo$family)

singular_Aripo= rep(NA,p)
resid_Aripo = matrix(nrow =p, ncol=n)
rownames(resid_Aripo) = rownames(norm_counts_Aripo_Keep)
colnames(resid_Aripo) = behave_Aripo$fish
#Construct residual matrix 
for(i in 1:p){
  y = unlist(log2(1+ norm_counts_Aripo_Keep[i,]))
  fit_lmer <- tryCatch(lmer(y ~ pop + rear 
                            + (1|week) + (1|family), data = behave_Aripo),
                       error = function(e){NA})
  if(is.na(fit_lmer)){next}
  singular_Aripo[i] = isSingular(fit_lmer)
  resid_Aripo[i,] = resid(fit_lmer)
}
#Quare
p = dim(norm_counts_Quare_Keep)[1]
n = dim(norm_counts_Quare_Keep)[2]

singular_Quare= rep(NA,p)
resid_Quare = matrix(nrow =p, ncol=n)
rownames(resid_Quare) = rownames(norm_counts_Quare_Keep)
colnames(resid_Quare) = behave_Quare$fish
#Construct residual matrix 
for(i in 1:p){
  if(i %% 100== 0) print(paste(i, "th iteration"))
  y = unlist(log2(1+ norm_counts_Quare_Keep[i,]))
  fit_lmer <- lmer(y ~ pop + rear 
                   + (1|week) + (1|family_new), data = behave_Quare)
  singular_Quare[i] = isSingular(fit_lmer)
  resid_Quare[i,] = resid(fit_lmer)
}

#Save residual from lmer
# write.csv(resid_Aripo, file= "../outputs/network/residual_aripo.csv", row.names = TRUE)
# write.csv(resid_Quare, file= "../outputs/network/residual_quare.csv", row.names = TRUE)