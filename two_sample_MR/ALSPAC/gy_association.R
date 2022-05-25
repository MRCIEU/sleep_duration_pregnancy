# Clear the work environment
rm(list = ls())

#input data
alspacgenotype<-read.csv("/path/to/the/genotype/data")
alspacphenotype<-read.csv("/path/to/the/phenotype/data")
pc<-read.table("/path/to/the/principal/components")
pc<-pc[,-1]
colnames(pc)<-c("aln","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
pc$aln<-substr(pc$aln,1,5)

#combine datasets
alspac<-merge(alspacgenotype,alspacphenotype,by="aln")
alspac<-merge(alspac,pc,by="aln")

#depression
alspac$depressAP_Diag2<-alspac$depressAP_Diag
alspac$depressAP_Diag2[alspac$d171a==1]<-NA
table(alspac$depressAP_Diag2,exclude=NULL)

#check n
table(alspac$miscarriage,exclude=NULL)
table(alspac$loss,exclude=NULL)
table(alspac$stillbirth,exclude=NULL)
table(alspac$gdm, exclude=NULL)
table(alspac$HDP,exclude=NULL)
table(alspac$depressAP_Diag,exclude=NULL)
table(alspac$preterm,exclude=NULL)
table(alspac$lbw,exclude=NULL)
table(alspac$hbw,exclude=NULL)
summary(alspac$kz030)

alspac$b032.x[alspac$b032.x<0]<-NA
summary(alspac$b032.x)
sd(alspac$b032.x,na.rm=TRUE)

#G-Y in ALSPAC
vars <- names(alspac[c(2:79)])
col.names <- c("miscarriage","miscarriageSE","stillbirth","stillbirthSE","loss","lossSE","gdm","gdmSE","hdp","hdpSE",
               "pdM","pdMSE","pt","ptSE","lbw","lbwSE","hbw","hbwSE","bw","bwSE","eaf")
covar <- matrix(, ncol=21, nrow=length(vars))
dimnames(covar) <- list(vars, col.names)

#binary outcomes
outcome<-c("miscarriage","stillbirth","loss","gdm","HDP","depressAP_Diag2","preterm","lbw","hbw")
i<-1
for (i in 1:9){
  p<-2*i-1
  q<-2*i
  j<-1
  for (j in 1:length(vars)){
    m<-glm(alspac[,outcome[i]] ~ ., data=alspac[,c(vars[j],"mz028b","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")],family = "binomial")
    covar[j,p] <- summary(m)$coefficients[2,1]
    covar[j,q] <- summary(m)$coefficients[2,2]
  }
}
#continuous coutcome
j<-1
for(j in 1:length(vars)) {
  n<-lm(alspac$kz030 ~ ., data=alspac[,c(vars[j],"mz028b","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")])
  covar[j,19] <- summary(n)$coefficients[2,1]
  covar[j,20] <- summary(n)$coefficients[2,2]
  covar[j,21] <- mean(alspac[[vars[j]]])/2
}  


write.csv(covar,file="/path/to/save/the/gy/associations")
#################################################################################################
#to calculate maternal gy associations adjusting for fetal genotype
# Clear the work environment
rm(list = ls())

#input data
alspacgenotype<-read.csv("/path/to/maternal/genotype")
alspacgenotypeF<-read.csv("/path/to/fetal/genotype")
alspacphenotype<-read.csv("/path/to/phenotype")
pc<-read.table("/path/to/principal/components")
pc<-pc[,-1]
colnames(pc)<-c("aln","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
pc$aln<-substr(pc$aln,1,5)


#combine datasets
alspac<-merge(alspacgenotype,alspacphenotype,by="aln")
alspac<-merge(alspac,pc,by="aln")


#add  fetal genotype
alspac$aln_use<-paste(alspac$aln,alspac$qlet,sep="")
colnames(alspacgenotypeF)[1]<-"aln_use"
alspac2<-merge(alspac,alspacgenotypeF,by="aln_use")

#depression
alspac2$depressAP_Diag2<-alspac2$depressAP_Diag
alspac2$depressAP_Diag2[alspac2$d171a==1]<-NA
table(alspac2$depressAP_Diag2,exclude=NULL)

#check n
table(alspac2$gdm, exclude=NULL)
table(alspac2$HDP,exclude=NULL)
table(alspac2$depressAP_Diag,exclude=NULL)
table(alspac2$preterm,exclude=NULL)
table(alspac2$lbw,exclude=NULL)
table(alspac2$hbw,exclude=NULL)


#G-Y in ALSPAC
vars1 <- names(alspac2[c(3:80)])
vars2 <- names(alspac2[c(183:260)])
col.names <- c("miscarriage","miscarriageSE","stillbirth","stillbirthSE","gdm","gdmSE","hdp","hdpSE","pdM","pdMSE","pt","ptSE","lbw","lbwSE","hbw","hbwSE","bw","bwSE","eaf_check")
covar <- matrix(, ncol=19, nrow=length(vars1))
dimnames(covar) <- list(vars1, col.names)
#binary outcomes
outcome<-c("miscarriage","stillbirth","gdm","HDP","depressAP_Diag2","preterm","lbw","hbw")
i<-1
for (i in 1:8){
  p<-2*i-1
  q<-2*i
  j<-1
  for (j in 1:length(vars1)){
    m<-glm(alspac2[,outcome[i]] ~ ., data=alspac2[,c(vars1[j],vars2[j],"mz028b","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")],family = "binomial")
    covar[j,p] <- summary(m)$coefficients[2,1]
    covar[j,q] <- summary(m)$coefficients[2,2]
  }
}

#bw, eaf
m<-1
for (m in 1:length(vars1)){
  n<-lm(alspac2$kz030 ~ ., data=alspac2[,c(vars1[m],vars2[m],"mz028b","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")])
  covar[m,17] <- summary(n)$coefficients[2,1]
  covar[m,18] <- summary(n)$coefficients[2,2]
  covar[m,19] <- mean(alspac2[[vars1[m]]])/2
}



write.csv(covar,file="/path/to/save/gy/associations/adjusting/for/fetal/genotype")
