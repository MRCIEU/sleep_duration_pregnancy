# Clear the work environment
rm(list = ls())

#package
#install.packages("devtools")
#library(devtools)
#Sys.setenv(GITHUB_PAT="/a/path/from/github")
#install_github("MRCIEU/TwoSampleMR")
require(TwoSampleMR)

#input data
gxA<-read.csv("/path/to/the/data")
gxB<-read.csv("/path/to/the/data")
#gyA<-read.csv("/path/to/the/data")
#gyB<-read.csv("/path/to/the/data")
eanea<-read.csv("/path/to/the/data")
snpgroup<-read.csv("/path/to/the/data")

colnames(eanea)[1]<-"X"
colnames(snpgroup)[1]<-"X"
eanea$ea<-substr(eanea$eanea,1,1)
eanea$nea<-substr(eanea$eanea,3,3)


#prepare G-X
prepareGX<-function(gx,group){
  gx<-merge(eanea,gxA,by="X")
  gx<-merge(gx,snpgroup,by="X")
  association<-gx[which(gx[group]==1),]
  association<-association[,c("X","ea","nea","EAF","duration_beta","duration_SE","duration_Pval")]
  colnames(association)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval")
  gxforMR<-format_data(association,type="exposure")
  return(gxforMR)
}

gxA55<-prepareGX(gxA,"grp6")
gxB55<-prepareGX(gxB,"grp6")
gxA43<-prepareGX(gxA,"grp1")
gxB43<-prepareGX(gxB,"grp1")

#prepare G-Y
#gyA<-merge(eanea,gyA,by="X")
#gyB<-merge(eanea,gyB,by="X")
#write.csv(gyA,"/path/to/the/data")
#write.csv(gyB,"/path/to/the/data")

prepareGY<-function(gxsnp,gyfile,outcome,outcome_SE){
  gy<-read_outcome_data(snps=gxsnp,filename=paste("/path/to/the/data",sep=""),
                        sep=",",snp_col="X",beta_col=outcome,se_col=outcome_SE,effect_allele_col="ea",other_allele_col="nea")
  return(gy)
}

#55 SNPs, 43 SNPs
outcome<-c("loss","stillbirthM","miscarriageM","gdmM","hdp","pdM","pt","lbw","hbw","bw")
outcomeSE<-c("lossSE","stillbirthMSE","miscarriageMSE","gdmMSE","hdpSE","pdMSE","ptSE","lbwSE","hbwSE","bwSE")
i<-1
for (i in 1:length(outcome)){
  assign(paste("gyA55_",outcome[i],sep=""),prepareGY(gxA55$SNP,"UKB_gyA2",outcome[i],outcomeSE[i]))
  assign(paste("gyB55_",outcome[i],sep=""),prepareGY(gxA55$SNP,"UKB_gyB2",outcome[i],outcomeSE[i]))
  assign(paste("gyA43_",outcome[i],sep=""),prepareGY(gxA43$SNP,"UKB_gyA2",outcome[i],outcomeSE[i]))
  assign(paste("gyB43_",outcome[i],sep=""),prepareGY(gxA43$SNP,"UKB_gyB2",outcome[i],outcomeSE[i]))
}

#MR estimates
MRresults<-function(gx,gy){
  gxgyforMR<-harmonise_data(exposure_dat=gx,outcome_dat=gy,action = 1)
  results<-mr(gxgyforMR,method_list=c("mr_ivw"))
  results$exposure<-deparse(substitute(gx))
  results$outcome<-deparse(substitute(gy))
  return(results)
}

#55 SNPs
#G-Y in B / G-X in A
a1<-MRresults(gxA55,gyB55_loss)
a2<-MRresults(gxA55,gyB55_stillbirthM)
a3<-MRresults(gxA55,gyB55_miscarriageM)
a4<-MRresults(gxA55,gyB55_gdmM)
a5<-MRresults(gxA55,gyB55_hdp)
a6<-MRresults(gxA55,gyB55_pdM)
a7<-MRresults(gxA55,gyB55_pt)
a8<-MRresults(gxA55,gyB55_lbw)
a9<-MRresults(gxA55,gyB55_hbw)
a10<-MRresults(gxA55,gyB55_bw)
#G-Y in A / G-X in B
b1<-MRresults(gxB55,gyA55_loss)
b2<-MRresults(gxB55,gyA55_stillbirthM)
b3<-MRresults(gxB55,gyA55_miscarriageM)
b4<-MRresults(gxB55,gyA55_gdmM)
b5<-MRresults(gxB55,gyA55_hdp)
b6<-MRresults(gxB55,gyA55_pdM)
b7<-MRresults(gxB55,gyA55_pt)
b8<-MRresults(gxB55,gyA55_lbw)
b9<-MRresults(gxB55,gyA55_hbw)
b10<-MRresults(gxB55,gyA55_bw)

#43 SNPs
#G-Y in B / G-X in A
c1<-MRresults(gxA43,gyB43_loss)
c2<-MRresults(gxA43,gyB43_stillbirthM)
c3<-MRresults(gxA43,gyB43_miscarriageM)
c4<-MRresults(gxA43,gyB43_gdmM)
c5<-MRresults(gxA43,gyB43_hdp)
c6<-MRresults(gxA43,gyB43_pdM)
c7<-MRresults(gxA43,gyB43_pt)
c8<-MRresults(gxA43,gyB43_lbw)
c9<-MRresults(gxA43,gyB43_hbw)
c10<-MRresults(gxA43,gyB43_bw)
#G-Y in A / G-X in B
d1<-MRresults(gxB43,gyA43_loss)
d2<-MRresults(gxB43,gyA43_stillbirthM)
d3<-MRresults(gxB43,gyA43_miscarriageM)
d4<-MRresults(gxB43,gyA43_gdmM)
d5<-MRresults(gxB43,gyA43_hdp)
d6<-MRresults(gxB43,gyA43_pdM)
d7<-MRresults(gxB43,gyA43_pt)
d8<-MRresults(gxB43,gyA43_lbw)
d9<-MRresults(gxB43,gyA43_hbw)
d10<-MRresults(gxB43,gyA43_bw)

#output results
AonB55<-rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
BonA55<-rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
AonB43<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
BonA43<-rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10)

output<-function(allresults,resultsname){
  allresults["lci"]<-allresults["b"]-1.96*allresults["se"]
  allresults["uci"]<-allresults["b"]+1.96*allresults["se"]
  allresults["OR"]<-round(exp(allresults["b"]),2)
  allresults["ORlci"]<-round(exp(allresults["lci"]),2)
  allresults["ORuci"]<-round(exp(allresults["uci"]),2)
  write.csv(allresults,paste("/path/to/the/data",resultsname,"_v2.csv",sep=""),row.names = FALSE)
}

output(AonB55,"AonB55")
output(BonA55,"BonA55")
output(AonB43,"AonB43")
output(BonA43,"BonA43")

###############################################################################################################################################
#meta-analyse MR estimates from AonB and BonA
# Clear the work environment
rm(list = ls())

#package
require(metafor)

#input data
AonB55<-read.csv("/path/to/the/data")
BonA55<-read.csv("/path/to/the/data")
AonB43<-read.csv("/path/to/the/data")
BonA43<-read.csv("/path/to/the/data")

#create empty dataframes
pooldata55<-data.frame("study"=c("AonB","BonA"),"beta"=c(NA,NA),"se"=c(NA,NA))
pooldata43<-data.frame("study"=c("AonB","BonA"),"beta"=c(NA,NA),"se"=c(NA,NA))
vars55 <- AonB55[,c("outcome","method")]
vars43 <- AonB43[,c("outcome","method")]
poolresults55 <- matrix(, ncol=3, nrow=10)
poolresults43 <- matrix(, ncol=3, nrow=10)
colnames(poolresults55)<-c("beta","se","pFORq")
colnames(poolresults43)<-c("beta","se","pFORq")
poolresults55 <- cbind(vars55,poolresults55)
poolresults43 <- cbind(vars43,poolresults43)

#meta-analysis
for (i in 1:10){
  #extract estimates
  pooldata55[1,2]<-AonB55[i,"b"]
  pooldata55[1,3]<-AonB55[i,"se"]
  pooldata55[2,2]<-BonA55[i,"b"]
  pooldata55[2,3]<-BonA55[i,"se"]
  pooldata43[1,2]<-AonB43[i,"b"]
  pooldata43[1,3]<-AonB43[i,"se"]
  pooldata43[2,2]<-BonA43[i,"b"]
  pooldata43[2,3]<-BonA43[i,"se"]
  #meta
  a<-rma(yi=beta, sei=se, slab=study, data=pooldata55,method="FE")
  b<-rma(yi=beta, sei=se, slab=study, data=pooldata43,method="FE")
  poolresults55[i,"beta"]<-a$beta
  poolresults55[i,"se"]<-a$se
  poolresults55[i,"pFORq"]<-a$QEp
  poolresults43[i,"beta"]<-b$beta
  poolresults43[i,"se"]<-b$se
  poolresults43[i,"pFORq"]<-b$QEp
}

#output results
output<-function(allresults,resultsname){
  allresults["lci"]<-allresults["beta"]-1.96*allresults["se"]
  allresults["uci"]<-allresults["beta"]+1.96*allresults["se"]
  allresults["OR"]<-round(exp(allresults["beta"]),2)
  allresults["ORlci"]<-round(exp(allresults["lci"]),2)
  allresults["ORuci"]<-round(exp(allresults["uci"]),2)
  write.csv(allresults,paste("/path/to/the/data",resultsname,"_v2.csv",sep=""),row.names = FALSE)
}

output(poolresults55,"pool55")
output(poolresults43,"pool43")
