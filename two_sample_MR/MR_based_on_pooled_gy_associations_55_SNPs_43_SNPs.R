#55 SNPs as IVs
# Clear the work environment
rm(list = ls())

#package
#install.packages("devtools")
#library(devtools)
#Sys.setenv(GITHUB_PAT="/a/path/from/github")
#install_github("MRCIEU/TwoSampleMR")
require(TwoSampleMR)

#input data
gx<-read.csv("/path/to/the/data")
gy<-read.csv("/path/to/the/data")
eanea<-read.csv("/path/to/the/data")
snpgroup<-read.csv("/path/to/the/data")

colnames(eanea)[1]<-"X"
colnames(snpgroup)[1]<-"X"
eanea$ea<-substr(eanea$eanea,1,1)
eanea$nea<-substr(eanea$eanea,3,3)

#prepare G-X
gx<-merge(eanea,gx,by="X")
gx<-merge(gx,snpgroup,by="X")
gx<-gx[which(gx$grp6==1),]
gx<-gx[,c("X","ea","nea","EAF","duration_beta","duration_SE","duration_Pval")]
colnames(gx)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval")
write.csv(gx,file="/path/to/the/data")
gxforMR<-format_data(gx,type="exposure")


#prepare G-Y
gy<-gy[-1]
colnames(gy)[1]<-"X"
gy<-merge(eanea,gy,by="X")
write.csv(gy,"/path/to/the/data")

prepareGY<-function(outcome,outcome_SE){
  gy<-read_outcome_data(snps=gx$SNP,filename="/path/to/the/data",sep=",",
                        snp_col="X",beta_col=outcome,se_col=outcome_SE,effect_allele_col="ea",other_allele_col="nea")
  return(gy)}

outcome<-c("loss","stillbirth","miscarriage","gdm","hdp","pd","pt","lbw","hbw","bw")
outcomeSE<-c("lossSE","stillbirthSE","miscarriageSE","gdmSE","hdpSE","pdSE","ptSE","lbwSE","hbwSE","bwSE")
i<-1
for (i in 1:length(outcome)){assign(paste(outcome[i],sep=""),prepareGY(outcome[i],outcomeSE[i]))}

#MR estimates
MRresults<-function(gx,gy){
  gxgyforMR<-harmonise_data(exposure_dat=gx,outcome_dat=gy,action = 1)
  results<-mr(gxgyforMR,
              method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  results$exposure<-deparse(substitute(gx))
  results$outcome<-deparse(substitute(gy))
  return(results)
}

#G-Y in B / G-X in A
a1<-MRresults(gxforMR,loss)
a2<-MRresults(gxforMR,stillbirth)
a3<-MRresults(gxforMR,miscarriage)
a4<-MRresults(gxforMR,gdm)
a5<-MRresults(gxforMR,hdp)
a6<-MRresults(gxforMR,pd)
a7<-MRresults(gxforMR,pt)
a8<-MRresults(gxforMR,lbw)
a9<-MRresults(gxforMR,hbw)
a10<-MRresults(gxforMR,bw)

#output results
tsmr<-rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
tsmr["lci"]<-tsmr["b"]-1.96*tsmr["se"]
tsmr["uci"]<-tsmr["b"]+1.96*tsmr["se"]
tsmr["OR"]<-round(exp(tsmr["b"]),2)
tsmr["ORlci"]<-round(exp(tsmr["lci"]),2)
tsmr["ORuci"]<-round(exp(tsmr["uci"]),2)

write.csv(tsmr,file="/path/to/save/the/results",row.names=FALSE)

############################################################################################################################
#heterogeniety
MRQ<-function(gx,gy){
  gxgyforQ<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  Qresults<-mr_heterogeneity(gxgyforQ,method_list=c("mr_ivw","mr_egger_regression"))
  Qresults$exposure<-deparse(substitute(gx))
  Qresults$outcome<-deparse(substitute(gy))
  return(Qresults)
}

b1<-MRQ(gxforMR,loss)
b2<-MRQ(gxforMR,stillbirth)
b3<-MRQ(gxforMR,miscarriage)
b4<-MRQ(gxforMR,gdm)
b5<-MRQ(gxforMR,hdp)
b6<-MRQ(gxforMR,pd)
b7<-MRQ(gxforMR,pt)
b8<-MRQ(gxforMR,lbw)
b9<-MRQ(gxforMR,hbw)
b10<-MRQ(gxforMR,bw)

tsmr_Q<-rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
write.csv(tsmr_Q,file="/path/to/save/the/results",row.names=FALSE)

############################################################################################################################
#Egger intercept
MRintercept<-function(gx,gy){
  gxgyforEgger<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  Eggerresults<-mr_pleiotropy_test(gxgyforEgger)
  Eggerresults$exposure<-deparse(substitute(gx))
  Eggerresults$outcome<-deparse(substitute(gy))
  return(Eggerresults)
}

c1<-MRintercept(gxforMR,loss)
c2<-MRintercept(gxforMR,stillbirth)
c3<-MRintercept(gxforMR,miscarriage)
c4<-MRintercept(gxforMR,gdm)
c5<-MRintercept(gxforMR,hdp)
c6<-MRintercept(gxforMR,pd)
c7<-MRintercept(gxforMR,pt)
c8<-MRintercept(gxforMR,lbw)
c9<-MRintercept(gxforMR,hbw)
c10<-MRintercept(gxforMR,bw)

tsmr_Egger<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
write.csv(tsmr_Egger,file="/path/to/save/the/results",row.names=FALSE)
###################################################################################################################################
#43 SNPs as IVs
# Clear the work environment
rm(list = ls())

#package
#install.packages("devtools")
#library(devtools)
#Sys.setenv(GITHUB_PAT="/a/path/from/github")
#install_github("MRCIEU/TwoSampleMR")
require(TwoSampleMR)

#input data
gx<-read.csv("/path/to/the/data")
gy<-read.csv("/path/to/the/data")
eanea<-read.csv("/path/to/the/data")
snpgroup<-read.csv("/path/to/the/data")

colnames(eanea)[1]<-"X"
colnames(snpgroup)[1]<-"X"
eanea$ea<-substr(eanea$eanea,1,1)
eanea$nea<-substr(eanea$eanea,3,3)

#prepare G-X
gx<-merge(eanea,gx,by="X")
gx<-merge(gx,snpgroup,by="X")
gx<-gx[which(gx$grp1==1),]
gx<-gx[,c("X","ea","nea","EAF","duration_beta","duration_SE","duration_Pval")]
colnames(gx)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval")
write.csv(gx,file="/path/to/the/data")
gxforMR<-format_data(gx,type="exposure")


#prepare G-Y
gy<-gy[-1]
colnames(gy)[1]<-"X"
gy<-merge(eanea,gy,by="X")
write.csv(gy,"/path/to/the/data")

prepareGY<-function(outcome,outcome_SE){
  gy<-read_outcome_data(snps=gx$SNP,filename="/path/to/the/data",sep=",",
                        snp_col="X",beta_col=outcome,se_col=outcome_SE,effect_allele_col="ea",other_allele_col="nea")
  return(gy)}

#43 SNPs
outcome<-c("loss","stillbirth","miscarriage","gdm","hdp","pd","pt","lbw","hbw","bw")
outcomeSE<-c("lossSE","stillbirthSE","miscarriageSE","gdmSE","hdpSE","pdSE","ptSE","lbwSE","hbwSE","bwSE")
i<-1
for (i in 1:length(outcome)){assign(paste(outcome[i],sep=""),prepareGY(outcome[i],outcomeSE[i]))}

#MR estimates
MRresults<-function(gx,gy){
  gxgyforMR<-harmonise_data(exposure_dat=gx,outcome_dat=gy,action = 1)
  results<-mr(gxgyforMR,
              method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  results$exposure<-deparse(substitute(gx))
  results$outcome<-deparse(substitute(gy))
  return(results)
}

#G-Y in B / G-X in A
a1<-MRresults(gxforMR,loss)
a2<-MRresults(gxforMR,stillbirth)
a3<-MRresults(gxforMR,miscarriage)
a4<-MRresults(gxforMR,gdm)
a5<-MRresults(gxforMR,hdp)
a6<-MRresults(gxforMR,pd)
a7<-MRresults(gxforMR,pt)
a8<-MRresults(gxforMR,lbw)
a9<-MRresults(gxforMR,hbw)
a10<-MRresults(gxforMR,bw)

#output results
tsmr<-rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
tsmr["lci"]<-tsmr["b"]-1.96*tsmr["se"]
tsmr["uci"]<-tsmr["b"]+1.96*tsmr["se"]
tsmr["OR"]<-round(exp(tsmr["b"]),2)
tsmr["ORlci"]<-round(exp(tsmr["lci"]),2)
tsmr["ORuci"]<-round(exp(tsmr["uci"]),2)

write.csv(tsmr,file="/path/to/save/the/results",row.names=FALSE)

############################################################################################################################
#heterogeniety
MRQ<-function(gx,gy){
  gxgyforQ<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  Qresults<-mr_heterogeneity(gxgyforQ,method_list=c("mr_ivw","mr_egger_regression"))
  Qresults$exposure<-deparse(substitute(gx))
  Qresults$outcome<-deparse(substitute(gy))
  return(Qresults)
}

b1<-MRQ(gxforMR,loss)
b2<-MRQ(gxforMR,stillbirth)
b3<-MRQ(gxforMR,miscarriage)
b4<-MRQ(gxforMR,gdm)
b5<-MRQ(gxforMR,hdp)
b6<-MRQ(gxforMR,pd)
b7<-MRQ(gxforMR,pt)
b8<-MRQ(gxforMR,lbw)
b9<-MRQ(gxforMR,hbw)
b10<-MRQ(gxforMR,bw)

tsmr_Q<-rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
write.csv(tsmr_Q,file="/path/to/save/the/results",row.names=FALSE)

############################################################################################################################
#Egger intercept
MRintercept<-function(gx,gy){
  gxgyforEgger<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  Eggerresults<-mr_pleiotropy_test(gxgyforEgger)
  Eggerresults$exposure<-deparse(substitute(gx))
  Eggerresults$outcome<-deparse(substitute(gy))
  return(Eggerresults)
}

c1<-MRintercept(gxforMR,loss)
c2<-MRintercept(gxforMR,stillbirth)
c3<-MRintercept(gxforMR,miscarriage)
c4<-MRintercept(gxforMR,gdm)
c5<-MRintercept(gxforMR,hdp)
c6<-MRintercept(gxforMR,pd)
c7<-MRintercept(gxforMR,pt)
c8<-MRintercept(gxforMR,lbw)
c9<-MRintercept(gxforMR,hbw)
c10<-MRintercept(gxforMR,bw)

tsmr_Egger<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)
write.csv(tsmr_Egger,file="/path/to/save/the/results",row.names=FALSE)
