# Clear the work environment
rm(list = ls())

#package
#install.packages("devtools")
#library(devtools)
#Sys.setenv(GITHUB_PAT="/a/path/provided/by/github")
#install_github("MRCIEU/TwoSampleMR")
require(TwoSampleMR)

#input data
gxA<-read.csv("/path/to/the/data")
gxB<-read.csv("/path/to/the/data")
gyA<-read.csv("/path/to/the/data")
gyB<-read.csv("/path/to/the/data")
eanea<-read.csv("/path/to/the/data")
snpgroup<-read.csv("/path/to/the/data")

colnames(eanea)[1]<-"X"
colnames(snpgroup)[1]<-"X"
eanea$ea<-substr(eanea$eanea,1,1)
eanea$nea<-substr(eanea$eanea,3,3)


#prepare G-X
prepareGX<-function(gx,group){
  gx<-merge(eanea,gx,by="X")
  gx<-merge(gx,snpgroup,by="X")
  association<-gx[which(gx[group]==1),]
  association<-association[,c("X","ea","nea","EAF","duration_beta","duration_SE","duration_Pval")]
  colnames(association)<-c("SNP","effect_allele","other_allele","eaf","beta","se","pval")
  gxforMR<-format_data(association,type="exposure")
  return(gxforMR)
}

gxA78<-prepareGX(gxA,"grp4")
gxB78<-prepareGX(gxB,"grp4")
write.csv(gxA78,"/path/to/the/data")
write.csv(gxB78,"/path/to/the/data")


#prepare G-Y
gyA<-merge(eanea,gyA,by="X")
gyB<-merge(eanea,gyB,by="X")
write.csv(gyA,"/path/to/the/data")
write.csv(gyB,"/path/to/the/data")

prepareGY<-function(gxsnp,gyfile,outcome,outcome_SE){
  gy<-read_outcome_data(snps=gxsnp,filename=paste("/path/to/the/data",gyfile,"_v2.csv",sep=""),
    sep=",",snp_col="X",beta_col=outcome,se_col=outcome_SE,effect_allele_col="ea",other_allele_col="nea")
  return(gy)
}

#78 SNPs
outcome<-c("loss","stillbirthM","stillbirthS","miscarriageM","miscarriageS","gdmM","gdmS","hdp","pdM","pdS","pt","lbw","hbw","bw")
outcomeSE<-c("lossSE","stillbirthMSE","stillbirthSSE","miscarriageMSE","miscarriageSSE","gdmMSE","gdmSSE","hdpSE","pdMSE","pdSSE","ptSE","lbwSE","hbwSE","bwSE")
i<-1
for (i in 1:length(outcome)){
  assign(paste("gyA78_",outcome[i],sep=""),prepareGY(gxA78$SNP,"UKB_gyA2",outcome[i],outcomeSE[i]))
}
i<-1
for (i in 1:length(outcome)){
  assign(paste("gyB78_",outcome[i],sep=""),prepareGY(gxA78$SNP,"UKB_gyB2",outcome[i],outcomeSE[i]))
}

#MR estimates
MRresults<-function(gx,gy){
  gxgyforMR<-harmonise_data(exposure_dat=gx,outcome_dat=gy,action = 1)
  results<-mr(gxgyforMR,method_list=c("mr_ivw","mr_weighted_median","mr_weighted_mode","mr_egger_regression"))
  results$exposure<-deparse(substitute(gx))
  results$outcome<-deparse(substitute(gy))
  return(results)
}


#G-Y in B / G-X in A
a1<-MRresults(gxA78,gyB78_loss)
a2<-MRresults(gxA78,gyB78_stillbirthM)
a3<-MRresults(gxA78,gyB78_stillbirthS)
a4<-MRresults(gxA78,gyB78_miscarriageM)
a5<-MRresults(gxA78,gyB78_miscarriageS)
a6<-MRresults(gxA78,gyB78_gdmM)
a7<-MRresults(gxA78,gyB78_gdmS)
a8<-MRresults(gxA78,gyB78_hdp)
a9<-MRresults(gxA78,gyB78_pdM)
a10<-MRresults(gxA78,gyB78_pdS)
a11<-MRresults(gxA78,gyB78_pt)
a12<-MRresults(gxA78,gyB78_lbw)
a13<-MRresults(gxA78,gyB78_hbw)
a14<-MRresults(gxA78,gyB78_bw)


#G-Y in A / G-X in B
b1<-MRresults(gxB78,gyA78_loss)
b2<-MRresults(gxB78,gyA78_stillbirthM)
b3<-MRresults(gxB78,gyA78_stillbirthS)
b4<-MRresults(gxB78,gyA78_miscarriageM)
b5<-MRresults(gxB78,gyA78_miscarriageS)
b6<-MRresults(gxB78,gyA78_gdmM)
b7<-MRresults(gxB78,gyA78_gdmS)
b8<-MRresults(gxB78,gyA78_hdp)
b9<-MRresults(gxB78,gyA78_pdM)
b10<-MRresults(gxB78,gyA78_pdS)
b11<-MRresults(gxB78,gyA78_pt)
b12<-MRresults(gxB78,gyA78_lbw)
b13<-MRresults(gxB78,gyA78_hbw)
b14<-MRresults(gxB78,gyA78_bw)


#output results
AonB<-rbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14)
BonA<-rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14)

output<-function(allresults,resultsname){
  allresults["lci"]<-allresults["b"]-1.96*allresults["se"]
  allresults["uci"]<-allresults["b"]+1.96*allresults["se"]
  allresults["OR"]<-round(exp(allresults["b"]),2)
  allresults["ORlci"]<-round(exp(allresults["lci"]),2)
  allresults["ORuci"]<-round(exp(allresults["uci"]),2)
  write.csv(allresults,paste("/path/to/save/the/results",resultsname,"_v2.csv",sep=""),row.names = FALSE)
}

output(AonB,"AonB")
output(BonA,"BonA")
############################################################################################################################
#heterogeniety
MRQ<-function(gx,gy){
  gxgyforQ<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  Qresults<-mr_heterogeneity(gxgyforQ,method_list=c("mr_ivw","mr_egger_regression"))
  Qresults$exposure<-deparse(substitute(gx))
  Qresults$outcome<-deparse(substitute(gy))
  return(Qresults)
}

#G-Y in B / G-X in A
c1<-MRQ(gxA78,gyB78_loss)
c2<-MRQ(gxA78,gyB78_stillbirthM)
c3<-MRQ(gxA78,gyB78_stillbirthS)
c4<-MRQ(gxA78,gyB78_miscarriageM)
c5<-MRQ(gxA78,gyB78_miscarriageS)
c6<-MRQ(gxA78,gyB78_gdmM)
c7<-MRQ(gxA78,gyB78_gdmS)
c8<-MRQ(gxA78,gyB78_hdp)
c9<-MRQ(gxA78,gyB78_pdM)
c10<-MRQ(gxA78,gyB78_pdS)
c11<-MRQ(gxA78,gyB78_pt)
c12<-MRQ(gxA78,gyB78_lbw)
c13<-MRQ(gxA78,gyB78_hbw)
c14<-MRQ(gxA78,gyB78_bw)

#G-Y in A / G-X in B
d1<-MRQ(gxB78,gyA78_loss)
d2<-MRQ(gxB78,gyA78_stillbirthM)
d3<-MRQ(gxB78,gyA78_stillbirthS)
d4<-MRQ(gxB78,gyA78_miscarriageM)
d5<-MRQ(gxB78,gyA78_miscarriageS)
d6<-MRQ(gxB78,gyA78_gdmM)
d7<-MRQ(gxB78,gyA78_gdmS)
d8<-MRQ(gxB78,gyA78_hdp)
d9<-MRQ(gxB78,gyA78_pdM)
d10<-MRQ(gxB78,gyA78_pdS)
d11<-MRQ(gxB78,gyA78_pt)
d12<-MRQ(gxB78,gyA78_lbw)
d13<-MRQ(gxB78,gyA78_hbw)
d14<-MRQ(gxB78,gyA78_bw)

#putput results
AonB_Q<-rbind(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)
BonA_Q<-rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14)

write.csv(AonB_Q,"/path/to/save/the/results",row.names = FALSE)
write.csv(BonA_Q,"/path/to/save/the/results",row.names = FALSE)
############################################################################################################################
#Egger intercept
MRintercept<-function(gx,gy){
  gxgyforEgger<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  Eggerresults<-mr_pleiotropy_test(gxgyforEgger)
  Eggerresults$exposure<-deparse(substitute(gx))
  Eggerresults$outcome<-deparse(substitute(gy))
  return(Eggerresults)
}

#G-Y in B / G-X in A
e1<-MRintercept(gxA78,gyB78_loss)
e2<-MRintercept(gxA78,gyB78_stillbirthM)
e3<-MRintercept(gxA78,gyB78_stillbirthS)
e4<-MRintercept(gxA78,gyB78_miscarriageM)
e5<-MRintercept(gxA78,gyB78_miscarriageS)
e6<-MRintercept(gxA78,gyB78_gdmM)
e7<-MRintercept(gxA78,gyB78_gdmS)
e8<-MRintercept(gxA78,gyB78_hdp)
e9<-MRintercept(gxA78,gyB78_pdM)
e10<-MRintercept(gxA78,gyB78_pdS)
e11<-MRintercept(gxA78,gyB78_pt)
e12<-MRintercept(gxA78,gyB78_lbw)
e13<-MRintercept(gxA78,gyB78_hbw)
e14<-MRintercept(gxA78,gyB78_bw)

#G-Y in A / G-X in B
f1<-MRintercept(gxB78,gyA78_loss)
f2<-MRintercept(gxB78,gyA78_stillbirthM)
f3<-MRintercept(gxB78,gyA78_stillbirthS)
f4<-MRintercept(gxB78,gyA78_miscarriageM)
f5<-MRintercept(gxB78,gyA78_miscarriageS)
f6<-MRintercept(gxB78,gyA78_gdmM)
f7<-MRintercept(gxB78,gyA78_gdmS)
f8<-MRintercept(gxB78,gyA78_hdp)
f9<-MRintercept(gxB78,gyA78_pdM)
f10<-MRintercept(gxB78,gyA78_pdS)
f11<-MRintercept(gxB78,gyA78_pt)
f12<-MRintercept(gxB78,gyA78_lbw)
f13<-MRintercept(gxB78,gyA78_hbw)
f14<-MRintercept(gxB78,gyA78_bw)

#putput results
AonB_Egger<-rbind(e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14)
BonA_Egger<-rbind(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14)

write.csv(AonB_Egger,"/path/to/save/the/results",row.names = FALSE)
write.csv(BonA_Egger,"/path/to/save/the/results",row.names = FALSE)

############################################################################################################################
#Leave-one-out
MRfig<-function(gx,gy){
  gxgyforfigure<-harmonise_data(exposure_dat = gx,outcome_dat = gy,action=1)
  res_loo<-mr_leaveoneout(gxgyforfigure)
  p3<-mr_leaveoneout_plot(res_loo)
  return(p3[[1]])
}


#G-Y in B / G-X in A
MRfig(gxA78,gyB78_loss)
MRfig(gxA78,gyB78_stillbirthM)
MRfig(gxA78,gyB78_miscarriageM)
MRfig(gxA78,gyB78_gdmM)
MRfig(gxA78,gyB78_hdp)
MRfig(gxA78,gyB78_pdM)
MRfig(gxA78,gyB78_pt)
MRfig(gxA78,gyB78_lbw)
MRfig(gxA78,gyB78_hbw)
MRfig(gxA78,gyB78_bw)

#G-Y in A / G-X in B
MRfig(gxB78,gyA78_loss)
MRfig(gxB78,gyA78_stillbirthM)
MRfig(gxB78,gyA78_miscarriageM)
MRfig(gxB78,gyA78_gdmM)
MRfig(gxB78,gyA78_hdp)
MRfig(gxB78,gyA78_pdM)
MRfig(gxB78,gyA78_pt)
MRfig(gxB78,gyA78_lbw)
MRfig(gxB78,gyA78_hbw)
MRfig(gxB78,gyA78_bw)

