# Clear the work environment
rm(list = ls())

#package
require(metafor)

nonlinear<-read.csv("/path/to/the/data",header = FALSE)
colnames(nonlinear)<-c("strata","outcome","beta1","se1","lower1","upper1","beta2","se2","estimate","se","lci","uci","xmean")

nonlinear2<-read.csv("/path/to/the/data",header = FALSE)
colnames(nonlinear2)<-c("strata","outcome","beta1","se1","lower1","upper1","beta2","se2","estimate","se","lci","uci","xmean")


#function
Qp<-function(outcomename){
  group3<-nonlinear[which(nonlinear$outcome==outcomename & nonlinear$strata%in%c("short","healthy","long")),]
  group5<-nonlinear[which(nonlinear$outcome==outcomename & nonlinear$strata%in%c("grp1","grp2","grp3","grp4","grp5")),]
  group3PLUS<-nonlinear2[which(nonlinear2$outcome==outcomename & nonlinear2$strata%in%c("short","healthy","long")),]
  a<-rma(yi=log(as.numeric(group3$estimate)),sei=group3$se)
  b<-rma(log(as.numeric(group3$estimate))~group3$xmean,vi=group3$se^2,method="FE")
  c<-rma(yi=log(as.numeric(group5$estimate)),sei=group5$se)
  d<-rma(log(as.numeric(group5$estimate))~group5$xmean,vi=group5$se^2,method="FE")
  e<-rma(yi=log(as.numeric(group3PLUS$estimate)),sei=group3PLUS$se)
  f<-rma(log(as.numeric(group3PLUS$estimate))~group3PLUS$xmean,vi=group3PLUS$se^2,method="FE")
  return(c(a$QEp,b$pval[2],c$QEp,d$pval[2],e$QEp,f$pval[2]))
}

Qp("stillbirth")
Qp("miscarriage")
Qp("gdm_define")
Qp("HDP_sDiag")
Qp("depress")
Qp("PTD_subsamp_obs")
Qp("lbw")
Qp("hbw")


#birthweight
bw3<-nonlinear[which(nonlinear$outcome=="bwchild" & nonlinear$strata%in%c("short","healthy","long")),]
bw5<-nonlinear[which(nonlinear$outcome=="bwchild" & nonlinear$strata%in%c("grp1","grp2","grp3","grp4","grp5")),]
bw3PLUS<-nonlinear2[which(nonlinear2$outcome=="bwchild" & nonlinear2$strata%in%c("short","healthy","long")),]
a1<-rma(yi=bw3$estimate,sei=bw3$se)
a2<-rma(bw3$estimate~bw3$xmean,vi=bw3$se^2,method="FE")
b1<-rma(yi=bw5$estimate,sei=bw5$se)
b2<-rma(bw5$estimate~bw5$xmean,vi=bw5$se^2,method="FE")
c1<-rma(yi=bw3PLUS$estimate,sei=bw3PLUS$se)
c2<-rma(bw3PLUS$estimate~bw3PLUS$xmean,vi=bw3PLUS$se^2,method="FE")
c(a1$QEp,a2$pval[2],b1$QEp,b2$pval[2],c1$QEp,c2$pval[2])


nonlinear$estimate_use<-round(nonlinear$estimate,2)
nonlinear$lci_use<-round(nonlinear$lci,2)
nonlinear$uci_use<-round(nonlinear$uci,2)

nonlinear2$estimate_use<-round(nonlinear2$estimate,2)
nonlinear2$lci_use<-round(nonlinear2$lci,2)
nonlinear2$uci_use<-round(nonlinear2$uci,2)


#GRS-SLEEP
#3 strata
beta1_3<-nonlinear$beta1[1:3]
se1_3<-nonlinear$se1[1:3]
xmean_3<-nonlinear$xmean[1:3]
results1_3<-rma(yi=beta1_3,sei=se1_3)
results2_3<-rma(beta1_3~xmean_3,vi=se1_3^2,method="FE")

results1_3$QEp
results2_3$pval[2]


#5 strata
beta1_5<-nonlinear$beta1[4:8]
se1_5<-nonlinear$se1[4:8]
xmean_5<-nonlinear$xmean[4:8]
results1_5<-rma(yi=beta1_5,sei=se1_5)
results2_5<-rma(beta1_5~xmean_5,vi=se1_5^2,method="FE")

results1_5$QEp
results2_5$pval[2]


#thirds (sensitivity analysis)
beta1<-nonlinear2$beta1[1:3]
se1<-nonlinear2$se1[1:3]
xmean<-nonlinear2$xmean[1:3]
plus1<-rma(yi=beta1,sei=se1)
plus2<-rma(beta1~xmean,vi=se1^2,method="FE")

plus1$QEp
plus2$pval[2]

