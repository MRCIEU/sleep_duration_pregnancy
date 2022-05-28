# Clear the work environment
rm(list = ls())

#package
require(metafor)

#input data
ukb<-read.csv("/path/to/the/data")
other<-read.csv("/path/to/the/data")


#merge results
ukb$outcome<-substring(ukb$outcome,7)
ukb<-ukb[which(ukb$outcome%in%c("loss","stillbirthM","miscarriageM","gdmM",
                                "hdp","pdM","pt","lbw","hbw","bw")),]
ukb$outcome<-other$outcome
other<-other[-c(1,2,4,6)]
ukbother<-merge(ukb,other,by=c("outcome","method"))

#create empty dataframes
pooldata<-data.frame("study"=c("UKB","Other"),"beta"=c(NA,NA),"se"=c(NA,NA))
vars <- ukbother[,c("outcome","method")]
poolresults <- matrix(, ncol=3, nrow=40)
colnames(poolresults)<-c("beta","se","pFORq")
poolresults <-cbind(vars,poolresults)


#meta-analysis
i<-1
for (i in 1:40){
  #extract estimates
  pooldata[1,2]<-ukbother[i,"beta"]
  pooldata[1,3]<-ukbother[i,"se.x"]
  pooldata[2,2]<-ukbother[i,"b"]
  pooldata[2,3]<-ukbother[i,"se.y"]
  #meta
  a<-rma(yi=beta, sei=se, slab=study, data=pooldata, method="FE")
  poolresults[i,"beta"]<-a$beta
  poolresults[i,"se"]<-a$se
  poolresults[i,"pFORq"]<-a$QEp
}

poolresults["lci"]<-poolresults["beta"]-1.96*poolresults["se"]
poolresults["uci"]<-poolresults["beta"]+1.96*poolresults["se"]
poolresults["OR"]<-round(exp(poolresults["beta"]),2)
poolresults["ORlci"]<-round(exp(poolresults["lci"]),2)
poolresults["ORuci"]<-round(exp(poolresults["uci"]),2)

write.csv(poolresults,file="/path/to/save/the/results",row.names=FALSE)

