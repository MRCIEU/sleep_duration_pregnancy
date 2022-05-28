# Clear the work environment
rm(list = ls())

#package
require(metafor)

#input data
AonB<-read.csv("/path/to/the/MR/results")
BonA<-read.csv("/path/to/the/MR/results")

#create empty dataframes
pooldata<-data.frame("study"=c("AonB","BonA"),"beta"=c(NA,NA),"se"=c(NA,NA))
vars <- AonB[,c("outcome","method")]
poolresults <- matrix(, ncol=4, nrow=56)
colnames(poolresults)<-c("beta","se","q","pFORq")
poolresults <-cbind(vars,poolresults)

#meta-analysis
  for (i in 1:56){
    #extract estimates
    pooldata[1,2]<-AonB[i,"b"]
    pooldata[1,3]<-AonB[i,"se"]
    pooldata[2,2]<-BonA[i,"b"]
    pooldata[2,3]<-BonA[i,"se"]
    #meta
    a<-rma(yi=beta, sei=se, slab=study, data=pooldata,method="FE")
    poolresults[i,"beta"]<-a$beta
    poolresults[i,"se"]<-a$se
    poolresults[i,"q"]<-a$QE
    poolresults[i,"pFORq"]<-a$QEp
  }

#output results
output<-function(allresults,resultsname){
  allresults["lci"]<-allresults["beta"]-1.96*allresults["se"]
  allresults["uci"]<-allresults["beta"]+1.96*allresults["se"]
  allresults["OR"]<-round(exp(allresults["beta"]),2)
  allresults["ORlci"]<-round(exp(allresults["lci"]),2)
  allresults["ORuci"]<-round(exp(allresults["uci"]),2)
  write.csv(allresults,paste("/path/to/save/the/results",resultsname,"_v2.csv",sep=""),row.names = FALSE)
}

output(poolresults,"pool")
