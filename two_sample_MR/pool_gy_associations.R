# Clear the work environment
rm(list = ls())

#package
require(metafor)


#input data
ALSPAC<-read.csv("/path/to/the/data")
bibresults<-read.csv("/path/to/the/data")
moba<-read.csv("/path/to/the/data")
fin<-read.csv("/path/to/the/data")


#clean data
ALSPAC<-ALSPAC[c(1:21)]
bibresults<-bibresults[c(1:21)]
moba<-moba[c(1,3:22)]


#clean gy in BIB
for (i in 1:77){
  if(bibresults$X[i]=="rs7115226"|bibresults$X[i]=="rs11621908"|bibresults$X[i]=="rs34556183"){
    bibresults$miscarriage[i]<-NA
    bibresults$miscarriageSE[i]<-NA
  }
}

for (i in 1:77){
  if(bibresults$X[i]=="rs34556183"){
    bibresults$stillbirth[i]<-NA
    bibresults$stillbirthSE[i]<-NA
  }
}

for (i in 1:77){
  if(bibresults$X[i]=="rs34556183"){
    bibresults$loss[i]<-NA
    bibresults$lossSE[i]<-NA
  }
}

for (i in 1:76){
  if(moba$X[i]=="rs7951019"){
    moba$stillbirth[i]<-NA
    moba$stillbirthSE[i]<-NA
  }
}

gy<-merge(ALSPAC,bibresults,by="X",all.x=TRUE)
gy<-merge(gy,moba,by="X",all.x=TRUE)
gy<-merge(gy,fin,by="X",all.x=TRUE)



#create empty dataframes
poolgy<-data.frame("study"=c("ALSPAC","BIB","MOBA","FIN"),"beta"=c(NA,NA,NA,NA),"se"=c(NA,NA,NA,NA))
vars <- gy[1]
poolresults <- matrix(, ncol=30, nrow=78)
poolresults <-cbind(vars,poolresults)
colnames(poolresults)<-c("snp","miscarriage","miscarriageSE","miscarriagepFORq","stillbirth","stillbirthSE","stillbirthpFORq",
                         "loss","lossSE","losspFORq","gdm","gdmSE","gdmpFORq","hdp","hdpSE","hdppFORq","pd","pdSE","pdpFORq",
                         "pt","ptSE","ptpFORq","lbw","lbwSE","lbwpFORq","hbw","hbwSE","hbwpFORq","bw","bwSE","bwpFORq")

#meta-analysis
i<-1
j<-1
for (i in 1:78){
  for (j in 1:10){
    a<-2*j
    b<-2*j+1
    c<-a+20
    d<-b+20
    e<-a+40
    f<-b+40
    g<-a+60
    h<-b+60
    o<-3*j-1
    p<-3*j
    q<-3*j+1
    #extract estimates
    poolgy[1,2]<-gy[i,a]
    poolgy[1,3]<-gy[i,b]
    poolgy[2,2]<-gy[i,c]
    poolgy[2,3]<-gy[i,d]
    poolgy[3,2]<-gy[i,e]
    poolgy[3,3]<-gy[i,f]
    poolgy[4,2]<-gy[i,g]
    poolgy[4,3]<-gy[i,h]
    #meta
    runpoolgy<-rma(yi=beta,sei=se,slab=study,data=poolgy,method="FE")
    poolresults[i,o]<-runpoolgy$beta
    poolresults[i,p]<-runpoolgy$se
    poolresults[i,q]<-runpoolgy$QEp
  }
}

write.csv(poolresults,file="/path/to/save/gy/associations")

