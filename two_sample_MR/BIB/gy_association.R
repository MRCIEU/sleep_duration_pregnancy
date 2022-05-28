# Clear the work environment
rm(list = ls())


#set environment vatiables
Sys.setenv(mydata="/path/to/my/data")


#input data
outcome<-read.csv(paste(Sys.getenv('mydata'),'bib_outcome_clean_random_pregnancy.csv',sep=''))
duration<-read.csv("/path/to/my/data")


#merge
bibdata<-merge(duration,outcome,by="MotherID")


#check n
table(bibdata$loss,exclude=NULL)
table(bibdata$stillbirth,exclude=NULL)
table(bibdata$miscarriage,exclude=NULL)
table(bibdata$depression,exclude=NULL)
table(bibdata$gdm,exclude=NULL)
table(bibdata$hdp,exclude=NULL)
table(bibdata$ptb,exclude=NULL)
table(bibdata$lbw,exclude=NULL)
table(bibdata$hbw,exclude=NULL)
summary(bibdata$eclbirthwt)
sd(bibdata$eclbirthwt,na.rm=TRUE)

summary(bibdata$eclregpart)
sd(bibdata$eclregpart,na.rm=TRUE)

#G-Y in BIB
vars <- names(bibdata[c(6:82)])
col.names <- c("miscarriage","miscarriageSE","stillbirth","stillbirthSE","loss","lossSE","gdm","gdmSE","hdp","hdpSE",
               "pdM","pdMSE","pt","ptSE","lbw","lbwSE","hbw","hbwSE","bw","bwSE","eaf")
covar <- matrix(, ncol=21, nrow=length(vars))
dimnames(covar) <- list(vars, col.names)

#binary outcomes
outcome<-c("miscarriage","stillbirth","loss","gdm","hdp","depression","ptb","lbw","hbw")
i<-1
for (i in 1:9){
  p<-2*i-1
  q<-2*i
  j<-1
  for (j in 1:length(vars)){
    m<-glm(bibdata[,outcome[i]] ~ ., data=bibdata[,c(vars[j],"agemy_mbqall","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")],
           family = "binomial")
    covar[j,p] <- summary(m)$coefficients[2,1]
    covar[j,q] <- summary(m)$coefficients[2,2]
  }
}


#continuous coutcome
j<-1
for(j in 1:length(vars)) {
  n<-lm(bibdata$eclbirthwt ~ ., data=bibdata[,c(vars[j],"agemy_mbqall","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
  covar[j,19] <- summary(n)$coefficients[2,1]
  covar[j,20] <- summary(n)$coefficients[2,2]
  covar[j,21] <- mean(bibdata[[vars[j]]])/2
}

write.csv(covar,file="/path/to/save/gy/association")
########################################################################################################################################
#gy associations adjusting for fetal genotype
#input data
outcome<-read.csv(paste(Sys.getenv('mydata'),'bib_outcome_clean_random_pregnancy.csv',sep=''))
duration<-read.csv(paste(Sys.getenv('mydata'),'bib_genotype78_v2_clean.csv',sep=''))
durationF<-read.csv(paste(Sys.getenv('mydata'),'bib_genotype78_v2_clean_fetal.csv',sep=''))


#merge
bibdata<-merge(duration,outcome,by="MotherID")
bibdata2<-merge(bibdata,durationF,by="ChildID")


#G-Y in BIB
vars1 <- names(bibdata2[c(7:83)])
vars2 <- names(bibdata2[c(168:244)])
col.names <- c("eaf_check","miscarriage","miscarriageSE","stillbirth","stillbirthSE","gdm","gdmSE","hdp","hdpSE","pdM","pdMSE","pt","ptSE","lbw","lbwSE","hbw","hbwSE","bw","bwSE")
covar <- matrix(, ncol=19, nrow=length(vars1))
dimnames(covar) <- list(vars1, col.names)


#binary outcomes
outcome<-c("miscarriage","stillbirth","gdm","hdp","depression","ptb","lbw","hbw")
for (i in 1:8){
  p<-2*i
  q<-2*i+1
  j<-1
  for (j in 1:length(vars1)){
    m<-glm(bibdata2[,outcome[i]] ~ ., data=bibdata2[,c(vars1[j],vars2[j],"agemy_mbqall","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")],family = "binomial")
    covar[j,p] <- summary(m)$coefficients[2,1]
    covar[j,q] <- summary(m)$coefficients[2,2]
  }
}


a<-1
for(a in 1:length(vars1)) {
  b<-lm(bibdata2$eclbirthwt ~ ., data=bibdata2[,c(vars1[a],vars2[a],"agemy_mbqall","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
  covar[a,18] <- summary(b)$coefficients[2,1]
  covar[a,19] <- summary(b)$coefficients[2,2]
  covar[a,1] <- mean(bibdata2[[vars1[a]]])/2
}

write.csv(covar,file="/path/to/save/gy/association")

