#clear the work environment
rm(list = ls())

#package
require(lmtest)

#input data
duration<-read.csv("/path/to/the/data")

#sleep duration
table(duration$CC1027,exclude=NULL)
#>=10 hours=1, 8-9 hours=2, 6-7 hours=3, 4-5hours=4, less than 4 hours=5, 4-5 & less than 4 hours=9
duration$sleep[duration$CC1027==1]<-1
duration$sleep[duration$CC1027==2]<-0
duration$sleep[duration$CC1027==3]<--1
duration$sleep[duration$CC1027==4|duration$CC1027==5|duration$CC1027==9]<--2
table(duration$sleep,exclude=NULL)

#exclude participants whose sleep duration is NA
duration<-duration[which(!is.na(duration$sleep)),]
duration$sleepL[duration$sleep==1]<-11
duration$sleepL[duration$sleep==0]<-8.5
duration$sleepL[duration$sleep==-1]<-6.5
duration$sleepL[duration$sleep==-2]<-3.5
duration$sleepL<-as.numeric(duration$sleepL)
duration$sleepNL<-as.factor(duration$sleep)

#add parity as a confounder
table(duration$PARITET_5,exclude = NULL)
duration_adj<-duration[which(!is.na(duration$BMI)&!is.na(duration$age)&!is.na(duration$smoke)&!is.na(duration$alcohol)&!is.na(duration$income)&!is.na(duration$AA1124)&!is.na(duration$PARITET_5)),]

#crude model (no changes)
crude<-function(outcomename){
  fit1=glm(duration[,outcomename]~sleepL,data=duration,family=binomial())
  sumx1=summary(fit1)
  logor1=sumx1$coefficient[2,"Estimate"]
  se1=sumx1$coefficient[2,"Std. Error"]
  p1=sumx1$coefficient[2,"Pr(>|z|)"]
  fit2=glm(duration[,outcomename]~relevel(factor(sleepNL),ref="0"),data=duration,family=binomial())
  sumx2=summary(fit2)
  logor_group1=sumx2$coefficient[2,"Estimate"]
  se_group1=sumx2$coefficient[2,"Std. Error"]
  logor_group2=sumx2$coefficient[3,"Estimate"]
  se_group2=sumx2$coefficient[3,"Std. Error"]
  logor_group3=sumx2$coefficient[4,"Estimate"]
  se_group3=sumx2$coefficient[4,"Std. Error"]
  a<-lrtest(fit1,fit2)
  b<-a$`Pr(>Chisq)`[2]
  write.table(cbind(outcomename,logor1,se1,p1,b),
              file="/path/to/save/the/results",
              append = TRUE,quote=FALSE,sep=",",col.names = FALSE,row.names =FALSE)
  write.table(cbind(outcomename,logor_group1,se_group1,logor_group2,se_group2,logor_group3,se_group3),
              file="/path/to/save/the/results",
              append = TRUE,quote=FALSE,sep=",",col.names = FALSE,row.names =FALSE)
  
}

crude("stillbirth_use")
crude("miscarriage_use")
crude("gest_diab")
crude("PIH_PE")
crude("perinatal_diag")
crude("ptb")
crude("lbw")
crude("hbw")

#BW
bwc1<-lm(duration[,"VEKT"]~sleepL,data=duration)
bwc2<-lm(duration[,"VEKT"]~relevel(factor(sleepNL),ref="0"),data=duration)
round(summary(bwc1)$coefficient[2,"Estimate"],2)
round(summary(bwc1)$coefficient[2,"Estimate"]-1.96*summary(bwc1)$coefficient[2,"Std. Error"],2)
round(summary(bwc1)$coefficient[2,"Estimate"]+1.96*summary(bwc1)$coefficient[2,"Std. Error"],2)
summary(bwc1)$coefficient[2,"Pr(>|t|)"]
round(summary(bwc2)$coefficient[2,"Estimate"],2)
round(summary(bwc2)$coefficient[2,"Estimate"]-1.96*summary(bwc2)$coefficient[2,"Std. Error"],2)
round(summary(bwc2)$coefficient[2,"Estimate"]+1.96*summary(bwc2)$coefficient[2,"Std. Error"],2)
round(summary(bwc2)$coefficient[3,"Estimate"],2)
round(summary(bwc2)$coefficient[3,"Estimate"]-1.96*summary(bwc2)$coefficient[3,"Std. Error"],2)
round(summary(bwc2)$coefficient[3,"Estimate"]+1.96*summary(bwc2)$coefficient[3,"Std. Error"],2)
round(summary(bwc2)$coefficient[4,"Estimate"],2)
round(summary(bwc2)$coefficient[4,"Estimate"]-1.96*summary(bwc2)$coefficient[4,"Std. Error"],2)
round(summary(bwc2)$coefficient[4,"Estimate"]+1.96*summary(bwc2)$coefficient[4,"Std. Error"],2)
lrtest(bwc1,bwc2)


#adjust model (add one more confounder)
adjusted<-function(outcomename){
  fit1=glm(duration[,outcomename]~sleepL+BMI+age+smoke+alcohol+income+AA1124+PARITET_5,data=duration,family=binomial())
  sumx1=summary(fit1)
  logor1=sumx1$coefficient[2,"Estimate"]
  se1=sumx1$coefficient[2,"Std. Error"]
  p1=sumx1$coefficient[2,"Pr(>|z|)"]
  fit2=glm(duration[,outcomename]~relevel(factor(sleepNL),ref="0")+BMI+age+smoke+alcohol+income+AA1124+PARITET_5,data=duration,family=binomial())
  sumx2=summary(fit2)
  logor_group1=sumx2$coefficient[2,"Estimate"]
  se_group1=sumx2$coefficient[2,"Std. Error"]
  logor_group2=sumx2$coefficient[3,"Estimate"]
  se_group2=sumx2$coefficient[3,"Std. Error"]
  logor_group3=sumx2$coefficient[4,"Estimate"]
  se_group3=sumx2$coefficient[4,"Std. Error"]
  c<-lrtest(fit1,fit2)
  d<-c$`Pr(>Chisq)`[2]
  write.table(cbind(outcomename,logor1,se1,p1,d),
              file="/path/to/save/the/results",
              append = TRUE,quote=FALSE,sep=",",col.names = FALSE,row.names =FALSE)
  write.table(cbind(outcomename,logor_group1,se_group1,logor_group2,se_group2,logor_group3,se_group3),
              file="/path/to/save/the/results",
              append = TRUE,quote=FALSE,sep=",",col.names = FALSE,row.names =FALSE)
}
adjusted("pregloss")
adjusted("stillbirth_use")
adjusted("miscarriage_use")
adjusted("gest_diab")
adjusted("PIH_PE")
adjusted("perinatal_diag")
adjusted("ptb")
adjusted("lbw")
adjusted("hbw")


#BW
bwc1<-lm(duration[,"VEKT"]~sleepL+BMI+age+smoke+alcohol+income+AA1124+PARITET_5,data=duration)
bwc2<-lm(duration[,"VEKT"]~relevel(factor(sleepNL),ref="0")+BMI+age+smoke+alcohol+income+AA1124+PARITET_5,data=duration)
round(summary(bwc1)$coefficient[2,"Estimate"],2)
round(summary(bwc1)$coefficient[2,"Estimate"]-1.96*summary(bwc1)$coefficient[2,"Std. Error"],2)
round(summary(bwc1)$coefficient[2,"Estimate"]+1.96*summary(bwc1)$coefficient[2,"Std. Error"],2)
summary(bwc1)$coefficient[2,"Pr(>|t|)"]
round(summary(bwc2)$coefficient[2,"Estimate"],2)
round(summary(bwc2)$coefficient[2,"Estimate"]-1.96*summary(bwc2)$coefficient[2,"Std. Error"],2)
round(summary(bwc2)$coefficient[2,"Estimate"]+1.96*summary(bwc2)$coefficient[2,"Std. Error"],2)
round(summary(bwc2)$coefficient[3,"Estimate"],2)
round(summary(bwc2)$coefficient[3,"Estimate"]-1.96*summary(bwc2)$coefficient[3,"Std. Error"],2)
round(summary(bwc2)$coefficient[3,"Estimate"]+1.96*summary(bwc2)$coefficient[3,"Std. Error"],2)
round(summary(bwc2)$coefficient[4,"Estimate"],2)
round(summary(bwc2)$coefficient[4,"Estimate"]-1.96*summary(bwc2)$coefficient[4,"Std. Error"],2)
round(summary(bwc2)$coefficient[4,"Estimate"]+1.96*summary(bwc2)$coefficient[4,"Std. Error"],2)
lrtest(bwc1,bwc2)


#results
crude_results1<-read.csv("/path/to/the/results",header=FALSE)
colnames(crude_results1)<-c("outcome","logor","se","p","likelihoodp")
crude_results1$OR<-round(exp(crude_results1$logor),2)
crude_results1$lci<-round(exp(crude_results1$logor-1.96*crude_results1$se),2)
crude_results1$uci<-round(exp(crude_results1$logor+1.96*crude_results1$se),2)

crude_results2<-read.csv("/path/to/the/results",header=FALSE)
colnames(crude_results2)<-c("outcome","logor1","se1","logor2","se2","logor3","se3")
crude_results2$OR1<-round(exp(crude_results2$logor1),2)
crude_results2$lci1<-round(exp(crude_results2$logor1-1.96*crude_results2$se1),2)
crude_results2$uci1<-round(exp(crude_results2$logor1+1.96*crude_results2$se1),2)
crude_results2$OR2<-round(exp(crude_results2$logor2),2)
crude_results2$lci2<-round(exp(crude_results2$logor2-1.96*crude_results2$se2),2)
crude_results2$uci2<-round(exp(crude_results2$logor2+1.96*crude_results2$se2),2)
crude_results2$OR3<-round(exp(crude_results2$logor3),2)
crude_results2$lci3<-round(exp(crude_results2$logor3-1.96*crude_results2$se3),2)
crude_results2$uci3<-round(exp(crude_results2$logor3+1.96*crude_results2$se3),2)

adjusted_results1<-read.csv("/path/to/the/results",header=FALSE)
colnames(adjusted_results1)<-c("outcome","logor","se","p","likelihoodp")
adjusted_results1$OR<-round(exp(adjusted_results1$logor),2)
adjusted_results1$lci<-round(exp(adjusted_results1$logor-1.96*adjusted_results1$se),2)
adjusted_results1$uci<-round(exp(adjusted_results1$logor+1.96*adjusted_results1$se),2)

adjusted_results2<-read.csv("/path/to/the/results",header=FALSE)
colnames(adjusted_results2)<-c("outcome","logor1","se1","logor2","se2","logor3","se3")
adjusted_results2$OR1<-round(exp(adjusted_results2$logor1),2)
adjusted_results2$lci1<-round(exp(adjusted_results2$logor1-1.96*adjusted_results2$se1),2)
adjusted_results2$uci1<-round(exp(adjusted_results2$logor1+1.96*adjusted_results2$se1),2)
adjusted_results2$OR2<-round(exp(adjusted_results2$logor2),2)
adjusted_results2$lci2<-round(exp(adjusted_results2$logor2-1.96*adjusted_results2$se2),2)
adjusted_results2$uci2<-round(exp(adjusted_results2$logor2+1.96*adjusted_results2$se2),2)
adjusted_results2$OR3<-round(exp(adjusted_results2$logor3),2)
adjusted_results2$lci3<-round(exp(adjusted_results2$logor3-1.96*adjusted_results2$se3),2)
adjusted_results2$uci3<-round(exp(adjusted_results2$logor3+1.96*adjusted_results2$se3),2)

