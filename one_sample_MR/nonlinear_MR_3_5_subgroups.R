# Clear the work environment
rm(list = ls())

#package
require(metafor)
require(ggplot2)
require(AER)

#prepare environment
#setwd("/path/to/the/data")
Sys.setenv(UKB="/path/to/the/data")
Sys.setenv(mydata="/path/to/the/data")

#input data
ukbg<-read.csv(paste(Sys.getenv('UKB'),'duration/duration_long.csv',sep=''))
ukbp<-read.csv(paste(Sys.getenv('UKB'),'ukbb-sleep/phenotype/outcome.csv',sep=''))
colnames(ukbp)[1]<-"app"
ukb<-merge(ukbg,ukbp,by="app")

#clean duration
table(ukb$duration,exclude=NULL)
ukb$duration[ukb$duration==-3|ukb$duration==-1|ukb$duration==1|ukb$duration>12]<-NA
ukb<-ukb[which(!is.na(ukb$duration)),]

##############################################################################################
#remove withdraw & never being pregnant
#new withdraw list
ukbwithdraw<-read.csv("/path/to/the/data")
colnames(ukbwithdraw)<-"app"
ukbwithdraw_check<-merge(ukbwithdraw,ukb,by="app")
#new useful data
ukb<-ukb[-which(ukb$app%in%ukbwithdraw$app),]

#remove live birth = 0 & pregnancy loss = 0
ukb<-ukb[-which(ukb$nbirths==0 & ukb$pregloss==0),]
#############################################################################################

#unweighted allele score for IV
snpgroup<-read.csv(paste(Sys.getenv('mydata'),'snpgroup.csv',sep=''))
snpgroup<-snpgroup[which(snpgroup$grp4==1),]
a<-as.character(snpgroup$SNP)
ukb$grs<-rowSums(ukb[,c(a)])

#split into 3/5 strata
b<-lm(duration~grs+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+
        pc21+pc22+pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=ukb)
summary(b)$coefficients[2,1]
summary(b)$coefficients[2,1]-1.96*summary(b)$coefficients[2,2]
summary(b)$coefficients[2,1]+1.96*summary(b)$coefficients[2,2]
mean0<-mean(b$fit)
ukb$duration0<-ukb$duration-b$fit+mean0
hist(ukb$duration0,breaks=100,xlim=c(1.5,12.5),ylim=c(0,30000),main="",xlab="residual sleep duration (hours/day)")

##############################################################################################################################
#3 strata
ukb$group3[ukb$duration0<7]<-1
ukb$group3[ukb$duration0>=7&ukb$duration0<9]<-2
ukb$group3[ukb$duration0>=9]<-3
table(ukb$group3,ukb$duration)
short<-ukb[which(ukb$group3==1),]
healthy<-ukb[which(ukb$group3==2),]
long<-ukb[which(ukb$group3==3),]
#5 strata
ukb$group5[ukb$duration0<6]<-1
ukb$group5[ukb$duration0>=6&ukb$duration0<7]<-2
ukb$group5[ukb$duration0>=7&ukb$duration0<8]<-3
ukb$group5[ukb$duration0>=8&ukb$duration0<9]<-4
ukb$group5[ukb$duration0>=9]<-5
table(ukb$group5,ukb$duration)
grp1<-ukb[which(ukb$group5==1),]
grp2<-ukb[which(ukb$group5==2),]
grp3<-ukb[which(ukb$group5==3),]
grp4<-ukb[which(ukb$group5==4),]
grp5<-ukb[which(ukb$group5==5),]
########################################################################################################################
#sample size
table(ukb$group5,ukb$stillbirth,exclude=NULL)
table(ukb$group5,ukb$miscarriage,exclude=NULL)
table(ukb$group5,ukb$gdm_define,exclude=NULL)
table(ukb$group5,ukb$HDP_sDiag,exclude=NULL)
table(ukb$group5,ukb$depress,exclude=NULL)
table(ukb$group5,ukb$PTD_subsamp_obs,exclude=NULL)
table(ukb$group5,ukb$lbw,exclude=NULL)
table(ukb$group5,ukb$hbw,exclude=NULL)
summary(grp1$bwchild)
summary(grp2$bwchild)
summary(grp3$bwchild)
summary(grp4$bwchild)
summary(grp5$bwchild)

table(ukb$group3,ukb$stillbirth,exclude=NULL)
table(ukb$group3,ukb$miscarriage,exclude=NULL)
table(ukb$group3,ukb$gdm_define,exclude=NULL)
table(ukb$group3,ukb$HDP_sDiag,exclude=NULL)
table(ukb$group3,ukb$depress,exclude=NULL)
table(ukb$group3,ukb$PTD_subsamp_obs,exclude=NULL)
table(ukb$group3,ukb$lbw,exclude=NULL)
table(ukb$group3,ukb$hbw,exclude=NULL)
summary(short$bwchild)
summary(healthy$bwchild)
summary(long$bwchild)
########################################################################################################################
#function
nonlinear<-function(database,outcomename){
  firststage<-lm(database[,"duration"] ~., data=database[,c("grs","array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
                                                            "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20",
                                                            "pc21","pc22","pc23","pc24","pc25","pc26","pc27","pc28","pc29","pc30",
                                                            "pc31","pc32","pc33","pc34","pc35","pc36","pc37","pc38","pc39","pc40")])

  secondstage<-glm(database[,outcomename] ~., data=database[,c("grs","array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
                                                               "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20",
                                                               "pc21","pc22","pc23","pc24","pc25","pc26","pc27","pc28","pc29","pc30",
                                                               "pc31","pc32","pc33","pc34","pc35","pc36","pc37","pc38","pc39","pc40")],family=binomial())
  logor1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  lower1<-logor1-1.96*se1
  upper1<-logor1+1.96*se1
  logor2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  logor<-logor2/logor1
  se<-se2/logor1
  lower = exp(logor-1.96*se)
  upper = exp(logor+1.96*se)
  xmean = mean(database[,"duration"])
  write.table(cbind(deparse(substitute(database)),outcomename,logor1,se1,lower1,upper1,logor2,se2,exp(logor),se,lower,upper,xmean),file=paste(Sys.getenv('mydata'),'UKB_nonlinear_ratio_v2.csv',sep=''), append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}


outcomelist<-c("pregloss","stillbirth","stillbirth_s","miscarriage","miscarriage_s","gdm_define","gdm_define_s","depress","depress_s","lbw","hbw","HDP_sDiag","PTD_subsamp_obs")

for (i in 1:length(outcomelist)){
  nonlinear(short,outcomelist[i])
  nonlinear(healthy,outcomelist[i])
  nonlinear(long,outcomelist[i])
  nonlinear(grp1,outcomelist[i])
  nonlinear(grp2,outcomelist[i])
  nonlinear(grp3,outcomelist[i])
  nonlinear(grp4,outcomelist[i])
  nonlinear(grp5,outcomelist[i])
}


#birthweight (continuous)
nonlinearC<-function(database){
  firststage<-lm(database[,"duration"] ~., data=database[,c("grs","array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
                                                            "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20",
                                                            "pc21","pc22","pc23","pc24","pc25","pc26","pc27","pc28","pc29","pc30",
                                                            "pc31","pc32","pc33","pc34","pc35","pc36","pc37","pc38","pc39","pc40")])

  secondstage<-lm(database[,"bwchild"] ~., data=database[,c("grs","array","age","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10",
                                                            "pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20",
                                                            "pc21","pc22","pc23","pc24","pc25","pc26","pc27","pc28","pc29","pc30",
                                                            "pc31","pc32","pc33","pc34","pc35","pc36","pc37","pc38","pc39","pc40")])

  beta1<-summary(firststage)$coefficient[2,"Estimate"]
  se1<-summary(firststage)$coefficient[2,"Std. Error"]
  lower1<-beta1-1.96*se1
  upper1<-beta1+1.96*se1
  beta2<-summary(secondstage)$coefficient[2,"Estimate"]
  se2<-summary(secondstage)$coefficient[2,"Std. Error"]
  beta<-beta2/beta1
  se<-se2/beta1
  lower = beta-1.96*se
  upper = beta+1.96*se
  xmean = mean(database[,"duration"])
  write.table(cbind(deparse(substitute(database)),"bwchild",beta1,se1,lower1,upper1,beta2,se2,beta,se,lower,upper,xmean),file=paste(Sys.getenv('mydata'),'UKB_nonlinear_ratio_v2.csv',sep=''), append=TRUE, quote=FALSE, sep=',',row.names=FALSE, col.names=FALSE)
}

nonlinearC(short)
nonlinearC(healthy)
nonlinearC(long)
nonlinearC(grp1)
nonlinearC(grp2)
nonlinearC(grp3)
nonlinearC(grp4)
nonlinearC(grp5)

