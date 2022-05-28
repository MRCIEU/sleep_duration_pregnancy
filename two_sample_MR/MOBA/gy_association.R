# Clear the work environment
rm(list = ls())

#package


#input data
moba17k_trios<-read.csv("/path/to/the/data")

#select mum
mumID<-read.table("/path/to/the/data",header=TRUE,sep=",")
mumID<-mumID[,c("M_ID_2552","SENTRIXID","Role","BATCH")]
moba17k_mum1<-merge(mumID,moba17k_trios,by="SENTRIXID")

#QC: relatedness exclusion
#qc1<-read.table("/path/to/the/data")
#qc2<-read.table("/path/to/the/data")
relatedness<-read.table("/path/to/the/data",header=TRUE)
colnames(relatedness)[1]<-"SENTRIXID"
moba17k_mum2<-merge(relatedness,moba17k_mum1,by="SENTRIXID",all.y=TRUE)

#setting 2: 14584 mums, only include unrelated mums
moba17k_mum3<-moba17k_mum2[which(moba17k_mum2$mothers_only_analysis==0),]

#phenotype data
outcome<-read.csv("/path/to/the/data")
moba17k_mum4<-merge(moba17k_mum3,outcome,by="M_ID_2552")

#PC
m12pc<-read.table("/path/to/the/data")
m24pc<-read.table("/path/to/the/data")
rotterdam1pc<-read.table("/path/to/the/data",header=TRUE)
#clean PC
colnames(m12pc)<-c("SENTRIXID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
colnames(m24pc)<-c("SENTRIXID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
rotterdam1pc<-rotterdam1pc[c(2:12)]
colnames(rotterdam1pc)[1]<-"SENTRIXID"
pc<-rbind(m12pc,m24pc)
pc<-rbind(pc,rotterdam1pc)

#add PC
moba17k_mum5<-merge(moba17k_mum4,pc,by="SENTRIXID")

####################################################################################################################################
#sex of children, 1=male, 2=female
table(moba17k_mum5$KJONN,exclude=NULL)
#age
moba17k_mum5$MORS_ALDER[moba17k_mum5$MORS_ALDER==917|moba17k_mum5$MORS_ALDER==945]<-NA
round(mean(moba17k_mum5$MORS_ALDER,na.rm=TRUE),1)
round(sd(moba17k_mum5$MORS_ALDER,na.rm=TRUE),1)
#BMI
round(mean(moba17k_mum5$BMI,na.rm=TRUE),1)
round(sd(moba17k_mum5$BMI,na.rm=TRUE),1)
#height
round(mean(moba17k_mum5$AA87,na.rm=TRUE),1)
round(sd(moba17k_mum5$AA87,na.rm=TRUE),1)
#education, 1=secondary, 2=high school, 3=college or above
moba17k_mum5$edu[moba17k_mum5$AA1124==1]<-1
moba17k_mum5$edu[moba17k_mum5$AA1124==2|moba17k_mum5$AA1124==3|moba17k_mum5$AA1124==4]<-2
moba17k_mum5$edu[moba17k_mum5$AA1124==5|moba17k_mum5$AA1124==6]<-3
table(moba17k_mum5$edu,exclude=NULL)
#income
moba17k_mum5$income<-(moba17k_mum5$incomeM+moba17k_mum5$incomeF)/2
moba17k_mum5$income[is.na(moba17k_mum5$AA1315)&(is.na(moba17k_mum5$AA1316)|moba17k_mum5$AA1316==8)]<-NA
round(mean(moba17k_mum5$income,na.rm=TRUE),1)
round(sd(moba17k_mum5$income,na.rm=TRUE),1)
#parity
round(mean(moba17k_mum5$PARITET_5,na.rm=TRUE),1)
round(sd(moba17k_mum5$PARITET_5,na.rm=TRUE),1)


#G-Y associations
####################################################################################################################################
vars <- names(moba17k_trios[c(2:77)])
col.names<-c("EAF","miscarriage","miscarriageSE","stillbirth","stillbirthSE","loss","lossSE","gdm","gdmSE","hdp","hdpSE","pd","pdSE","pt","ptSE","lbw","lbwSE","hbw","hbwSE","bw","bwSE")
duration<- matrix(, ncol=21, nrow=length(vars))
dimnames(duration) <- list(vars, col.names)

#function
outcomelist<-c("miscarriage_use","stillbirth_use","pregloss","gest_diab","PIH_PE","perinatal_diag","ptb","lbw","hbw")
rungy<-function(mydata,myresults,varlist){
  i<-1
  for (i in 1:9){
    p<-2*i
    q<-2*i+1
    j<-1
    for (j in 1:length(varlist)){
      m<-glm(mydata[,outcomelist[i]] ~ ., data=mydata[,c(varlist[j],"BATCH","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")],family=binomial())
      myresults[j,p]<-round(summary(m)$coefficients[2,"Estimate"],3)
      myresults[j,q]<-round(summary(m)$coefficients[2,"Std. Error"],3)
    }
  }
  return(myresults)
}

#duration
duration<-rungy(moba17k_mum5, duration, vars)

i<-1
for(i in 1:length(vars)) {
  duration[i,1] <- round(mean(moba17k_mum5[[vars[i]]])/2,3)
  m<-lm(VEKT~moba17k_mum5[[vars[i]]]+BATCH+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=moba17k_mum5)
  duration[i,20] <- round(summary(m)$coefficients[2,"Estimate"],3)
  duration[i,21] <- round(summary(m)$coefficients[2,"Std. Error"],3)
}


#output summary results
write.csv(duration,"/path/to/the/data")
