# Clear the work environment
rm(list = ls())


#input allele information and GY associations
iv<-read.csv("/path/to/data/eanea.csv")
miscarriage<-read.csv("/path/to/data/FinnGen/miscarriage.csv")
gdm<-read.csv("/path/to/data/FinnGen/gdm.csv")
hdp<-read.csv("/path/to/data/FinnGen/hdp.csv")
ptb<-read.csv("/path/to/data/FinnGen/ptb.csv")


#merge data
#EA, NEA
iv$EA<-substr(iv$eanea,1,1)
iv$NEA<-substr(iv$eanea,3,3)
colnames(iv)[1]<-"rsids"
iv<-iv[,c(1,3,4)]

miscarriage1<-miscarriage[,c(1:10)]
gdm1<-gdm[,c(1:10)]
hdp1<-hdp[,c(1:10)]
ptb1<-ptb[,c(1:10)]
miscarriage2<-merge(iv,miscarriage1,by="rsids")
gdm2<-merge(iv,gdm1,by="rsids")
hdp2<-merge(iv,hdp1,by="rsids")
ptb2<-merge(iv,ptb1,by="rsids")

#check EA, NEA
#0=no need to change, 1=need to recode
i<-1
for (i in 1:86){
  if(miscarriage2$alt[i]==miscarriage2$EA[i]){miscarriage2$check[i]<-0}else{miscarriage2$check[i]<-1}
  if(gdm2$alt[i]==gdm2$EA[i]){gdm2$check[i]<-0}else{gdm2$check[i]<-1}
  if(hdp2$alt[i]==hdp2$EA[i]){hdp2$check[i]<-0}else{hdp2$check[i]<-1}
  if(ptb2$alt[i]==ptb2$EA[i]){ptb2$check[i]<-0}else{ptb2$check[i]<-1}
}


#recode estimate & allele frequency
m<-1
for (m in 1:86){
  if(miscarriage2$check[m]==0){miscarriage2$eaf[m]<-miscarriage2$maf[m]
  miscarriage2$betaUSE[m]<-miscarriage2$beta[m]}else{miscarriage2$eaf[m]<-1-miscarriage2$maf[m]
  miscarriage2$betaUSE[m]<-0-miscarriage2$beta[m]}
  if(gdm2$check[m]==0){gdm2$eaf[m]<-gdm2$maf[m]
  gdm2$betaUSE[m]<-gdm2$beta[m]}else{gdm2$eaf[m]<-1-gdm2$maf[m]
  gdm2$betaUSE[m]<-0-gdm2$beta[m]}
  if(hdp2$check[m]==0){hdp2$eaf[m]<-hdp2$maf[m]
  hdp2$betaUSE[m]<-hdp2$beta[m]}else{hdp2$eaf[m]<-1-hdp2$maf[m]
  hdp2$betaUSE[m]<-0-hdp2$beta[m]}
  if(ptb2$check[m]==0){ptb2$eaf[m]<-ptb2$maf[m]
  ptb2$betaUSE[m]<-ptb2$beta[m]}else{ptb2$eaf[m]<-1-ptb2$maf[m]
  ptb2$betaUSE[m]<-0-ptb2$beta[m]}}


#output G-Y associations
write.csv(miscarriage2,"/path/to/data/FinnGen/miscarriage_use.csv")
write.csv(gdm2,"/path/to/data/FinnGen/gdm_use.csv")
write.csv(hdp2,"/path/to/data/FinnGen/hdp_use.csv")
write.csv(ptb2,"/path/to/data/FinnGen/ptb_use.csv")


miscarriage3<-miscarriage2[,c(1,15,11)]
gdm3<-gdm2[,c(1,15,11)]
hdp3<-hdp2[,c(1,15,11)]
ptb3<-ptb2[,c(1,15,11)]

colnames(miscarriage3)<-c("X","Fmiscarriage","FmiscarriageSE")
colnames(gdm3)<-c("X","Fgdm","FgdmSE")
colnames(hdp3)<-c("X","Fhdp","FhdpSE")
colnames(ptb3)<-c("X","Fptb","FptbSE")


pool<-miscarriage3
pool$Fstillbirth<-NA
pool$FstillbirthSE<-NA
pool$Floss<-NA
pool$FlossSE<-NA
pool<-merge(pool,gdm3,by="X")
pool<-merge(pool,hdp3,by="X")
pool$FpdM<-NA
pool$FpdMSE<-NA
pool<-merge(pool,ptb3,by="X")
pool$Flbw<-NA
pool$FlbwSE<-NA
pool$Fhbw<-NA
pool$FhbwSE<-NA
pool$Fbw<-NA
pool$FbwSE<-NA


write.csv(pool,"/path/to/data/FinnGen/gy_use.csv",row.names = FALSE)
