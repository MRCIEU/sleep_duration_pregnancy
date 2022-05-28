# Clear the work environment
rm(list = ls())

#input data
stime<-read.csv("/path/to/the/data")
head(stime)

#merge with pregnancy outcomes
outcome<-read.csv("/path/to/the/data")
summary(outcome)
colnames(outcome)[1]<-"app"
stime<-merge(stime,outcome,by="app")

#clean duration
table(stime$duration,exclude=NULL)
stime$duration[stime$duration==-3|stime$duration==-1|stime$duration==1|stime$duration>12]<-NA

#remove withdraw & never being pregnant
#new withdraw list
ukbwithdraw<-read.csv("/path/to/the/data")
colnames(ukbwithdraw)<-"app"
ukbwithdraw_check<-merge(ukbwithdraw,stime,by="app")
#new useful data
stime<-stime[-which(stime$app%in%ukbwithdraw$app),]

#remove live birth = 0 & pregnancy loss = 0
stime<-stime[-which(stime$nbirths==0 & stime$pregloss==0),]

#Part 1. G-X associations in all UKB women
####################################################################################################################################
vars <- names(stime[c(3:93)])
col.names <- c("EAF","duration_beta","duration_SE","duration_Pval")
gxall <- matrix(, ncol=4, nrow=length(vars))
dimnames(gxall) <- list(vars, col.names)

i<-1
for(i in 1:length(vars)) {
  gxall[i,1] <- mean(stime[[vars[i]]])/2
  m1<-lm(duration~stime[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
           pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stime)
  gxall[i,2] <- summary(m1)$coefficients[2,"Estimate"]
  gxall[i,3] <- summary(m1)$coefficients[2,"Std. Error"]
  gxall[i,4] <- summary(m1)$coefficients[2,"Pr(>|t|)"]
}

write.csv(gxall,"/path/to/save/gx/associations")


#Part 2. split into sample A and sample B
####################################################################################################################################
set.seed(1234)
n<-nrow(stime)
stime$sample<-rbinom(n, 1, 0.5)
table(stime$sample,exclude=NULL)

#split sample
stimeA<-stime[which(stime$sample==0), ]
stimeB<-stime[which(stime$sample==1), ]


#Part 3. G-X associations in sample A and sample B
####################################################################################################################################
gxA <- matrix(, ncol=4, nrow=length(vars))
dimnames(gxA) <- list(vars, col.names)

gxB <- matrix(, ncol=4, nrow=length(vars))
dimnames(gxB) <- list(vars, col.names)

i<-1
for(i in 1:length(vars)) {
  gxA[i,1] <- mean(stimeA[[vars[i]]])/2
  m1<-lm(duration~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
           pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA)
  gxA[i,2] <- summary(m1)$coefficients[2,"Estimate"]
  gxA[i,3] <- summary(m1)$coefficients[2,"Std. Error"]
  gxA[i,4] <- summary(m1)$coefficients[2,"Pr(>|t|)"]
}

i<-1
for(i in 1:length(vars)) {
  gxB[i,1] <- mean(stimeB[[vars[i]]])/2
  m1<-lm(duration~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
           pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB)
  gxB[i,2] <- summary(m1)$coefficients[2,"Estimate"]
  gxB[i,3] <- summary(m1)$coefficients[2,"Std. Error"]
  gxB[i,4] <- summary(m1)$coefficients[2,"Pr(>|t|)"]
}

write.csv(gxA,"/path/to/save/gx/associations")
write.csv(gxB,"/path/to/save/gx/associations")


#Part 4. G-Y associations in sample A and sample B
####################################################################################################################################
col.names2 <- c("loss","lossSE","stillbirthM","stillbirthMSE","stillbirthS","stillbirthSSE","miscarriageM","miscarriageMSE",
               "miscarriageS","miscarriageSSE","gdmM","gdmMSE","gdmS","gdmSSE","pdM","pdMSE","pdS","pdSSE","lbw","lbwSE",
               "hbw","hbwSE","bw","bwSE","hdp","hdpSE","pt","ptSE")
gyA <- matrix(, ncol=28, nrow=length(vars))
dimnames(gyA) <- list(vars, col.names2)

gyB <- matrix(, ncol=28, nrow=length(vars))
dimnames(gyB) <- list(vars, col.names2)

i<-1
for(i in 1:length(vars)) {
  m1<-glm(pregloss~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
           pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m2<-glm(stillbirth~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m3<-glm(stillbirth_s~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m4<-glm(miscarriage~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m5<-glm(miscarriage_s~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m6<-glm(gdm_define~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m7<-glm(gdm_define_s~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m8<-glm(depress~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m9<-glm(depress_s~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m10<-glm(lbw~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m11<-glm(hbw~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m12<-lm(bwchild~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA)
  m13<-glm(HDP_sDiag~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  m14<-glm(PTD_subsamp_obs~stimeA[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeA,family=binomial())
  gyA[i,1] <- summary(m1)$coefficients[2,"Estimate"]
  gyA[i,2] <- summary(m1)$coefficients[2,"Std. Error"]
  gyA[i,3] <- summary(m2)$coefficients[2,"Estimate"]
  gyA[i,4] <- summary(m2)$coefficients[2,"Std. Error"]
  gyA[i,5] <- summary(m3)$coefficients[2,"Estimate"]
  gyA[i,6] <- summary(m3)$coefficients[2,"Std. Error"]
  gyA[i,7] <- summary(m4)$coefficients[2,"Estimate"]
  gyA[i,8] <- summary(m4)$coefficients[2,"Std. Error"]
  gyA[i,9] <- summary(m5)$coefficients[2,"Estimate"]
  gyA[i,10] <- summary(m5)$coefficients[2,"Std. Error"]
  gyA[i,11] <- summary(m6)$coefficients[2,"Estimate"]
  gyA[i,12] <- summary(m6)$coefficients[2,"Std. Error"]
  gyA[i,13] <- summary(m7)$coefficients[2,"Estimate"]
  gyA[i,14] <- summary(m7)$coefficients[2,"Std. Error"]
  gyA[i,15] <- summary(m8)$coefficients[2,"Estimate"]
  gyA[i,16] <- summary(m8)$coefficients[2,"Std. Error"]
  gyA[i,17] <- summary(m9)$coefficients[2,"Estimate"]
  gyA[i,18] <- summary(m9)$coefficients[2,"Std. Error"]
  gyA[i,19] <- summary(m10)$coefficients[2,"Estimate"]
  gyA[i,20] <- summary(m10)$coefficients[2,"Std. Error"]
  gyA[i,21] <- summary(m11)$coefficients[2,"Estimate"]
  gyA[i,22] <- summary(m11)$coefficients[2,"Std. Error"]
  gyA[i,23] <- summary(m12)$coefficients[2,"Estimate"]
  gyA[i,24] <- summary(m12)$coefficients[2,"Std. Error"]
  gyA[i,25] <- summary(m13)$coefficients[2,"Estimate"]
  gyA[i,26] <- summary(m13)$coefficients[2,"Std. Error"]
  gyA[i,27] <- summary(m14)$coefficients[2,"Estimate"]
  gyA[i,28] <- summary(m14)$coefficients[2,"Std. Error"]
}

i<-1
for(i in 1:length(vars)) {
  m1<-glm(pregloss~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m2<-glm(stillbirth~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m3<-glm(stillbirth_s~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m4<-glm(miscarriage~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m5<-glm(miscarriage_s~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m6<-glm(gdm_define~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m7<-glm(gdm_define_s~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m8<-glm(depress~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m9<-glm(depress_s~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m10<-glm(lbw~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
             pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m11<-glm(hbw~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
             pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m12<-lm(bwchild~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
            pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB)
  m13<-glm(HDP_sDiag~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
             pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  m14<-glm(PTD_subsamp_obs~stimeB[[vars[i]]]+array+age+pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14+pc15+pc16+pc17+pc18+pc19+pc20+pc21+pc22+
             pc23+pc24+pc25+pc26+pc27+pc28+pc29+pc30+pc31+pc32+pc33+pc34+pc35+pc36+pc37+pc38+pc39+pc40,data=stimeB,family=binomial())
  gyB[i,1] <- summary(m1)$coefficients[2,"Estimate"]
  gyB[i,2] <- summary(m1)$coefficients[2,"Std. Error"]
  gyB[i,3] <- summary(m2)$coefficients[2,"Estimate"]
  gyB[i,4] <- summary(m2)$coefficients[2,"Std. Error"]
  gyB[i,5] <- summary(m3)$coefficients[2,"Estimate"]
  gyB[i,6] <- summary(m3)$coefficients[2,"Std. Error"]
  gyB[i,7] <- summary(m4)$coefficients[2,"Estimate"]
  gyB[i,8] <- summary(m4)$coefficients[2,"Std. Error"]
  gyB[i,9] <- summary(m5)$coefficients[2,"Estimate"]
  gyB[i,10] <- summary(m5)$coefficients[2,"Std. Error"]
  gyB[i,11] <- summary(m6)$coefficients[2,"Estimate"]
  gyB[i,12] <- summary(m6)$coefficients[2,"Std. Error"]
  gyB[i,13] <- summary(m7)$coefficients[2,"Estimate"]
  gyB[i,14] <- summary(m7)$coefficients[2,"Std. Error"]
  gyB[i,15] <- summary(m8)$coefficients[2,"Estimate"]
  gyB[i,16] <- summary(m8)$coefficients[2,"Std. Error"]
  gyB[i,17] <- summary(m9)$coefficients[2,"Estimate"]
  gyB[i,18] <- summary(m9)$coefficients[2,"Std. Error"]
  gyB[i,19] <- summary(m10)$coefficients[2,"Estimate"]
  gyB[i,20] <- summary(m10)$coefficients[2,"Std. Error"]
  gyB[i,21] <- summary(m11)$coefficients[2,"Estimate"]
  gyB[i,22] <- summary(m11)$coefficients[2,"Std. Error"]
  gyB[i,23] <- summary(m12)$coefficients[2,"Estimate"]
  gyB[i,24] <- summary(m12)$coefficients[2,"Std. Error"]
  gyB[i,25] <- summary(m13)$coefficients[2,"Estimate"]
  gyB[i,26] <- summary(m13)$coefficients[2,"Std. Error"]
  gyB[i,27] <- summary(m14)$coefficients[2,"Estimate"]
  gyB[i,28] <- summary(m14)$coefficients[2,"Std. Error"]
}

write.csv(gyA,"/path/to/save/gy/associations")
write.csv(gyB,"/path/to/save/gy/associations")

