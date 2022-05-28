# Clear the work environment
rm(list = ls())

#package
require(data.table)
require(reshape2)

#input data
duration<-fread("/path/to/the/data")
information<-duration[,list(chromosome,SNPID,rsid,position,alleleA,alleleB)]
write.csv(information,"/path/to/save/information")
duration<-duration[, c("chromosome","SNPID","position","alleleA","alleleB"):=NULL]


#reshape data
duration_long<-dcast(melt(as.matrix(duration)), Var2~paste0('r', Var1), value.var='value')
colnames(duration_long) <- duration_long[1,]
duration_long <- duration_long[-1, ]
colnames(duration_long)[1]<-"ieu"

#linker
link<-read.csv("/path/to/the/data")

#covariates
covariate<-read.table("path/to/the/data",header = TRUE)
covariate<-covariate[c(2:4)]
colnames(covariate)<-c("app","sex","array")
covariate<-merge(link,covariate,by="app")
covariate$sex<-as.numeric(covariate$sex)
female<-covariate[ which(covariate$sex==1),]
female_duration<-merge(duration_long,female,by="ieu")


pc40<-read.table("/path/to/the/data",header = TRUE)
pc40<-pc40[c(2:42)]
colnames(pc40)<-c("app","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10","pc11","pc12","pc13","pc14","pc15","pc16","pc17","pc18","pc19","pc20","pc21","pc22","pc23","pc24","pc25","pc26","pc27","pc28","pc29","pc30","pc31","pc32","pc33","pc34","pc35","pc36","pc37","pc38","pc39","pc40")
female_duration<-merge(female_duration,pc40,by="app")


#relatedness
minimal<-read.table("/path/to/the/data")
high_related<-read.table("/path/to/the/data")
related<-rbind(minimal,high_related)
colnames(related)<-c("app","app2")
exclude_related<-merge(female_duration,related,by="app",all=TRUE)
female_use<-exclude_related[ which(is.na(exclude_related$app2)),]
female_use<-female_use[c(1:135)]
colnames(female_use)[46]<-"rs1776776"


#check EA, NEA
information<-read.csv("/path/to/the/data")
eanea<-read.csv("/path/to/the/data")
#EA, NEA
eanea$EA<-substr(eanea$eanea,1,1)
eanea$NEA<-substr(eanea$eanea,3,3)
information<-merge(information,eanea,by="rsid")
#0=no need to change, 1=need to recode
i<-1
for (i in 1:91){
  if(information$alleleB[i]==information$EA[i]){
    information$check[i]<-0
  }else{
    information$check[i]<-1
  }
}
write.csv(information,"/path/to/save/the/information")

#recode
female_use$rs1073160<-2-as.numeric(female_use$rs1073160)
female_use$rs10761674<-2-as.numeric(female_use$rs10761674)
female_use$rs11155606<-2-as.numeric(female_use$rs11155606)
female_use$rs11190970<-2-as.numeric(female_use$rs11190970)
female_use$rs112230981<-2-as.numeric(female_use$rs112230981)
female_use$rs113113059<-2-as.numeric(female_use$rs113113059)
female_use$rs11602180<-2-as.numeric(female_use$rs11602180)
female_use$rs11614986<-2-as.numeric(female_use$rs11614986)
female_use$rs11621908<-2-as.numeric(female_use$rs11621908)
female_use$rs12246842<-2-as.numeric(female_use$rs12246842)
female_use$rs12569901<-2-as.numeric(female_use$rs12569901)
female_use$rs12607679<-2-as.numeric(female_use$rs12607679)
female_use$rs12611523<-2-as.numeric(female_use$rs12611523)
female_use$rs1263056<-2-as.numeric(female_use$rs1263056)
female_use$rs13109404<-2-as.numeric(female_use$rs13109404)
female_use$rs17427571<-2-as.numeric(female_use$rs17427571)
female_use$rs17732997<-2-as.numeric(female_use$rs17732997)
female_use$rs1776776<-2-as.numeric(female_use$rs1776776)
female_use$rs180769<-2-as.numeric(female_use$rs180769)
female_use$rs1939455<-2-as.numeric(female_use$rs1939455)
female_use$rs1991556<-2-as.numeric(female_use$rs1991556)
female_use$rs2072727<-2-as.numeric(female_use$rs2072727)
female_use$rs2079070<-2-as.numeric(female_use$rs2079070)
female_use$rs2192528<-2-as.numeric(female_use$rs2192528)
female_use$rs308604<-2-as.numeric(female_use$rs308604)
female_use$rs3095508<-2-as.numeric(female_use$rs3095508)
female_use$rs34354917<-2-as.numeric(female_use$rs34354917)
female_use$rs34556183<-2-as.numeric(female_use$rs34556183)
female_use$rs365663<-2-as.numeric(female_use$rs365663)
female_use$rs374153<-2-as.numeric(female_use$rs374153)
female_use$rs3788337<-2-as.numeric(female_use$rs3788337)
female_use$rs460692<-2-as.numeric(female_use$rs460692)
female_use$rs55658675<-2-as.numeric(female_use$rs55658675)
female_use$rs62120041<-2-as.numeric(female_use$rs62120041)
female_use$rs6575005<-2-as.numeric(female_use$rs6575005)
female_use$rs73219758<-2-as.numeric(female_use$rs73219758)
female_use$rs7503199<-2-as.numeric(female_use$rs7503199)
female_use$rs7616632<-2-as.numeric(female_use$rs7616632)
female_use$rs7644809<-2-as.numeric(female_use$rs7644809)
female_use$rs7806045<-2-as.numeric(female_use$rs7806045)
female_use$rs7915425<-2-as.numeric(female_use$rs7915425)
female_use$rs8038326<-2-as.numeric(female_use$rs8038326)
female_use$rs8050478<-2-as.numeric(female_use$rs8050478)
female_use$rs8074498<-2-as.numeric(female_use$rs8074498)
female_use$rs915416<-2-as.numeric(female_use$rs915416)
female_use$rs9382445<-2-as.numeric(female_use$rs9382445)
female_use$rs9903973<-2-as.numeric(female_use$rs9903973)
female_use$rs9940646<-2-as.numeric(female_use$rs9940646)


write.csv(female_use,"/path/to/save/the/data",row.names = FALSE)

