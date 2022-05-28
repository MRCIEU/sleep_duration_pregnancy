# Clear the work environment
rm(list = ls())


#set environment vatiables
Sys.setenv(BiB="/path/to/the/data")
Sys.setenv(mydata="/path/to/the/data")

#input data
gene<-read.csv(paste(Sys.getenv('mydata'),'bib_genotype78_v2.csv',sep=''))
mumID<-read.csv(paste(Sys.getenv('BiB'),'data.mother.linkage.csv',sep=''))

whiteID<-read.table(paste(Sys.getenv('BiB'),'data.white.txt',sep=''))
#pID<-read.table(paste(Sys.getenv('BiB'),'data.pakist.txt',sep=''))
#deminf<-read.csv(paste(Sys.getenv('BiB'),'Mother_Baseline_Questionnaire_v7_deminf_Data.csv',sep=''))
#qc_sex<-read.table(paste(Sys.getenv('BiB'),'data.sexcheckfails_dropped.txt',sep=''))
degreefail<-read.table(paste(Sys.getenv('BiB'),'data.firstdegfails_mumch.txt',sep=''),header=TRUE)
#qc_duplicated<-read.table(paste(Sys.getenv('BiB'),'data.dups_127ofthemdropped.txt',sep=''),header=TRUE)
pc20<-read.table(paste(Sys.getenv('mydata'),'directly_genotyped_released_2016-09-06_pca.eigenvec',sep=''),header=TRUE)

#select mum
mumID<-mumID[c("PregnancyID","SentrixID","PersonID","ParticipantType","HasReplicate")]
colnames(mumID)<-c("PregnancyID","gwassentrixps","MotherID","ParticipantType","HasReplicate")
gene_mum<-merge(mumID,gene,by="gwassentrixps")
#gene_mum_check<-gene_mum[with(gene_mum, order(-gene_mum$HasReplicate, gene_mum$MotherID)),]

#white European only
colnames(whiteID)<-c("gwassentrixps","whiteEuropean")
gene_White<-merge(gene_mum,whiteID,by="gwassentrixps")
#colnames(pID)<-c("gwassentrixps","whiteEuropean")
#gene_ancestry<-merge(gene_mum,pID,by="gwassentrixps",all.x=TRUE)
#gene_ancestry<-merge(gene_ancestry,whiteID,by="gwassentrixps",all.x=TRUE)
#deminf<-deminf[c("MotherID","deminfeth3gpcomb")]
#gene_mum2<-merge(gene_mum,deminf,by="MotherID")
#table(gene_mum2$deminfeth3gpcomb,exclude=NULL)
#gene_mum2<-gene_mum2[which(gene_mum2$deminfeth3gpcomb=="White British"),]

#QC
#sex mismatched=0
#colnames(qc_sex)<-"gwassentrixps"
#qc_sex$wrongsex<-1
#gene_mum_sex<-merge(gene_White,qc_sex,by="gwassentrixps",all=TRUE)

#first degree fail mum-child
mumdelete<-as.data.frame(unique(degreefail[,"gwassentrix_mum"]))
colnames(mumdelete)<-"gwassentrixps"
mumdelete$mumdelete<-1
gene_mum_degree<-merge(gene_White,mumdelete,by="gwassentrixps",all.x=TRUE)
gene_mum2<-gene_mum_degree[which(is.na(gene_mum_degree$mumdelete)),]

#duplicated
gene_mum3<-gene_mum2[with(gene_mum2, order(gene_mum2$MotherID,gene_mum2$PregnancyID)),]
duplicatedmum<-as.data.frame(unique(gene_mum3[,"MotherID"]))


#Add PCs
colnames(pc20)[1]<-"gwassentrixps"
pc20<-pc20[c(1,3:22)]
gene_mum4<-merge(gene_mum3,pc20,by="gwassentrixps")

#output
write.csv(gene_mum4,file=paste(Sys.getenv('mydata'),'bib_genotype78_v2_clean.csv',sep=''),row.names=FALSE)

####################################################################################################################
#for fetal genotype
#input data
gene<-read.csv(paste(Sys.getenv('mydata'),'bib_genotype78_v2.csv',sep=''))
childID<-read.csv(paste(Sys.getenv('BiB'),'data.child.linkage.csv',sep=''))


#select child
childID<-childID[c("PregnancyID","SentrixID","PersonID","ParticipantType","HasReplicate")]
colnames(childID)<-c("PregnancyID","gwassentrixps","ChildID","ParticipantType_C","HasReplicate_C")
gene_child<-merge(childID,gene,by="gwassentrixps")
#gene_child_check<-gene_child[with(gene_child, order(-gene_child$HasReplicate, gene_child$PregnancyID)),]


write.csv(gene_child,file=paste(Sys.getenv('mydata'),'bib_genotype78_v2_clean_fetal.csv',sep=''),row.names=FALSE)
