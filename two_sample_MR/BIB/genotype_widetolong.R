# Clear the work environment
rm(list = ls())

#set environment vatiables
Sys.setenv(mydata="/path/to/the/data")

#package
require(data.table)

#input data
gene<-fread("/path/to/the/data")
snp1<-read.table("/path/to/the/data")

#reshape data
information<-gene[,list(chromosome,SNPID,rsid,position,alleleA,alleleB)]
gene_use<-gene[, c("chromosome","SNPID","position","alleleA","alleleB"):=NULL]
gene_use_long<-dcast(melt(as.matrix(gene_use)), Var2~paste0('r', Var1), value.var='value')
colnames(gene_use_long) <- gene_use_long[1,]
gene_use_long <- gene_use_long[-1, ]
colnames(gene_use_long)[1]<-"gwassentrixps"

#EA, NEA
i<-2
for (i in 2:77){gene_use_long[,i]<-as.numeric(gene_use_long[,i])}
gene_use_long$rs10761674<-2-gene_use_long$rs10761674
gene_use_long$rs11190970<-2-gene_use_long$rs11190970
gene_use_long$rs112230981<-2-gene_use_long$rs112230981
gene_use_long$rs113113059<-2-gene_use_long$rs113113059
gene_use_long$rs11602180<-2-gene_use_long$rs11602180
gene_use_long$rs11614986<-2-gene_use_long$rs11614986
gene_use_long$rs11621908<-2-gene_use_long$rs11621908
gene_use_long$rs12246842<-2-gene_use_long$rs12246842
gene_use_long$rs12607679<-2-gene_use_long$rs12607679
gene_use_long$rs12611523<-2-gene_use_long$rs12611523
gene_use_long$rs1263056<-2-gene_use_long$rs1263056
gene_use_long$rs13109404<-2-gene_use_long$rs13109404
gene_use_long$rs17427571<-2-gene_use_long$rs17427571
gene_use_long$rs17732997<-2-gene_use_long$rs17732997
gene_use_long$rs1776776<-2-gene_use_long$rs1776776
gene_use_long$rs180769<-2-gene_use_long$rs180769
gene_use_long$rs1939455<-2-gene_use_long$rs1939455
gene_use_long$rs1991556<-2-gene_use_long$rs1991556
gene_use_long$rs2072727<-2-gene_use_long$rs2072727
gene_use_long$rs2079070<-2-gene_use_long$rs2079070
gene_use_long$rs2192528<-2-gene_use_long$rs2192528
gene_use_long$rs3095508<-2-gene_use_long$rs3095508
gene_use_long$rs34354917<-2-gene_use_long$rs34354917
#gene_use_long$rs34556183<-2-gene_use_long$rs34556183
gene_use_long$rs365663<-2-gene_use_long$rs365663
gene_use_long$rs374153<-2-gene_use_long$rs374153
gene_use_long$rs460692<-2-gene_use_long$rs460692
gene_use_long$rs55658675<-2-gene_use_long$rs55658675
gene_use_long$rs62120041<-2-gene_use_long$rs62120041
gene_use_long$rs6575005<-2-gene_use_long$rs6575005
gene_use_long$rs73219758<-2-gene_use_long$rs73219758
gene_use_long$rs7503199<-2-gene_use_long$rs7503199
gene_use_long$rs7616632<-2-gene_use_long$rs7616632
gene_use_long$rs7644809<-2-gene_use_long$rs7644809
gene_use_long$rs7806045<-2-gene_use_long$rs7806045
gene_use_long$rs7915425<-2-gene_use_long$rs7915425
gene_use_long$rs8038326<-2-gene_use_long$rs8038326
gene_use_long$rs8050478<-2-gene_use_long$rs8050478
gene_use_long$rs915416<-2-gene_use_long$rs915416
gene_use_long$rs9382445<-2-gene_use_long$rs9382445
gene_use_long$rs9903973<-2-gene_use_long$rs9903973
gene_use_long$rs9940646<-2-gene_use_long$rs9940646

#add rs34556183
snp1$gwassentrixps<-paste(snp1$V1,"_",snp1$V2,sep="")
snp1$V9<-paste(snp1$V7,snp1$V8, sep="")
snp1$rs34556183[snp1$V9=="AA"]<-2
snp1$rs34556183[snp1$V9=="AG"|snp1$V9=="GA"]<-1
snp1$rs34556183[snp1$V9=="GG"]<-0
snp1<-snp1[c("gwassentrixps","rs34556183")]
gene_use_long2<-merge(snp1,gene_use_long,by="gwassentrixps")

write.csv(gene_use_long2,file=paste(Sys.getenv('mydata'),'bib_genotype78_v2.csv',sep=''),row.names=FALSE)
