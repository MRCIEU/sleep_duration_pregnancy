# Clear the work environment
rm(list = ls())

#package
require(data.table)

#input genetic data
gene<-fread("/path/to/the/data")
information<-gene[,list(chromosome,SNPID,rsid,position,alleleA,alleleB)]
gene_use<-gene[, c("chromosome","SNPID","position","alleleA","alleleB"):=NULL]
gene_use$rsid<-c("rs915416","rs269054","rs61796569","rs12567114","rs62120041","rs374153","rs75539574","rs72804080","rs7556815","rs12611523","rs4128364",
                 "rs4538155","rs11885663","rs10173260","rs112230981","rs17732997","rs7644809","rs13088093","rs7616632","rs2192528","rs17427571","rs35531607",
                 "rs13109404","rs365663","rs460692","rs56372231","rs180769","rs11567976","rs151014368","rs34556183","rs80193650","rs113113059","rs9382445",
                 "rs2231265","rs9345234","rs34731055","rs2079070","rs7806045","rs330088","rs73219758","rs10973207","rs1776776","rs12246842","rs10761674",
                 "rs11190970","rs7915425","rs1517572","rs4592416","rs11602180","rs174560","rs12791153","rs1553132","rs1939455","rs7115226","rs1263056",
                 "rs7951019","rs1057703","rs34354917","rs11614986","rs4767550","rs6575005","rs10483350","rs61985058","rs55658675","rs11621908","rs8038326",
                 "rs3095508","rs11643715","rs9940646","rs8050478","rs7503199","rs205024","rs2139261","rs1991556","rs9903973","rs12607679","rs10421649","rs2072727")

#reshape data
gene_use_long<-dcast(melt(as.matrix(gene_use)), Var2~paste0('r', Var1), value.var='value')
colnames(gene_use_long) <- gene_use_long[1,]
gene_use_long <- gene_use_long[-1, ]
colnames(gene_use_long)[1]<-"aln"

#select mum
gene_use_long <- gene_use_long[ which(substring(gene_use_long$aln,6,6)=="M"), ]
mum<-read.table("/path/to/the/data")
colnames(mum)[1]<-"aln"
mum<-mum["aln"]
gene_use_correct <- merge(gene_use_long,mum,by="aln")
gene_use_correct$aln<-substr(gene_use_correct$aln,1,5)
i<-2
for (i in 2:79){
  gene_use_correct[,i]<-as.numeric(gene_use_correct[,i])
}

#EA, NEA (same as UKB)
gene_use_correct$rs10761674<-2-gene_use_correct$rs10761674
gene_use_correct$rs11190970<-2-gene_use_correct$rs11190970
gene_use_correct$rs112230981<-2-gene_use_correct$rs112230981
gene_use_correct$rs113113059<-2-gene_use_correct$rs113113059
gene_use_correct$rs11602180<-2-gene_use_correct$rs11602180
gene_use_correct$rs11614986<-2-gene_use_correct$rs11614986
gene_use_correct$rs11621908<-2-gene_use_correct$rs11621908
gene_use_correct$rs12246842<-2-gene_use_correct$rs12246842
gene_use_correct$rs12607679<-2-gene_use_correct$rs12607679
gene_use_correct$rs12611523<-2-gene_use_correct$rs12611523
gene_use_correct$rs1263056<-2-gene_use_correct$rs1263056
gene_use_correct$rs13109404<-2-gene_use_correct$rs13109404
gene_use_correct$rs17427571<-2-gene_use_correct$rs17427571
gene_use_correct$rs17732997<-2-gene_use_correct$rs17732997
gene_use_correct$rs1776776<-2-gene_use_correct$rs1776776
gene_use_correct$rs180769<-2-gene_use_correct$rs180769
gene_use_correct$rs1939455<-2-gene_use_correct$rs1939455
gene_use_correct$rs1991556<-2-gene_use_correct$rs1991556
gene_use_correct$rs2072727<-2-gene_use_correct$rs2072727
gene_use_correct$rs2079070<-2-gene_use_correct$rs2079070
gene_use_correct$rs2192528<-2-gene_use_correct$rs2192528
gene_use_correct$rs3095508<-2-gene_use_correct$rs3095508
gene_use_correct$rs34354917<-2-gene_use_correct$rs34354917
gene_use_correct$rs34556183<-2-gene_use_correct$rs34556183
gene_use_correct$rs365663<-2-gene_use_correct$rs365663
gene_use_correct$rs374153<-2-gene_use_correct$rs374153
gene_use_correct$rs460692<-2-gene_use_correct$rs460692
gene_use_correct$rs55658675<-2-gene_use_correct$rs55658675
gene_use_correct$rs62120041<-2-gene_use_correct$rs62120041
gene_use_correct$rs6575005<-2-gene_use_correct$rs6575005
gene_use_correct$rs73219758<-2-gene_use_correct$rs73219758
gene_use_correct$rs7503199<-2-gene_use_correct$rs7503199
gene_use_correct$rs7616632<-2-gene_use_correct$rs7616632
gene_use_correct$rs7644809<-2-gene_use_correct$rs7644809
gene_use_correct$rs7806045<-2-gene_use_correct$rs7806045
gene_use_correct$rs7915425<-2-gene_use_correct$rs7915425
gene_use_correct$rs8038326<-2-gene_use_correct$rs8038326
gene_use_correct$rs8050478<-2-gene_use_correct$rs8050478
gene_use_correct$rs915416<-2-gene_use_correct$rs915416
gene_use_correct$rs9382445<-2-gene_use_correct$rs9382445
gene_use_correct$rs9903973<-2-gene_use_correct$rs9903973
gene_use_correct$rs9940646<-2-gene_use_correct$rs9940646

write.csv(gene_use_correct,"/path/to/the/data",row.names=FALSE)
####################################################################################
#to clean child's genotype data, scripts are the same except for selecting ALSPAC children
#select child
gene_use_long <- gene_use_long[ which(substring(gene_use_long$aln,6,6)!="M"), ]
child<-read.table("/path/to/the/data")
colnames(child)[1]<-"aln"
child<-child["aln"]
gene_use_correct <- merge(gene_use_long,child,by="aln")
