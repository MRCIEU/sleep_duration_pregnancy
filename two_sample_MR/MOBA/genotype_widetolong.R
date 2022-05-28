# Clear the work environment
rm(list = ls())

#package
require(data.table)
require(reshape2)

#input data
Sys.setenv(duration17k_harvest12m="/path/to/the/data")
Sys.setenv(duration17k_harvest24m="/path/to/the/data")
Sys.setenv(duration17k_rotterdam1="/path/to/the/data")
i<-1
for(i in 1:22){
  harvest12m = paste("harvest12m_", i, sep="")
  harvest24m = paste("harvest24m_", i, sep="")
  rotterdam1 = paste("rotterdam1_", i, sep="")
  assign(harvest12m, fread(paste(Sys.getenv('duration17k_harvest12m'),harvest12m,".dosage", sep="")))
  assign(harvest24m, fread(paste(Sys.getenv('duration17k_harvest24m'),harvest24m,".dosage", sep="")))
  assign(rotterdam1, fread(paste(Sys.getenv('duration17k_rotterdam1'),rotterdam1,".dosage", sep="")))}
harvest12_wide<-rbindlist(list(harvest12m_1,harvest12m_2,harvest12m_3,harvest12m_4,harvest12m_5,harvest12m_6,harvest12m_7,
                               harvest12m_8,harvest12m_9,harvest12m_10,harvest12m_11,harvest12m_12,harvest12m_13,harvest12m_14,
                               harvest12m_15,harvest12m_16,harvest12m_17,harvest12m_18,harvest12m_19,harvest12m_20,harvest12m_21,harvest12m_22))
harvest24_wide<-rbindlist(list(harvest24m_1,harvest24m_2,harvest24m_3,harvest24m_4,harvest24m_5,harvest24m_6,harvest24m_7,
                               harvest24m_8,harvest24m_9,harvest24m_10,harvest24m_11,harvest24m_12,harvest24m_13,harvest24m_14,
                               harvest24m_15,harvest24m_16,harvest24m_17,harvest24m_18,harvest24m_19,harvest24m_20,harvest24m_21,harvest24m_22))
rotterdam1_wide<-rbindlist(list(rotterdam1_1,rotterdam1_2,rotterdam1_3,rotterdam1_4,rotterdam1_5,rotterdam1_6,rotterdam1_7,
                                rotterdam1_8,rotterdam1_9,rotterdam1_10,rotterdam1_11,rotterdam1_12,rotterdam1_13,rotterdam1_14,
                                rotterdam1_15,rotterdam1_16,rotterdam1_17,rotterdam1_18,rotterdam1_19,rotterdam1_20,rotterdam1_21,rotterdam1_22))
moba17k_wide<-merge(harvest12_wide,harvest24_wide,by=c("chromosome","SNPID","rsid","position","alleleA","alleleB"))
moba17k_wide<-merge(moba17k_wide,rotterdam1_wide,by=c("chromosome","SNPID","rsid","position","alleleA","alleleB"))

#SNP information
information<-moba17k_wide[,list(chromosome,SNPID,rsid,position,alleleA,alleleB)]

#reshape data
moba17k_wide<-moba17k_wide[, c("chromosome","SNPID","position","alleleA","alleleB"):=NULL]
moba17k_long<-dcast(melt(as.matrix(moba17k_wide)), Var2~paste0('r', Var1), value.var='value')
colnames(moba17k_long) <- moba17k_long[1,]
moba17k_long <- moba17k_long[-1, ]
colnames(moba17k_long)[1]<-"SENTRIXID"

#EA, NEA (same as UKB)
i<-2
for (i in 2:77){
  moba17k_long[,i]<-as.numeric(moba17k_long[,i])
}
moba17k_long$rs10761674<-2-moba17k_long$rs10761674
moba17k_long$rs11190970<-2-moba17k_long$rs11190970
moba17k_long$rs112230981<-2-moba17k_long$rs112230981
moba17k_long$rs113113059<-2-moba17k_long$rs113113059
moba17k_long$rs11602180<-2-moba17k_long$rs11602180
moba17k_long$rs11614986<-2-moba17k_long$rs11614986
moba17k_long$rs11621908<-2-moba17k_long$rs11621908
moba17k_long$rs12246842<-2-moba17k_long$rs12246842
moba17k_long$rs12607679<-2-moba17k_long$rs12607679
moba17k_long$rs12611523<-2-moba17k_long$rs12611523
moba17k_long$rs1263056<-2-moba17k_long$rs1263056
moba17k_long$rs13109404<-2-moba17k_long$rs13109404
moba17k_long$rs17427571<-2-moba17k_long$rs17427571
moba17k_long$rs17732997<-2-moba17k_long$rs17732997
moba17k_long$rs1776776<-2-moba17k_long$rs1776776
moba17k_long$rs180769<-2-moba17k_long$rs180769
moba17k_long$rs1939455<-2-moba17k_long$rs1939455
moba17k_long$rs1991556<-2-moba17k_long$rs1991556
moba17k_long$rs2072727<-2-moba17k_long$rs2072727
moba17k_long$rs2079070<-2-moba17k_long$rs2079070
moba17k_long$rs2192528<-2-moba17k_long$rs2192528
moba17k_long$rs3095508<-2-moba17k_long$rs3095508
moba17k_long$rs34354917<-2-moba17k_long$rs34354917
#moba17k_long$rs34556183<-2-moba17k_long$rs34556183
moba17k_long$rs365663<-2-moba17k_long$rs365663
moba17k_long$rs374153<-2-moba17k_long$rs374153
moba17k_long$rs460692<-2-moba17k_long$rs460692
moba17k_long$rs55658675<-2-moba17k_long$rs55658675
moba17k_long$rs62120041<-2-moba17k_long$rs62120041
moba17k_long$rs6575005<-2-moba17k_long$rs6575005
moba17k_long$rs73219758<-2-moba17k_long$rs73219758
moba17k_long$rs7503199<-2-moba17k_long$rs7503199
moba17k_long$rs7616632<-2-moba17k_long$rs7616632
moba17k_long$rs7644809<-2-moba17k_long$rs7644809
moba17k_long$rs7806045<-2-moba17k_long$rs7806045
moba17k_long$rs7915425<-2-moba17k_long$rs7915425
moba17k_long$rs8038326<-2-moba17k_long$rs8038326
moba17k_long$rs8050478<-2-moba17k_long$rs8050478
moba17k_long$rs915416<-2-moba17k_long$rs915416
moba17k_long$rs9382445<-2-moba17k_long$rs9382445
moba17k_long$rs9903973<-2-moba17k_long$rs9903973
moba17k_long$rs9940646<-2-moba17k_long$rs9940646

write.csv(moba17k_long,"/path/to/the/data",row.names=FALSE)
