#for miscarriage
# Clear the work environment
rm(list = ls())


#package
require(data.table)


#environment variable
Sys.setenv(Mydata="/path/to/my/data")


#input data
miscarriage<-fread(paste(Sys.getenv('Mydata'),'summary_stats_finngen_R5_O15_ABORT_SPONTAN',sep=''), sep='\t')


#extract 91 SNPs
miscarriage_use<-miscarriage[ which(miscarriage$rsids%in%c("rs7556815","rs75539574","rs12607679","rs915416","rs9940646","rs13109404","rs8050478","rs56372231","rs13088093","rs2079070","rs34556183",
                                                           "rs3095508","rs34731055","rs73219758","rs10973207","rs2139261","rs4592416","rs365663","rs1517572","rs7915425","rs330088",
                                                           "rs8038326","rs460692","rs9382445","rs4767550","rs11885663","rs1991556","rs1057703","rs4128364","rs61796569","rs10483350",
                                                           "rs7115226","rs269054","rs112230981","rs11602180","rs2192528","rs12246842","rs205024","rs12567114","rs7616632","rs6575005",
                                                           "rs1776776","rs11621908","rs10421649","rs2072727","rs113113059","rs374153","rs151014368","rs62120041","rs7503199","rs1939455",
                                                           "rs7951019","rs17732997","rs61985058","rs17427571","rs7806045","rs35531607","rs7644809","rs9345234","rs12791153","rs1263056",
                                                           "rs55658675","rs11567976","rs180769","rs1553132","rs9903973","rs11614986","rs2231265","rs174560","rs10173260","rs72804080",
                                                           "rs12611523","rs11643715","rs4538155","rs34354917","rs80193650","rs10761674","rs11190970","rs4988235","rs11957190","rs12569901",
                                                           "rs1073160","rs8074498","rs4780834","rs3788337","rs308604","rs11155606","rs11135570","rs10273733","rs2279681","rs7683893")), ]



write.csv(miscarriage_use,"/path/to/save/gy/associations",row.names=FALSE)
#####################################################################################################
#for gestational diabetes
# Clear the work environment
rm(list = ls())


#package
require(data.table)


#environment variable
Sys.setenv(Mydata="/path/to/my/data")


#input data
gdm<-fread(paste(Sys.getenv('Mydata'),'summary_stats_finngen_R5_GEST_DIABETES',sep=''), sep='\t')


#extract 91 SNPs
gdm_use<-gdm[ which(gdm$rsids%in%c("rs7556815","rs75539574","rs12607679","rs915416","rs9940646","rs13109404","rs8050478","rs56372231","rs13088093","rs2079070","rs34556183",
                                   "rs3095508","rs34731055","rs73219758","rs10973207","rs2139261","rs4592416","rs365663","rs1517572","rs7915425","rs330088",
                                   "rs8038326","rs460692","rs9382445","rs4767550","rs11885663","rs1991556","rs1057703","rs4128364","rs61796569","rs10483350",
                                   "rs7115226","rs269054","rs112230981","rs11602180","rs2192528","rs12246842","rs205024","rs12567114","rs7616632","rs6575005",
                                   "rs1776776","rs11621908","rs10421649","rs2072727","rs113113059","rs374153","rs151014368","rs62120041","rs7503199","rs1939455",
                                   "rs7951019","rs17732997","rs61985058","rs17427571","rs7806045","rs35531607","rs7644809","rs9345234","rs12791153","rs1263056",
                                   "rs55658675","rs11567976","rs180769","rs1553132","rs9903973","rs11614986","rs2231265","rs174560","rs10173260","rs72804080",
                                   "rs12611523","rs11643715","rs4538155","rs34354917","rs80193650","rs10761674","rs11190970","rs4988235","rs11957190","rs12569901",
                                   "rs1073160","rs8074498","rs4780834","rs3788337","rs308604","rs11155606","rs11135570","rs10273733","rs2279681","rs7683893")), ]



write.csv(gdm_use,"/path/to/save/gy/associations",row.names=FALSE)
#####################################################################################################
#for hypertensive disorders of pregnancy
# Clear the work environment
rm(list = ls())


#package
require(data.table)


#environment variable
Sys.setenv(Mydata="/path/to/my/data")


#input data
hdp<-fread(paste(Sys.getenv('Mydata'),'summary_stats_finngen_R5_O15_GESTAT_HYPERT',sep=''), sep='\t')


#extract 91 SNPs
hdp_use<-hdp[ which(hdp$rsids%in%c("rs7556815","rs75539574","rs12607679","rs915416","rs9940646","rs13109404","rs8050478","rs56372231","rs13088093","rs2079070","rs34556183",
                                   "rs3095508","rs34731055","rs73219758","rs10973207","rs2139261","rs4592416","rs365663","rs1517572","rs7915425","rs330088",
                                   "rs8038326","rs460692","rs9382445","rs4767550","rs11885663","rs1991556","rs1057703","rs4128364","rs61796569","rs10483350",
                                   "rs7115226","rs269054","rs112230981","rs11602180","rs2192528","rs12246842","rs205024","rs12567114","rs7616632","rs6575005",
                                   "rs1776776","rs11621908","rs10421649","rs2072727","rs113113059","rs374153","rs151014368","rs62120041","rs7503199","rs1939455",
                                   "rs7951019","rs17732997","rs61985058","rs17427571","rs7806045","rs35531607","rs7644809","rs9345234","rs12791153","rs1263056",
                                   "rs55658675","rs11567976","rs180769","rs1553132","rs9903973","rs11614986","rs2231265","rs174560","rs10173260","rs72804080",
                                   "rs12611523","rs11643715","rs4538155","rs34354917","rs80193650","rs10761674","rs11190970","rs4988235","rs11957190","rs12569901",
                                   "rs1073160","rs8074498","rs4780834","rs3788337","rs308604","rs11155606","rs11135570","rs10273733","rs2279681","rs7683893")), ]



write.csv(hdp_use,"/path/to/save/gy/associations",row.names=FALSE)
#####################################################################################################
#for preterm birth
# Clear the work environment
rm(list = ls())


#package
require(data.table)


#environment variable
Sys.setenv(Mydata="/path/to/my/data")


#input data
pt<-fread(paste(Sys.getenv('Mydata'),'summary_stats_finngen_R5_O15_PRETERM',sep=''), sep='\t')


#extract 81 SNPs
pt_use<-pt[ which(pt$rsids%in%c("rs7556815","rs75539574","rs12607679","rs915416","rs9940646","rs13109404","rs8050478","rs56372231","rs13088093","rs2079070","rs34556183",
                                "rs3095508","rs34731055","rs73219758","rs10973207","rs2139261","rs4592416","rs365663","rs1517572","rs7915425","rs330088",
                                "rs8038326","rs460692","rs9382445","rs4767550","rs11885663","rs1991556","rs1057703","rs4128364","rs61796569","rs10483350",
                                "rs7115226","rs269054","rs112230981","rs11602180","rs2192528","rs12246842","rs205024","rs12567114","rs7616632","rs6575005",
                                "rs1776776","rs11621908","rs10421649","rs2072727","rs113113059","rs374153","rs151014368","rs62120041","rs7503199","rs1939455",
                                "rs7951019","rs17732997","rs61985058","rs17427571","rs7806045","rs35531607","rs7644809","rs9345234","rs12791153","rs1263056",
                                "rs55658675","rs11567976","rs180769","rs1553132","rs9903973","rs11614986","rs2231265","rs174560","rs10173260","rs72804080",
                                "rs12611523","rs11643715","rs4538155","rs34354917","rs80193650","rs10761674","rs11190970","rs4988235","rs11957190","rs12569901",
                                "rs1073160","rs8074498","rs4780834","rs3788337","rs308604","rs11155606","rs11135570","rs10273733","rs2279681","rs7683893")), ]



write.csv(pt_use,"/path/to/save/gy/associations",row.names=FALSE)


