module add apps/qctool-2.0 
### IMPUTED BGEN DATA
bgen="/path_to_the_data/data_#.bgen"

### NAME OUTPUT FILE (INDIVIDUAL DATA FOR DOSAGES FOR EACH VARIANT)
output="/path_to_the_data/alspac78.dosage"

### SAMPLE ID FILE
scp newblue1:/path_to_the_data/data.sample .
sample=data.sample

### LIST OF RSIDs
rsid="/path_to_the_data/UKB_78SNPs_forALSPAC.txt"

### NAME LOG FILE
log="/path_to_the_data/log.txt"

### SNP STATS
stats="/path_to_the_data/stats.txt"

qctool -g ${bgen} \
		-og ${output} \
		-s ${sample} \
		-incl-rsids ${rsid} \
		-log ${log} 
