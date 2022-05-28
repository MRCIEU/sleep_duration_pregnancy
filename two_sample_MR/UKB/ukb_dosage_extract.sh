module add apps/qctool/2.0rc4 

############################################################################################		

# 1) Set your working directory

# 2) Define your input and output files
### IMPUTED BGEN DATA
bgen="/path/to/the/data"
# Note: If the input filename contains a # character, e.g. example_#.gen this is treated as a chromosomal wildcard and will match all (human) chromosomes

### NAME OUTPUT FILE (INDIVIDUAL DATA FOR DOSAGES FOR EACH VARIANT)
output="/path/to/the/data"

### SAMPLE ID FILE
sample="/path/to/the/data"

### LIST OF RSIDs
rsid="/path/to/the/data"

### NAME LOG FILE
log="/path/to/save/the/log"

############################################################################################		

# 3) Run QCTOOL to extract dosages
qctool -g ${bgen} \
		-og ${output} \
		-s ${sample} \
		-incl-rsids ${rsid} \
		-log ${log} 
