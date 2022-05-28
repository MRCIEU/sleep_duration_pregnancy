module add apps/qctool-2.0 

############################################################################################		

# 1) Set your working directory
export WORK_DIR=/local/$PBS_JOBID
mkdir -p /local/$PBS_JOBID
cd $WORK_DIR

# 2) Define your input and output files
### IMPUTED BGEN DATA
vcffile="/path/to/the/data"
# Note: If the input filename contains a # character, e.g. example_#.gen this is treated as a chromosomal wildcard and will match all (human) chromosomes

### NAME OUTPUT FILE (INDIVIDUAL DATA FOR DOSAGES FOR EACH VARIANT)
output="/path/to/the/data"

### LIST OF RSIDs
rsid="/path/to/the/data"

### NAME LOG FILE
log="/path/to/the/data"

### SNP STATS
stats="/path/to/the/data"

############################################################################################		

# 3) Run QCTOOL to extract dosages
qctool -g ${vcffile} \
       -filetype vcf \
	   -vcf-genotype-field GP \
	   -incl-rsids ${rsid} \
		-og ${output} \
		-ofiletype dosage \
		-log ${log} 
    
    
