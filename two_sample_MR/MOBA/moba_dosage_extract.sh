#data on the three batches are extracted separately using the code below

source /cluster/bin/jobsetup
module purge
set -o errexit

##Copy input files to the work directory:
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do 
cp /path/to/the/data $SCRATCH
done
cp /path/to/the/SNPs/list $SCRATCH

##Mark outfiles for automatic copying to $SLURM_SUBMIT_DIR:
for m in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
cleanup "cp $SCRATCH/harvest12m_$m.dosage /path/to/the/results"
done

##Then type your code below
cd $SCRATCH
for n in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
qctool_v2.0.6 -g $n.vcf.gz -incl-rsids snp78.txt -filetype vcf -og harvest12m_$n.dosage -ofiletype dosage
done
