#!bin/sh
#$ -S /bin/sh
#$ -l s_vmem=200G
#$ -l mem_req=200G
#$ -pe def_slot 1
#$ -o /home/xzeng/project/Ciona/log/Analysis
#$ -e /home/xzeng/project/Ciona/log/Analysis
#$ -N Clark2019neuron

time=$(date "+%Y-%m-%d-%H:%M:%S")
echo "This Job was Ran on ${time}"
echo "Integrating samples in $1"

module load R/4.1.0

OUTPUT_PATH=${1%/*}
BASENAME=$(basename $1 .rds)
OUTPUT_FILE=${OUTPUT_PATH}/${BASENAME}_integrated.rds
echo "The output file name will be ${OUTPUT_FILE}"

Rscript /home/xzeng/project/Ciona/codes/Analysis/seurat_integration.r $1 ${OUTPUT_FILE}

echo "task completed!"

#qsub -q mjobs.q,ljobs.q rscript_integration.sh /home/xzeng/project/Ciona/report/rds/Clark2019neuron_raw.rds
