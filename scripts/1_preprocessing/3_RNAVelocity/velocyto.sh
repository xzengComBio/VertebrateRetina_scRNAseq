#!bin/sh
#$ -S /bin/sh
#$ -l s_vmem=10G
#$ -l mem_req=10G
#$ -pe def_slot 4
#$ -o /home/xzeng/project/Ciona/log/velocyto
#$ -e /home/xzeng/project/Ciona/log/velocyto
#$ -cwd

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate velocyto

echo "*****************************************"
echo "Running velocyto...."
echo "	"
echo "SAMPLE PATH: ${1}"
echo "GTF FILE PATH: ${2}"
echo "MASK FILE PATH ${3}"
echo "  "
echo "*****************************************"

velocyto run10x	 $1  \
$2 \
-m $3 

echo "Job finish!"

