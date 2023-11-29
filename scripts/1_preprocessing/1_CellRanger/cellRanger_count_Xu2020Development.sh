#!/bin/bash

cat /home/xzeng/project/Ciona/raw_data/Xu2020Development/sample.txt | while read sample
do
	
	#Internal script name
	SHORT="cellranger_${sample}"
		
	#Construct shell file
	echo "Creating script cellranger_${sample}"
	cat > .${SHORT}.sh <<EOF
#!bin/sh
#$ -l s_vmem=8G,mem_req=8G
#$ -pe def_slot 8
#$ -o /home/xzeng/project/Ciona/log/cellranger/Xu2020Development/
#$ -e /home/xzeng/project/Ciona/log/cellranger/Xu2020Development/
#$ -cwd
#$ -N ${SHORT}

echo ""
echo "  ******************** Job starts ********************"
echo "    $(date)"
echo ""
echo "  ****************** Shirokane info ******************"
echo "  **  User: Xin Zeng                                **"
echo "  **  Project name: Ciona                           **"
echo "  **  Job name: ${SHORT}                            **"
echo "  ****************************************************"
echo ""

## List current modules for reproducibility
module list

export PATH=/home/xzeng/tools/cellranger-6.1.2:$PATH

echo "  running cellranger count for ${sample}"

time cellranger count --id=Xu2020Development_${sample} \
		      --fastqs=/home/xzeng/project/Ciona/raw_data/Xu2020Development/merge_fastaq \
		      --include-introns \
 		      --sample=${sample} \
		      --transcriptome=/home/xzeng/data/ref_data/Danio_rerio.GRCz11.107/DanioGRCz11_cellranger \
		      --localcores=8 \
		      --localmem=64 \
		      --nosecondary \

date
echo "  ******************** Job ends ********************"
EOF

    call="qsub -q mjobs.q,ljobs.q .${SHORT}.sh"
    echo $call
    $call
done

