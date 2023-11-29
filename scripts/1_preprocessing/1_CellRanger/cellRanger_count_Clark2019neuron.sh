#!/bin/bash




cat /home/xzeng/project/Ciona/raw_data/Clark2019neuron/sample.txt | while read sample
do
	
	#Internal script name
	SHORT="cellranger_${sample}"
		
	#Construct shell file
	echo "Creating script cellranger_${sample}"
	cat > .${SHORT}.sh <<EOF
#!bin/sh
#$ -l s_vmem=8G,mem_req=8G
#$ -pe def_slot 8
#$ -o /home/xzeng/project/Ciona/log/cellranger/Clark2019neuron
#$ -e /home/xzeng/project/Ciona/log/cellranger/Clark2019neuron
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

time cellranger count --id=Clark2019neuron_${sample} \
		      --fastqs=/home/xzeng/project/Ciona/raw_data/Clark2019neuron/merge_fastaq \
 		      --sample=${sample} \
 		      --include-introns \
		      --transcriptome=/home/xzeng/data/ref_data/refdata-gex-mm10-2020-A \
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

