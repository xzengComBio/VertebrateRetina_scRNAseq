#!/bin/bash

# Define MASKFILE as a GTF file located in a specific path, likely a reference file for masking
MASKFILE="/home/xzeng/data/ref_data/Danio_rerio.GRCz11.107/GRCz11_rmsk.gtf"
GTFFILE="/home/xzeng/data/ref_data/Danio_rerio.GRCz11.107/DanioGRCz11_cellranger/genes/genes_forVelocyto.gtf"

# Read the text file containing sample names, and then iterate over each sample using a while loop
cat /home/xzeng/project/Ciona/raw_data/Xu2020Development/sample.txt | while read sample
do
	SAMPLEPATH="/home/xzeng/project/Ciona/raw_data/Xu2020Development/Xu2020Development_${sample}"
	JOB_NAME="zebrafish_${sample}"
	
	# Pass three arguments to velocyto.sh: the sample's path, the gene GTF file, and the MASK file
	qsub -N ${JOB_NAME} velocyto.sh ${SAMPLEPATH} ${GTFFILE} ${MASKFILE}

done
