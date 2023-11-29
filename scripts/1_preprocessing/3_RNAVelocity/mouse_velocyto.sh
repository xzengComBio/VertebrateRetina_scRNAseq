#!/bin/bash

# Define MASKFILE as a GTF file located in a specific path, likely a reference file for masking
MASKFILE="/home/xzeng/data/ref_data/refdata-gex-mm10-2020-A/mm10_rmsk.gtf"
GTFFILE="/home/xzeng/data/ref_data/refdata-gex-mm10-2020-A/genes/genes.gtf"

# Read the text file containing sample names, and then iterate over each sample using a while loop
cat /home/xzeng/project/Ciona/raw_data/Clark2019neuron/sample.txt | while read sample
do
	SAMPLEPATH="/home/xzeng/project/Ciona/raw_data/Clark2019neuron/Clark2019neuron_${sample}"
	JOB_NAME="mouse_${sample}"
	
	# Pass three arguments to velocyto.sh: the sample's path, the gene GTF file, and the MASK file
	qsub -q mjobs.q,ljobs.q -N ${JOB_NAME} velocyto.sh ${SAMPLEPATH} ${GTFFILE} ${MASKFILE}

done