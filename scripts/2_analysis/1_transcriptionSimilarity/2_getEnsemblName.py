#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from gtfparse import read_gtf
import argparse
import os

def getOutputFilePath(file_path):

	dir_name = os.path.dirname(file_path)
	output_file = dir_name + '/geneSymbolENSEMBL.txt'

	return output_file


def getgeneSymbolENSEMBL(file_path):

	df = read_gtf(file_path)
	genes = df[df['feature'] == 'gene']
	genes_filtered = genes[['gene_id', 'gene_name']]

	return genes_filtered


def main():

	parser = argparse.ArgumentParser(
        description="Generate (Gene Symbol)-(ENSEMBL ID) pair from gtf file")
	parser.add_argument(
		"-i",
		"--input",
		dest="gtf",
		action='store',
		type=str,
		required=True,
		help="Input GTF file.")

	args = parser.parse_args()
	output_file =  getOutputFilePath(args.gtf)

	gene_filtered = getgeneSymbolENSEMBL(args.gtf)
	gene_filtered.to_csv(output_file, sep = '\t', header=False, index=False)
	print('The output file are saved in ' + output_file)
	print('Job finish!')

if __name__ == "__main__":
    main()