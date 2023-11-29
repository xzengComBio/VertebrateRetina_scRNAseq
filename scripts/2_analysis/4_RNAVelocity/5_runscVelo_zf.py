import warnings
import scvelo as scv

import pandas as pd
import numpy as np
import os
import sys

scv.logging.print_version()
warnings.filterwarnings('ignore')

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  

def myAddMetadata(meta_path, loom_object):

	'''
	add metadata from seurat object to loom object
	
	input:
	meta_path: the file path of metadata
	loom_object: the loom object that was generated from velocity.

	output:
	anndata object that contain metadata information     
	'''

	#read metadata from seurat object
	barcodes = pd.read_csv(os.path.join(meta_path, 'cellID_obs.csv'))
	cell_umap = pd.read_csv(os.path.join(meta_path, "cell_embeddings.csv"),
							header=0, 
							names=["Cell ID", "UMAP_1", "UMAP_2"])
	cell_metadata = pd.read_csv(os.path.join(meta_path, "metadata.csv"),
								header=0, 
								#names=["Cell ID", "cluster", 'cellType', 'phase', 's_score', 'g2m_score'])
								names=["Cell ID", "cluster", 'cellType', 'module', 'phase', 's_score', 'g2m_score']) #48hpf


	#filter out cells
	cells_filtered = loom_object[np.isin(loom_object.obs.index, barcodes)]

	#get the filtered cell names            
	cell_index = pd.DataFrame(cells_filtered.obs.index)
	cell_index = cell_index.rename(columns = {0:'Cell ID'})

	# reordered metadata
	umap_ordered = cell_index.merge(cell_umap, on = "Cell ID")
	cell_metadata_reordered = cell_index.merge(cell_metadata, on = "Cell ID")

	umap_ordered = umap_ordered.iloc[:,1:]
	clusters_ordered = cell_metadata_reordered[['cluster']]
	cell_type_ordered = cell_metadata_reordered[['cellType']]
   # new_cell_type_ordered = cell_metadata_reordered[['module']]
	phase_ordered = cell_metadata_reordered[['phase']]
	s_score = cell_metadata_reordered[['s_score']]
	g2m_score = cell_metadata_reordered[['g2m_score']]

	# add metadata to the object
	cells_filtered.obsm['X_umap'] = umap_ordered.values
	cells_filtered.uns['clusters'] = clusters_ordered.values
	cells_filtered.obs['celltype'] = cell_type_ordered.values
	#cells_filtered.obs['module'] = new_cell_type_ordered.values
	cells_filtered.obs['phase'] = phase_ordered.values
	cells_filtered.obs['S_score'] = s_score.values
	cells_filtered.obs['G2M_score'] = g2m_score.values
	adata = cells_filtered
	adata.var_names_make_unique()

	return adata

def filter_low_unspliced_cells(loom_data):
	
	num_spliced_per_cell = loom_data.layers['spliced'].todense().sum(axis = 1)
	num_unspliced_per_cell = loom_data.layers['unspliced'].todense().sum(axis = 1)

	unspliced_proportion = num_unspliced_per_cell/(num_unspliced_per_cell + num_spliced_per_cell)

	cell_selected = loom_data.obs[unspliced_proportion > 0.07].index.to_list()
	selected_loom_data = loom_data[cell_selected].copy()

	return  selected_loom_data


def main():

	loom_data = scv.read('../results/zebrafish/merge/merge_flitered.h5ad', cache=False)

	print(loom_data.obs_names)
	meta_path = '../Xu2020Development/merge/'
	adata = myAddMetadata(meta_path, loom_data)
	print(adata)

	gene_retain_file = pd.read_table('/Users/xzeng/Desktop/Research/ciona/SCENIC/results/R_downstream/zebrafish/regulon_for_analysis.txt', sep= '\t', header=0)
	gene_retain = gene_retain_file['tf_name'].to_list()
	gene_retain = list(set(gene_retain))
	print(gene_retain)

	# preprocessing	
	scv.pp.filter_genes(adata, min_shared_cells = 10,retain_genes = gene_retain)
	scv.pp.normalize_per_cell(adata)
	scv.pp.filter_genes_dispersion(adata, n_top_genes=2000, retain_genes = gene_retain)
	scv.pp.log1p(adata)
	
	scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

	scv.tl.recover_dynamics(adata, n_jobs=12)
	scv.tl.velocity(adata,mode='dynamical')
	scv.tl.velocity_graph(adata, n_jobs=12)

	print('WRITING ADATA...')
	adata.write('../results/zebrafish/merge/merge.h5ad')

	print('job finish!')

if __name__ == "__main__":
	main()









