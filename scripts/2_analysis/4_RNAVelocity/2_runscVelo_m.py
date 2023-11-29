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

    # for merged samples
    cell_metadata = pd.read_csv(os.path.join(meta_path, "metadata.csv"),
                                header=0, 
                                names=["Cell ID", "cellType", 'stage', 'phase', 's_score', 'g2m_score'])

    #filter out cells
    cells_filtered = loom_object[np.isin(loom_object.obs.index, barcodes)]

    #get the filtered cell names            
    cell_index = pd.DataFrame(cells_filtered.obs.index)
    cell_index.rename(columns={'CellID' : 'Cell ID'}, inplace=True)

    # reordered metadata
    umap_ordered = cell_index.merge(cell_umap, on = "Cell ID")
    cell_metadata_reordered = cell_index.merge(cell_metadata, on = "Cell ID")

    umap_ordered = umap_ordered.iloc[:,1:]
    cell_type_ordered = cell_metadata_reordered[['cellType']]
    stage_order =  cell_metadata_reordered[['stage']]
    phase_ordered = cell_metadata_reordered[['phase']]
    s_score = cell_metadata_reordered[['s_score']]
    g2m_score = cell_metadata_reordered[['g2m_score']]

    # add metadata to the object
    cells_filtered.obsm['X_umap'] = umap_ordered.values
    cells_filtered.obs['stage'] = stage_order.values[:, 0]
    cells_filtered.obs['celltype'] = cell_type_ordered.values[:, 0]

    cells_filtered.obs['phase'] = phase_ordered.values[:, 0]
    cells_filtered.obs['S_score'] = s_score.values[:, 0]
    cells_filtered.obs['G2M_score'] = g2m_score.values[:, 0]
    adata = cells_filtered
    adata.var_names_make_unique()

    return adata

def main():

    loom_data = scv.read('../Clark2019neuron/Clark2019neuron_merge.loom', cache=True)

    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('x',''))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E11:','E11_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E12:','E12_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E14_rep1:','E14R1_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E14_rep2:','E14R2_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E16:','E16_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E18_rep2:','E18R2_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_E18_rep3:','E18R3_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_P0:','P0_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_P2_rep3:','P2R2_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_P5:','P5_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_P8_rep1:','P8R1_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_P8_rep2:','P8R2_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Clark2019neuron_P14:','P14_'))
    
    meta_path = '../Clark2019neuron/merge/'
    print('-----------------------------------')
    print("Adding metadata to the loom project")
    print('-----------------------------------')
    adata = myAddMetadata(meta_path, loom_data)

    gene_retain_file = pd.read_table('/Users/xzeng/Desktop/Research/ciona/SCENIC/results/R_downstream/mouse/regulon_for_analysis.txt', sep= '\t', header=0)
    gene_retain = gene_retain_file['tf_name'].to_list()
    gene_retain = list(set(gene_retain))

    # preprocessing 
    scv.pp.filter_genes(adata, min_shared_counts=10,retain_genes = gene_retain)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000, retain_genes = gene_retain)
    scv.pp.log1p(adata)

    # Run scVelo
    print('-----------------------------------')
    print("Runing scVelo...")
    print('-----------------------------------')
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    scv.tl.recover_dynamics(adata, n_jobs=12)
    scv.tl.velocity(adata,mode='dynamical')
    scv.tl.velocity_graph(adata, n_jobs=12)

    print('-----------------------------------')
    print('Writing data...')
    print('-----------------------------------')

    adata.write('../Clark2019neuron/Clark2019neuron_merge_processed1.h5ad')

if __name__ == "__main__":
    main()




