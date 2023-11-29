import scvelo as scv
import pandas as pd
import numpy as np

scv.logging.print_version()

def filter_low_unspliced_cells(loom_data):
    
    num_spliced_per_cell = loom_data.layers['spliced'].todense().sum(axis = 1)
    num_unspliced_per_cell = loom_data.layers['unspliced'].todense().sum(axis = 1)

    unspliced_proportion = num_unspliced_per_cell/(num_unspliced_per_cell + num_spliced_per_cell)

    cell_selected = loom_data.obs[unspliced_proportion > 0.07].index.to_list()
    selected_loom_data = loom_data[cell_selected].copy()

    return  selected_loom_data,cell_selected

def main():

    loom_data = scv.read('../Xu2020Development/Xu2020Development_merge.loom', cache=False)
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('x',''))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_24hpf_possorted:','24hpf_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_36hpf_possorted:','36hpf_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_48hpf_possorted:','48hpf1_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_48hpf_rep_possorted:','48hpf2_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_72hpf_possorted:','72hpf1_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_72hpf_rep_possorted:','72hpf2_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_14dpf_possorted:','14dpf1_'))
    loom_data.obs = loom_data.obs.rename(index = lambda x: x.replace('Xu2020Development_14dpf_rep_possorted:','14dpf2_'))
    loom_data.obs_names.name = 'Cell ID'
    loom_data,selected_cell_names = filter_low_unspliced_cells(loom_data)


    loom_data.write('../results/zebrafish/merge/merge_flitered.h5ad')

    f = open('../results/zebrafish/merge/filtered_cells.txt', 'w')
    for x in selected_cell_names:
        f.write(x + '\n')
    f.close()
    print('Job filished!')

if __name__ == "__main__":
    main()
