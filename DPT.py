#Diffusion Pseudotime

#Imports
import scanpy as sc
import os
import numpy as np
import pandas as pd

wdpath = os.getcwd()
def DPT(adata, writer):
    group_names = adata.cluster
    interested_genes = adata.markerGenes
    
    book = writer.book        
    worksheet = book.add_worksheet("DPT")  #Create PAGA page in report
    
    paths = []
    print('View paga and umap_paga .png files in figures for path analysis.')
    i=1
    while True:
        #Enter i path name
        if i==1:
            name_input = input('Enter 1st Path Name: ')
        else:
            name_input = input(f'Enter Path {i} Name (or press enter to end): ')
            
        if not name_input:
            break
        
        j=1
        path = []
        while True:
            path_input = input(f'Enter step {j} in path {name_input} (or press enter to end path):' )
            if not path_input:
                break
            elif path_input not in set(list(adata.obs[group_names[0]])):
                print(f"Invalid Input \n Valid Path Names: {set(list(adata.obs[group_names[0]]))}")
            else:
                path.append(path_input)
                j=j+1 
        
        paths.append((name_input,path))
        i=i+1
        
    for i in range(len(paths)):
        root_cell = paths[i][1][0]
        adata.uns['iroot'] = np.flatnonzero(adata.obs[group_names[0]]  == root_cell)[0]
        sc.tl.dpt(adata)

        sc.tl.draw_graph(adata, init_pos='umap',layout='fr')
        sc.pl.draw_graph(adata, color=[group_names[0], 'dpt_pseudotime'],
                        legend_loc='on data',title=[f'Path {paths[i][0]}', 'dpt_pseudotime'],
                        show=False, save=f'_{paths[i][0]}.png')

        adata.obs['distance'] = adata.obs['dpt_pseudotime']
        
        tup = paths[i]
        lst = [tup[0]] + tup[1]
        df = pd.DataFrame(lst).T
        df.to_excel(writer, sheet_name='DPT',startrow=25*i, startcol=0)
        
        sc.pl.paga_path(adata, nodes=paths[i][1], keys=interested_genes,
                        title=f'Path {paths[i][0]}', show=False, save=f'_{paths[i][0]}.png')
        worksheet.insert_image(f'A{25*i+5}', wdpath+f'\\figures\\draw_graph_fr_{paths[i][0]}.png')
        worksheet.insert_image(f'T{25*i+5}', wdpath+f'\\figures\\paga_path_{paths[i][0]}.png')
        
    return(adata, writer)