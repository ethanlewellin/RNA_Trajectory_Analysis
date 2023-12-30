#Partition-based graph abstraction

#Imports
import pandas as pd
import scanpy as sc
import os

path=os.getcwd()

def PAGA(adata, writer):
    interested_genes = adata.markerGenes
    group_names = adata.cluster
    sc.tl.paga(adata, groups=group_names[0])

    pos=pd.DataFrame(adata.obsm['X_umap'],index=adata.obs_names)
    pos['group']=adata.obs[adata.uns['paga']['groups']]
    pos=pos.groupby('group').mean()

    book = writer.book        
    worksheet = book.add_worksheet("PAGA")  #Create PAGA page in report
    worksheet.write(0, 0, "Partition-based graph abstraction")

    # Loop over the chunks and create separate plots for each chunk
    i=0
    for chunk in [group_names[0]]+interested_genes:
            sc.pl.paga(adata, color=chunk,
                    pos=pos.values,
                    title=chunk,
                    show=False,
                    save=f'_{chunk}.png')
            sc.pl.umap(adata,
                        color=chunk,
                        title=chunk,
                        show = False,
                        save = f'_paga_{chunk}.png')
            
            worksheet.insert_image(f'A{20*i+8}', path+f'\\figures\\paga_{chunk}.png')
            worksheet.insert_image(f'M{20*i+8}', path+f'\\figures\\umap_paga_{chunk}.png')
            i=i+1
    return(adata, writer)


