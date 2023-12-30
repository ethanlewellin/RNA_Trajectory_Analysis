#Cluster Analysis

#Imports
import warnings
import os
from click import group
import scanpy as sc
import pandas as pd
warnings.filterwarnings("ignore")
path = os.getcwd()

def ClusterAnalysis(adata, writer):
    
    group_names = None
    while group_names == None:
        group_names = [input('What is the name of your cluster column: ')]
        if group_names[0] not in adata.obs.columns:
            print(f'Column name not in data. \n Valid Column Names: {list(adata.obs.columns)}')
            group_names = None
    
    adata.cluster = group_names
    # Visualize UMAP
    sc.pl.pca(adata,  title='PCA Map', show=False, save='.png')
    sc.pl.umap(adata, color=group_names,legend_loc='on data',
            title='Clusters', add_outline=True, show=False, save='_clusters.png')

    sc.tl.rank_genes_groups(adata, groupby=group_names[0], method='t-test')

    # Extract top N genes for each cluster
    top_n_genes = 10
    marker_genes_df = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(top_n_genes)

    # Get potential cell types for each cluster
    potential_cell_types = {}

    for cluster_label in marker_genes_df.columns:
        top_genes = marker_genes_df[cluster_label]
        
        # Use your own logic to map genes to potential cell types
        # For simplicity, we'll use the gene names as potential cell types
        potential_cell_types[cluster_label] = top_genes.tolist()

    # Print potential cell types for each cluster and
    #Ask user if they want to assign cluster names
    group_name_change=False
    for cluster_label, annotations in potential_cell_types.items():
        print(f'Cluster {cluster_label}: Potential Cell Types: {annotations}')
        user_input = input(f"Assign a name for Cluster {cluster_label} (or press Enter to skip): ")
        
        if user_input:
            group_name_change=True
            adata.obs[group_names[0]][adata.obs.leiden == cluster_label] = user_input
            marker_genes_df.rename(columns={cluster_label:user_input}, inplace=True)
    
    df = marker_genes_df.T
    df.to_excel(writer, sheet_name='Cluster Analysis',startrow=1, startcol=0)
    
    worksheet = writer.sheets['Cluster Analysis']
    worksheet.write(0, 0, "Potential cell types for each cluster ")
    worksheet.insert_image(f'A{marker_genes_df.shape[1]+5}', path+'\\figures\\umap_clusters.png')        

    if group_name_change:
        sc.pl.umap(adata, color=group_names,legend_loc='on data',
            title='Clusters', add_outline=True, show=False, save='_clusters_named.png')
        worksheet.insert_image(f'P{marker_genes_df.shape[1]+5}', path+'\\figures\\umap_clusters_named.png')  
        
    worksheet.insert_image(f'A{marker_genes_df.shape[1]+25}', path+'\\figures\\pca.png')        
    
    # Visualize top-ranked genes for each cluster
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False, save='.png')
    worksheet.insert_image(f'A{marker_genes_df.shape[1]+40}', path+'\\figures\\rank_genes_groups.png')
    
    unique_group_genes = {}
    for cluster_label, genes in potential_cell_types.items():
        # Check if the gene is present only in the current group
        unique_genes = set(genes)
        for other_group, other_genes in potential_cell_types.items():
            if other_group != cluster_label:
                unique_genes -= set(other_genes)
        
        # Add the unique gene(s) to the dictionary
        unique_group_genes[cluster_label] = list(unique_genes)

    # Print or use unique_group_genes as needed
    for group, genes in unique_group_genes.items():
        print(f'Group {group}: Unique Genes: {genes}')
    
    df = pd.DataFrame.from_dict(unique_group_genes, orient='index')
    df.to_excel(writer, sheet_name='Cluster Analysis',startrow=1, startcol=marker_genes_df.shape[0]+5)
    worksheet.write(0, marker_genes_df.shape[0]+5, "Unique Genes Per Group")

    marker_genes = []
    i=1
    while True:
        user_input=None
        if i==1:
            user_input = input('Select 1st Marker Gene (enter "all" for all unique genes or press enter to have 1 unique gene from each cluster as marker gene): ')
        else:
            user_input = input(f'Select Marker Gene {i} (or press enter to end): ')
        
        if not user_input:
            break
        elif user_input=='all':
            break
        elif user_input not in list(adata.raw.var.index): #and user_input not in [gene for genes in potential_cell_types.values() for gene in genes]:
            print('Invalid Gene')
        else:
            marker_genes.append(user_input)
            i=i+1

    #Ask for manually generated marker genes or automatically generated
    if len(marker_genes)==0:
        if user_input=='all':
            marker_genes=[gene for genes in unique_group_genes.values() for gene in genes]
        else:
            num_marker_genes_per_group = 1
            marker_genes=[gene for genes in unique_group_genes.values() for gene in genes[0:num_marker_genes_per_group]]
        
    adata.markerGenes = marker_genes
    return(adata, writer)
