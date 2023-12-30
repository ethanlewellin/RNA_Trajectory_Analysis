#Create Clusters Based off of Adata

#Imports
from fileinput import filename
from sqlite3 import Row
import scanpy as sc
import os 
path = os.getcwd()

def CreateClusters(adata, writer):
    
    while True:
        user_input = ('Is your Data Pre-Clustered (Y/n): ')
        if user_input == 'Y' or user_input == 'n':
            break
        else:
            print('Invalid Input')
            
    if user_input == 'Y':
        return(adata,writer)
    
    resolutions = [0.1,0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    clusterMethod = None

    book = writer.book        
    worksheet = book.add_worksheet("CreateClusters")  #Create EDA page in report
    worksheet.write_row(0, 0, "Creation of Clusters")
    
    while clusterMethod != 'leiden' and clusterMethod != 'louvain':
        clusterMethod = input("Select clustering method ('leiden' or 'louvain'): ")
        if clusterMethod != 'leiden' and clusterMethod != 'louvain':
            print('Invalid Choice')
    
    if clusterMethod=='leiden':
        for resolution in resolutions:
            sc.tl.leiden(adata, resolution=resolution)
            adata.obs[f'leiden_{resolution}'] = adata.obs['leiden']
    
        sc.pl.umap(adata, color=[f'leiden_{resolution}' for resolution in resolutions],
                   legend_loc='on data', add_outline=True, save='_resolution.png')
        worksheet.insert_image('A5', path+f'\\figures\\umap_resolution.png')        
        columns_to_drop = [f'leiden_{resolution}' for resolution in resolutions]
            
    elif clusterMethod=='louvain':
        for resolution in resolutions:
            sc.tl.louvain(adata, resolution=resolution)
            adata.obs[f'louvain_{resolution}'] = adata.obs['louvain']
    
        sc.pl.umap(adata, color=[f'louvain_{resolution}' for resolution in resolutions],
                    legend_loc='on data', add_outline=True, save=f'_resolution.png')
        worksheet.insert_image('A5', path+'\\figures\\umap_resolution.png')        
        columns_to_drop = [f'louvain_{resolution}' for resolution in resolutions]
            
    adata.obs.drop(columns=columns_to_drop, inplace=True)

    best_resolution = 1
    while True:
        user_input = input("Which Resolution Do you want to Continue with? (Default=1) (View umap_resolution.png in figure folder): ")
        try:
            float(user_input)
            break
        except:
            print('Invalid Input, must enter float')

    best_resolution = float(user_input)
    if clusterMethod=='leiden':
        sc.tl.leiden(adata, resolution=best_resolution)
    elif clusterMethod=='louvain':
        sc.tl.louvain(adata, resolution=best_resolution)
        
    return(adata, writer)