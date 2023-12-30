#Principal Component Analysis

#Imports
import warnings
from numpy import save
import pandas as pd
import scanpy as sc
import os
def PCA(adata, writer):
    warnings.filterwarnings("ignore",category=pd.errors.PerformanceWarning)
    warnings.filterwarnings("ignore",category=UserWarning)
    warnings.filterwarnings("ignore",category=FutureWarning)

    df = adata.var
    df.to_excel(writer, sheet_name='Processed Variables')
    
    df = adata.obs
    df.to_excel(writer, sheet_name='Processed Observations')

    warnings.filterwarnings("ignore",category=pd.errors.PerformanceWarning)
    warnings.filterwarnings("ignore",category=UserWarning)
    warnings.filterwarnings("ignore",category=FutureWarning)
    path = os.getcwd()

    book = writer.book        
    worksheet = book.add_worksheet("PCA")  #Create EDA page in report
    worksheet.write(0, 0, "Principal Component Analysis")
    
    sc.tl.pca(adata, svd_solver='randomized')
    sc.pl.pca_variance_ratio(adata,log=True, n_pcs=50, save='.png')
    worksheet.insert_image('A5', path+'\\figures\\pca_variance_ratio.png')        

    N_pcs = 50
    N_pcs = input("How many principal components would you like to use (default 50): ")
    
    sc.pp.neighbors(adata,n_pcs=int(N_pcs))
    sc.tl.umap(adata)
    sc.tl.diffmap(adata)
    
    return(adata,writer)