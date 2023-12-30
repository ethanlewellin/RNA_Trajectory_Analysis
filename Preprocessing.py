#Preprocess the data

#Imports 
import scanpy as sc
import os
import contextlib
import warnings
import xlsxwriter

def Preprocessing(adata, writer):
    
    while True:
        user_input = ('Is your Data Preprocessed (Y/n): ')
        if user_input == 'Y' or user_input == 'n':
            break
        else:
            print('Invalid Input')
            
    if user_input == 'Y':
        return(adata,writer)
    
    #EDA Visualizations
    warnings.filterwarnings("ignore")
    path = os.getcwd()
    while True:
        mitochondrial_input = input("What is the mitochondrial prefix for this data ex(mt-,Mt-,MT-): ")
        if sum(adata.var_names.str.startswith(mitochondrial_input)) == 0:
            user_confirm = input('0 mitochondrial genes found. Do you want to re-enter prefix? (Y/N): ')
            if user_confirm != 'Y' and user_confirm != 'N':
                print('Invalid Input')
            elif user_confirm == 'N':
                    break
        else:
            break
    
    book = writer.book        
    worksheet = book.add_worksheet("EDA")  #Create EDA page in report
    worksheet.write(0, 0, "Exploratory Data Analysis")

    with contextlib.redirect_stdout(None):
        sc.pl.highest_expr_genes(adata, n_top=20, show=False, save='.png') #highest expressed genes
        worksheet.insert_image('A5', path+'\\figures\\highest_expr_genes.png')        
        
        adata.var[mitochondrial_input] = adata.var_names.str.startswith(mitochondrial_input)  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=[mitochondrial_input], percent_top=None, log1p=False, inplace=True)
        
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', f'pct_counts_{mitochondrial_input}'],
                    jitter=0.4, multi_panel=True, show=False, save='.png') #violin plot
        worksheet.insert_image('A30',path+'\\figures\\violin.png')
        
        sc.pl.scatter(adata, x='total_counts', y=f'pct_counts_{mitochondrial_input}',
                    show=False, save=f'_pct_counts_{mitochondrial_input}_.png')
        worksheet.insert_image('A55',path+f'\\figures\\scatter_pct_counts_{mitochondrial_input}_.png')

        
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                    show=False, save='_n_genes_by_counts.png')
        worksheet.insert_image('A75',path+'\\figures\\scatter_n_genes_by_counts.png')

    df = adata.var
    df.to_excel(writer, sheet_name='Raw Variables')
    
    df = adata.obs
    df.to_excel(writer, sheet_name='Raw Observations')

    print('View Figuers Folder to View EDA and Determine Preprocessing Values')

    limGeneCounts = limCellCounts = numGenes = dataScale = None
    limMTPct = False

    while True:
        limGeneCounts = input("Set minimum number of genes expressed required for a cell to pass filtering (or press Enter to skip): ")
        if not limGeneCounts:
            break
        else:
            try:
                int(limGeneCounts)
                break
            except:
                print('Invalid input, integer required')
                
    while True:
        limCellCounts = input("Set minimum number of counts required for a cell to pass filtering (or press Enter to skip): ")
        if not limCellCounts:
            break
        else:
            try:
                int(limCellCounts)
                break
            except:
                print('Invalid input, integer required')

    while True:
        limMTPct = input("Set percentage of mitochondrial counts for a gene to be excluded (1=1%) (or press Enter to skip): ")
        if not limMTPct:
            break
        else:
            try:
                float(limMTPct)
                break
            except:
                print('Invalid input, float required')
                
    while True:
        numGenes = input("Set number of genes to keep in data (selects the top n highly variable genes) (or press Enter to skip): ")
        if not numGenes:
            break
        else:
            try:
                int(numGenes)
                break
            except:
                print('Invalid input, integer required')
                
    while True:
        dataScale = input("Set Max value for data scaling (or press Enter to skip): ")
        if not dataScale:
            break
        else:
            try:
                float(dataScale)
                break
            except:
                print('Invalid input, float required')
                
    with contextlib.redirect_stdout(None):
        
        if limGeneCounts:
            sc.pp.filter_cells(adata, min_genes=int(limGeneCounts))
        if limCellCounts:
            sc.pp.filter_cells(adata, min_counts=int(limCellCounts))
        if limMTPct:
            adata = adata[adata.obs[f'pct_counts_{mitochondrial_input}'] < float(limMTPct), :] #Get rid of genes with high mitochondrial counts
        
        sc.pp.normalize_total(adata, target_sum=1e4) #Normalize counts per cell.
        sc.pp.log1p(adata) #Logarithmize the data matrix.
        
        if numGenes:
            sc.pp.highly_variable_genes(adata, n_top_genes=int(numGenes)) #Select only highly variable genes (change number)
            sc.pl.highly_variable_genes(adata, show=False, save='.png')
            
            worksheet.insert_image(row=100,col=0, filename=path+'\\figures\\filter_genes_dispersion.png')
            adata = adata[:, adata.var.highly_variable]

        if dataScale:
            sc.pp.scale(adata, max_value=float(dataScale))        
    return(adata,writer)