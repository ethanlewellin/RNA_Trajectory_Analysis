#imports
#import loompy
from tkinter.filedialog import askopenfilename
import pathlib
import xlsxwriter
import scanpy as sc
import warnings
import pandas as pd

def Read_Data():
    file = askopenfilename()
    #file = "C:\\Users\\ethan\\Downloads\\dentate_gyrus_A_10X_V1.loom"
    file_extension = pathlib.Path(file).suffix

    ### Open File
    if file_extension == ".loom":
        warnings.filterwarnings("ignore")
        adata = sc.read_loom(filename = file # type: ignore
                        ,var_names='Gene') 
        adata.var_names_make_unique()
        warnings.filterwarnings("default")
    else:
        raise Exception('Data Type not accepted, currently only takes .loom')
        
    adata.X = adata.X.astype('float64')  # type: ignore # convert to get higher percision
    writer = pd.ExcelWriter(f'{pathlib.Path(file).stem}.xlsx', engine = 'xlsxwriter')
    writer.path = str(pathlib.Path(file).parent)+f'\\{pathlib.Path(file).stem}.xlsx' # type: ignore
    return(adata,writer)

