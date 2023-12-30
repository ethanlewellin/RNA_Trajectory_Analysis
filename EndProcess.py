#End Process and Close Report

def EndProcess(adata, writer):
    
    df = adata.var
    df.to_excel(writer, sheet_name='Final Variables')
    
    df = adata.obs
    df.to_excel(writer, sheet_name='Final Observations')
    writer.close()
    
    return