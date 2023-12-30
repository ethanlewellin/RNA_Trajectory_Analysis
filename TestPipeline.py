from EndProcess import EndProcess
from DPT import DPT
from PAGA import PAGA
from ClusterAnalysis import ClusterAnalysis
from CreateClusters import CreateClusters
from PCA import PCA
from Read_Data import Read_Data
from Preprocessing import Preprocessing

adata, results = Read_Data()
adata, results = Preprocessing(adata=adata,writer=results)
adata, results = PCA(adata=adata, writer=results)
adata, results = CreateClusters(adata=adata, writer=results)
adata, results = ClusterAnalysis(adata=adata, writer=results)
adata, results = PAGA(adata=adata, writer=results)
adata, results = DPT(adata=adata, writer=results)
EndProcess(adata,results)