#load libraries
from scipy.stats import ttest_ind,mannwhitneyu
import scanpy as sc
import pandas as pd 
from sklearn import metrics 
import numpy as np

"""
Function to perform a grid search between the hyper-parameters of cluster analysis
"""
        
def find_bh(
        ras_adata,                         #AnnData object used to store the RAS values
        resolutions=[1],                   #default range for resolution of Leiden algorithm
        n_pcs=[10],                        #default range for the number of PCs used to compute the neighboorhood of a cell 
        n_neighbors=[5],                   #default range for the number of neighboorhoods used to compute the neighboorhood of a cell 
        names_of_groud_truth=[]            #compute the Normalized Mutual Information of clustering using any element of this list as ground truth
        ):

    ras_adata_clustering=ras_adata.copy()
    
    if names_of_groud_truth:
        for name in names_of_groud_truth:
            if name not in ras_adata_clustering.obs.columns:
                print("feature not present!")
                return -1
    
    #%%pca analysis
    sc.tl.pca(ras_adata_clustering, svd_solver='arpack',n_comps=max(n_pcs))

    cluster_values_nmi=dict()
    cluster_values_sil=list()
    
    for name in names_of_groud_truth:
        cluster_values_nmi[name]=list()
    
    num_cluster=list()
    
    res_values=list()
    pcs_values=list()
    neigh_values=list()

    for npc in n_pcs:          
        for neigh in n_neighbors:
            adata=ras_adata_clustering.copy()

            sc.pp.neighbors(adata, n_neighbors=neigh, n_pcs=npc)   
            
            for res in resolutions:

                sc.tl.leiden(adata,resolution=res,key_added = "leiden"+str(res))
    
                res_values.append(res)
                pcs_values.append(npc)
                neigh_values.append(neigh)
                if len(set(adata.obs["leiden"+str(res)].values))>1:
                    cluster_values_sil.append(metrics.silhouette_score(adata.obsm['X_pca'][:,0:npc],adata.obs["leiden"+str(res)],metric='euclidean'))
                else:
                    cluster_values_sil.append(None)
                
                for name in names_of_groud_truth:
                    cluster_values_nmi[name].append(metrics.normalized_mutual_info_score(adata.obs[name],adata.obs["leiden"+str(res)]))
                num_cluster.append(len(set(adata.obs["leiden"+str(res)].values)))                   

    df=pd.DataFrame()

    df["res"]=res_values
    df["pcs_values"]=pcs_values
    df["neigh_values"]=neigh_values
    df["num_cluster"]=num_cluster

    df["cluster_values_sil"]=cluster_values_sil
    for name in names_of_groud_truth:
        df["cluster_values_nmi_"+str(name)]=cluster_values_nmi[name]

    return df

"""
Wrapper to function rank_genes_groups of Scanpy
"""
def rank_reactions_groups(*args, **kwargs):
    sc.tl.rank_genes_groups(*args, **kwargs)

"""
Function used to perform t-test between two lists
"""
def ttest(list1,list2):

    [val,pval]=ttest_ind(list1,list2)
    
    return pval

"""
Function used to perform MannWhitney Utest between two lists
"""
def mtest(list1,list2):

    if  np.sum(list1)==0 and np.sum(list2)==0:
        pval=1
        return pval
    else:
        [val,pval]=mannwhitneyu(list1,list2)
        return pval

    