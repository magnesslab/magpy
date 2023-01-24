# Import libraries
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

from scipy import stats
from sklearn_extra.cluster import KMedoids
from statsmodels.discrete.discrete_model import NegativeBinomial

def fit_nbinom(y, method='nb2', verbose=True):
    '''
    Fits a negative binomial distribution to 1d input data (ndarray)
    Returns a tuple of (size,prob) for the fitted distribution
    '''
    x = np.ones(y.shape)
    res = NegativeBinomial(y, x, loglike_method=method,disp=verbose).fit()

    #Convert alpha and mu parameters to size and prob
    mu = np.exp(res.params[0])
    if method =='nb2': size = 1. / res.params[1]
    elif method == 'nb1': size = 1. / res.params[1] * mu
    else: raise Exception("Method must be nb1 or nb2") 
    prob = size / (size + mu)

    #Generate summary
    if verbose: 
        print(f"mu = {mu}")
        print(f"size = {size}")
        print(f"prob = {prob}")
        print()
    
    return (size,prob)


def evaluate_hashtags(adata, n_clusters=None, plot=True, plot_type='tsne', write_pfx = None):
    #Assign figure directory for outputs
    figdir = sc.settings.figdir
    
    #Subset out hashtag expression data and generate a new AnnData object
    adata.X = adata.X.todense()
    adata.raw = adata

    #Normalize to the geometric mean and log transform
    adata.X += 1 #Prevent errors with dividing by 0 or log(0)
    adata.X = np.log(adata.X / stats.gmean(adata.X))
    
    #Perform k-mediods clustering
    adata.obs['cluster'] = KMedoids(n_clusters=n_clusters or adata.shape[1], random_state=0).fit(adata.X).labels_
    adata.obs['cluster_str'] = adata.obs['cluster'].astype('str')
    
    cluster_map = {}
    for cluster in range(n_clusters or adata.shape[1]+1):
        raw_subset = adata.raw.to_adata()[adata.obs['cluster']==cluster,:]
        adata.var[f'cluster{cluster}_mean'] = raw_subset.X.mean(axis=0)
        cluster_map[cluster] = adata.var[f'cluster{cluster}_mean'].idxmax()
    
    #Do stuff n things
    adata.obs['num_positive'] = 0
    raw_adata = adata.raw.to_adata()
    if plot: fig,axes = plt.subplots(1,adata.shape[1],figsize=(30,5))
    
    #Fit neg binomial distribution for each hashtag
    print("Fitting noise distribution models...")
    for i,hashtag in enumerate(adata.var_names):
        #Get values for all negative clusters for given hashtag
        neg_clusters = [cluster for cluster in cluster_map if cluster_map[cluster] != hashtag]
        subset = raw_adata[adata.obs['cluster'].isin(neg_clusters),hashtag]
        counts = subset.X.ravel()
        
        #Remove outliers (top 1%)
        subset = subset[counts <= np.percentile(counts,99),:].copy()
        counts = subset.X.ravel()
        
        #Run fitting algorithm
        size,prob = fit_nbinom(counts,verbose=False)
        
        #Plot results
        if plot:
            x_vals = np.arange(1,int(counts.max()),1)
            y_vals = stats.nbinom.pmf(x_vals, size, prob)
            axes[i].hist(counts,bins=50,density=True)
            axes[i].plot(x_vals,y_vals,'r-')
            axes[i].set_title(hashtag)
            p_cutoff = stats.nbinom.isf(0.01,size, prob)
        
            axes[i].plot((p_cutoff,p_cutoff),(0,y_vals.max()),'b--')

        #Calculate p values for every cell
        adata.obs[hashtag+'_pval'] = stats.nbinom.sf(raw_adata[:,hashtag].X,size,prob).ravel()
        
        #Assign +/- based on pval, increment total num_positive for each significant hashtag
        adata.obs[hashtag+'_positive'] = adata.obs[hashtag+'_pval'] < 0.01
        adata.obs['num_positive'] += adata.obs[hashtag+'_positive']
        
    if plot:
        fig.savefig(f"{figdir}/{write_pfx}_negbinom_cutoffs.svg", dpi = 300, format = 'svg')
        plt.show()
 

    pd.options.display.float_format = "{:,.2f}".format
    print(adata.var)
    if plot:
        sc.pp.pca(adata)
        fig, axes = plt.subplots(1,2,figsize=(14,7))
        
        if plot_type == 'umap':
            sc.pp.neighbors(adata,use_rep='X')
            sc.tl.umap(adata)
            sc.pl.umap(adata,color='cluster_str',ax=axes[0],show=False)
        elif plot_type == 'tsne':
            sc.tl.tsne(adata,use_rep='X')
            sc.pl.tsne(adata,color='cluster_str',ax=axes[0],show=False)
        else: 
            raise Exception("Invalid plot_type - choose umap or tsne")

    #Count number of singlets/doublets
    unique, counts = np.unique(adata.obs['num_positive'].to_numpy(), return_counts=True)
    print("Multiplets - ",dict(zip(unique, counts)))
    
    #Find column with minimum p_value for remaining cells
    pval_cols = [col for col in adata.obs if '_pval' in col]
    Y = adata.obs[pval_cols].to_numpy()
    adata.obs['idxmin'] = Y.argmin(axis=1)
    adata.obs['hash_label'] = adata.obs['idxmin'].map({i:col[:-5] for i,col in enumerate(pval_cols)})
    
    #Label negative cells
    adata.obs['hash_label'] = adata.obs['hash_label'].where(adata.obs['num_positive']==1,'Negative')
    
    #Label multiplets
    adata.obs['hash_label'] = adata.obs['hash_label'].where(adata.obs['num_positive']<2,'Multiplet')
    
    #Display counts for each hashtag
    unique, counts = np.unique(adata.obs['hash_label'].to_numpy(), return_counts=True)
    print("Tags - ",dict(zip(unique, counts)))
    
    #Subset out multiplets
    initial = adata.shape[0]
    subset = adata[adata.obs['num_positive'] < 2,:].copy()
    final = subset.shape[0]
    print(f"{initial - final} cells identified as multiplets. {final}/{initial} cells remaining.")
    
    if plot: 
        if plot_type == 'tsne':
            sc.pl.tsne(subset,color='hash_label',ax=axes[1], show = False)
        elif plot_type == 'umap':
            sc.pl.umap(subset,color='hash_label',ax=axes[1], show = False)
        plt.show()
        fig.savefig(f"{figdir}/{write_pfx}_hashtag_clusters.svg", dpi = 300, format = 'svg')
    print()
    
    return adata