##### FUNCTION BACKUP #####
# #DEPRECATED - Testing should be done with general linear modeling
# #Perform Benjamini-Hochberg multiple-testing correction on an ndarray of p-values
# def bh_correct(pvals):
# 	ranked_pvals = rankdata(pvals)
# 	adjusted_pvals = (pvals * len(pvals)) / ranked_pvals
# 	adjusted_pvals[adjusted_pvals>1] = 1
# 	return(adjusted_pvals)

	
#############################
##### PLOT SELECT GENES #####
#############################

# expressed_genes = adata.var_names
# genes_to_plot = search(adata,"myh")

# for gene in genes_to_plot:
#     if gene not in expressed_genes:
#         print(f"Sorry, {gene} is not present in this dataset")
#     else:
#         sc.pl.umap(adata, color=gene, alpha=0.95)


##############################
##### PLOT PC COMPONENTS #####
##############################

# #Create a dict of PC loading values {'PC1':pd.Series}, and sort to descending order
# col_names = ["PC"+str(n) for n in range(1,11)]
# pc_loadings = pd.DataFrame(adata.varm['PCs'][:,:10], index=adata.var.index, columns=col_names)
# PC_dict = {key:pc_loadings[key].abs().sort_values(ascending=False) for key in col_names}

# for gene in PC_dict["PC2"].index.tolist()[:11]:
#     sc.pl.umap(adata, color=gene, alpha=0.95)


############################################
##### T-TEST CLUSTERS FOR SIGNIFICANCE #####
############################################

# #Subset cells based on louvain cluser into dict of Anndata objects
# clusters = adata.obs['louvain'].unique()
# cluster_dict = {}
# anticluster_dict = {}

# #Build Anndata views for each cluster and anti-cluster
# for cluster in clusters:
#     cluster_adata = adata[adata.obs['louvain'] == cluster,:]
#     cluster_dict[cluster] = cluster_adata
    
#     anti_cluster = adata[adata.obs['louvain'] != cluster,:]
#     anticluster_dict[cluster] = anti_cluster

# #Run t-testing
# for cluster in clusters:
#     #Compare each cluster to its anti-cluster
#     print(f"Comparing cluster {str(cluster)} to its anti-cluster...")
#     _t,pvals = stats.ttest_ind(cluster_dict[cluster].X, anticluster_dict[cluster].X, equal_var=False)
#     adata.var[f'clust_{str(cluster)}'] = bh_correct(pvals)
    
#     #Perform cluster x cluster pairwise comparisons
#     for other_clust in clusters:
#         if cluster == other_clust: pass
#         else:
#             print(f"Comparing clusters {str(cluster)} and {str(other_clust)}...")
#             _t,pvals = stats.ttest_ind(cluster_dict[cluster].X, cluster_dict[other_clust].X, equal_var=False)
#             adata.var[f'clust_{str(cluster)}v{str(other_clust)}'] = bh_correct(pvals)
#     print()
    
# sorted_pvals = adata.var['clust_2v4'].sort_values()


##############################################
##### Process and cluster paired samples #####
##############################################

# expt_path += "/paired_expts"
# sc.settings.figdir = expt_path + "/figures"
# adata = load(expt_path, 2)
# for hash1 in range(1,5):
#     for hash2 in range(hash1+1,6):
#         print(f"Current samples: {hash1} + {hash2}")
#         hash_data = adata[(adata.obs["Hash_id"] == hash1) | (adata.obs["Hash_id"] == hash2)].copy()
#         hash_data.var["total_expression"] = hash_data.X.sum(axis=0).A1
#         hash_data = hash_data[:,hash_data.var["total_expression"]>0]
#         hash_data = process(expt_path, data=hash_data, fig_pfx=f"_paired_{hash1}{hash2}_", write_file=f"processed_paired_{hash1}{hash2}.h5ad")
#         hash_data = cluster(expt_path, data=hash_data, fig_pfx=f"_paired_{hash1}{hash2}_", write_file=f"clustered_paired_{hash1}{hash2}.h5ad")
#         sc.pl.umap(hash_data, color='Hash_id', save=f'_paired_{hash1}{hash2}_hash_id.pdf', legend_loc='right margin', alpha=0.95)
#         print()       


###################################################
##### Process and cluster individual hashtags #####
###################################################

# for hashtag in range(1,6):
#     print(f"Current sample: {hashtag}")
#     hash_data = adata[adata.obs["Hash_id"] == hashtag].copy()
#     hash_data.var["total_expression"] = hash_data.X.sum(axis=0).A1
#     hash_data = hash_data[:,hash_data.var["total_expression"]>0]
#     hash_data = process(expt_path, data=hash_data, fig_pfx=f"Expt{hashtag}_", write_file=f"processed_hashtag_{hashtag}.h5ad")
#     hash_data = cluster(expt_path, data=hash_data, fig_pfx=f"Expt{hashtag}_", write_file=f"clustered_hashtag_{hashtag}.h5ad")
#     print()