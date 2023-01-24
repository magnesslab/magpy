files = """

# ================================================== #
# ============= FILES & DIRECTORIES ================ #
# ================================================== #

# Default pipeline file names

raw_file = "filtered_feature_bc_matrix.h5"
adata_file = "annotated_adata.h5ad"
pp_file = "preprocessed_adata.h5ad"
pca_file = "processed_adata.h5ad"
cluster_file = "clustered_adata.h5ad"
"""

autosave = True
autoshow = True
figdir = '/proj/lost+found/figures'
file_format_figs = '.svg'
plot_prefix = 'fig_'

raw_file = "filtered_feature_bc_matrix.h5"
adata_file = "annotated_adata.h5ad"
pp_file = "preprocessed_adata.h5ad"
merged_file = "merged_adata.h5ad"
pca_file = "processed_adata.h5ad"
cluster_file = "clustered_adata.h5ad"
loom_file = "velocyto.loom"
merged_loom_file = "merged_velocyto.loom"
clustered_velocity_file = "clustered_velo_adata.h5ad"
cluster_subset_file = "cluster_subset_adata.h5ad"


# --------------------------------------------------- #



preprocessing = """
# ================================================== #
# ================ PREPROCESSING =================== #
# ================================================== #

# Quality Filtering
# Set values based on histograms in annotation step
max_percent_mito = 50
max_percent_ribo = 40
max_percent_hb = 40
min_genes_per_cell = 2000
min_counts_per_cell = 5000

# Doublet Filtering
max_genes_per_cell = 10000
max_counts_per_cell = 80000

# Gene Filtering
min_cells_per_gene = 3
"""

max_percent_mito = 50
# max_percent_ribo = 40
# max_percent_hb = 40
min_genes_per_cell = 2000
min_counts_per_cell = 5000
max_genes_per_cell = 10000
max_counts_per_cell = 80000
min_cells_per_gene = 3



# --------------------------------------------------- #



processing = """
# ================================================== #
# ================== PROCESSING ==================== #
# ================================================== #

# While scaling, clip extreme outliers
max_scaled_value = 10

# Parameters for determining highly variable genes
# Values are in normalized deviations from the mean
min_mean = 0.0125
max_mean = 6
min_disp = 0.2
"""
min_mean = 0.0125
max_mean = 6
min_disp = 0.2
max_scaled_value = 10

# --------------------------------------------------- #



clustering = """
# ================================================== #
# ================== CLUSTERING ==================== #
# ================================================== #

# k nearest neighbor
num_neighbors = 50
num_pcs = 10

#leiden clustering
leiden_resolution = 0.6

# UMAP settings
min_dist = 0.3 
maxiter = None 
spread = 1
gamma = 1
"""
num_neighbors = 10
num_pcs = 50
leiden_resolution = 0.6
min_dist = 0.3 
maxiter = None 
spread = 1
gamma = 1


# --------------------------------------------------- #


marker_genes = '''
# ================================================== #
# ================= MARKER Genes =================== #
# ================================================== #
'''

min_in_group_fraction = 0.4
min_fold_change = 2
max_out_group_fraction = 0.2

# --------------------------------------------------- #