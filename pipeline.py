import os
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import magpy.settings as settings
from magpy.functions import *
from shutil import copyfile

#########################
##### Main Pipeline #####
#########################

sep = '-----------------------------------------------------------------------------------------------------'

## This reads in the meta-data table for each sample and holds it in memory
def load_metadata(sample_ID):
	
	annotation_dict = dict()
	table_path = '/proj/lab_data/single_cell_meta_data_table.tsv'
	
	for line in open(table_path, 'r'):
		elem = str.split(line.rstrip())
		if elem:
			if elem[0] not in annotation_dict:
				annotation_dict[elem[0]] = elem[1:]
		
	metadata = annotation_dict
	return(metadata[sample_ID])

#Make a copy of the raw data file into an experiment directory
def copy_data(expt_path='', sample_ID=None, read_file=None, write_file=None,):
	
	if read_file is None: read_file = settings.raw_file
	if write_file is None: write_file = settings.raw_file
	
	if sample_ID is not None:
		metadata=load_metadata(sample_ID)
		raw_path = str('/proj/lab_data/'+metadata[1]).rstrip(read_file)
		new_path = expt_path+metadata[0]
	
	read_path = os.path.join(raw_path,read_file)
	write_path = os.path.join(new_path,write_file)

	print(f"Copying raw data from {read_path}...")

	if not os.path.exists(raw_path):
		raise Exception("The specified raw data directory does not exist.")

	if not os.path.exists(read_path):
		raise Exception(f"The specifed directory does not contain a file named {read_file}")

	if not os.path.exists(write_path):
		print("The specified experiment directory does not exist.")
		print(f"Creating a new directory at {expt_path}.")
		os.makedirs(new_path)

	copyfile(read_path,write_path)
	print(f"Raw data copied to {write_path}.\n")
	
	return(new_path)

#Load an h5 file containing 10X raw filtered output and add some annotations
def annotate(expt_path='', sample_ID=None, data=None, hashed=False, save=True, read_file=None, write_file=None, show=True):
	print('Running initial filtering for:', sample_ID,'...')
	
	#Generate full paths for reading/writing
	if read_file is None: read_file = settings.raw_file
	if write_file is None: write_file = settings.adata_file
	raw_path = os.path.join(expt_path,read_file)
	adata_path = os.path.join(expt_path,write_file)	
	
	if sample_ID is not None:
		metadata=load_metadata(sample_ID)
		if metadata[4]=='mouse':
			mito_prefix = 'mt-'
			tissue = 'mm' 
		elif metadata[4]=='human':
			mito_prefix = 'MT-'
			tissue = 'hs'
		if metadata[6]=='True': 
			hashed = True

	print(sep,'\n')
		
	if data is None: 
		print(f"Reading data from {raw_path}")
		adata = sc.read_10x_h5(raw_path,gex_only=False)		
	else: adata = data
	print(f'Data contains {len(adata.obs_names)} cells and {len(adata.var_names)} genes.') 

	adata.var_names_make_unique()
	adata.obs_names_make_unique()

	# Annotate hashtag information
	if hashed:
		genes = []
		hashtags = []
		for item in adata.var_names.tolist():
			if 'hash' in item.lower(): hashtags.append(item)
			elif 'totalseq' in item.lower(): hashtags.append(item)
			else: genes.append(item)
		if len(hashtags) == 0: raise Exception("No hashtags found.")
		adata.obs['Hash_id'] = adata[:,hashtags].X.argmax(axis=1).A1 + 1
		for tag in hashtags:
			adata.obs[tag] = adata[:,tag].to_df()
		adata = adata[:,genes]
	
	#Initial pre-filtering based on minimum thresholds
	sc.pp.filter_cells(adata, min_genes=500)
	sc.pp.filter_genes(adata, min_cells=1)
	print(f'{len(adata.obs_names)} cells and {len(adata.var_names)} genes kept.')
	
	#Add annotation for total genes per cell
	num_counts = np.sum(adata.X, axis=1).A1
	adata.obs['n_counts'] = num_counts
	
	#Add annotations from metadata_table to adata objects for each sample
	if sample_ID is not None:
		metadata_fields = load_metadata('ID')
		for count,field in enumerate(metadata_fields[4:],4):
			meta_name = field
			meta_value = metadata[count]
			adata.obs[meta_name] = meta_value
	
	#Add Annotation for what percent of counts are mitochondrial genes
	print("Annotating mitochondral genes...")
	mito_genes = [gene for gene in adata.var_names if 'mt-' in gene.lower()]
	adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / num_counts
	
	#Annotate cell-cycle genes, call phase for each cell
	print("Annotating cell-cycle genes...")
	dir_path = os.path.dirname(os.path.realpath(__file__))
	with open(os.path.join(dir_path,'regev_lab_cell_cycle_genes.txt'),'r') as file:
		mouse_cell_cycle_genes = []
		human_cell_cycle_genes = []
		for gene in file:
			gene = gene.strip()
			mouse_cell_cycle_genes.append(gene.lower().capitalize())
			human_cell_cycle_genes.append(gene)
	s_genes = mouse_cell_cycle_genes[:43] + human_cell_cycle_genes[:43]
	g2m_genes = mouse_cell_cycle_genes[43:] + human_cell_cycle_genes[43:]

	#Annotate cell phase
	print("Determining cell phases...")
	adata.var['s_phase'] = adata.var_names.isin(s_genes)
	adata.var['g2m_phase'] = adata.var_names.isin(g2m_genes)
	sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
	
	#Summarize the data
	#summarize_adata(expt_path, adata, show, file_name=expt_path+"/raw_data_summary.png")
	
	#Save annotated data to file
	if save:
		print(f"Saving annotated data to {adata_path}")
		adata.write(adata_path)
	print("Annotation complete.\n")
	
	return(adata)

#Take annotated data; Filter and normalize
def preprocess(expt_path='', sample_ID=None, data=None, save=True, read_file=None, write_file=None, show=True):
	#Generate full paths for reading/writing
	if read_file is None: read_file = settings.adata_file
	if write_file is None: write_file = settings.pp_file
	adata_path = os.path.join(expt_path,read_file)
	pp_path = os.path.join(expt_path,write_file)
	
	if sample_ID is not None:
		metadata=load_metadata(sample_ID)
	
	if data is None:
		print(f"Reading adata from {adata_path}")
		adata = sc.read_h5ad(adata_path)
	else: adata = data
	print("Filtering data...")
	
	start_cells = len(adata.obs_names)
	start_genes = len(adata.var_names)
	
	#Filter out cells with a high mitochondrial ratio
	adata = adata[adata.obs['percent_mito'] < settings.max_percent_mito, :]

	#Filter out cells with low signal
	adata = adata[adata.obs['n_genes'] > settings.min_genes_per_cell, :]  
	adata = adata[adata.obs['n_counts'] > settings.min_counts_per_cell, :]
	
	#Filter out doublets, indicated by high counts
	adata = adata[adata.obs['n_genes'] < settings.max_genes_per_cell, :]  
	adata = adata[adata.obs['n_counts'] < settings.max_counts_per_cell, :]

	#Refilter genes to get rid of genes that are only in a tiny number of cells
	sc.pp.filter_genes(adata, min_cells=settings.min_cells_per_gene)

	#Normalize values to percent of total read counts in that cell
	sc.pp.normalize_total(adata,target_sum=1e6)
	
	#Summarize the data
	print(f'{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells in final dataset.')
	#summarize_adata(expt_path, adata, show, file_name="preprocessed_data_summary.png") 
	
	#Save preprocessed data to file
	if save:
		print(f"Saving preprocessed data to {pp_path}")
		adata.write(pp_path)
	print("Preprocessing complete.\n")
	
	return(adata)	

def combine_experiments(sample_list='', expt_path='', data=None, save=True, read_file=None, write_file=None):
	
	if read_file is None: read_file = settings.pp_file
	if write_file is None: write_file = settings.merged_file
	merged_path = os.path.join(expt_path,write_file)
	
	adata_list = []
	min_obs = float('inf')
	
	for sample_ID in sample_list:
		metadata=load_metadata(sample_ID)
		new_path = expt_path+metadata[0]
		pp_path = os.path.join(new_path,read_file)
		adata_list.append(sc.read_h5ad(pp_path))
	
	for adata in adata_list:
			if min_obs > adata.n_obs:
				min_obs = adata.n_obs
	
	for item in range(len(adata_list)):
		adata_list[item] = sc.pp.subsample(adata_list[item],n_obs=min_obs,copy=True)

	print('Combining all samples into one adata object and subsampling to lowest number of filtered cells')
	print('\nNumber of cells per sample included in final combined adata object: ',min_obs)
	adata = adata_list[0].concatenate(adata_list[1:],join='outer')
	print(adata)
	
	if save:
		print(f"Saving preprocessed data to {merged_path}")
		adata.write(merged_path)
	print("Sample combining complete.\n")	

#Log-transform, regress out confounding variables, scale, and run PCA
def process(expt_path='', data=None, save=True, fig_pfx="", read_file=None, write_file=None, show=True, merged=False):
	#Generate full read path
	
	if read_file is None and merged is False: read_file = settings.pp_file
	elif read_file is None and merged is True: 
		read_file = settings.merged_file
	if write_file is None: write_file = settings.pca_file
	pp_path = os.path.join(expt_path,read_file) 
	pca_path = os.path.join(expt_path,write_file)
	
	if data is None:
		print(f"Reading adata from {pp_path}")
		adata = sc.read_h5ad(pp_path)
	else: adata = data
	print("Processing data...") 
	#Log transform, then save raw data for later plotting and differential expression testing
	print("Log-tranforming data...")
	sc.pp.log1p(adata)
	adata.raw = adata

	#Regress out effects of total reads per cell and percentage mitochondrial genes
	print("Regressing out confounding variables...")
	sc.pp.regress_out(adata, ['n_counts','percent_mito']) 

	#Annotate genes that do not show enough variance
	print("Determining highly variable genes...")
	sc.pp.highly_variable_genes(adata, min_mean=settings.min_mean, max_mean=settings.max_mean, min_disp=settings.min_disp)
	
	#Scale each gene to unit variance. Clip extreme outliers
	print("Scaling data...")
	sc.pp.scale(adata, max_value=settings.max_scaled_value)
	
	#Run PCA to compute the default number of components
	print("Computing principal components...")
	sc.tl.pca(adata, svd_solver='arpack',use_highly_variable=True)

	#Rank genes according to contributions to PCs.
	sc.pl.pca_loadings(adata, show=False, components=list(range(1,21)), save=f'{fig_pfx}_PCA-loadings.pdf')
	
	#Draw the PCA elbow plot to determine which PCs to use
	sc.pl.pca_variance_ratio(adata, log=True, n_pcs=100, save=f'{fig_pfx}_elbowPlot.pdf', show=False)
	
	#Save processed data to a file
	if save:
		print(f"Saving processed data to {pca_path}")
		adata.write(pca_path)
	print("Processing complete.\n")
	
	return(adata)

#Calculate nearest neighbors, clustering, and dimesionality reduction
def cluster(expt_path='', data=None, save=True, fig_pfx="", read_file=None, write_file=None, show=True):
	#Generate full read path
	if read_file is None: read_file = settings.pca_file
	if write_file is None: write_file = settings.cluster_file
	pca_path = os.path.join(expt_path,read_file)
	cluster_path = os.path.join(expt_path,write_file)
		
	if data is None:
		print(f"Reading adata from {pca_path}")
		adata = sc.read_h5ad(pca_path)
	else: adata = data
	print("Clustering data...")
	
	#Compute nearest-neighbors graph
	sc.pp.neighbors(adata, n_neighbors=settings.num_neighbors, n_pcs=settings.num_pcs)

	#Calculate cell clusters via leiden algorithm
	print("Calculating Leiden clusters...")
	sc.tl.leiden(adata, resolution=settings.leiden_resolution)

	#Run UMAP Dim reduction
	print("Running UMAP dimesionality reduction...")
	sc.tl.umap(adata, min_dist=settings.min_dist,maxiter=settings.maxiter,spread=settings.spread,gamma=settings.gamma)
	
	#Plot some basic umap QCs
	sc.pl.umap(adata, color='phase', save=f'_{fig_pfx}cell_cycle.pdf', legend_loc='right margin', alpha=0.95, show=show)
	sc.pl.umap(adata, color='leiden', save=f'_{fig_pfx}leiden_clusters.pdf', legend_loc='right margin', alpha=0.95, show=show)
	# sc.pl.umap(adata, color='treatment', save=f'_{fig_pfx}treatment.pdf', legend_loc='right margin', alpha=0.95, show=show)
	
	#Save processed data to a file
	if save:
		print(f"Saving clustered data to {cluster_path}")
		adata.write(cluster_path)
	print("Clustering complete.\n")
	
	return(adata)
