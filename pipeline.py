import os
import scanpy as sc
import magpy.settings as settings
from shutil import copyfile

# import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# from magpy.functions import *

#########################
##### Main Pipeline #####
#########################

sep = '-----------------------------------------------------------------------------------------------------'

#Reads in the meta-data table and returns meta-data for the speficied sample
def load_metadata(sample_ID=None):

	annotation_dict = dict()
	table_path = '/data/lab_data/single_cell_meta_data_table.tsv'

	if sample_ID:
		for line in open(table_path, 'r'):
			elem = str.split(line.rstrip())
			if elem:
				if elem[0] not in annotation_dict:
					annotation_dict[elem[0]] = elem[1:]

		metadata = annotation_dict
		return(metadata[sample_ID])

	else:
		df = pd.read_csv(table_path,sep="\t",index_col=0)
		return df

#Make a copy of the raw data file into an experiment directory
def copy_data(raw_path='',expt_path='',sample_ID=None, read_file=None, write_file=None):

	if read_file is None: read_file = settings.raw_file
	if write_file is None: write_file = settings.raw_file

	if sample_ID is not None:
		metadata = load_metadata(sample_ID)
		raw_path = str('/data/lab_data/'+metadata[1]).rstrip(read_file)
		new_path = expt_path+metadata[0]

	elif raw_path == '':
		raise Exception("Please provide either a sample ID or read path.")

	else:
		new_path = expt_path

	read_path = os.path.join(raw_path,read_file)
	write_path = os.path.join(new_path,write_file)

	print(f"Copying raw data from {read_path}...")

	if not os.path.exists(raw_path):
		raise Exception("The specified raw data directory does not exist.")

	if not os.path.isfile(read_path):
		raise Exception(f"The specifed directory does not contain a file named {read_file}")

	if not os.path.exists(new_path):
		print("The specified experiment directory does not exist.")
		print(f"Creating a new directory at {expt_path}.")
		os.makedirs(new_path)

	copyfile(read_path,write_path)
	print(f'\nRaw data copied to: {write_path}.\n')

	return(new_path)

#Load an h5 file containing 10X raw filtered output and add some annotations
def annotate(expt_path='', sample_ID=None, data=None, hash_clusters=False, hashtag='totalseq', save=True, min_filter = True, read_file=None, write_file=None, write_pfx = None, show=True, species=None):
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
	elif species is None: raise Exception("Please specify either a sample_ID or species of \'mouse\' or \'human\'")
	elif species=='mouse':
		mito_prefix = 'mt-'
		tissue = 'mm'
	elif species=='human':
		mito_prefix = 'MT-'
		tissue = 'hs'
	else: raise Exception("Species name invalid, please enter either \'mouse\' or \'human\'")


	print(sep,'\n')

	if data is None: 
		print(f"Reading data from {raw_path}")
		adata = sc.read_10x_h5(raw_path,gex_only=False)
		print(adata)
	else: adata = data
	print(f'Data contains {len(adata.obs_names)} cells and {len(adata.var_names)} genes.') 

	adata.var_names_make_unique()
	adata.obs_names_make_unique()

	# Annotate hashtag information
	if hash_clusters:
		from magpy.evaluate_hashtags import evaluate_hashtags
		print("Annotating data with hashtag information...")
		genes = []
		hashtags = []

		#Find all gene names containing the hashtag search term
		for item in adata.var_names.tolist():
			if hashtag.lower() in item.lower(): hashtags.append(item)
			else: genes.append(item)
		if len(hashtags) == 0: raise Exception(f"No hashtags found matching {hashtag}.")

		#Run hashtag demultiplexing algorithm
		hdata = adata[:,hashtags].copy()
		hdata = evaluate_hashtags(hdata, n_clusters=hash_clusters, write_pfx=write_pfx)
		adata.obs['hash_label'] = hdata.obs['hash_label']
		adata.obs['hash_cluster'] = hdata.obs['cluster']

		#Subset out hashtag info from gene expression matrix
		print("Subsetting out detected multiplets")
		adata = adata[hdata.obs['num_positive']<2,genes]
		
	#Initial pre-filtering based on absolute minimum thresholds
	if min_filter is True:
		sc.pp.filter_cells(adata, min_genes=200)
		sc.pp.filter_genes(adata, min_cells=3)
		print(f'{len(adata.obs_names)} cells and {len(adata.var_names)} genes kept based on minimum QC.')
	
	#Add annotations from metadata_table to adata objects for each sample
	if sample_ID is not None:
		metadata_fields = load_metadata('ID')
		for count,field in enumerate(metadata_fields[4:],4):
			meta_name = field
			meta_value = metadata[count]
			adata.obs[meta_name] = meta_value
	
	#Calculate and annotate qc metrics
	print(adata)
	adata.var['mt'] = adata.var_names.str.startswith(mito_prefix)  # annotate the group of mitochondrial genes as 'mt'
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	
	#Save annotated data to file
	if save:
		print(f"Saving annotated data to {adata_path}")
		adata.write(adata_path)
	print("Annotation complete.\n")
	
	return(adata)

#Take annotated data and filter based on QC thresholds
def preprocess(expt_path='', sample_ID=None, data=None, save=True, read_file=None, write_file=None, show=True, normalize=True, groupby = None):
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
	print("\nFiltering data...")
	
	print(adata)

	start_cells = len(adata.obs_names)
	start_genes = len(adata.var_names)	

	if show: 
		sc.pl.violin(adata, ['total_counts','n_genes_by_counts','pct_counts_mt'],jitter=0.4, multi_panel=True)

	#Filter out cells with low UMI counts
	print('\nFiltering out cells with counts below threshold: ', settings.min_counts_per_cell)
	adata = adata[adata.obs['total_counts'] > settings.min_counts_per_cell, :]
	if show: print(f'{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells pass filtering.')

	#Filter out doublets, indicated by abnormally high UMI counts
	print('\nFiltering out likely doublets with number of counts above: ', settings.max_counts_per_cell)
	adata = adata[adata.obs['total_counts'] < settings.max_counts_per_cell, :]
	if show: print(f'{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells pass filtering.')

	#Filter out cells with low transcriptome complexity (num unique genes)
	print('\nFiltering out cells with genes below threshold: ', settings.min_genes_per_cell)
	adata = adata[adata.obs['n_genes_by_counts'] > settings.min_genes_per_cell, :]
	if show: print(f'{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells pass filtering.')

	#Filter out doublets, indicated by abnormally high number of unique genes
	print('\nFiltering out likely doublets with number of genes above: ', settings.max_genes_per_cell)
	adata = adata[adata.obs['n_genes_by_counts'] < settings.max_genes_per_cell, :]
	if show: print(f'{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells pass filtering.')
	
	if show: 
		sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color = groupby)
		sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color = groupby)

	#Filter out cells with a high mitochondrial ratio
	print('\nFiltering out cells with more than ', settings.max_percent_mito,'% mitochondrial gene content')
	adata = adata[adata.obs['pct_counts_mt'] < settings.max_percent_mito, :]
	if show: print(f'{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells pass filtering.')

	if show:
		sc.pl.violin(adata, ['total_counts','n_genes_by_counts','pct_counts_mt'],jitter=0.4, multi_panel=True)
		sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', color = groupby)
		sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color = groupby)
	
	#Summarize the data
	print(f'\n{len(adata.var_names)}/{start_genes} genes and {len(adata.obs_names)}/{start_cells} cells in final dataset.')
	
	#Save preprocessed data to file
	if save:
		print(f"Saving preprocessed data to {pp_path}")
		adata.write(pp_path)
	print("Preprocessing complete.\n")
	
	return(adata)	

#Log-transform, regress out confounding variables, scale, and run PCA
def process(expt_path='', data=None, save=True, fig_pfx="", read_file=None, write_file=None, show=True, 
	merged=False, annotate_cell_cycle=True, regress_cell_cycle=False, integrate_cycling_cells=False,
	calc_hvg = True, save_pca = False, n_comps = 50):
	
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

	#Normalize values to percent of total read counts in that cell
	print("Normalizing per-cell counts to dataset median...")
	sc.pp.normalize_total(adata)
	#print(f"Normalized to {adata.X.todense()[0].sum()} counts per cell.")

	#Log transform, then save raw data for later plotting and differential expression testing
	print("Log-tranforming data...")
	sc.pp.log1p(adata)
	adata.raw = adata
	
	#Regress out effects of total reads per cell and percentage mitochondrial genes
	print("Regressing out effects of cell read quality...")
	sc.pp.regress_out(adata,['total_counts','pct_counts_mt']) 

	#Annotate genes that show high variability in dataset
	if calc_hvg:
		print("Determining highly variable genes...")
		sc.pp.highly_variable_genes(adata, min_mean=settings.min_mean, max_mean=settings.max_mean, min_disp=settings.min_disp)
	
	#Scale each gene to unit variance. Clip extreme outliers
	print("Scaling data...")
	sc.pp.scale(adata, max_value=settings.max_scaled_value)
	
	#Scrublet testing
	# scrub = scr.Scrublet(adata.raw.X)
	# adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets()
	# scrub.plot_histogram()
	# sum(adata.obs['predicted_doublets'])
	# adata.obs['doublet_info'] = adata.obs["predicted_doublets"].astype(str)
	# sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, groupby = 'doublet_info', rotation=45)

	#Load cell cycle genes from regev data set
	if annotate_cell_cycle or integrate_cycling_cells:
		print("Finding cell-cycle genes...")
		dir_path = os.path.dirname(os.path.realpath(__file__))
		with open(os.path.join(dir_path,'regev_lab_cell_cycle_genes.txt'),'r') as file:
			human_cell_cycle_genes = [gene.strip() for gene in file]
		mouse_cell_cycle_genes = [gene.lower().capitalize() for gene in human_cell_cycle_genes]

		#Determine if working with mouse or human genes
		s_genes = mouse_cell_cycle_genes[:43] + human_cell_cycle_genes[:43]
		s_genes = [gene for gene in s_genes if gene in adata.var_names]

		g2m_genes = mouse_cell_cycle_genes[43:] + human_cell_cycle_genes[43:]
		g2m_genes = [gene for gene in g2m_genes if gene in adata.var_names]

		#Annotate cell-cycle genes and cell phase
		print("Determining cell phases...")
		adata.var['s_phase'] = adata.var_names.isin(s_genes)
		adata.var['g2m_phase'] = adata.var_names.isin(g2m_genes)
		sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

	#Regress out cell cycle effects if desired
	if regress_cell_cycle:
		print("Regressing out effects of cell cycle...")
		sc.pp.regress_out(adata,['S_score', 'G2M_score'])
		sc.pp.scale(adata)

	#Re-calculate highly variable genes with S and G2M phase cells temporarily removed
	if integrate_cycling_cells:
		print("Re-calculating highly variable genes using only G1 cells")
		raw_adata = adata.raw.to_adata()
		raw_data = raw_adata[raw_adata.obs['phase'] == 'G1',:]
		df = sc.pp.highly_variable_genes(raw_adata, inplace=False, min_mean=settings.min_mean, max_mean=settings.max_mean, min_disp=settings.min_disp)
		for col in df:
			adata.obs[col] = df[col]
	
	#Run PCA as initial dimension reduction
	print("Computing principal components...")
	if not n_comps:	
		print(n_comps)
		n_comps = settings.num_pcs
		print(n_comps)
	print(n_comps)
	sc.tl.pca(adata, svd_solver='arpack',use_highly_variable=True, n_comps = n_comps)

	#Rank genes according to contributions to PCs.
	if save_pca: 
		sc.pl.pca_loadings(adata, show=show, components=list(range(1,11)), save=f'{fig_pfx}_PCA-loadings.pdf')
		sc.pl.pca_variance_ratio(adata, log=True, n_pcs=100, save=f'{fig_pfx}_elbowPlot.pdf', show=show)
	else: 
		sc.pl.pca_loadings(adata, show=show, components=list(range(1,11)), save=False, )
		sc.pl.pca_variance_ratio(adata, log=True, n_pcs=100, save=False, show=show)
	
	#Draw the PCA elbow plot to determine which PCs to use
	
	
	#Save processed data to a file
	if save:
		print(f"Saving processed data to {pca_path}")
		adata.write(pca_path)
	print("Processing complete.\n")
	
	return(adata)

#Calculate nearest neighbors, clustering, and dimesionality reduction
def cluster(expt_path='', data=None, save=True, fig_pfx="", read_file=None, write_file=None, show=True, neighbors_key=None, harmonize=False, calc_leiden = True, **kwargs):
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
	
	if harmonize:
		print("Running harmonypy integration of principal components...")
		sc.external.pp.harmony_integrate(adata,key='batch', **kwargs)

	#Compute nearest-neighbors graph
	if harmonize: neighbors_key = "X_pca_harmony"
	elif neighbors_key is None: neighbors_key = 'X_pca'
	sc.pp.neighbors(adata, n_neighbors=settings.num_neighbors, n_pcs=settings.num_pcs, use_rep = neighbors_key)

	#Calculate cell clusters via leiden algorithm
	if calc_leiden:
		print("Calculating Leiden clusters...")
		sc.tl.leiden(adata, resolution=settings.leiden_resolution)

	#Run PAGA to generate initial positions
	sc.tl.paga(adata, groups='leiden')
	sc.pl.paga(adata, show=False)
	
	#Run UMAP Dim reduction
	print("Running UMAP dimesionality reduction...")
	sc.tl.umap(adata, init_pos='paga', min_dist=settings.min_dist,maxiter=settings.maxiter,spread=settings.spread,gamma=settings.gamma)
	
	#Plot some basic umap QCs
	sc.pl.umap(adata, color='phase', save=f'_{fig_pfx}cell_cycle.pdf', legend_loc='right margin', alpha=0.95, show=show)
	sc.pl.umap(adata, color='leiden', save=f'_{fig_pfx}leiden_clusters.pdf', legend_loc='on data', alpha=0.95, show=show)
	# sc.pl.umap(adata, color='treatment', save=f'_{fig_pfx}treatment.pdf', legend_loc='right margin', alpha=0.95, show=show)
	
	#Save processed data to a file
	if save:
		print(f"Saving clustered data to {cluster_path}")
		adata.write(cluster_path)
	print("Clustering complete.\n")
		
	return(adata)

def recluster(adata, save=False, expt_path=None, write_file=None, recalc_hvg=True, hvg_num =None, recalc_pca=True, recalc_leiden=True, reembed=True, use_raw = True, regress_qc_params = False, 
	harmonize=True,  annotate_cell_cycle=True, regress_cell_cycle=False, num_pcs=None, num_neighbors=None, leiden_res=1.0, rescale_data = True, **kwargs):

	#Generate full write path
	if save:
		if expt_path is None: raise Exception('expt_path not specified.')
		if write_file is None: raise Exception('write_file not specified.')
		write_path = os.path.join(expt_path,write_file)
	
	#Prevents error from recomuputing PAGA
	if 'leiden_colors' in adata.uns: del adata.uns['leiden_colors']
	if num_pcs is None: num_pcs = settings.num_pcs
	
	#Pull out raw depth-normalized data and copy all data layers
	if use_raw:
		raw_adata = adata.raw.to_adata()
		for layer in adata.layers:
			raw_adata.layers[layer] = adata.layers[layer]
	else:
		raw_adata = adata.copy()
		
	#Regress out effects of total reads per cell and percentage mitochondrial genes
	if regress_qc_params:
		print("Regressing out effects of cell read quality...")
		sc.pp.regress_out(adata,['total_counts','pct_counts_mt']) 

	#Annotate genes that show high variability in dataset
	if recalc_hvg:
		print("Determining highly variable genes...")
		sc.pp.highly_variable_genes(raw_adata, min_mean=settings.min_mean, 
			max_mean=settings.max_mean, min_disp=settings.min_disp, batch_key='donor', n_top_genes = hvg_num)
		print("Highly variable genes: %d"%sum(raw_adata.var.highly_variable))
	else: raw_adata.var['highly_variable'] = adata.var['highly_variable']
	
	#Scale each gene to unit variance.
	if rescale_data:
		print("Scaling data...")
		sc.pp.scale(raw_adata, max_value=settings.max_scaled_value)
	
	#Annotate and/or regress cell cycle information
	if annotate_cell_cycle or regress_cell_cycle:
		raw_adata.obs['phase'] = adata.obs['phase']
		raw_adata.obs['S_score'] = adata.obs['S_score']
		raw_adata.obs['G2M_score'] = adata.obs['G2M_score']
	
		#Regress out cell cycle effects if desired
		if regress_cell_cycle:
			print("Regressing out effects of cell cycle...")
			sc.pp.regress_out(raw_adata,['S_score', 'G2M_score'])
			sc.pp.scale(raw_adata)
	
	#Run PCA as initial dimension reduction
	if recalc_pca: 
		print("Computing principal components...")
		sc.tl.pca(raw_adata, svd_solver='arpack', use_highly_variable=True, n_comps=num_pcs)

		#Display data on PCAs
		#sc.pl.pca_loadings(raw_adata, components=list(range(1,21)))
		sc.pl.pca_variance_ratio(raw_adata, log=True, n_pcs=num_pcs)
	
	#Integrate principal component data with HarmonyPy
	if harmonize:
		print("Running harmonypy integration of principal components...")
		sc.external.pp.harmony_integrate(raw_adata,key='batch',**kwargs)
		raw_adata.obsm['X_pca'] = raw_adata.obsm['X_pca_harmony']
	
	#Compute nearest neighbors
	print("Computing neighborhood graph...")
	if num_neighbors is None: num_neighbors = settings.num_neighbors
	sc.pp.neighbors(raw_adata, n_neighbors=num_neighbors, n_pcs=num_pcs)

	#Re-calculate cell clusters via leiden algorithm
	if recalc_leiden:
		print("Calculating Leiden clusters...")
		sc.tl.leiden(raw_adata, resolution=leiden_res)

	#Re-run PAGA and UMAP Dim reduction
	if reembed:
		print("Running UMAP dimensionality reduction...")
		sc.tl.paga(raw_adata, groups='leiden')
		sc.pl.paga(raw_adata)
		sc.tl.umap(raw_adata, init_pos='paga', min_dist=settings.min_dist, 
			maxiter=settings.maxiter, spread=settings.spread, gamma=settings.gamma)
	
	#Plot some basic umap QCs
	sc.pl.umap(raw_adata, color='lineage', legend_loc='right margin', alpha=0.95)
	sc.pl.umap(raw_adata, color='leiden', legend_loc='right margin', alpha=0.95)
	
	
	#Save processed data to a file
	if save:
		print(f"Saving reclustered data to {write_path}")
		raw_adata.write(write_path)

	print("Clustering complete.\n")
	return(raw_adata)
