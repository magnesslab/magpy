import os
import numpy as np
import scanpy as sc
import seaborn as sns
import magpy as mp
import loompy as lpy
import magpy.settings as settings
import matplotlib.pyplot as plt
import matplotlib as mpl
import scvelo as scv
import pandas as pd
from shutil import copyfile


def copy_loom(expt_path='', sample_ID=None, loom_file=None, read_file=None, write_file=None):
	
	if loom_file is None: loom_file = '.loom'
	if read_file is None: read_file = settings.raw_file
	if write_file is None: write_file = settings.loom_file
	
	if sample_ID is not None:
		metadata=mp.pipeline.load_metadata(sample_ID)
		raw_path = str('/proj/lab_data/'+metadata[1]).rstrip(read_file)
		new_path = expt_path+metadata[0]
		loom_file = raw_path.split('/')[5]+'.loom'
		velo_path = str(raw_path.rstrip('/outs')+'/velocyto/'+loom_file)
		read_path = str(raw_path.rstrip('/outs')+'/velocyto/')
	
	loom_path = os.path.join(velo_path)
	read_path = os.path.join(read_path)
	write_path = os.path.join(new_path,write_file)
	
	if not os.path.exists(read_path):
		raise Exception("The specified velocyto data directory does not exist.")

	if not os.path.exists(read_path):
		raise Exception(f"The specifed directory does not contain a file named {_file}")

	if not os.path.exists(new_path):
		print("The specified experiment directory does not exist.")
		print(f"Creating a new directory at {expt_path}.")
		os.makedirs(new_path)

	copyfile(loom_path,write_path)
	
	return()

def combine_loom(sample_list='', expt_path='', save=True, read_file=None, write_file=None):
	
	if read_file is None: read_file = settings.loom_file
	if write_file is None: write_file = settings.merged_loom_file
	merged_path = os.path.join(expt_path,write_file)
	
	ldata_list = []
	
	for sample_ID in sample_list:
		metadata=mp.pipeline.load_metadata(sample_ID)
		new_path = expt_path+metadata[0]
		loom_path = os.path.join(new_path,read_file)
		ldata_list.append(loom_path)
		print(loom_path)
	
	#print(ldata_list)
	
	print('Combining all loom files into one merged loom file')
	
	#files = ["/proj/lab_data/10x_data/m-011-012_SPEM_1mo_10x_output/011_ctl/velocyto/ctl.loom" ,"/proj/lab_data/10x_data/m-009-010_SPEM_4mo_10x_output/010_ctl/velocyto/010_ctl.loom"]

	lpy.combine(ldata_list, merged_path)
	
	print("Loom file merging complete.\n")

def process_velocity(expt_path='', data=None, save=True, read_file=None, write_file=None, loom_file=None, merged=False):
		
	#Build read path
	if loom_file is None and merged is False: loom_file = settings.loom_file
	elif loom_file is None and merged is True: 
		loom_file = settings.merged_loom_file
	if read_file is None: read_file = settings.cluster_file
	if write_file is None: write_file = settings.clustered_velocity_file

	adata_path = os.path.join(expt_path,read_file)
	save_path = os.path.join(expt_path,write_file)

	print(save_path)
	if data is None:
		print(f"Reading adata from {adata_path}")
		adata = sc.read_h5ad(adata_path)
	else: adata = data

	read_path = os.path.join(expt_path, loom_file)

	#Read data
	print(f"Reading loom file from {read_path}")
	
	ldata = scv.read(read_path, cache = True)
	adata = scv.utils.merge(adata,ldata)
	
	print("Annotating pre-processed adata object with RNA velocity data ...")
	print(adata)
	
	scv.pp.filter_and_normalize(adata, log = False)
	scv.pp.moments(adata)
	
	scv.tl.velocity(adata)
	scv.tl.velocity_graph(adata)
	scv.tl.velocity_embedding(adata)

	if save:
		print(f"Saving updated preprocessed data to {save_path}")
		adata.write(save_path)
	print("RNA velocity addition complete.\n")
	
	return(adata)

#Generate a panel of summary figures for preprocessing QC
def summarize_adata(expt_path, adata, show=False, file_name="summary.png"):
	save_path = os.path.join(expt_path,file_name)
	
	fig,axes = plt.subplots(3,2,figsize=(10,15))

	# axes[0,0].hist(adata.obs['n_genes'],bins=50)
	# axes[0,0].set_title("n_genes_per_cell")

	# axes[0,1].hist(adata.obs['n_counts'],bins=50)
	# axes[0,1].set_title("n_counts_per_cell")

	# axes[1,0].hist(adata.obs['percent_mito'],bins=50)
	# axes[1,0].set_title("percent_mito_per_cell")
	
	# axes[1,1].hist(adata.var['n_cells'],bins=50)
	# axes[1,1].set_title("n_cells_per_gene")
	
	sc.pl.violin(adata, 'n_genes', jitter=0.4, show = False, ax=axes[0,0])
	axes[0,0].set_title("n_genes_per_cell")
	sc.pl.violin(adata, 'n_counts', jitter=0.4, show = False, ax=axes[0,1])
	axes[0,1].set_title("n_counts_per_cell")
	sc.pl.violin(adata, 'percent_mito', jitter=0.4, show = False, ax=axes[1,0])
	axes[1,0].set_title("percent_mito_per_cell")
	sc.pl.scatter(adata, x='n_counts', y='percent_mito', show = False, ax=axes[1,1])
	axes[1,1].set_title("percent_mito_per_cell")
	sc.pl.scatter(adata, x='n_genes', y='percent_mito', show = False, ax=axes[2,0])
	axes[2,0].set_title("n_genes v percent_mito")
	
	plt.savefig(save_path)
	if show ==True:
		plt.show()

def create_gene_expression_list():
	expressed_dict = dict()
	
	if not draw_dpt:
		for gene in adata.raw.var_names.values.tolist():
			if gene not in expressed_dict:
				expressed_dict[str(gene)] = 1
	else:
		for gene in adata.raw.var_names.values.tolist():
			if gene not in expressed_dict:
				expressed_dict[str(gene)] = 1
	
	return expressed_dict

def calc_gene_scores():
	expressed_dict = create_gene_expression_list()
	
	genes_to_score= []
	unscored_genes = []
	key_list = []
	
	#print(adata.var_names.tolist())
	
	with open(score_list,'r') as file:
		for line in file:
			if ('=' in line):
				plot = line.strip().split(' = ')
				label = plot[0]
				genes = plot[1]
				list = genes.split(',')
				for gene in list:
					if gene in expressed_dict:
						genes_to_score.append(gene.rstrip())
					else:
						unscored_genes.append(gene.strip())
				if showProgress:
					print('Generating', label, 'scores . . . ')
					print('\nGenes included in scoring: ', label, ' - ', str(genes_to_score))
					print('\nGenes excluded from scoring: ', label, ' - ', str(unscored_genes),'\n') 
					print('Plotting gene scores with UMAP embedding . . . ')
					print('\nPlot being saved to: ',figdir+file_tag+label+'-scored_gene_expression.pdf' )
				name = '_'+file_tag+label+'_score_plot'
				sc.tl.score_genes(adata, gene_list = genes_to_score, score_name = label, use_raw = True)
				sc.pl.dotplot(adata, var_names = genes_to_score, groupby=dot_groups, save = '_'+label+'_scores', color_map = my_dot_cmap, log = True, show = False, mean_only_expressed=True, dot_min=0.2, dot_max=1)
				key_list = [cluster_key,label]
				sc.pl.umap(adata, color = key_list, save = '_'+label+'-scores', show = False, size = 25, use_raw = True, legend_loc = 'on data')
			genes_to_score = []
			unplotted_score = []

def draw_scatter_plots(adata):
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito','percent_Rp'], jitter=0.4, multi_panel=True, save = '_'+file_tag+'postFiltering1_plot', show = False)
	sc.pl.scatter(adata, x='n_counts', y='percent_mito', color = dot_groups, show = False, save = '_'+file_tag+'post_counts_v_mito')
	sc.pl.scatter(adata, x='n_counts', y='n_genes', color = dot_groups, show = False, save = '_'+file_tag+'post_counts_v_genes')
	sc.pl.scatter(adata, x='n_genes', y='percent_mito', color = dot_groups, show = False, save = '_'+file_tag+'post_genes_v_mito')
	sc.pl.scatter(adata, x='n_counts', y='percent_Rp', color = dot_groups, show = False, save = '_'+file_tag+'post_counts_v_Rp')
	sc.pl.scatter(adata, x='n_genes', y='percent_Rp', color = dot_groups, show = False, save = '_'+file_tag+'post_genes_v_Rp')
	sc.pl.scatter(adata, x='percent_mito', y='percent_Rp', color = dot_groups, show = False, save = '_'+file_tag+'post_mito_v_Rp')
	
	sc.pl.scatter(adata, x = 'SPEM', y = 'Proliferation', color = cluster_key, use_raw = True, save = '_SPEM_v_Prolif')
	sc.pl.scatter(adata, x = 'Proliferation', y = 'SPEM', color = cluster_key, use_raw = True, save = '_Prolif_v_SPEM')
	sc.pl.scatter(adata, x = 'Muc6', y = 'Proliferation', color = cluster_key, use_raw = True, save = '_Muc6_v_Prolif')
	sc.pl.scatter(adata, x = 'Muc5ac', y = 'Proliferation', color = cluster_key, use_raw = True, save = '_Muc5ac_v_Prolif')
	sc.pl.scatter(adata, x = 'Gif', y = 'Proliferation', color = cluster_key, use_raw = True, save = '_SPEM_v_Prolif')
	sc.pl.scatter(adata, x = 'Muc5ac', y = 'Dclk1', color = 'sampleName', use_raw = True, save = '_Muc5ac_v_Lgr5')
	sc.pl.scatter(adata, x = 'Muc5ac', y = 'Tnfrsf19', color = cluster_key, use_raw = True, save = '_Muc5ac_v_Troy')
	sc.pl.scatter(adata, x = 'Lgr5', y = 'Tnfrsf19', color = cluster_key, use_raw = True, save = '_Lgr5_v_Troy')
	sc.pl.scatter(adata, x = 'Lgr5', y = 'Bhlha15', color = cluster_key, use_raw = True, save = '_Lgr5_v_Mist1')
	sc.pl.scatter(adata, x = 'Bhlha15', y = 'Tnfrsf19', color = cluster_key, use_raw = True, save = '_Mist1_v_Troy')

def subset_by_cluster(subset_cluster_list, expt_path='', data=None, save=True, fig_pfx="", read_file=None, write_file=None):

	if read_file is None: read_file = settings.cluster_file
	if write_file is None: write_file = settings.cluster_subset_file
	cluster_path = os.path.join(expt_path,read_file)
	subset_path = os.path.join(expt_path,write_file)

	print('\nSubsetting dataset by Leiden Clustering')

	
	if data is None:
		print(f"Reading adata from {cluster_path}")
		adata = sc.read_h5ad(cluster_path)
	else: adata = data
	print("Clustering data...")
	print(adata)
	
	#print(adata.obs_names.tolist())
	
	adata = adata[adata.obs['leiden'].isin(subset_cluster_list)]
	
	#print(adata)
	
	#print(adata.obs_names.tolist())
	#print(adata.obs['tuft'])
	
	#Save subset data to a file
	
	if save:
		print(f"Saving clustered data to {subset_path}")
		adata.write(subset_path)

	print("Clustering complete.\n")
	print(adata)
	return(adata)
	
def in_silico_FACS(adata,cell_IDs,gene,cutoff,criteria,form,sort_count,gate_logic):

	#print('\nPulling Cell IDs from index for gene:', gene)
	df = pd.DataFrame(adata[:, gene].X,columns=[gene],index=adata.obs_names)
	all_cell_IDs = adata.obs_names.tolist()
	
	if sort_count == 1:
		cell_total = len(all_cell_IDs)
	else:
		cell_total = len(cell_IDs)
		
	if form == 'percent':
		cutoff = float(cutoff)
		score_cutoff = df.quantile(cutoff)
		if criteria == 'greater':
			df1 = df[df > score_cutoff]
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for',gene,': ',len(sorted_cell_IDs))
			print('Cell Total:', cell_total,'\n\n')
		if criteria == 'less':
			df1 = df[df < score_cutoff]
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for',gene,': ',len(sorted_cell_IDs))
			print('Cell Total:', cell_total,'\n\n')
	if (form == 'number' and cutoff == 'mean'):
		score_cutoff = df.mean()[0]
		if criteria == 'greater':
			df1 = df[df > score_cutoff]
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for',gene,': ',df[df < df.mean()[0]].count())
			print('Cell Total:', cell_total,'\n\n')
		elif criteria == 'less':
			df1 = df[df < score_cutoff]
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for ',gene,': ',len(sorted_cell_IDs))
			print('Cell Total:', cell_total,'\n\n')
	elif form == 'number':
		score_cutoff = float(cutoff)
		if criteria == 'greater':
			df1 = df[df > score_cutoff]
			#(df1)
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for',gene+': [',len(sorted_cell_IDs),']')
			#print('Cell Total:', cell_total,'\n\n')
		if criteria == 'less':
			df1 = df[df < score_cutoff]
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for',gene,': ',len(sorted_cell_IDs))
			print('Cell Total:', cell_total,'\n\n')
		elif criteria =='exists':
			df1 = df
			sorted_cell_IDs = df1.loc[df1[gene].notnull()].index.tolist()
			print('Cells that pass filtering threshold for',gene,': ',len(sorted_cell_IDs))
			print('Cell Total:', cell_total,'\n\n')
	elif not sorted_cell_IDs:
		sorted_cell_IDs = adata.obs_names.tolist()
		print('No sorting for ', gene, ' occurred')
	if cell_IDs:
		if gate_logic == 'OR':
			sorted_cell_IDs = set(cell_IDs).union(sorted_cell_IDs)
		if gate_logic == 'AND':
			sorted_cell_IDs = set(cell_IDs).intersection(sorted_cell_IDs)
	else:
		if gate_logic == 'AND':
			sorted_cell_IDs = set(all_cell_IDs).intersection(sorted_cell_IDs)
	
	#print('Pulling out [',len(sorted_cell_IDs),'] cells that meet sorting criteria for',gene,' \n')
	sort_count=sort_count+1 
	#print(sort_count)
	return(sorted_cell_IDs,sort_count)

def cluster_markers(expt_path='', data=None, min_in_group_fraction=None, min_fold_change=None, max_out_group_fraction=None):
	
	if data is None:
		print(f"Reading adata from {merged_path}")
		adata = sc.read_h5ad(merged_path)
	else: adata = data
	
	if min_in_group_fraction is None: min_in_group_fraction = settings.min_in_group_fraction
	if min_fold_change is None: min_fold_change = settings.min_fold_change
	if max_out_group_fraction is None: max_out_group_fraction = settings.max_out_group_fraction
	
	print(min_in_group_fraction)
	print(min_fold_change)
	print(max_out_group_fraction)
	
	sc.tl.rank_genes_groups(data, 'leiden', method='wilcoxon')
	
	sc.tl.filter_rank_genes_groups(data, groupby='leiden', log=True, key_added='rank_genes_groups_filtered', min_in_group_fraction=min_in_group_fraction, min_fold_change=min_fold_change, max_out_group_fraction=max_out_group_fraction)
	
	sc.tl.dendrogram(data, groupby='leiden')
	
	sc.pl.rank_genes_groups(data, key_added='rank_genes_groups_filtered', n_genes = 10, save = True, show = True)
	sc.pl.rank_genes_groups_dotplot(data, key='rank_genes_groups_filtered', n_genes = 10)
	
	#sc.pl.heatmap(adata, chief_genes, groupby = cluster_key, log = True, use_raw = True, show = False, standard_scale = 'var', swap_axes = True)
	#sc.pl.rank_genes_groups_heatmap(adata, chief_genes, key = 'rank_genes_groups_filtered', groupby = cluster_key, save = '_'+file_tag+'DE_genes', use_raw=True, show = False)

	#write_marker_file(expt_path, data)

def write_marker_file(expt_path='', data=None):
	print('\nParsing marker genes and saving to file. . . \n')

	file_out = expt_path+'ranked_genes.csv'
	
	if data is None:
		print(f"Reading adata from {cluster_path}")
		adata = sc.read_h5ad(cluster_path)
	else: adata = data
	print("Clustering data...")
	
	marker_dict = adata.uns['rank_genes_groups']
	unique_list = [] 
	for x in adata.obs['leiden'].values.tolist(): 
		if x not in unique_list: 
			unique_list.append(str(x))
	
	outFile = open(file_out, 'w')
	outFile.write('logFC,gene,score,pval,padj,leiden_cluster\n')
	
	parsed_dict = dict()
	
	for item in marker_dict:
		if type(marker_dict[item]) is not dict:
			cluster_data = []
			for subitem in marker_dict[item]:
				cluster_data.append(subitem.tolist())
				if str(item) not in parsed_dict:
					parsed_dict[str(item)] = cluster_data
				
	for cluster in range(0, len(unique_list)):
		for num in range(0, ranked_genes):
			line_out = []
			for marker_value in marker_dict:
				if type(marker_dict[marker_value]) is not dict:
					line_out.append(str(marker_dict[marker_value][num].tolist()[cluster]))
			line_out.append(str(cluster))
			outFile.write(','.join(line_out) + '\n')

	outFile.close()

def make_cmap(colors, position=None, bit=False):

	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit('position length must be the same as colors')
		elif position[0] != 0 or position[-1] != 1:
			sys.exit('position must start with 0 and end with 1')
	if bit:
		for i in range(len(colors)):
			colors[i] = (bit_rgb[colors[i][0]],
						 bit_rgb[colors[i][1]],
						 bit_rgb[colors[i][2]])
	cdict = {'red':[], 'green':[], 'blue':[]}
	for pos, color in zip(position, colors):
		cdict['red'].append((pos, color[0], color[0]))
		cdict['green'].append((pos, color[1], color[1]))
		cdict['blue'].append((pos, color[2], color[2]))

	cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	return cmap
	
def scRNA_velocity(expt_path='',data=None):
	
	if read_file is None: read_file = settings.cluster_file
	if write_file is None: write_file = settings.tuft_cluster_subset_file
	cluster_path = os.path.join(expt_path,read_file)
	
	
	
	
	if data is None:
		print(f"Reading adata from {cluster_path}")
		adata = sc.read_h5ad(cluster_path)
	else: adata = data
	print("Clustering data...")
	
	adata = scv.read(d+results_file, cache = True)
	ldata = scv.read(loom_file, cache = True)
	adata = scv.utils.merge(adata,ldata)
	
	scv.pp.filter_and_normalize(adata)
	scv.pp.moments(adata)
	
	scv.tl.velocity(adata, mode = 'stochastic')
	scv.tl.velocity_graph(adata)
	
	scv.pl.velocity_embedding(adata, basis='umap', save = '_'+file_tag+'_velocity_embedding')
	scv.pl.velocity_embedding_grid(adata, basis='umap', save = '_'+file_tag+'_velocity_embedding_grid')
	scv.pl.velocity_embedding_stream(adata, basis='umap', save = '_'+file_tag+'_velocity_embedding_stream')