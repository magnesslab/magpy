import sys
import bbknn
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pathlib import Path
#import mjc_functions as mjc
import matplotlib.pyplot as plt
#import scanorama
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 3			 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.settings.set_figure_params(dpi_save=1200, dpi=1200)
sns.set(style="white", color_codes=True)

greatestPalette = ["#f44336","#265bcf","#36b3ba","#ffeb3b","#e91e63","#00cc92","#4caf50","#ffb65c","#9c27b0","#03a9f4","#43d61f","#ff9800","#673ab7","#cddc39","#81bdc5","#ff5722"]
palestPalette = ["#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb","#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff","#edf3ba","#ffccbd"]

#############################################################################
##         Flag to toggle rerunning the analysis or just plotting          ##
#############################################################################

rerun = False

#############################################################################
##                          Plot section flags                             ##
#############################################################################

do_ec_analysis = False
plot_raw = False
plot_normalized = False

redraw_featureplots = True
redraw_umaps = True

run_marker_analysis = False

expression_cutoff = 0.01  # Sets the threshold for min expression to be shown in plots

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
mistorage_mount_point = '/mnt/black/scRNA-seq/'

#############################################################################
## Change this to list the sampleIDs of the samples you want to see        ##
#############################################################################
sample_list = ['m-007']


#############################################################################
## Change this to contain all genes you want to see expression for         ##
#############################################################################
genes_of_interest = ['EPCAM','VIM','CDX2','CHGA','KRAS','EGFR','ERBB2','ERBB3','ERBB4','PTPRC','HES1',
'PDX1','SMO','PAX6','WNT5A','DLL1','GATA6','RFX3','NKX6-1','NEUROD1',
'CDK6','NKX2-2','RFX6','HNF4A','INSM1','BMP4','BMP6','BMP5','MEN1','PDPK1','SIDT2']

epi_cell_type_genes = ['EPCAM','CDX2','LGR5', 'OLFM4', 'CHGA', 'LYZ', 'MUC2', 'MUC5', 'VIL1', 'DPP4', 'FABP1', 'ALPI', 'SI', 'MGAM', 'SPDEF', 'SHH', 'IHH', 'PDGFA', 'SOX9','AXIN2','LGR5','ASCL2','CMYC','CD44','CYCLIND','TCF1','LEF1','PPARD','CJUN','MMP7','MSI1','LRP5','LRP6']

marker_genes = ['RGS5','APOD','UGT2B7','GNG2','SLC17A4','PITX2','TRPA1','PIK3C2G','ADH6',
'ECM1','MAP3K20','UACA','RGS5','LGALS3','PIGR','TSPAN8','S100A14','FAT1','ANXA2','PTPRC',
'SRGN','MYL12A','ARHGDIB','CD52','ETS1','CTSS','TYROBP','COTL1','AIF1','LYZ','TYMP','CD79A',
'MEF2C','IGKC','CYBA','IGHA1','VIM','MUC2','CLCA1','FCGBP','ST6GALNAC1','ITLN1','SPINK4']

#############################################################################
##        Adjust these parameters to get a nicer looking UMAP plot         ##
#############################################################################
# UMAP arguments
num_neighbors_use = 50
num_pcs_use = 8
umap_spread = 1
umap_min_dist = 0.3
maxiter = None
umap_gamma=1

dot_size = 15

# Louvain arguments
louv_res = 0.6

# PAGA arguments
size=20
paga_layout='fr'
threshold=0.005
node_size_scale=3

#############################################################################
##               DO NOT CHANGE ANYTHING BEYOND THIS POINT                  ##
#############################################################################




## Location to output the anndata h5ad files
raw_data_file = ''.join(['./data/Data_', '_'.join(sample_list), '.scanpy.raw.h5ad'])  # the file that will store the raw combined data
results_file = ''.join(['./data/Data_', '_'.join(sample_list), '.processed.h5ad'])  # the file that will store the analysis results
filtered_data_file = ''.join(['./data/Data_', '_'.join(sample_list), '.scanpy.filtered.h5ad'])  # the file that will store the raw combined data
endo_results_file = ''.join(['./data/Data_', '_'.join(sample_list), '.endothelial.processed.h5ad'])  # the file that will store the analysis results
neuro_results_file = ''.join(['./data/Data_', '_'.join(sample_list), '.neuronal.processed.h5ad'])  # the file that will store the analysis results

## Define function to generate a color gradient from a defined starting and ending color
def make_cmap(colors, position=None, bit=False):
	'''
	make_cmap takes a list of tuples which contain RGB values. The RGB
	values may either be in 8-bit [0 to 255] (in which bit must be set to
	True when called) or arithmetic [0 to 1] (default). make_cmap returns
	a cmap with equally spaced colors.
	Arrange your tuples so that the first color is the lowest value for the
	colorbar and the last is the highest.
	position contains values from 0 to 1 to dictate the location of each color.
	'''
	import matplotlib as mpl
	import numpy as np
	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit("position length must be the same as colors")
		elif position[0] != 0 or position[-1] != 1:
			sys.exit("position must start with 0 and end with 1")
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

annotation_dict = dict()

for line in open(''.join([mistorage_mount_point, 'single_cell_meta_data_table.tsv']), 'r'):
	#print(line)
	elem = str.split(line.rstrip())
	#print(elem)
	if elem:
		if elem[0] not in annotation_dict:
			annotation_dict[elem[0]] = elem[1:]

def Create_Scanpy_Anndata(mistorage_mount_point, sampleID):
	metadata_list = annotation_dict[sampleID][1:]
	newAdata = sc.read_10x_h5(''.join([mistorage_mount_point, annotation_dict[sampleID][0]]), gex_only=False)
	## Set gene names to be unique since there seem to be duplicate names from Cellranger
	newAdata.var_names_make_unique()
	## Add metadata for each sample to the observation (cells) annotations in the Anndata objects
	print('\nAdding Metadata for sample',sampleID,'\n')
	for field in metadata_list:
		field_list = str.split(field, ':')
		meta_name = field_list[0]
		meta_value = field_list[1]
		newAdata.obs[meta_name] = meta_value
	return(newAdata)

# function to get unique values 
def unique(list1): 
  
    # intilize a null list 
    unique_list = [] 
      
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return(unique_list)

## Pasring function to get marker gene data similar to Seurat's marker output
def write_marker_file(adata, file_out='markers_output.csv', n_genes=100):
	print('Parsing markers...')
	marker_dict = adata.uns['rank_genes_groups']
	unique_list = [] 
	for x in adata.obs['louvain'].values.tolist(): 
		if x not in unique_list: 
			unique_list.append(str(x))
	
	outFile = open(file_out, 'w')
	outFile.write('logFC,gene,score,pval,padj,cluster\n')
	
	#i = 0
	
	parsed_dict = dict()
	
	for item in marker_dict:
		if type(marker_dict[item]) is not dict:
			cluster_data = []
			for subitem in marker_dict[item]:
				cluster_data.append(subitem.tolist())
			
			if str(item) not in parsed_dict:
				parsed_dict[str(item)] = cluster_data
	for cluster in range(0, len(unique_list)):
		for i in range(0, n_genes):
			line_out = []
			for marker_value in marker_dict:
				if type(marker_dict[marker_value]) is not dict:
					line_out.append(str(marker_dict[marker_value][i].tolist()[cluster]))
			line_out.append(str(cluster))
			outFile.write(','.join(line_out) + '\n')
	
	print('Saving marker data to:', file_out)
	outFile.close()

## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
#feature_colors = [(230,230,230), (35,35,142), (255,127,0)]
#position=[0, expression_cutoff, 1]
#my_feature_cmap = make_cmap(feature_colors, position=position, bit=True)
#dot_colors = [(230,230,230), (153,0,0), (255,145,0)]
#my_dot_cmap = make_cmap(dot_colors, position=position, bit=True)

#feature_colors = [(210,210,210), (210,210,210), (0,51,102), (255,141,41)]
#position=[0, 0.019999, 0.02, 1]
#my_feature_cmap = make_cmap(feature_colors, position=position, bit=True)
#dot_position=[0, 0.019999, 0.02, 1]
#dot_colors = [(210,210,210), (210,210,210), (0,51,102), (255,141,41)]
#my_dot_cmap = make_cmap(dot_colors, position=dot_position, bit=True)

feature_colors = [(210,210,210), (210,210,210), (230,230,250), (218,112,214), (75,0,130)]
position=[0, 0.019999, 0.02, 0.55, 1]
my_feature_cmap = make_cmap(feature_colors, position=position, bit=True)
dot_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
my_dot_cmap = make_cmap(dot_colors, position=position, bit=True)

hypox_colors = [(210,210,210), (210,210,210), (54,38,144), (54,38,144), (54,38,144)]
hypox_position=[0, 0.000001, 0.0000011,0.00000111, 1]
hypox_cmap = make_cmap(hypox_colors, position=hypox_position, bit=True)

grey_colors = [(210,210,210), (210,210,210)]
grey_position=[0, 1]
grey_cmap = make_cmap(grey_colors, position=grey_position, bit=True)




## Read the raw Cellranger filtered data matrices into new Anndata objects
def runBasicAnalysis():
	
	first_adata = False
	
	if Path(raw_data_file).is_file():
		print(''.join(['Data_', '_'.join(sample_list), '.scanpy.raw.h5ad']), 'found, using this existing raw data file\n')
		adata = sc.read_h5ad(raw_data_file)
	else:
		print('\nNo existing h5ad raw data file found, reading in 10x h5 data for each sample\n')
		for sample in sample_list:
			if not first_adata:
				adata = Create_Scanpy_Anndata(mistorage_mount_point, sample)
				first_adata = True
			else:
				adata = adata.concatenate(Create_Scanpy_Anndata(mistorage_mount_point, sample))
		
		## Make cell names unique by adding _1, _2, _3 sequentially to each duplicated 10x barcode/name
		adata.obs_names_make_unique()
		## Write the raw combined dataset to disk so you can skip combining samples next time
		print('\nSaving raw combined sample data to', raw_data_file, '\n')
		adata.write(raw_data_file)
		open('./data/Prefiltered_gene_list.txt', 'w').write('\n'.join(adata.var_names.tolist()))
	
	## Basic filtering to get rid of useless cells and unexpressed genes
	
	sc.pp.filter_cells(adata, min_genes=1000)
	sc.pp.filter_genes(adata, min_cells=5)
	
	print('\nDoing initial filtering...\nKeeping', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
	mito_genes = adata.var_names.str.startswith('MT-')
	# Calculate the percent of genes derived from mito vs genome
	# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
	adata.obs['percent_mito'] = np.sum(
		adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	# add the total counts per cell as observations-annotation to adata
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1
	
	cell_cycle_genes = [x.strip() for x in open('/mnt/black/scRNA-seq/regev_lab_cell_cycle_genes.txt')]
	
	s_genes = cell_cycle_genes[:43]
	g2m_genes = cell_cycle_genes[43:]
	cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
	
	sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
	
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
				 jitter=0.4, multi_panel=True, save = '_preFiltering_plot.pdf', show = False)
	
	## Actually do the filtering.
	
	adata = adata[adata.obs['n_genes'] > 1000, :]   # Keep cells with more than 1000 genes
	adata = adata[adata.obs['n_genes'] < 10000, :]   # Keep cells with less than 5000 genes to remove most doublets
	adata = adata[adata.obs['n_counts'] < 80000, :] # Keep cells with less than 15000 UMIs to catch a few remaining doublets
	adata = adata[adata.obs['percent_mito'] < 0.8, :]   # Keep cells with less than 0.1 mito/genomic gene ratio
	sc.pp.filter_genes(adata, min_cells=5)	# Refilter genes to get rid of genes that are only in a tiny number of cells
	
	print('\nDoing final filtering...\nKeeping', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
	open('./data/Final_filtered_gene_list.txt', 'w').write('\n'.join(adata.var_names.tolist()))
	
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
				 jitter=0.4, multi_panel=True, save = '_postFiltering_plot.pdf', show = False)
	
	## Normalize the expression matrix to 10,000 reads per cell, so that counts become comparable among cells.
	# This corrects for differences in sequencing depth between cells and samples
	
	#sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.normalize_total(adata)
	
	## Log transform the data.
	
	sc.pp.log1p(adata)
	
	## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
	# We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
	
	adata.write(filtered_data_file)
	adata.raw = adata
	
	## Identify highly-variable genes based on dispersion relative to expression level.
	
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=6, min_disp=0.2)
	
	## Filter the genes to remove non-variable genes since they are uninformative
	
	adata = adata[:, adata.var['highly_variable']]
	
	## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
	
	sc.pp.regress_out(adata, ['n_counts'])
	
	## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
	
	sc.pp.scale(adata, max_value=10)
	
	## Run PCA to compute the default number of components
	
	sc.tl.pca(adata, svd_solver='arpack')
	
	## Rank genes according to contributions to PCs.
	
	sc.pl.pca_loadings(adata, show=False, components=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], save='_PCA-loadings.pdf')
	
	## Draw the PCA elbow plot to determine which PCs to use
	sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 100, save = '_elbowPlot.pdf', show = False)
	
	## Compute nearest-neighbors
	
	sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, n_pcs=num_pcs_use)
	
	## fix batch differences based on XX/XY
	#bbknn.bbknn(adata, batch_key='sampleName', n_pcs=50, neighbors_within_batch=3, copy=False)
	
	## Calculate cell clusters via Louvain algorithm
	
	sc.tl.louvain(adata, resolution = louv_res)
	
	## Run UMAP Dim reduction
	
	sc.tl.umap(adata, min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma)
	
	## Save the final clustered and dim-reduced data file as h5ad
	print('\nSaving processed Anndata data to', results_file, '\n')
	adata.write(results_file)
	return(adata)


if rerun:
	adata = runBasicAnalysis()
	adata.raw = sc.read_h5ad(filtered_data_file)
	print('Rerunning filtering and normalization...')
elif not rerun:
	print('Rerun toggled to "off"...\nMoving directly to subclustering and plotting\n')
	
	figure_dir = './figures/sample_1'
	sc.settings.figdir = figure_dir
	results_adata = sc.read_h5ad(results_file)
	filtered_adata = sc.read_h5ad(filtered_data_file)
	
	filtered_adata.obs['louvain'] = results_adata.obs['louvain']
	
	filtered_adata[filtered_adata.obs['louvain'].isin(['5', '1', '3'])].write('./data/sample1.anndata.h5ad')
	
	adata = sc.read_h5ad('./data/sample1.anndata.h5ad')
	adata.raw = sc.read_h5ad('./data/sample1.anndata.h5ad')
	
	df = adata.to_df()
	
	hypox_gene_list = ['LGR5','OLFM4','ASCL2','SMOC2','MKI67','PCNA','MUC2','LYZ','ALPI','HOPX']
	
	for gene in hypox_gene_list:
		df1 = df[gene]
		print(gene, df1[df1 > 0].count())
	
	'''
	df1 = df['LGR5']
	print('LGR5', df1[df1 > 0].count())
	
	df1 = df['OLFM4']
	print('OLFM4', df1[df1 > 0].count())
	
	df1 = df['ASCL2']
	print('ASCL2', df1[df1 > 0].count())
	
	df1 = df['SMOC2']
	print('SMOC2', df1[df1 > 0].count())
	
	df1 = df['MKI67']
	print('MKI67', df1[df1 > 0].count())
	
	df1 = df['PCNA']
	print('PCNA', df1[df1 > 0].count())
	'''
	
	print('Total', df['EPCAM'].count())
	
	expressed_dict = dict()
	
	for gene in adata.var_names.values.tolist():
		if gene not in expressed_dict:
			expressed_dict[str(gene)] = 1
			#print(gene)
	
	il_stem_cell = ['LGR5','OLFM4','ASCL2','IL10RA','IL10RB','IL11RA','IL12RB1','IL12RB2','IL13RA1','IL13RA2','IL15RA','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL17REL','IL18R1','IL1R1','IL1R2','IL1RAP','IL1RAPL1','IL1RAPL2','IL1RL1','IL1RL2','IL1RN','IL20RA','IL20RB','IL21R','IL22RA1','IL22RA2','IL23R','IL27RA','IL2RA','IL2RB','IL2RG','IL31RA','IL36RN','IL3RA','IL4R','IL5RA','IL6R','IL6RP1','IL7R','IL9R','IL9RP2','IL9RP3','IL9RP4']
	
	genes_to_plot = []
	for gene in il_stem_cell:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	df2 = pd.concat([df['LGR5'], df['OLFM4'],df['ASCL2'],df['SMOC2'],df['MKI67'],df['PCNA'],df['IL10RA'],df['IL10RB'],df['IL11RA'],df['IL12RB2'],df['IL13RA1'],df['IL15RA'],df['IL17RA'],df['IL17RB'],df['IL17RC'],df['IL17RD'],df['IL17RE'],df['IL18R1'],df['IL1R1'],df['IL1R2'],df['IL1RAP'],df['IL1RL2'],df['IL1RN'],df['IL20RA'],df['IL20RB'],df['IL22RA1'],df['IL23R'],df['IL27RA'],df['IL2RG'],df['IL3RA'],df['IL4R'],df['IL6R']], axis=1)
	
	print(df2)
	
	df2.to_csv(path_or_buf='MarkerGenes_dataframe.csv', sep=',')
	
	print('LGR5 and OLFM4', df2[df2['LGR5'] > 0][df2['OLFM4'] > 0].count())
	print('OLFM4 and PCNA', df2[df2['OLFM4'] > 0][df2['PCNA'] > 0].count())
	print('LGR5 and PCNA', df2[df2['LGR5'] > 0][df2['PCNA'] > 0].count())
	print('LGR5 and MKI67', df2[df2['LGR5'] > 0][df2['MKI67'] > 0].count())
	print('PCNA or MKI67', df2[(df2.PCNA > 0) | (df2.MKI67 > 0)].count())
	print('LGR5 or OLFM4 or ASCL2', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)].count())
	
	print('---------------------------------\n')
	
	print('IL10RA', df2[df2['IL10RA'] > 0].count())
	print('IL10RB', df2[df2['IL10RB'] > 0].count())
	print('IL11RA', df2[df2['IL11RA'] > 0].count())
	print('IL12RB2', df2[df2['IL12RB2'] > 0].count())
	print('IL13RA1', df2[df2['IL13RA1'] > 0].count())
	print('IL15RA', df2[df2['IL15RA'] > 0].count())
	print('IL17RA', df2[df2['IL17RA'] > 0].count())
	print('IL17RB', df2[df2['IL17RB'] > 0].count())
	print('IL17RC', df2[df2['IL17RC'] > 0].count())
	print('IL17RD', df2[df2['IL17RD'] > 0].count())
	print('IL17RE', df2[df2['IL17RE'] > 0].count())
	print('IL18R1', df2[df2['IL18R1'] > 0].count())
	print('IL1R1', df2[df2['IL1R1'] > 0].count())
	print('IL1R2', df2[df2['IL1R2'] > 0].count())
	print('IL1RAP', df2[df2['IL1RAP'] > 0].count())
	print('IL1RL2', df2[df2['IL1RL2'] > 0].count())
	print('IL1RN', df2[df2['IL1RN'] > 0].count())
	print('IL20RA', df2[df2['IL20RA'] > 0].count())
	print('IL20RB', df2[df2['IL20RB'] > 0].count())
	print('IL22RA1', df2[df2['IL22RA1'] > 0].count())
	print('IL23R', df2[df2['IL23R'] > 0].count())
	print('IL27RA', df2[df2['IL27RA'] > 0].count())
	print('IL2RG', df2[df2['IL2RG'] > 0].count())
	print('IL3RA', df2[df2['IL3RA'] > 0].count())
	print('IL4R', df2[df2['IL4R'] > 0].count())
	print('IL6R', df2[df2['IL6R'] > 0].count())
	
	print('---------------------------------\n')
	
	print('ASCL2 and IL10RA', df2[df2['ASCL2'] > 0][df2['IL10RA'] > 0].count())
	print('ASCL2 and IL10RB', df2[df2['ASCL2'] > 0][df2['IL10RB'] > 0].count())
	print('ASCL2 and IL11RA', df2[df2['ASCL2'] > 0][df2['IL11RA'] > 0].count())
	print('ASCL2 and IL12RB2', df2[df2['ASCL2'] > 0][df2['IL12RB2'] > 0].count())
	print('ASCL2 and IL13RA1', df2[df2['ASCL2'] > 0][df2['IL13RA1'] > 0].count())
	print('ASCL2 and IL15RA', df2[df2['ASCL2'] > 0][df2['IL15RA'] > 0].count())
	print('ASCL2 and IL17RA', df2[df2['ASCL2'] > 0][df2['IL17RA'] > 0].count())
	print('ASCL2 and IL17RB', df2[df2['ASCL2'] > 0][df2['IL17RB'] > 0].count())
	print('ASCL2 and IL17RC', df2[df2['ASCL2'] > 0][df2['IL17RC'] > 0].count())
	print('ASCL2 and IL17RD', df2[df2['ASCL2'] > 0][df2['IL17RD'] > 0].count())
	print('ASCL2 and IL17RE', df2[df2['ASCL2'] > 0][df2['IL17RE'] > 0].count())
	print('ASCL2 and IL18R1', df2[df2['ASCL2'] > 0][df2['IL18R1'] > 0].count())
	print('ASCL2 and IL1R1', df2[df2['ASCL2'] > 0][df2['IL1R1'] > 0].count())
	print('ASCL2 and IL1R2', df2[df2['ASCL2'] > 0][df2['IL1R2'] > 0].count())
	print('ASCL2 and IL1RAP', df2[df2['ASCL2'] > 0][df2['IL1RAP'] > 0].count())
	print('ASCL2 and IL1RL2', df2[df2['ASCL2'] > 0][df2['IL1RL2'] > 0].count())
	print('ASCL2 and IL1RN', df2[df2['ASCL2'] > 0][df2['IL1RN'] > 0].count())
	print('ASCL2 and IL20RA', df2[df2['ASCL2'] > 0][df2['IL20RA'] > 0].count())
	print('ASCL2 and IL20RB', df2[df2['ASCL2'] > 0][df2['IL20RB'] > 0].count())
	print('ASCL2 and IL22RA1', df2[df2['ASCL2'] > 0][df2['IL22RA1'] > 0].count())
	print('ASCL2 and IL23R', df2[df2['ASCL2'] > 0][df2['IL23R'] > 0].count())
	print('ASCL2 and IL27RA', df2[df2['ASCL2'] > 0][df2['IL27RA'] > 0].count())
	print('ASCL2 and IL2RG', df2[df2['ASCL2'] > 0][df2['IL2RG'] > 0].count())
	print('ASCL2 and IL3RA', df2[df2['ASCL2'] > 0][df2['IL3RA'] > 0].count())
	print('ASCL2 and IL4R', df2[df2['ASCL2'] > 0][df2['IL4R'] > 0].count())
	print('ASCL2 and IL6R', df2[df2['ASCL2'] > 0][df2['IL6R'] > 0].count())
	
	print('---------------------------------\n')
	
	print('OLFM4 and IL10RA', df2[df2['OLFM4'] > 0][df2['IL10RA'] > 0].count())
	print('OLFM4 and IL10RB', df2[df2['OLFM4'] > 0][df2['IL10RB'] > 0].count())
	print('OLFM4 and IL11RA', df2[df2['OLFM4'] > 0][df2['IL11RA'] > 0].count())
	print('OLFM4 and IL12RB2', df2[df2['OLFM4'] > 0][df2['IL12RB2'] > 0].count())
	print('OLFM4 and IL13RA1', df2[df2['OLFM4'] > 0][df2['IL13RA1'] > 0].count())
	print('OLFM4 and IL15RA', df2[df2['OLFM4'] > 0][df2['IL15RA'] > 0].count())
	print('OLFM4 and IL17RA', df2[df2['OLFM4'] > 0][df2['IL17RA'] > 0].count())
	print('OLFM4 and IL17RB', df2[df2['OLFM4'] > 0][df2['IL17RB'] > 0].count())
	print('OLFM4 and IL17RC', df2[df2['OLFM4'] > 0][df2['IL17RC'] > 0].count())
	print('OLFM4 and IL17RD', df2[df2['OLFM4'] > 0][df2['IL17RD'] > 0].count())
	print('OLFM4 and IL17RE', df2[df2['OLFM4'] > 0][df2['IL17RE'] > 0].count())
	print('OLFM4 and IL18R1', df2[df2['OLFM4'] > 0][df2['IL18R1'] > 0].count())
	print('OLFM4 and IL1R1', df2[df2['OLFM4'] > 0][df2['IL1R1'] > 0].count())
	print('OLFM4 and IL1R2', df2[df2['OLFM4'] > 0][df2['IL1R2'] > 0].count())
	print('OLFM4 and IL1RAP', df2[df2['OLFM4'] > 0][df2['IL1RAP'] > 0].count())
	print('OLFM4 and IL1RL2', df2[df2['OLFM4'] > 0][df2['IL1RL2'] > 0].count())
	print('OLFM4 and IL1RN', df2[df2['OLFM4'] > 0][df2['IL1RN'] > 0].count())
	print('OLFM4 and IL20RA', df2[df2['OLFM4'] > 0][df2['IL20RA'] > 0].count())
	print('OLFM4 and IL20RB', df2[df2['OLFM4'] > 0][df2['IL20RB'] > 0].count())
	print('OLFM4 and IL22RA1', df2[df2['OLFM4'] > 0][df2['IL22RA1'] > 0].count())
	print('OLFM4 and IL23R', df2[df2['OLFM4'] > 0][df2['IL23R'] > 0].count())
	print('OLFM4 and IL27RA', df2[df2['OLFM4'] > 0][df2['IL27RA'] > 0].count())
	print('OLFM4 and IL2RG', df2[df2['OLFM4'] > 0][df2['IL2RG'] > 0].count())
	print('OLFM4 and IL3RA', df2[df2['OLFM4'] > 0][df2['IL3RA'] > 0].count())
	print('OLFM4 and IL4R', df2[df2['OLFM4'] > 0][df2['IL4R'] > 0].count())
	print('OLFM4 and IL6R', df2[df2['OLFM4'] > 0][df2['IL6R'] > 0].count())
	
	print('---------------------------------\n')
	
	print('LGR5 and IL10RA', df2[df2['LGR5'] > 0][df2['IL10RA'] > 0].count())
	print('LGR5 and IL10RB', df2[df2['LGR5'] > 0][df2['IL10RB'] > 0].count())
	print('LGR5 and IL11RA', df2[df2['LGR5'] > 0][df2['IL11RA'] > 0].count())
	print('LGR5 and IL12RB2', df2[df2['LGR5'] > 0][df2['IL12RB2'] > 0].count())
	print('LGR5 and IL13RA1', df2[df2['LGR5'] > 0][df2['IL13RA1'] > 0].count())
	print('LGR5 and IL15RA', df2[df2['LGR5'] > 0][df2['IL15RA'] > 0].count())
	print('LGR5 and IL17RA', df2[df2['LGR5'] > 0][df2['IL17RA'] > 0].count())
	print('LGR5 and IL17RB', df2[df2['LGR5'] > 0][df2['IL17RB'] > 0].count())
	print('LGR5 and IL17RC', df2[df2['LGR5'] > 0][df2['IL17RC'] > 0].count())
	print('LGR5 and IL17RD', df2[df2['LGR5'] > 0][df2['IL17RD'] > 0].count())
	print('LGR5 and IL17RE', df2[df2['LGR5'] > 0][df2['IL17RE'] > 0].count())
	print('LGR5 and IL18R1', df2[df2['LGR5'] > 0][df2['IL18R1'] > 0].count())
	print('LGR5 and IL1R1', df2[df2['LGR5'] > 0][df2['IL1R1'] > 0].count())
	print('LGR5 and IL1R2', df2[df2['LGR5'] > 0][df2['IL1R2'] > 0].count())
	print('LGR5 and IL1RAP', df2[df2['LGR5'] > 0][df2['IL1RAP'] > 0].count())
	print('LGR5 and IL1RL2', df2[df2['LGR5'] > 0][df2['IL1RL2'] > 0].count())
	print('LGR5 and IL1RN', df2[df2['LGR5'] > 0][df2['IL1RN'] > 0].count())
	print('LGR5 and IL20RA', df2[df2['LGR5'] > 0][df2['IL20RA'] > 0].count())
	print('LGR5 and IL20RB', df2[df2['LGR5'] > 0][df2['IL20RB'] > 0].count())
	print('LGR5 and IL22RA1', df2[df2['LGR5'] > 0][df2['IL22RA1'] > 0].count())
	print('LGR5 and IL23R', df2[df2['LGR5'] > 0][df2['IL23R'] > 0].count())
	print('LGR5 and IL27RA', df2[df2['LGR5'] > 0][df2['IL27RA'] > 0].count())
	print('LGR5 and IL2RG', df2[df2['LGR5'] > 0][df2['IL2RG'] > 0].count())
	print('LGR5 and IL3RA', df2[df2['LGR5'] > 0][df2['IL3RA'] > 0].count())
	print('LGR5 and IL4R', df2[df2['LGR5'] > 0][df2['IL4R'] > 0].count())
	print('LGR5 and IL6R', df2[df2['LGR5'] > 0][df2['IL6R'] > 0].count())
	
	print('---------------------------------\n')
	
	print('LGR5 or OLFM4 or ASCL2 and IL10RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL10RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL10RB', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL10RB'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL11RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL11RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL12RB2', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL12RB2'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL13RA1', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL13RA1'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL15RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL15RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL17RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL17RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL17RB', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL17RB'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL17RC', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL17RC'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL17RD', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL17RD'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL17RE', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL17RE'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL18R1', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL18R1'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL1R1', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL1R1'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL1R2', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL1R2'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL1RAP', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL1RAP'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL1RL2', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL1RL2'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL1RN', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL1RN'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL20RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL20RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL20RB', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL20RB'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL22RA1', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL22RA1'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL23R', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL23R'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL27RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL27RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL2RG', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL2RG'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL3RA', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL3RA'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL4R', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL4R'] > 0].count())
	print('LGR5 or OLFM4 or ASCL2 and IL6R', df2[(df2.LGR5 > 0) | (df2.OLFM4 > 0) | (df2.ASCL2 > 0)][df2['IL6R'] > 0].count())
	
	print('---------------------------------\n')
	
	#sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.normalize_total(adata)
	
	## Log transform the data.
	
	sc.pp.log1p(adata)
	adata.raw = adata
	
	## Identify highly-variable genes based on dispersion relative to expression level.
	
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=6, min_disp=0.2)
	
	## Filter the genes to remove non-variable genes since they are uninformative
	
	adata = adata[:, adata.var['highly_variable']]
	
	## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
	
	sc.pp.regress_out(adata, ['n_counts'])
	
	## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
	sc.pp.scale(adata, max_value=10)
	
	## Run PCA to compute the default number of components
	sc.tl.pca(adata, svd_solver='arpack')
	
	## Rank genes according to contributions to PCs.
	sc.pl.pca_loadings(adata, show=False, components=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], save='_PCA-loadings.pdf')
	
	## Draw the PCA elbow plot to determine which PCs to use
	sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 100, save = '_elbowPlot.pdf', show = False)
	
	## Compute nearest-neighbors
	sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, n_pcs=num_pcs_use)
	
	## fix batch differences based on XX/XY
	#bbknn.bbknn(adata, batch_key='sampleName', n_pcs=50, neighbors_within_batch=3, copy=False)
	
	## Calculate cell clusters via Louvain algorithm
	sc.tl.louvain(adata, resolution = louv_res)
	
	## Run UMAP Dim reduction
	sc.tl.umap(adata, min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma)
	
	sc.pl.umap(adata, color='phase', save = '_cell_cycle.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color=hypox_gene_list, save = '_figure7_features.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size*3, cmap = hypox_cmap, alpha = 0.85)
	sc.pl.umap(adata, color=hypox_gene_list, save = '_figure7-grey_features.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size*3, cmap = grey_cmap, alpha = 0.85)
	sc.pl.umap(adata, color=['ACE2'], save = '_ACE2.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size*3, cmap = hypox_cmap, alpha = 0.85)
	sc.pl.umap(adata, color=['ACE2'], save = '_ACE2_binary.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size*3, cmap = my_feature_cmap, alpha = 0.85)
	
	expressed_dict = dict()
	
	for gene in adata.raw.var_names.values.tolist():
		if gene not in expressed_dict:
			expressed_dict[str(gene)] = 1
			#print(gene)
	
	genes_to_plot = []
	jbs_genes = ['Hashtag_1','Hashtag_2','Hashtag_3','Hashtag_4','Hashtag_5','Hashtag_6','MKI67','PCNA','LGR5','OLFM4','HOPX','MSI1','TERT','SOX9','RB1','RBL2','CCND1','CDK3','CDK4','FOXO3','EZH1','CYCS','KRT19','KRT20','CDKN2A','CDKN1B','H3C8','CYP2F1','GHR','TCP11L2','H2BC4','OGN','GBP2','PCDHA1','VNN1','RAB27B','SAA3P','H2BC21','PTGFR','DDX58','H2AC18','H2BC21','KDM7A','KDM5B','KDM6B','JARID2','PDCD4','SELENBP1','TOB1','IDH1','DDIT4','NFE2L2','YPEL2','ABCC3','CDKN1A']
	
	for gene in jbs_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_JBsGenes_featureplots_sample1.pdf', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	sc.pl.umap(adata, color=genes_to_plot, save = '_JBsGenes_featureplots_sample1.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	
	rISC = ['HOPX', 'BMI1', 'LRIG1']
	Paneth = ['LYZ', 'DEFA4', 'DEFA5', 'DEFA6', 'DLL1', 'DLL4', 'WNT3', 'MMP7']
	Goblet = ['MUC1', 'MUC2', 'MUC4', 'MUC13', 'MUC5B', 'TFF3', 'CLCA1']
	EE = ['CHGB', 'SST', 'CCK', 'TPH']
	Tuft = ['POU2F3', 'DCLK1', 'PTGS1', 'GNAT3', 'CHAT', 'GFI1B', 'TRPM5']
	AE = ['SI', 'ALPI']
	carbohydrate_metabolism = ['FUCA1', 'GLO1', 'HUM2DD', 'ALDOB', 'PFKP', 'ACLY']
	fatty_acid_metabolism = ['APOA1', 'APOA4', 'APOBEC1', 'ADIPOR2', 'FACL5']
	protein_metabolism = ['DPP4', 'RAB17']
	additional_genes = ['MKI67','PCNA','SOX4','SOX9','SOX17','TMPRSS2','LGR5','ASCL2','OLFM4','SMOC2','BMI1','HOPX','LRIG1','LYZ','DEFA5','DLL1','DLL4','WNT3','MMP7','MUC2','MUC4','MUC13','MUC5B','CHGA','CHGB','SST','TPH','CCK','POU2F3','TRMP5','SI','APOA1','APOA3','DPP4','RAB17']
	il_stem_cell = ['LGR5','OLFM4','ASCL2','IL10RA','IL10RB','IL11RA','IL12RB1','IL12RB2','IL13RA1','IL13RA2','IL15RA','IL17RA','IL17RB','IL17RC','IL17RD','IL17RE','IL17REL','IL18R1','IL1R1','IL1R2','IL1RAP','IL1RAPL1','IL1RAPL2','IL1RL1','IL1RL2','IL1RN','IL20RA','IL20RB','IL21R','IL22RA1','IL22RA2','IL23R','IL27RA','IL2RA','IL2RB','IL2RG','IL31RA','IL36RN','IL3RA','IL4R','IL5RA','IL6R','IL6RP1','IL7R','IL9R','IL9RP2','IL9RP3','IL9RP4']
	extra_genes = ['ALPI','APOA4','FUCA1','ALDOB','DCLK1','KRT20','RAB4']
	
	genes_to_plot = []
	for gene in rISC:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_rSIC_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in Paneth:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_Paneth_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in Goblet:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_Goblet_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in EE:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_EE_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in Tuft:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_Tuft_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in AE:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_AE_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in carbohydrate_metabolism:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_carbohydrate_metabolism_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in fatty_acid_metabolism:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_fatty_acid_metabolism_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in protein_metabolism:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_protein_metabolism_sample1.pdf', show = False, cmap = hypox_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in additional_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_additional_genes_sample1.pdf', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	sc.pl.umap(adata, color=genes_to_plot, save = '_additional_genes_sample1_gray.pdf', show = False, cmap = grey_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in il_stem_cell:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_il_stem_cell_sample1.pdf', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	for gene in extra_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	sc.pl.umap(adata, color=genes_to_plot, save = '_extra_genes_sample1.pdf', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)

	
	
	
	

orig_adata = sc.read_h5ad(raw_data_file)








if redraw_umaps:
	adata = sc.read_h5ad(results_file)
	figure_dir = './figures'
	sc.settings.figdir = figure_dir
	
	print('\nRedrawing the umap plots...\n---------------------------\n')
	sc.pl.umap(adata, color='louvain', save = '_clusterIdentity.pdf', show = False, legend_loc = 'on data', edges = True, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color='louvain', save = '_clusterIdentity_noEdge.pdf', show = False, legend_loc = 'on data', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	#sc.pl.umap(adata, color=['louvain', 'age'], save = '_clusterIdentity_age.pdf', show = False, legend_loc = 'right margin', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	#sc.pl.umap(adata, color='age', save = '_age.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	#sc.pl.umap(adata, color='tissue', save = '_tissue.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	#sc.pl.umap(adata, color='sex', save = '_sex.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	#sc.pl.umap(adata, color='gel', save = '_hydrogel.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='phase', save = '_cell_cycle.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='sampleName', save = '_sample.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color=['Hashtag_1','Hashtag_2','Hashtag_3','Hashtag_4','Hashtag_5','Hashtag_6'], save = '_antibody_features.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, cmap = my_feature_cmap, alpha = 0.95)
	sc.pl.umap(adata, color=['Hashtag_1','Hashtag_2','Hashtag_3','Hashtag_4','Hashtag_5','Hashtag_6','ACE2'], save = '_antibody_features.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, cmap = my_feature_cmap, alpha = 0.95)
	
	## Run PAGA to get nicer cell starting positions for UMAP
	sc.tl.paga(adata, groups='louvain')
	sc.pl.paga(adata, color='louvain', save='_paga-init-pattern.pdf', show=False, threshold=threshold, node_size_scale=node_size_scale, node_size_power=0.9, layout=paga_layout)
	
	adata.write(results_file)







if redraw_featureplots:
	print('\nRedrawing the feature plots...\n---------------------------\n')
	
	adata = sc.read_h5ad(results_file)
	figure_dir = './figures'
	sc.settings.figdir = figure_dir
	
	## Do FeaturePlots for select genes	
	expressed_dict = dict()
	
	for gene in adata.raw.var_names.values.tolist():
		if gene not in expressed_dict:
			expressed_dict[str(gene)] = 1
			#print(gene)
	
	genes_to_plot = []
	
	for gene in genes_of_interest:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_featureplots.png', show = False, cmap = my_feature_cmap, size = 25, use_raw = True)
	
	genes_to_plot = []
	scotts_genes = ['LGR4','LGR5','LGR6','OLFM4','ASCL2','SMOC2','MKI67','PCNA','LYZ','REG3A','DEFA1','DEFA2','DEFA3','DEFA4','DEFA5','DEFA6','CD24','CD44','DLL1','DLL4','TFF2','SOX4','SOX9','HES1','HES5','ATOH1','DCLK1','IL17RB','POU2F3','HOPX','BMI1','CLU1','CHGA','CHGB','TAC1','CCK','SCT','SST','GAST','GHRL','MLN','TPH1','THP2','PYY','GIP','ALPI','KRT20','MUC2','AXIN1']
	pharma_list = 	['CYP3A4','CYP3A5','CYP3A7','CYP2C8','CYP2C9','CYP2C19','CYP1A1','CYP2D6','CYP2E1','CYP2J2','UGT1A1','UGT1A3','UGT1A4','UGT1A6','UGT1A7','UGT1A8','UGT1A9','UGT1A10','SULT1A1','SULT1A3','SULT1B1','SULT2A1','SULT1E1','MDR1','BCRP','MRP1','MRP2','MRP3','MRP4','PEPT1','OCTN1','OCTN2','POU2F1','SLC22A2','POU5F1','PMAT','OAT1','OAT2','OAT3','OATP1A2','OATP1B3','OATP2B1','MCT1','SLC51A','SLC51B','MRP3','HNF4A','CD64','CRP','TNF','IL1B','IL6','DAO','GAST','HAMP','NFKB1','NOD2','MUC3','HLA-DRB1','HLA-DQB1','HLA-DQA1','HLA-DPB1','HLA-DRA','JAK1','JAK2','JAK3','STAT1','STAT2','STAT3','STAT4','STAT5A','STAT5B','STAT6','S1PR1','DHODH','TYK2','HIF1A','PXR','FGF19','FXR','NR3C1','CDX2','MUC2','MALRD1','ASBT','KLB','FGFR4','SMVT','CNT1','CNT2','ENT1','ENT2','CYB5R3','POR','RFK','CYB5A','CYB5R1','CYB5R2','CYB5R3','CYB5R4','MET','AHR','NR1I3','ARNT','SLC10A2','SLC5A6','SLC28A1','SLC28A2','SLC29A1','SLC29A2','ABCB1','ABCG2','ABCC1','ABCC2','ABCC3','ABCC4','SLC15A','SLC22A4','SLC22A5','SLC22A1','SLC22A2','SLC22A3','SLC29','SLC22A6','SLC22A7','SLC22A8','SLCO1A2','SLCO1B3','SLCO2B1','SLC16A','ABCC3']
	
	for gene in scotts_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_ScottsGenes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	
	
	genes_to_plot = []
	jbs_genes = ['Hashtag_1','Hashtag_2','Hashtag_3','Hashtag_4','Hashtag_5','Hashtag_6','MKI67','PCNA','LGR5','OLFM4','HOPX','MSI1','TERT','SOX9','RB1','RBL2','CCND1','CDK3','CDK4','FOXO3','EZH1','CYCS','KRT19','KRT20','CDKN2A','CDKN1B','H3C8','CYP2F1','GHR','TCP11L2','H2BC4','OGN','GBP2','PCDHA1','VNN1','RAB27B','SAA3P','H2BC21','PTGFR','DDX58','H2AC18','H2BC21','KDM7A','KDM5B','KDM6B','JARID2','PDCD4','SELENBP1','TOB1','IDH1','DDIT4','NFE2L2','YPEL2','ABCC3','CDKN1A']
	
	for gene in jbs_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_JBsGenes_featureplots.pdf', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	sc.pl.umap(adata, color=genes_to_plot, save = '_JBsGenes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	
	
	
	genes_to_plot = []
	
	for gene in pharma_list:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_pharma_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	
	for gene in epi_cell_type_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_epi_cell_type_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	genes_to_plot = []
	
	for gene in marker_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_marker_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)
	
	cytokine_receptors = ['EPOR','IL17RB','IL12B','CSF2RB','IL9R','IL23R','LIFR','IL12RB1','IFNAR1','IL17RE','CXCR5','CX3CR1','IL20RA','CCR7','CCR1','CSF3R','ACKR4','IL17RC','IL1RAPL2','CXCR3','GPR75','IL5RA','IL1R2','IL17RD','IL31RA','GHR','IL10RA','IL6ST','CD44','IL22RA2','EBI3','CCR9','IL21R','CRLF2','GPR17','OSMR','IL12RB2','IFNGR2','IFNGR1','CRLF1','GPR35','H0Y3Z8','IL1R1','IL22RA1','F3','ACKR3','CXCR4','IL10RB','CXCR2','CXCR1','IL7R','CMKLR1','IL17RA','FLT3','GFRA2','CCRL2','IL27RA','GFRAL','ACKR2','CXCR6','LEPR','PRLR','IL20RB','CD4','IL1RL2','IL2RA','CCR8','CCR6','CCR5','CCR4','CCR3','CNTFR','IL2RG','IL3RA','IFNAR2','IL2RB','IFNLR1','GFRA3','CD74','GFRA4','XCR1','CCR10','IL4R','IL1RL1','IL17REL','IL15RA','MPL','CCR2','IL18R1','IL13RA1','IL11RA','IL13RA2','IL1RAP','CSF2RA','interleukin23-receptor_human','il12-receptor_human','IL18RAP','FZD4','GFRA1','IL6R']

	
	genes_to_plot = []
	
	for gene in cytokine_receptors:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_cytokine_receptor_genes_featureplots.png', show = False, cmap = my_feature_cmap, size = dot_size*3, use_raw = True)

print('Starting gene expression plotting...\n------------------------------------\n')


print('Checking for expression of genes of interest\n')
expressed_dict = dict()

for gene in adata.raw.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1

genes_to_plot = []
	
for gene in genes_of_interest:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Found cells expressing', ' '.join(genes_to_plot), '\n')


if run_marker_analysis:
	print("\nAll done with general workflow... now finding marker genes.\n")
	## Find marker genes via Wilxocon test based on Louvain cluster assignment
	# Create a simple plot to show the top 25 most significant markers for each cluster
	
	sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
	
	write_marker_file(adata, file_out='markers_output.csv', n_genes=100)
	
	sc.tl.filter_rank_genes_groups(adata, groupby='louvain', use_raw=True, log=True, key_added='rank_genes_groups_filtered', min_in_group_fraction=0.25, min_fold_change=1.5, max_out_group_fraction=0.8)
	
	sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered', n_genes=30, sharey=False, save = '_markerPlots.pdf', show = False)
	sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered',  n_genes=20, save = '_markerDotPlots.pdf', color_map=my_dot_cmap, show = False, mean_only_expressed=True, dot_min=0.2, dot_max=1, standard_scale='var')
	

'''
sc.tl.diffmap(adata)

sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, use_rep='X_diffmap')

sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color='louvain', save='_paga-init-pattern_denoised.pdf', show=False, threshold=threshold, node_size_scale=node_size_scale, node_size_power=0.9, layout=paga_layout)

sc.tl.umap(adata, init_pos='paga', min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma)

sc.pl.umap(adata, color='louvain', save = '_clusterIdentity_denoised.pdf', show = False, legend_loc = 'on data', edges = True, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
sc.pl.umap(adata, color='louvain', save = '_clusterIdentity_noEdge_denoised.pdf', show = False, legend_loc = 'on data', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = dot_size, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
sc.pl.umap(adata, color='phase', save = '_cell_cycle_denoised.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
sc.pl.umap(adata, color='sampleName', save = '_sample_denoised.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, palette = greatestPalette, alpha = 0.95)
sc.pl.umap(adata, color=['Hashtag_1','Hashtag_2','Hashtag_3','Hashtag_4','Hashtag_5','Hashtag_6'], save = '_antibody_features_denoised.pdf', show = False, legend_loc = 'right margin', edges = False, size = dot_size, cmap = my_feature_cmap, alpha = 0.95)



for gene in adata.raw.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1
		#print(gene)

genes_to_plot = []

for gene in genes_of_interest:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Plotting genes:', ' '.join(genes_to_plot),'\n')

sc.pl.umap(adata, color=genes_to_plot, save = '_featureplots_denoised.png', show = False, cmap = my_feature_cmap, size = 25, use_raw = True)



sc.external.exporting.spring_project(adata, project_dir='spring_expt', embedding_method='umap', subplot_name=None, cell_groupings='louvain', custom_color_tracks=None, total_counts_key='n_counts', overwrite=False)
'''



print('\nDone with entire script execution')










