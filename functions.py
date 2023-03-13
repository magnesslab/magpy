import numpy as np
import pandas as pd
import matplotlib as mpl
import sys

# import os
# import scanpy as sc
# import seaborn as sns
# import magpy as mp
# import magpy.settings as settings
# import matplotlib.pyplot as plt
# from shutil import copyfile
# import loompy as lpy
# import scvelo as scv
# from matplotlib import rcParams
# from matplotlib.colors import ColorConverter, is_color_like
# from pandas import unique, isnull, Index


# Sort cells in preparation for heatmap plotting
# Can sort by any gene name or .obs annotations
def heatmap_sort(adata, groupby='lineage', sortby=None, ascending=True):
    
    if groupby not in adata.obs_keys():
        print(f'groupby parameter [{groupby}] is not present in adata.obs annotations. Please try again with one of these keys: \n\n',adata.obs_keys())
        return None
    if sortby in adata.obs_keys():
        df = adata.obs[[groupby,sortby]]
    elif sortby in adata.var_names.to_list():
        df = pd.DataFrame({groupby:adata.obs[groupby],sortby:adata[:,sortby].X.flatten()}, index = adata.obs_names)
    else:
        print(f'sortby parameter [{sortby}] is not present in adata.obs annotations or adata.var names. \nPlease try again with one of those values.')
        return None
    
    cells = []
    for category in adata.obs[groupby].cat.categories:
        cells += df[df[groupby]==category].sort_values(sortby,ascending=ascending).index.tolist()
    
    return adata[cells,:]

#Divide all cells by the mean of the highest cluster in groupby
def scale_adata(adata, groupby='lineage', layer='raw_scaled'):
    clusters = adata.obs[groupby].unique()
    cmeans = np.zeros((len(clusters),adata.shape[1]))
    for i,cluster in enumerate(clusters):
        subset = adata[adata.obs[groupby]==cluster]
        cmeans[i] = subset.layers['raw_normalized'].mean(axis=0).A1
    cmax = cmeans.max(axis=0)
    adata.layers[layer] = adata.layers['raw_normalized'] / cmax[None,:]
    return adata

#Returns gene_list, sorted by max % expression in groups specified by 'groupby'
def sort_list_by_pctexp(adata, gene_list, groupby='lineage', ascending=False):
    subset = adata[:,gene_list]
    clusters = subset.obs[groupby].unique()
    pctexp = np.zeros((len(clusters),subset.shape[1]))
    for i,cluster in enumerate(clusters):
        subset2 = subset[subset.obs[groupby]==cluster]
        pctexp[i] = (subset2.layers['raw_normalized']>0).sum(axis=0).A1 / subset.shape[0]
    subset.var['pctmax'] = pctexp.max(axis=0)
    sorted_list = subset.var['pctmax'].sort_values(ascending=ascending).index.tolist()
    return sorted_list

#Returns gene_list, sorted by mean expression in groups specified by 'groupby'
def sort_list_by_mean(adata, gene_list, groupby='lineage',ascending=False):
    subset = adata[:,gene_list]
    clusters = subset.obs[groupby].unique()
    cmeans = np.zeros((len(clusters),subset.shape[1]))
    for i,cluster in enumerate(clusters):
        subset2 = subset[subset.obs[groupby]==cluster]
        cmeans[i] = subset2.layers['raw_normalized'].mean(axis=0).A1
    subset.var['cmeans'] = cmeans.max(axis=0)
    sorted_list = subset.var['cmeans'].sort_values(ascending=ascending).index.tolist()
    return sorted_list

#Generates and (optionally) saves a dataframe of 'cluster x gene' mean expression values
def means_to_df(adata, gene_list, save=None, groupby='lineage'):
    subset = adata[:,gene_list]
    clusters = subset.obs[groupby].unique()
    cmeans = {}
    for i,cluster in enumerate(clusters):
        subset2 = subset[subset.obs[groupby]==cluster]
        cmeans[cluster] = subset2.layers['raw_normalized'].mean(axis=0).A1
    df = pd.DataFrame(cmeans,index=gene_list)
    if save: df.to_csv(save)
    return df

def sc_v_bulk(sc_data = None, bulk_data = None, v = None, gene_list_title = None):

    ## initialize dataframe
    df = pd.DataFrame(index = v)
    print(f'studying similarity of in vitro and in vivo data for {gene_list_title}')
    
    ## subset by genes and or cell types if not already completed
    # subset1 = adata[adata.obs.lineage=='C_lateACC',v]
    subset2 = adata[:,v]
#     print(subset2)
    lineage_list = [lin for lin in subset2.obs.lineage.unique()]
    lineage_order = ['C_ISC','C_earlyACC','C_lateACC']
    lineage_list = sorted(lineage_list, key = lineage_order.index)
    
    ## subset bulk RNAseq data to the same gene list 
    bulk_subset1 = bulk_data[:,v]
    
    ## calculate per cluster means for sc data
    for i,cluster in enumerate(lineage_list):
        subset2.varm[f'{cluster}_mean'] = np.zeros((subset2.n_vars, 4))
        subset3 = subset2[subset2.obs.lineage==cluster]
        subset2.varm[f'{cluster}_mean'][:,i] = subset3.X.mean(axis = 0).A1.flatten()
        # print(subset.varm['means'])
    
    ## add mean expression values for gene list of interest from single cell to initialized dataframe
    for i, cluster in enumerate(lineage_list):
        df[f'{cluster}_mean'] = subset2.varm[f'{cluster}_mean'][:,i]
    
    ## for differentiation, order the bulk RNAseq data based on the duration of differetiation
    day_list = [a for a in bulk_subset1.obs['diff'].unique()]
    day_order = ['DC_EM1', 'DC_EM4', 'DC_DM2', 'DC_DM4_bat1']
    day_list = sorted(day_list, key = day_order.index)

    ## calculate per differentiation time point mean expression values and add these 
    for i,day in enumerate(day_list):
        bulk_subset1.varm[f'{day}_mean'] = np.zeros((bulk_subset1.n_vars, 6)) #gene, day
        bulk_subset2 = bulk_subset1[bulk_subset1.obs['diff']==day]
        bulk_subset1.varm[f'{day}_mean'][:,i] = bulk_subset2.X.mean(axis=0)

    for i,day in enumerate(day_list):
        df[f'{day}_mean'] = bulk_subset1.varm[f'{day}_mean'][:,i]
    
    #### how to draw each day on a separate plot in matplotlib####
    
    ## this variable is used for highlighting certain linear regression plots for differentiation days (bulk) vs clusters (single cell)
    # best_fits = [(5,1),(4,1),(3,2),(2,2),(1,3),(0,2)]
    
    ## initialize matplotlib objects to hold the array of linear regressions
#     fig, axes = plt.subplots(6,4,figsize = [20,15], sharey = True)
    
#     color = 'bo'
#     slope = 1
    
#     # print(df)
    
#     for i,cluster in enumerate(lineage_list):
#         col = i%6
#         for j,day in enumerate(day_list):
#             row = j%6
#             for coords in best_fits:
#                 (x_coord, y_coord) = coords
#                 if i == y_coord and j == x_coord: 
#                     color = 'ro'
#                     # print(f'{cluster} ---- {day} ----- {color}')
#             axes[row,col].plot(df[f'{cluster}_mean'], df[f'{day}_mean'], color)
#             x_min, x_max = axes[row,col].get_xlim()
#             c = x_min
#             y_min, y_max = c, c + slope*(x_max-x_min)
#             axes[row,col].plot([x_min, x_max], [y_min, y_max], color = 'black', ls = '--')
#             axes[row,col].set_xlim([x_min, x_max])
#             if col ==0:
#                 axes[row,col].set(xlabel = f'{cluster}', ylabel = f'{day}')
#             elif col >0:
#                 axes[row,col].set(xlabel = f'{cluster}')
#             color = 'bo'

    #### how to draw each day on one figure ####

    # colors = ['bo','ro','go','yo','ko','mo']

    # fig,axis = plt.subplots()
    # for j, cluster in enumerate(lineage_list):
    #     for i,day in enumerate(day_list):
    #         axis.plot(df[cluster], df[day], colors[i])
    #         axis.legend(day_cols)
    
    ## initialize ditionary to hold results for residual calculation
    results = {}
    
    for i,cluster in enumerate(lineage_list):
        results[cluster] = []
        for j,day in enumerate(day_list):

            # axis = sns.residplot(x = df['sc_means'], y = df[day], data = df, label = day[:-8], lowess = True,
            #                      robust = True)
            # axis.legend()
            # axis.set(xlabel = 'sc_mean', ylabel = 'bulk_mean')

            #define response variable
            y = df[f'{day}_mean']

            #define predictor variables
            x = df[f'{cluster}_mean']

            #add constant to predictor variables
            x = sm.add_constant(x)

            #fit linear regression model
            model = sm.OLS(y, x).fit()

            #view model summary, calculate sum of squares
#             print(sep)
#             print(day)
            results[cluster].append(((df[f'{day}_mean'] - df[f'{cluster}_mean'])**2).sum())

#             print(model.summary())

#             print('Calculated residual sum of sqares: ',model.ssr)
    
    df2 = pd.DataFrame(results, index = day_list)
    df2 = df2[['C_ISC','C_earlyACC','C_lateACC']]
    
    # ax = plt.matshow(df2)
    # ax.set_xlabel('in_vivo data')
    # print(df2)
    
    return df2

## Function for pulling out cells with specific gene expression patterns
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

## useful for making custom colormaps if you want that
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

