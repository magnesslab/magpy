import os
import math
import scanpy as sc
import numpy as np
import pandas as pd
import diffxpy.api as de
from scipy.stats import rankdata

### FUNCTIONS FOR BATCHED TESTING ###

#Perform Benjamini-Hochberg multiple-testing correction on an ndarray of p-values
def _bh_correct(pvals):
    ranked_pvals = rankdata(pvals)
    adjusted_pvals = (pvals * len(pvals)) / ranked_pvals
    adjusted_pvals[adjusted_pvals>1] = 1
    return(adjusted_pvals)

#Limit compute requirements of differential expression by splitting into batches
def _batch_test(adata, formula, test_param, batch_size):
    #Determine number of batches
    num_batches = math.ceil(adata.shape[1]/batch_size)
    
    results = []
    for i in range(num_batches):
        #Slice the data to the current batch
        print(f"Current batch: {i+1} of {num_batches}")
        if (i+1)*batch_size > adata.shape[1]: subset = adata[:, i*batch_size:].copy()
        else: subset = adata[:, i*batch_size:(i+1)*batch_size].copy()
            
        #Run wald test on data subset
        test = de.test.wald(data=subset, formula_loc=formula, factor_loc_totest=test_param)
        results.append(test.summary())
        
    #Concatenate results for all batch tests
    results_df = pd.concat(results).set_index('gene')
    
    #Recalculate FDR-corrected q-value based on length of entire data set
    results_df['qval'] = _bh_correct(results_df['pval'].to_numpy())
    return  results_df

#Run a single set of tests on a dataset
def _run_test(prefix, adata, formula, test_param, in_group, min_cells, batch_size):

    #Calculate gene expression in each test category
    cdict = {}
    adata.varm['cmeans'] = np.zeros((adata.n_vars, len(adata.obs[test_param].unique())))  
    for i,category in enumerate(adata.obs[test_param].unique()):
        cdict[category] = i
        subset = adata[adata.obs[test_param]==category,:]
        adata.varm['cmeans'][:,i] = subset.X.mean(axis=0)
        
    #Check number of cells in each test category
    num_cells = {}
    for category in adata.obs[test_param].unique():
        num_cells[category] = (adata.obs[test_param] == category).sum()
            
    #Generate empty results dataframe
    results = pd.DataFrame(index=adata.var_names)
    
    #If the test param is boolean and in_group is not set, make True the in group
    if (adata.obs[test_param].dtype == bool) and (in_group is None): in_group = True

    ### CASE 1 - TEST PARAMETER IS NOT BOOL, IN_GROUP NOT SPECIFIED ###
    if in_group is None:
        #Check for adequate cell numbers
        if adata.n_obs < (min_cells * len(adata.obs[test_param].unique())):
            print(f"Not enough cells in {prefix} to perform DE.")
            print("To run DE with less cells, specify min_cells in the function call.")
            return None
        
        #Find mean expression and num cells for each category
        for category in adata.obs[test_param].unique():
            results[f"{prefix}_{category}_mean"] = adata.varm['cmeans'][:,cdict[category]]
            results[f"{prefix}_{category}_n-cells"] = num_cells[category]
            
    ### CASE 2 - TEST PARAMETER IS BINARY, IN_GROUP IS SPECIFIED ###
    elif len(adata.obs[test_param].unique()) == 2:
        #Check for adequate cell numbers
        for key in num_cells:
            if num_cells[key] < min_cells:
                print(f"Not enough cells in {prefix} to perform DE.")
                print("To run DE with less cells, specify min_cells in the function call.")
                return None
        
        #Calculate means and fold changes
        results[f'{prefix}_in_mean'] = adata.varm['cmeans'][:,cdict[in_group]]
        results[f'{prefix}_out_mean'] = adata.varm['cmeans'][:,1-cdict[in_group]]
        results[f'{prefix}_log2fc'] = np.log2(results[f'{prefix}_in_mean'] / results[f'{prefix}_out_mean'])
        
        #Count number of cells
        results[f'{prefix}_in_n-cells'] = num_cells[in_group]
        results[f'{prefix}_out_n-cells'] = adata.n_obs - num_cells[in_group]
    
    ### CASE 3 - TEST PARAMETER IS NOT BINARY, IN_GROUP IS SPECIFIED ###
    else:
        #Check for adequate cell numbers
        if (num_cells[in_group] < min_cells) or ((adata.n_obs - num_cells[in_group]) < min_cells):
            print(f"Not enough cells in {prefix} to perform DE.")
            print("To run DE with less cells, specify min_cells in the function call.")
            return None
            
        #Calculate means and fold changes
        in_mean = adata.varm['cmeans'][:,cdict[in_group]]
        fc = np.log2(in_mean[:,None] / adata.varm['cmeans'])
        results[f'{prefix}_in_mean'] = in_mean
        results[f'{prefix}_min_log2fc'] = fc.min(axis=1, where=(fc != 0), initial=100)
        results[f'{prefix}_max_log2fc'] = fc.max(axis=1, where=(fc != 0), initial=-100)

        #Count number of cells
        results[f'{prefix}_in_n-cells'] = num_cells[in_group]
        results[f'{prefix}_out_n-cells'] = adata.n_obs - num_cells[in_group]
        
        #After calculating fold-changes, make this into a one v rest test
        adata.obs['beta'] = adata.obs[test_param] == in_group
        formula = formula.replace(test_param,'beta')
        test_param = 'beta'
    
    #Run DE tests
    print(f"\nFitting models for {prefix} dataset...")
    test_results = _batch_test(adata, formula, test_param, batch_size)
    results[f'{prefix}_pval'] = test_results['pval']
    results[f'{prefix}_qval'] = test_results['qval']
    print(f"Fitting of {prefix} dataset complete.")
    
    return results
    
#UI for running tests
def de_test(adata, test_param, batch_param='donor', 
            in_group=None, test_each=None, test_order = None, 
            min_expression=0.1, min_cells=10, batch_size=250):
    
    if in_group:
        if in_group not in adata.obs[test_param].unique():
            raise Exception("in_group should be a category in test_param")
    
    #Calculate gene expression in each test category
    adata.varm['cmeans'] = np.zeros((adata.n_vars, len(adata.obs[test_param].unique())))  
    for i,category in enumerate(adata.obs[test_param].unique()):
        subset = adata[adata.obs[test_param]==category,:]
        adata.varm['cmeans'][:,i] = subset.X.mean(axis=0)
        
    #Subset to only genes with at least min_expression in one one test category
    adata = adata[:,np.any((adata.varm['cmeans'] > min_expression), axis=1)]
    print(f"Fitting distributions for {adata.n_vars} genes.")
    
    #Run tests on combined dataset
    results_list = []
    formula = f"~ 1 + {batch_param} + {test_param}"
    results = _run_test('combined', adata, formula, test_param, in_group, min_cells, batch_size)
    if results is None: return None #Exit out if not enough cells in combined dataset
    results_list.append(results)
    
    #Run tests on individual datasets
    if test_each:
        if not test_order: test_order = adata.obs[test_each].unique()
        if test_each == batch_param: formula = f"~ 1 + {test_param}"
        for category in test_order:
            #Subset the data
            subset = adata[adata.obs[test_each]==category,:]
            results = _run_test(category, subset, formula, test_param, in_group, min_cells, batch_size)
            if results is not None: results_list.append(results)
    
    #Concatenate results
    results = pd.concat(results_list, axis=1)
    
    #Order results columns so qvals are first
    results = results[[c for c in results if 'qval' in c]+[c for c in results if 'qval' not in c]]
    
    return results