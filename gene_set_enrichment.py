import os
import magpy as mp
import scanpy as sc
from glob import glob

mp_path = os.path.dirname(mp.__file__)
gsea_list_path = mp_path + "/gene_lists/msigdb.v7.2.symbols.gmt"
user_list_path = mp_path + "/gene_lists"


def load_gsea_gene_lists():
    """
    Returns a dict of {'list_name':[genes]} of all GSEA gene lists
    """
    gene_lists = {}
    with open(gsea_list_path) as file:
        for line in file:
            split_line = line.strip().split("\t")
            gene_lists[split_line[0]] = split_line[2:]
    return gene_lists


def load_user_gene_lists():
    """
    Returns a dict of {'list_name':[genes]} of all user-defined gene lists
    """
    files = glob(user_list_path + "/*.gene_list")
    gene_lists = {}
    for file in files:
        with open(file) as f:
            gene_list = [line.strip() for line in f]
        file_name = os.path.split(file)[-1][:-10]
        gene_lists[file_name] = gene_list
    return gene_lists


def load_all_gene_lists():
    """
    Returns a dict of {'list_name':[genes]} of all saved gene lists
    """
    gsea_gene_lists = load_gsea_gene_lists()
    user_gene_lists = load_user_gene_lists()
    gsea_gene_lists.update(user_gene_lists)
    return gsea_gene_lists


def load_gene_list(search_term):
    """
    Accepts a string or list of strings. Case-sensitive.
    Returns a dict of {'list_name':[genes]} from one or more specified GSEA gene lists
    """
    gene_lists = load_all_gene_lists()
    
    if type(search_term) == str:
        return {search_term:gene_lists[search_term]}
    
    if type(search_term) == list:
        return {list_name:gene_list for list_name,gene_list in gene_lists.items() if list_name in search_term}


def load_gene_lists(search_term):
    """
    Alias for load_gene_list
    """
    gene_lists = load_gene_list(search_term)
    return gene_lists


def save_gene_list(gene_list):
    """
    Accepts a dict of {'list_name':[genes]}
    Saves each gene list to a file, where each key becomes the list name
    """
    for key in gene_list:
        save_path = user_list_path + f"/{key}.gene_list"
        with open(save_path,"w") as f:
            f.writelines([gene+"\n" for gene in gene_list[key]])
    return None


def save_gene_lists(gene_list):
    """
    Alias for save_gene_list
    """
    save_gene_list(gene_list)
    return None


def search_gene_lists(search_term):
    """
    Accepts a string. Case-insensitive, name can be partial match.
    Returns list of all available gene lists containing search term in their name
    """ 
    gene_lists = []

    user_lists = load_user_gene_lists()
    for key in user_lists:
        if search_term.lower() in key:
            print(key)
            print()
            gene_lists.append(key)

    with open(gsea_list_path) as file:
        for line in file:
            split_line = line.strip().split("\t")
            if search_term.lower() in split_line[0].lower():
                print(f"{split_line[0]}")
                print(f"{split_line[1]}")
                print()
                gene_lists.append(split_line[0])
    return gene_lists
    

def search_gsea_by_gene(search_term):
    """
    Accepts a string. Case-sensitive, gene must be exact match.
    Returns list of all available gene lists containing the specified gene
    """
    gene_lists = []
    with open(gsea_list_path) as file:
        for line in file:
            split_line = line.strip().split("\t")
            if search_term in split_line[2:]:
                print(f"{split_line[0]}")
                print(f"{split_line[1]}")
                print()
                gene_lists.append(split_line[0])
    return gene_lists


def filter_genes(adata, genes):
    """
    Accepts an AnnData object and one of:
        Dict of {'list_name':[genes]}
        List of genes

    Returns a dict, filters all gene lists to only contain genes in the given adata object.
    """
    if type(genes) == dict:
        new_dict = {}
        for list_name,gene_list in genes.items():
            filtered_list = [item for item in gene_list if item in adata.var_names.tolist()]
            missing_list = [item for item in gene_list if item not in filtered_list]
            if len(missing_list) > 0:
                print(f"The following genes were filtered out for {list_name}:")
                print(missing_list)
                print()
            new_dict[list_name] = filtered_list
        return new_dict

    elif type(genes) == list:
        filtered_list = [item for item in genes if item in adata.var_names.tolist()]
        missing_list = [item for item in genes if item not in filtered_list]
        print(f"The following genes were filtered out:")
        print(missing_list)
        print()
        return filtered_list


def score_genes(adata, gene_list, plot=True, inplace=True, normalize=False):
    """
    Scores adata objects for the specified gene list(s)

    gene_list can be:
        Dict of {'list_name':[genes]}
        String specifying name of a saved gene list
        List of gene list names

    plot: Should the scored genes be displayed on UMAP plots
    inplace: Should annotation be added to the given adata object. If false a new adata object is returned.
    normalize: Should the score results be normalized to 1. 
    """ 

    if (type(gene_list) == str) or (type(gene_list) == list):
        gene_dict = load_gene_list(gene_list)
        
    elif type(gene_list) == dict:
        gene_dict = gene_list
    
    else:
        raise Exception("gene_list must be a string, list, or dict")

    if not inplace:
        adata = adata.copy()
    
    gene_dict = filter_genes(adata, gene_dict)
    for key in gene_dict:
        sc.tl.score_genes(adata, gene_list=gene_dict[key], score_name=key)
        if normalize:
            adata.obs[key] = adata.obs[key] / adata.obs[key].max()
 
    if plot:
        sc.pl.umap(adata, color=[key for key in gene_dict], ncols=3)
                
    return adata