import os
import scanpy as sc
import magpy.settings as settings
# import psutil

#Alias for sc.read_h5ad and sc.read_10x_h5, file_choice should be an integer or file name
def load(expt_path, file_choice, anndata=True, hashed=False):
	#Build read path
	if type(file_choice) == int:
		files = [settings.raw_file, settings.adata_file, settings.pp_file, settings.pca_file, settings.cluster_file, settings.merged_file, settings.clustered_velocity_file]
		read_path = os.path.join(expt_path, files[file_choice])
	else: read_path = os.path.join(expt_path, file_choice)

	#Read data
	print(f"Reading data from {read_path}")
	if anndata: adata = sc.read_h5ad(read_path)
	else: adata = sc.read_10x_h5(read_path, gex_only=(not hashed)) 
	print()

	return (adata)

#Alias for adata.write, to allow integer-based saving
def save(adata, expt_path, file_choice):
	if type(file_choice) == int:
		files = [settings.raw_file, settings.adata_file, settings.pp_file, settings.pca_file, settings.cluster_file, settings.merged_file, settings.clustered_velocity_file]
		write_path = os.path.join(expt_path, files[file_choice])
	else: write_path = os.path.join(expt_path, file_choice)

	print(f"Saving annotated data to {write_path}")
	adata.write(write_path)

#TODO - Add functionality for searching for a list of genes
#Search adata.var_names for terms containing search_term
def search(adata, search_term):
	
	if isinstance(search_term,str) is True:
		results = [gene for gene in adata.var_names if search_term.lower() in gene.lower()]
		if len(results) == 0: print(f"Sorry, searching for '{search_term}' yielded no results.")
		else:
			print(f"Search for '{search_term}' yielded {len(results)} result(s):")
			for result in results:
				print(result)
			print()
		return(results)
	else:
		gene_list=[]
		not_found=[]
		for item in search_term:
			results = [gene for gene in adata.var_names if item.lower() == gene.lower()]
			if len(results) == 0: 
				not_found.append(item)
			else:
				for result in results:
					gene_list.append(result)
		print('Genes not found: ',len(not_found))
		print('Genes found: ',len(gene_list))
		print('Gene hits: ',gene_list)
		return(gene_list)