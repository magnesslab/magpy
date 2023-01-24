#Archive
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

#Archive
def combine_loom(sample_list='', expt_path='', save=True, read_file=None, write_file=None, loom_files = None):
	
	if read_file is None: read_file = settings.loom_file
	if write_file is None: write_file = settings.merged_loom_file
	merged_path = os.path.join(expt_path,write_file)
	
	ldata_list = []
	
	if loom_files is None:
		for sample_ID in sample_list:
			metadata=mp.pipeline.load_metadata(sample_ID)
			new_path = expt_path+metadata[0]
			loom_path = os.path.join(new_path,read_file)
			ldata_list.append(loom_path)
			print(loom_path)
	elif type(loom_files) == list:
		for file in loom_files:
			print(file)
			ldata_list.append(file)
	else:
		print('looom files not identified! check your directories and make sure files are in the proper places')
	#print(ldata_list)
	
	print('Combining all loom files into one merged loom file')
	
	lpy.combine(ldata_list, merged_path)
	
	print("Loom file merging complete.\n")

#Archive
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
