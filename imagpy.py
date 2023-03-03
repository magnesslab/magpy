import os
import pandas as pd

def fit_gmm(data, file_path = None, file_name = None, num_components = None, ln_transform = True, columns = None, labels = None):

	#check for correct input variables#
	if data is None:
		if file_path and file_name:
			if '.csv' in file_name: 
				data = pd.read_csv(file_path+file_name)
			else:
				raise Exception('input file must be a csv')
		if not file_path and file_name:
			raise Exception('Input data or file location must be specified')
	
	if not num_components:
		raise Exception('Number of components must be specified using the num_components variable')
	else: n_components = num_components
	display('Showing input data')
	display(data)
	
	# if not labels:
	# 	raise Exception('labels of gmm clusters must be provided using the labels variable in list format.')
	# elif len(labels) != n_components:
	# 	raise Exception(f'Number of specified labels ({len(labels)}) is not equal to the number of components ({n_components}) to fit')

	if not columns: 
		raise Exception('Columns to generate gmm model must be specified using the columns variable in list format')

	if ln_transform: 
		data = np.log(data)
		display('Displaying log-transformed data')
		display(data)
		
	label_map = {}
	
	##begin fitting gmm on specified list of columns##
	
	if columns: 
		for column in columns:
			display(f'Histogram for {column}') 
			plt.hist(data[f'{column}'], density = True, bins = 100)
			plt.show()
			gmm= mixture.GaussianMixture(n_components=n_components, covariance_type='full', n_init = 100, max_iter = 10000).fit(data[f'{column}'].to_numpy().reshape(-1,1))
			data_labels = gmm.predict(data[f'{column}'].to_numpy().reshape(-1,1))
			data[f'{column}_label'] = data_labels
			data_means = gmm.means_
			data_means = data_means.ravel()
			data_sort = data_means.argsort().tolist()
			display(data)
			if n_components ==4:
				max_sort_pos = data_sort[3]
				min_sort_pos = data_sort[0]
				sublow_sort_pos = data_sort[1]
				low_sort_pos = data_sort[2]
				
				label_map = {
				max_sort_pos:'High',
				min_sort_pos:'Neg',
				low_sort_pos:'Mid',
				sublow_sort_pos:'Sublow'
				}
				
				print('Total cells in condition :',len(data[f'{column}']))
				d0 = pd.DataFrame(data[data[f'{column}_label']== 0], columns = ['predicted_label',f'{column}']) 
				print('Number cells in first cluster :',len(d0))
				d1 = pd.DataFrame(data[data[f'{column}_label']== 1], columns = ['predicted label',f'{column}']) 
				print('Number cells in second cluster :',len(d1))
				d2 = pd.DataFrame(data[data[f'{column}_label']== 2], columns = ['predicted_label',f'{column}'])
				print('Number cells in third cluster :',len(d2))
				d3 = pd.DataFrame(data[data[f'{column}_label']== 3], columns = ['predicted_label',f'{column}'])
				print('Number cells in fourth cluster :',len(d3))
				plt.hist(d0, density = True, bins = 50)
				plt.show()
				plt.hist(d1, density = True, bins = 50)
				plt.show()
				plt.hist(d2, density = True, bins = 50)
				plt.show()
				plt.hist(d3, density = True, bins = 50)
				plt.show()
			
			elif n_components ==3:
				max_sort_pos = data_sort[2]
				low_sort_pos = data_sort[0]
				mid_sort_pos = data_sort[1]
				
				label_map = {
				max_sort_pos:'High',
				low_sort_pos:'Low',
				mid_sort_pos:'Mid',
				}
				
				print('Total cells in condition :',len(data[f'{column}']))
				d0 = pd.DataFrame(data[data[f'{column}_label']== 0], columns = ['predicted_label',f'{column}']) 
				print('Number cells in first cluster :',len(d0))
				d1 = pd.DataFrame(data[data[f'{column}_label']== 1], columns = ['predicted label',f'{column}']) 
				print('Number cells in second cluster :',len(d1))
				d2 = pd.DataFrame(data[data[f'{column}_label']== 2], columns = ['predicted_label',f'{column}'])
				print('Number cells in third cluster :',len(d2))
				plt.hist(d0, density = True, bins = 50)
				plt.show()
				plt.hist(d1, density = True, bins = 50)
				plt.show()
				plt.hist(d2, density = True, bins = 50)
				plt.show()
				
			data[f'{column}_label']=data[f'{column}_label'].map(label_map)
			
	display(data)