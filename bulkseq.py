from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import string
import sys
import random

#############################################################################
##                          Define run parameters                          ##
#############################################################################

## Initialize global variables, dictionaries and empty lists for data handling 
geneValuesDict = dict()
clusterGenesDict = dict()
dataMatrix = []
symbolList = []
sse = []
sil = []

## load input files and convert files into formats usable for downstream processing ##
def process(expt_path = '' ,read_file = None, df = None):
		
	readpath = os.path.join(expt_path,read_file)
	
	inFile = open(readpath, 'r')

	## Read TPM data matrix for all samples
	for count,line in enumerate(inFile):
		elem = str.split(line.rstrip())
		if count == 0:
			header = 'geneSymbol\tcluster\t'+line
			samplelist = line
		else:
			floatList = [float(i) for i in elem[1:]]
			if sum(floatList) >= 10:
				if elem[0] not in geneValuesDict:
					geneValuesDict[elem[0]] = floatList
				symbolList.append(elem[0])
				dataMatrix.append(floatList)

	## Convert data matrix (list of lists) into NumPy array
	dataArray = np.array(dataMatrix)

	## print(out array to prompt before and after normalization)
	## Normalize array by calculating z-score across samples for each gene (axis 1 = <-->)
	#print('\nRaw TMM normalized TCM Data array')
	#print(dataArray)
	zArray = stats.zscore(dataArray, axis=1, ddof=df)
	## Alternative normalization method using log transformation instead of z-score (commented out now)
	#zArray = np.log2(dataArray)
	#print('\nZ-score by row/gene normalized data array')
	#print(zArray)
	
	return(zArray,header,samplelist)


## function returns WSS score for k values from 1 to kmax
def calculate_WSS(zArray, expt_path = '',write_file = None):

	if write_file is None: write_file = 'unclaimed_WSS_plot.pdf'

	plotFileName = os.path.join(expt_path+'WSS_plot.pdf')
	for k in range(1, kmax+1):
		kmeans = KMeans(n_clusters = k).fit(zArray)
		centroids = kmeans.cluster_centers_
		pred_clusters = kmeans.predict(zArray)
		curr_sse = 0
 
		# calculate square of Euclidean distance of each point from its cluster center and add to current WSS
		for i in range(len(zArray)):
			curr_center = centroids[pred_clusters[i]]
			curr_sse += (zArray[i, 0] - curr_center[0]) ** 2 + (zArray[i, 1] - curr_center[1]) ** 2
		sse.append(curr_sse)

	print('Number of kmeans clusters: ',list(range(1,kmax)),'\n')
	print('WSS scores: ',sse)
	plt.title('WSS score by number of kMeans clusters')

	x = list(range(1,kmax+1))
	y = sse
	plt.xlabel('Number of kmeans clusters')
	plt.ylabel('WSS Score')
	plt.plot(x, y, fmt, alpha=0.5)
	plt.xticks(x)
	plt.ylabel('WSS Score')
	plt.savefig(plotFileName)
	plt.clf()	
	print('\nThe resulting plots are saved in', plotFileName)


##	The Silhouette Score measures how similar a point is to its own cluster (cohesion) compared to other clusters (separation)
def calculate_Silhouette(zArray):


	plotFileName = file_tag+'silhouette_plot_lfc-2.pdf'
# dissimilarity would not be defined for a single cluster, thus, minimum number of clusters should be 2
	for k in range(2, kmax+1):
		kmeans = KMeans(n_clusters = k).fit(zArray)
		labels = kmeans.labels_
		sil.append(silhouette_score(zArray, labels, metric = 'euclidean'))
	
	print('Number of kmeans clusters: ',list(range(1,kmax)),'\n')
	print('Silhouette scores: ',sil)
	plt.title('Silhouette score by number of kMeans clusters')
	x = list(range(1,kmax))
	y = sil
	plt.xlabel('Number of kmeans clusters')
	plt.ylabel('Silhouette Score')
	plt.plot(x, y, fmt, alpha=0.5)
	plt.xticks(x)
	plt.ylabel('Silhouette Score')
	plt.savefig(plotFileName)
	plt.clf()	
	print('\nThe resulting plots are saved in', plotFileName)

## Generate the actual gene clusters via kmeans then write the cluster assignments to a file
def generate_kmeans(zArray, expt_path='',df=None,header='',staggerAmt = 0.15,samples = None):

	sampleArray = np.array(samples)
	print(sampleArray)

	for k in range(0, 3):
		outPutFileName = expt_path+'kmeansOutPython_k-' + str(k) + '.tsv'
		plotTitle = 'kMeans z-score Plots By Cluster: k=' + str(k)
		plotFileName = expt_path+'kMeans_dotPlot_z-scorePlots_k-' + str(k) + '.pdf'
		
		clusterGenesDict = dict()
		
		outFile = open(outPutFileName, 'w')
		
		## Generate the actual gene clusters via kmeans then write the cluster assignments to a file
		print('\nCalculating', k, 'clusters with', df, 'degrees of freedom')
		kmeans = KMeans(n_clusters=k, init='k-means++', random_state=0, verbose=0).fit(zArray)
		zScoredListOfLists = zArray.tolist()
		
		outFile.write(header)
		for i in range(len(kmeans.labels_)):
			if kmeans.labels_[i] not in clusterGenesDict:
				clusterGenesDict[kmeans.labels_[i]] = [symbolList[i]]
			elif kmeans.labels_[i] in clusterGenesDict:
				clusterGenesDict[kmeans.labels_[i]].append(symbolList[i])
			lineOut = str(symbolList[i]) + '\t' + str(kmeans.labels_[i]) + '\t' + '\t'.join(str(zScoredListOfLists[i]).strip('[]').split(',')) + '\n'
			outFile.write(lineOut)
		
		print('\nDone!\n\nResults are written in', outPutFileName)
		
		## Close the cluster results file
		outFile.close()
				
		print('\nNow drawing plots for each cluster')
		
		## Draw the dot-plot for each cluster
		
		print('Genes by kmeans cluster: ',clusterGenesDict.keys())	
		
		#plt.boxplot(clusterArray)
		
		i = 1
		
		plt.figure(figsize=(20,k*4))
		plt.title(plotTitle)
		
		for clusterNumber in clusterGenesDict:
			print(clusterNumber)
			currentGeneList = clusterGenesDict[clusterNumber]
			#print(currentGeneList)
			currentPlotTitle = 'Cluster ' + str(clusterNumber) + ' (n=' + str(len(currentGeneList)) + ')'
			plt.subplot(k, 0, i)
			for gene in currentGeneList:
				tpmVals = geneValuesDict[gene]
				tpmArray = np.array(tpmVals)
				y = stats.zscore(tpmArray, ddof=df)
				## Setup the x positions, 1 for each sample, then add in random stagger so that they aren't all on the same overlapping x-position
				gene[x] = [random.uniform(1-staggerAmt,1+staggerAmt), random.uniform(2-staggerAmt,2+staggerAmt), random.uniform(3-staggerAmt,3+staggerAmt), random.uniform(4-staggerAmt,4+staggerAmt), random.uniform(5-staggerAmt,5+staggerAmt), random.uniform(6-staggerAmt,6+staggerAmt)]
				#plt.plot(x, y, 'o', alpha=0.5)
				plt.boxplot(y)
			plt.xticks([1, 2, 3, 4, 5, 6], ('100ng-1','100ng-2','100ng-3','0.1ng-1','0.1ng-2','0.1ng-3'))
			plt.ylabel(currentPlotTitle)
			#plt.xlim(0.5, 18.5)
			#i = i + 1
		plt.savefig(plotFileName)
		plt.clf()	
		print('\nThe resulting plots are saved in', plotFileName)
