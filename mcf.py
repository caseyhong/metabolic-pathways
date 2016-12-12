import json 
from gpr_mapping import * #aggregate, parse_gpr_mapping 
from geneid_mapping import * #prepareEntrez,ensemblToEntrez
from diff_paths import * #difference(x,y), createNetwork
from make_hypergraph import * #make_hypergraph 
from igraph import * 
from colon_preprocessor import * #colon_translator
from parse_kegg_genes import * 

def createGeneIdMapping(): 
	entrezDict = prepareEntrez()
	ensemblEntrezDict = ensemblToEntrez(entrezDict)
	return ensemblEntrezDict

def getPatientData(fileName,pid,ensemblEntrezDict): 
	"""
	Returns the patient's expression data specified by pid 
	INPUTS: 
		fileName - file specifying expression data matrix 
		pid - patient id specifying which patient from file 
		ensemblEntrezDict - created from createGeneIdMapping(), ensembl id to entrez id dictionary 
	OUTPUTS: 
		patientName - patient (sample) name 
		data - patient's exprsesion val dictionary {Entrez Gene Id: Expression Value}
	"""
	data = {}
	if 'colon' in fileName: 
		colon_dict = colon_translator()

	with open(fileName,'r+') as expData: 
		header = expData.readline().split('\t')
		patientName = header[pid]
		numPatients = len(header) - 1 
		if pid > numPatients: 
			print 'Invalid patient index!'
			return 

		for row in expData: 
			row = row.split('\t')
			gene = row[0]
			if 'colon' in fileName and gene in colon_dict.keys(): 
				gene = colon_dict[gene.replace('\"','')]
			exp = row[pid]

			# Convert from Ensembl -> Entrez Gene Id 
			if gene in ensemblEntrezDict.keys(): 
				geneids = ensemblEntrezDict[gene]
				for g in geneids: #Several Gene Id Transcript Versions 
					data[g] = exp 

	return patientName,data


	#List of patients , rxnid, value 
def patientRxnMappings(RNASeqFileName,patientIds,recon1RxnMappings,ensemblEntrezDict): 
	patientToRxns = {}

	for pid in patientIds: 
		patientName, pdata = getPatientData(RNASeqFileName,pid,ensemblEntrezDict)
		patientToRxns[patientName] = {}

		for rxnId in recon1RxnMappings.keys(): 
			rxn = recon1RxnMappings[rxnId]
			rxnVal = aggregate(rxn,pdata,pdata.keys())
			patientToRxns[patientName][rxnId] = rxnVal

	return patientToRxns

def generateDigraph(patientToRxns): 

	#Initialize Hypergraph for all patients for tumor type 
	hg = make_hypergraph()

	#Add edge weights for each patient 
	for patient in patientToRxns.keys(): 
		rxn_expVal = patientToRxns[patient]
		hg = add_to_HG(hg,rxn_expVal)

	#Invert the weights 
	dg,inverted_weights = invert_weights(hg)
	#print inverted_weights
	
	#Create final digraph
	#dg = convert_to_DG(hg)
	#dg = convert_to_DG(hg,inverted_weights)
	return dg 

def mcf(lung_file,lung_ids,colon_file,colon_ids):
	"""
	mcf returns the final differential expressed graph for a particular cancer file 
	INPUTS: 
		lung_file - filename to parse for lung data 
		lung_ids - index values of which samples to take from file 
		colon_file - filename to parse for colon data 
		colon_ids - index values of which samples to take from file 
	"""
	ensemblEntrezDict = createGeneIdMapping()
	recon1RxnMappings = parse_gpr_mapping('./RECON1.json')
	rxnMetaboliteMapping = get_metabolite_associations('./RECON1.json')
	metabolite_pairs = get_metabolite_pairs(rxnMetaboliteMapping)

	# For Lung Patients 
	print 'Processing Lung DG'
	lung_patient_rxns = patientRxnMappings(lung_file,lung_ids,recon1RxnMappings,ensemblEntrezDict)
	lung_dg = generateDigraph(lung_patient_rxns)

	# For Colon Patients 
	print 'Processing Colon DG'
	colon_patient_rxns = patientRxnMappings(colon_file,colon_ids,recon1RxnMappings,ensemblEntrezDict)
	colon_dg = generateDigraph(colon_patient_rxns)

	# Perform graph subtraction 
	print 'Performing Graph Subtraction'
	diff_g = difference(lung_dg,colon_dg)
	#print diff_g
	return diff_g

 
# source - a list containing the source vertex IDs which should be included in the result. If None, all vertices will be considered.
# target - a list containing the target vertex IDs which should be included in the result. If None, all vertices will be considered.
# weights - a list containing the edge weights. It can also be an attribute name (edge weights are retrieved from the given attribute) or None (all edges have equal weight).
# mode - the type of shortest paths to be used for the calculation in directed graphs. OUT means only outgoing, IN means only incoming paths. ALL means to consider the directed graph as an undirected one.
# Returns:
# the shortest path lengths for given vertices in a matrix

def generate_features(diff_graph): 
	print 'Generating features for differential graph'
	src_target = parse()
	seeds = [] 
	for s in src_target: 
		if len(src_target[s]) >= 3: 
			seeds.append(s)

	for seed in seeds: 
		targets = src_target[seed]
		shortest_paths = run_dijkstra(diff_graph,seed,targets)
		print "For seed: "+ seed +"+, shortest path is: "
		print shortest_paths
		print 'targets: '+ ','.join(targets)
		print '###########################################'


def run_dijkstra(graph,source,targets):
	print 'Running Bellman Ford on Differential Graph'
	x = graph.shortest_paths_dijkstra(source=source, target=targets, weights='weight', mode=OUT)
	return x 



if __name__ == '__main__':
	lung_file = '../../../../../Desktop/RNASeq_Files/GSE81089_FPKM_cufflinks_nslc.tsv'
	colon_file = '../../../../../Desktop/RNASeq_Files/GSE41258_series_matrix_colon.txt'
	colon_30 = [26, 28, 30, 31, 33, 35, 37, 39, 40, 42, 43, 45, 47, 48, 49, 51, 54, 57, 59, 61, 62, 64, 67, 69, 73, 75, 77, 79, 81, 84]
	colon_50 = [26, 28, 30, 31, 33, 35, 37, 39, 40, 42, 43, 45, 47, 48, 49, 51, 54, 57, 59, 61, 62, 64, 67, 69, 73, 75, 77, 79, 81, 84, 86, 87, 89, 90, 91, 93, 94, 96, 97, 98, 99, 102, 104, 107, 108, 109, 110, 112, 114, 115]
	lung_30 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
	lung_50 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 45, 46, 47, 48, 49, 50, 51, 52]
	diff_g = mcf(lung_file,colon_30,colon_file,lung_30)
	generate_features(diff_g)
	#x = run_dijkstra(diff_g)




















