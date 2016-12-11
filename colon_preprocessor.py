def colon_translator(): 
	translator = {}
	x = open('../../../../../Desktop/RNASeq_Files/colon_biomart.txt','r+')
	header = x.readline()
	for line in x: 
		ensembl,probe = line.split(',')
		probe = probe.replace('\n','')
		translator[probe] = ensembl
	return translator