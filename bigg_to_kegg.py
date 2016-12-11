import numpy as np
import pandas as pd
import ast

def bigg_to_kegg(filename):
	df = pd.read_table(filename)
	bigg_to_kegg = {}
	for index, row in df.iterrows():
		links = ast.literal_eval(row['database_links'])
		if 'KEGG Reaction' in links.keys():
			kegg_id = links['KEGG Reaction'][0]['id']
			print links['KEGG Reaction'][0]['link']
			bigg_to_kegg[row['bigg_id']] = kegg_id
	return bigg_to_kegg

if __name__ == '__main__':
	bigg_to_kegg('./bigg_models_reactions.txt')