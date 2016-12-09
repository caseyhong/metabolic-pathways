import json
import re
from collections import defaultdict
from collections import Counter
from ast import literal_eval
import pyparsing

thecontent = pyparsing.Word(pyparsing.alphanums) | '+' | '-'
parens = pyparsing.nestedExpr('(', ')', content=thecontent)

def get_gene_set(filename):
	relevant_genes = set()
	bools = ['or', 'and']
	with open(filename) as json_file:
		data = json.load(json_file)
		for rxn in data["reactions"]:
			if len(rxn["gene_reaction_rule"]) > 0:
				rule = rxn["gene_reaction_rule"]
				words = rule.split()
				genes = [w for w in words if w not in bools]
				for g in genes:
					if g not in relevant_genes:
						relevant_genes.add(g)
	return relevant_genes

def parse_gpr_mapping(filename):
	parsed_mappings = []
	with open(filename) as json_file:
		data = json.load(json_file)
		for rxn in data["reactions"]:
			if len(rxn["gene_reaction_rule"]) > 0:
				rule = rxn["gene_reaction_rule"]
				words = rule.split()
				hashed_rule = []
				for word in words:
					if word[0] != '(' and word[-1] != ')':
						h = word.replace('_','ZZZZZZ')
						hashed_rule.append(h.encode('ascii'))
					elif word[0] == '(':
						h = word.replace('_','ZZZZZZ')
						hashed_rule.append(h.encode('ascii'))
					else:
						h = word.replace('_','ZZZZZZ')
						hashed_rule.append(h.encode('ascii'))
						hashed_rule.append(')')
				hashed = ''.join(hashed_rule)
				rule_clean = '(' + hashed + ')'
				rule_cleaner = rule_clean.replace('and', '+')
				rule_cleanest = rule_cleaner.replace('or', '-')
				res = parens.parseString(rule_cleanest)
				parsed = res.asList()
				parsed_mappings.append(parsed[0])
	# print parsed_mappings
	return parsed_mappings

def aggregate(mapping, expression, genes):
	s = 0
	for i in xrange(len(mapping)):
		if i == 0:
			if mapping[i] in genes:
				s += expression[mapping[i]]
			else:
				s += aggregate(mapping[i], expression, genes)

		if i%2==0 and i > 0:
			next = mapping[i]
			operator = mapping[i-1]

			if next in genes:
				if operator == '+':
					s += min(s, expression[next])
				else:
					s += max(s, expression[next])
			else:
				if operator == '+':
					s += min(s, aggregate(next, expression, genes))
				else:
					s += max(s, aggregate(next, expression, genes))
	return s


if __name__ == '__main__':
	mappings = parse_gpr_mapping('./RECON1.json')
	print aggregate(mappings[0], expression, genes)
	# print genes_of_interest

