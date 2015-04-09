#!/usr/bin/env python3

"""
Script to generate entrez gene IDs from HGNC gene symbols in the format
NCBIGene:1234
"""

import mygene
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', type=str, help='Path to input file')
parser.add_argument('--output', '-o', type=str, help='Path to output file')

args = parser.parse_args()
mg = mygene.MyGeneInfo()
symbol_list = []
output_file = open(args.output, 'w')

# Get symbols from mapping file

with open(args.input, 'rt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter="\t")
    for row in reader:
        (gene_id, gene_symbol) = row
        symbol_list.append(gene_symbol)

print("Querying mygene")
mygene_results = mg.querymany(symbol_list, scopes='symbol', fields='entrezgene', species='human', returnall=True)

mygene_output = mygene_results['out']
mygene_dict = {}

for result in mygene_output:
    if 'entrezgene' in result:
        mygene_dict[result['query']] = "NCBIGene:{0}".format(result['entrezgene'])

with open(args.input, 'rt') as tsvfile:
    reader = csv.reader(tsvfile, delimiter="\t")
    for row in reader:
        (gene_id, gene_symbol) = row
        if gene_symbol in mygene_dict:
            row.append(mygene_dict[gene_symbol])
            line = "\t".join(row)
            line += '\n'
            output_file.write(line)
        else:
            line = "\t".join(row)
            line += '\n'
            output_file.write(line)

output_file.close()