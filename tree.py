#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

import os, sys, time
import argparse
import dendropy

def get_parser():

	parser = argparse.ArgumentParser(description='Create tree from mash matrix')
	
	parser.add_argument('-t', action="store", dest='MATRIX', 
				type=str, required=True, help='mash matrix (REQUIRED)')
	
	parser.add_argument('-k', action="store", dest='OUTPUT',
				type=str, default='output', help='output newick name (default:output)')    
		
	parser.add_argument('--NJ', dest='NJ', action='store_true',
				help='use neighbour joining algorithm (default:UPGMA)', default=False)

	return parser

def make_nj_tree(mash_matrix):

	#Création de l'objet matrice à partir de la matrice au format tsv  obtenu avec Mash
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(mash_matrix), delimiter="\t")
	#Création de l'arbre avec la méthode neighbour-Joining
	tree = pdm.nj_tree()

	return tree

def make_upgma_tree(mash_matrix):

	#Création de l'objet matrice à partir de la matrice au format tsv  obtenu avec Mash
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(mash_matrix), delimiter="\t")
	#Création de l'arbre avec la méthode upgma
	tree = pdm.upgma_tree()

	return tree


def make_reroot_tree(tree):
	#Enracinement de l'arbre
	tree.reroot_at_midpoint(update_bipartitions=False, suppress_unifurcations=False)

	return tree

def write_reroot_tree(reroot_tree, output_file_name):
	output_file = open(output_file_name, 'w')
	output_file.write(reroot_tree.as_string("newick"))

	output_file.close()

	os.system("sed -i 's/\[&R\] //g' " + output_file_name)


#main function    
def main():

##################### gets arguments #####################
	parser=get_parser()

	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
    
	Arguments=parser.parse_args()

	if Arguments.NJ :
		tree = make_nj_tree(Arguments.MATRIX)
	else :
		tree = make_upgma_tree(Arguments.MATRIX)    

	reroot_tree = make_reroot_tree(tree)
	
	output_file = write_reroot_tree(reroot_tree, Arguments.OUTPUT + '.nwk')

if __name__ == "__main__":
    main()
