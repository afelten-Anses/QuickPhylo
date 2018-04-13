#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys, time
import argparse
import copy
import dendropy
from decimal import Decimal
from dendropy.calculate import treecompare


def get_parser():
	#Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description='compute distance between trees ')

	parser.add_argument('-i', action="store", dest='nwk', 
				type=str, required=True, nargs='+', help='newick files, minimum 2 (REQUIRED)')

	parser.add_argument('-o', action="store", dest='output', 
						type=str, default='output', help='output tsv name (default:output)')

	parser.add_argument('--U', dest='URFD', action='store_true',
						help='compute unweighted robinson foulds distance (default:weighted robinson foulds distance)',
						default=False)

	return parser

def get_weighted_robinson_foulds_distance(input_file, tree):
	
	tns = dendropy.TaxonNamespace()

	tree_1 = dendropy.Tree.get(
    		path=input_file,
    		schema="newick",
    		taxon_namespace=tns)

	tree_2 = dendropy.Tree.get(
    		path=tree,
    		schema="newick",
    		taxon_namespace=tns)

	distance = treecompare.weighted_robinson_foulds_distance(tree_1, tree_2)

	return distance


def get_unweighted_robinson_foulds_distance(input_file, tree):
	
	tns = dendropy.TaxonNamespace()

	tree_1 = dendropy.Tree.get(
    		path=input_file,
    		schema="newick",
    		taxon_namespace=tns)

	tree_2 = dendropy.Tree.get(
    		path=tree,
    		schema="newick",
    		taxon_namespace=tns)

	distance = treecompare.symmetric_difference(tree_1, tree_2)

	return distance


def make_distance_matrix(input_files, distance):

	dico_result = {}

	for input_file in input_files :
		trees = copy.copy(input_files)
		trees.remove(input_file)
		
		dico_dist = {}
		tree_name_1 = input_file.split('/')[-1].split('.')[0]
		
		for tree in trees :
			tree_name_2 = tree.split('/')[-1].split('.')[0]

			if distance :
				dico_dist[tree_name_2] = get_unweighted_robinson_foulds_distance(input_file, tree)
			else :
				dico_dist[tree_name_2] = get_weighted_robinson_foulds_distance(input_file, tree)

		dico_result[tree_name_1] = dico_dist

	return dico_result


def write_distance_matrix(dico_result, output):

	tree_names = dico_result.keys()

	matrix_file = open(output + '.tsv', 'w')

	matrix_file.write('\t' + '\t'.join(tree_names) + '\n')

	for tree_name in dico_result :
		matrix_file.write(tree_name)

		for name in tree_names :
			if name not in dico_result[tree_name] : # if name == tree_name
				matrix_file.write('\t0')

			else:
				matrix_file.write('\t' + str(dico_result[tree_name][name]))

		matrix_file.write('\n')		

	matrix_file.close()


def main():

	##################### gets arguments #####################

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	# mettre tout les arguments dans la variable Argument
	Arguments=parser.parse_args()

	matrix = make_distance_matrix(Arguments.nwk, Arguments.URFD)
	matrix = write_distance_matrix(matrix, Arguments.output)

    

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	            		           		