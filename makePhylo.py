#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys, time
import argparse
import dendropy



def get_parser():


	parser = argparse.ArgumentParser(description='Create tree from mash matrix')

	parser.add_argument('-t', action="store", dest='MATRIX', 
						type=str, required=True, help='mash matrix (REQUIRED)')

	parser.add_argument('-o', action="store", dest='OUTPUT', 
						type=str, default='output', help='output newick name (default:output)')

	parser.add_argument('--NJ', dest='NJ', action='store_true', 
		help='use neighbour joining algorithm (default:UPGMA)', default=False)

	return parser



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
		#Création de l'objet matrice à partir de la matrice au format tsv  obtenu avec Mash
		pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(Arguments.MATRIX),delimiter="\t")
		#Création de l'arbre avec la méthode neighbour-Joining
		tree = pdm.nj_tree()
	else :	
		#Création de l'objet matrice à partir de la matrice au format tsv  obtenu avec Mash
		pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(Arguments.MATRIX),delimiter="\t")
		#Création de l'arbre avec la méthode upgma
		tree = pdm.upgma_tree()
		
	#Enracinement de l'arbre
	tree.reroot_at_midpoint(update_bipartitions=False)
	output_file = open(Arguments.OUTPUT + ".nwk",'w')
	output_file.write(tree.as_string("newick"))
	output_file.close()

	'''
	command = "sumtrees.py -F newick --root-target-at-midpoint \
	--suppress-annotations --decimals=0 --percentages " + newick + \
	" > " + Arguments.OUTPUT + ".nwk"
	os.system(command)
	'''


if __name__ == "__main__":
	main()	            		


