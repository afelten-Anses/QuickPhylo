#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys, time
import argparse
import copy
from decimal import Decimal
import dendropy


##################################
#####  Arguments definition  #####
##################################


def get_parser():
	#Fonction permettant de pourvoir demander des arguments

	parser = argparse.ArgumentParser(description='Create mash matrix and rerooted taxonomic tree')

	#parser.add_argument('-f', action="store", dest='FASTA', 
	#					type=str, required=True, nargs='+', help='FASTA files, more than 2 (REQUIRED)')

	parser.add_argument('-i', action="store", dest='input', 
						type=str, required=True, help='TXT file containing path to reads or/and assembly, more than 2 (REQUIRED)')

	parser.add_argument('-o', action="store", dest='output', 
						type=str, default='output', help='output tsv name (default:output)')

	parser.add_argument('-T', action="store", dest='nbThreads', 
						type=int, default='1', help='maximum number of threads to use (default:1)')

	parser.add_argument('--S', action='store_true', dest='sketch',
						help='suppress sketch files (default:True)', default=True)
 
	#parser.add_argument('-t', action="store", dest='MATRIX', 
	#			type=str, required=True, help='mash matrix (REQUIRED)')
	
	parser.add_argument('-k', action="store", dest='OUTPUT',
						type=str, default='output', help='output newick name (default:output)')    
		
	parser.add_argument('--NJ', dest='NJ', action='store_true',
						help='use neighbour joining algorithm (default:UPGMA)', default=False)

	return parser


###########################################
#####  Create mash matrix (functions) #####
###########################################


def make_sketch_files(input_file, nbThreads):
	#Fonction qui sketch chaque fichier que se soient des assemblages ou des reads
	intFile = open (input_file, 'r')
	paths = intFile.readlines()
	intFile.close()

	os.system("mkdir sketch")
	os.system("mkdir concat_reads")

	sketch_files = []

	for path in paths :

		path = path.rstrip()
		path = path.split("\t")
		if len(path) == 1 :

			sketch_file_name = "sketch/" + path[0].split('/')[-1].split('.')[0].replace("_assembly",'_sketch')

			os.system("mash sketch -p " + str(nbThreads) + " -o " + sketch_file_name + " " + path[0])

			sketch_files.append(sketch_file_name + ".msh")

		else :
	
			cat_file_name = "concat_reads/" + path[0].split('/')[-1].replace("_R1",'_cat_reads')

			os.system("cat " + path[0] + " " + path[1] + " > " + cat_file_name)

			sketch_file_name = "sketch/" + cat_file_name.split('/')[-1].split('.')[0].replace("_cat_reads",'_sketch')

			os.system("mash sketch -p " + str(nbThreads) + " -r -o " + sketch_file_name + " " + cat_file_name)
			
			sketch_files.append(sketch_file_name + ".msh")

	return sketch_files
	

def make_dist_tab(fasta_querie, fasta_targets, nbThreads):
	# lance mash pour comparer fasta_querie avec tout les assemblages de le liste fasta_targets
	# la sortie de mash est écrite dans outputFile_name

	outputFile_name = ''.join(fasta_querie.split("_")) + "_dist.tsv"
	outputFile = open(outputFile_name,'w')

	os.system("mash dist -t -p " + str(nbThreads) + " " + fasta_querie + \
		" " + " ".join(fasta_targets) + " > " + outputFile_name)

	outputFile.close()

	return outputFile_name


def mash_dist_loop(sketch_files, nbThreads):
	# fonction qui lance make_dist_tab pour tous les génomes
	# on récupére une liste de nom de fichiers comprenant toutes les sortie de mash 

	distTable_files = []

	for sketch_file in sketch_files :

		sketch_targets = copy.copy(sketch_files)
		sketch_targets.remove(sketch_file)

		distTable_files.append(make_dist_tab(sketch_file, sketch_targets, nbThreads))

	return 	distTable_files


def make_dist_matrix(distTable_files):
	# ouvre tout les fichiers résultats de mash et stock les informations dans un dictionnaire
	# ce distionnaire contient en clef [1] les id des génomes et en valeur un autre dictionnaire
	# ce second dictionnaire contient en clef tout les id des génomes sauf ceux de la clef [1]
	# et en valeur la distance mash 
	
	dico_result = {}

	for distTable_file in distTable_files :

		distFile = open(distTable_file, 'r')
		lines = distFile.readlines()
		distFile.close()
		os.system("rm " + distTable_file)


		for line in lines :

			line = line.rstrip() # supprime retour chariot

			if line[0]=="#" :
				seqId = line.split('\t')[1].split('/')[-1].replace("_assembly.fasta",'')
				seqId = seqId.replace("_cat_reads.fastq.gz",'')
				dicoSeq = {}

			else :	
				seqName = line.split('\t')[0].split('/')[-1].replace("_assembly.fasta",'')
				seqName = seqName.replace("_cat_reads.fastq.gz",'')
				dist = line.split('\t')[1]
				dicoSeq[seqName] = dist	

		dico_result[seqId] = dicoSeq		


	return dico_result




def write_dist_matrix(dicoDist, matrixFileName):
	# "transpose le dictionnaire dicoDist en une matrice dans le fichier matrixFileName"
	
	seqNames = dicoDist.keys()
	new_dicoDist = {}

	matrixFile = open(matrixFileName, 'w')

	# écriture du header de la matrice
	matrixFile.write('\t' + '\t'.join(seqNames) + '\n')

	for seqName in dicoDist : 

		matrixFile.write(seqName)
		new_dicoDist[seqName] = []

		for name in seqNames :
			# si name == seqName alors distance = 0
			if name not in dicoDist[seqName] : # if name == seqName
				matrixFile.write('\t0')
				new_dicoDist[seqName].append(0)
			else:
				matrixFile.write('\t' + str(Decimal(dicoDist[seqName][name])))	
				
				new_dicoDist[seqName].append(float(Decimal(dicoDist[seqName][name])))

		matrixFile.write('\n')		

	matrixFile.close()
	return new_dicoDist


########################################################
#####  Create rerooted taxonomic tree (functions)  #####
########################################################

def make_nj_tree(mash_matrix):
	# Création de l'objet matrice à partir de la matrice au format tsv  obtenu avec Mash
	# Création de l'arbre avec la méthode neighbour-Joining

	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(mash_matrix), delimiter="\t")
	tree = pdm.nj_tree()

	return tree

def make_upgma_tree(mash_matrix):
	#Création de l'objet matrice à partir de la matrice au format tsv  obtenu avec Mash
	#Création de l'arbre avec la méthode upgma

	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=open(mash_matrix), delimiter="\t")
	tree = pdm.upgma_tree()

	return tree


def make_reroot_tree(tree):
	#Enracinement de l'arbre

	tree.reroot_at_midpoint(update_bipartitions=False, suppress_unifurcations=False)

	return tree

def write_reroot_tree(reroot_tree, tree_file_name):
	# Transpose l'arbre dans le fichier tree_file_name

	tree_file = open(tree_file_name, 'w')
	tree_file.write(reroot_tree.as_string("newick"))
	tree_file.close()

	os.system("sed -i 's/\[&R\] //g' " + tree_file_name)


###########################
#####  Main function  #####
###########################


def main():

	##################### gets arguments #####################

	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	# mettre tout les arguments dans la variable Argument
	Arguments=parser.parse_args()

    #####################  Create mash matrix  #####################

	sketch_files_name = make_sketch_files(Arguments.input, Arguments.nbThreads)
	distTable_files = mash_dist_loop(sketch_files_name, Arguments.nbThreads)
	distMatrix = make_dist_matrix(distTable_files)
	distMatrix = write_dist_matrix(distMatrix, Arguments.output + '.tsv')

	#####################  Create rerooted taxonomic tree  #####################

	if Arguments.NJ :
		tree = make_nj_tree(Arguments.output + '.tsv')
	else :
		tree = make_upgma_tree(Arguments.output + '.tsv')

	reroot_tree = make_reroot_tree(tree)
	
	tree_file = write_reroot_tree(reroot_tree, Arguments.OUTPUT + '.nwk')

# lancer la fonction main()  au lancement du script
if __name__ == "__main__":
	main()	            		           		