#!/usr/bin/python
# -*- coding: iso-8859-1 -*-
import os, sys, time
import argparse
import copy
from decimal import Decimal
#import numpy as np
#from scipy.cluster import hierarchy 
#import matplotlib.pyplot as plt


def get_parser():
	"""
	Parse arguments
	@return arguments list
	@rtype: parser object
	"""

	parser = argparse.ArgumentParser(description='Create NJ tree from ARTwork assemblies')

	parser.add_argument('-f', action="store", dest='FASTA', 
						type=str, required=True, nargs='+', help='FASTA files, more than 2 (REQUIRED)')

	parser.add_argument('-o', action="store", dest='output', 
						type=str, default='output', help='output name (default:output)')

	parser.add_argument('-T', action="store", dest='nbThreads', 
						type=int, default='1', help='maximum number of threads to use (default:1)')

	return parser



def make_dist_tab(fasta_querie, fasta_targets, nbThreads):

	outputFile_name = ''.join(fasta_querie.split("_")) + "_dit.tsv"
	outputFile = open(outputFile_name,'w')

	os.system("mash dist -t -p " + str(nbThreads) + " " + fasta_querie + \
		" " + " ".join(fasta_targets) + " > " + outputFile_name)

	outputFile.close()

	return outputFile_name


def mash_dist_loop(fasta_files, nbThreads):


	distTable_files = []

	for fasta_file in fasta_files :

		fasta_targets = copy.copy(fasta_files)
		fasta_targets.remove(fasta_file)

		distTable_files.append(make_dist_tab(fasta_file, fasta_targets, nbThreads))

	return 	distTable_files


def make_dist_matrix(distTable_files):
	
	dico_result = {}

	for distTable_file in distTable_files :

		distFile = open(distTable_file, 'r')
		lines = distFile.readlines()
		distFile.close()
		os.system("rm " + distTable_file)


		for line in lines :

			line = line.rstrip()

			if line[0]=="#" :
				seqId = line.split('\t')[1].replace("_assembly.fasta",'')
				dicoSeq = {}

			else :	
				seqName = line.split('\t')[0].replace("_assembly.fasta",'')
				dist = line.split('\t')[1]
				dicoSeq[seqName] = dist	

		dico_result[seqId] = dicoSeq		


	return dico_result			


def write_dist_matrix(dicoDist, matrixFileName):
	
	seqNames = dicoDist.keys()
	new_dicoDist = {}

	matrixFile = open(matrixFileName, 'w')

	matrixFile.write('\t' + '\t'.join(seqNames) + '\n')

	for seqName in dicoDist : 

		matrixFile.write(seqName)
		new_dicoDist[seqName] = []

		for name in seqNames :
			if name not in dicoDist[seqName] :
				matrixFile.write('\t0')
				new_dicoDist[seqName].append(0)
			else:
				matrixFile.write('\t' + str(Decimal(dicoDist[seqName][name])))	
				new_dicoDist[seqName].append(float(Decimal(dicoDist[seqName][name])))

		matrixFile.write('\n')		

	matrixFile.close()
	return new_dicoDist

'''
def matrix_to_array(matrix):

	matrix_list = []

	for element in matrix :
		matrix_list.append(matrix[element])

	array = np.array(matrix_list)
	return array
'''





#main function	
def main():

	##################### gets arguments #####################
	parser=get_parser()
	
	#print parser.help if no arguments
	if len(sys.argv)==1:
		parser.print_help()
		sys.exit(1)
	
	Arguments=parser.parse_args()

	distTable_files = mash_dist_loop(Arguments.FASTA, Arguments.nbThreads)
	distMatrix = make_dist_matrix(distTable_files)
	distMatrix = write_dist_matrix(distMatrix, Arguments.output + '.tsv')

if __name__ == "__main__":
	main()	            		


