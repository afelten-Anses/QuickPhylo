README
======
Authors: Pauline Barbet, Arnaud Felten

Afiliation: [Food Safety Laboratory - ANSES Maisons Alfort (France)](https://www.anses.fr/en/content/laboratory-food-safety-maisons-alfort-and-boulogne-sur-mer)

You can find the latest version of the tool at [https://github.com/afelten-Anses/QuickPhylo](https://github.com/afelten-Anses/QuickPhylo)


Workflow
========
This workflow aims to faslty build a root tree with Mash and dendropy. To make the root tree, two methods are available : UPGMA and nNeighbourg-Joining.

The script takes in input, assembly, reads (compressed or not) and sketch files. The Mash sketch function, sketch each fasta and fastq files into msh files.

Then, the Mash dist function create a distance matrix based on jaccard index at tsv format with previous sketch files and files already sketch.

Finally dendropy produce in output a newick file with the root tree. 

![](workflow.JPG?raw=true "script workflow")

Dependencies
============

The script has been developped with python 2.7 (tested with 2.7.12)

## Expernal dependencies

* [Mash](https://github.com/marbl/Mash/blob/master/INSTALL.txt) tested with 2.0
* [Dendropy](https://www.dendropy.org/) tested with 4.3.0


Parameters
==========

## Parameters

* -i : tsv file containing paths to reads, asssemly or/and sketched files, more than 2 (REQUIRED)
* -o : output tsv matrix name  (default:output)
* -T : maximum number of threads to use (default:1)
* -k : k-mer size for sketching (default:21)
* -s : sketch size = number of k-mer (default:1000)
* -e : output newick tree name (default:output)

## Options

* --S : supress sketch files (defaut:False)
* --NJ : use neighbour joinning algorithm (defaut:UPGMA)

Test
====

After installing Mash and dendropy you can test the script whith the command lines :

	mkdir test
	cd test
	python MashtoDist.py -T nbThreads -i input.tsv -o matrix_name -e tree_name

