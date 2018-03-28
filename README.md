QuickPhylo
==========

QuickPhylo create trees with UPGMA or NJ methods from assemblies, reads or sketch files with Mash and Dendropy module from python.

Workflow
========
This workflow aims to build a root tree with mash and dendropy. To make the root tree, two methods are available : UPGMA and nNeighbourg-Joining . The script can takes assembly, reads (compressed or not) or sketch files in input. Mash sketch, sketch each fasta and fastq files into msh files. Then mash dist create a distance matrix based on jaccard index at tsv format with previous sketch files and files already sketch. Finally dendropy produce in output a newick file with the root tree. 

![](workflow.jpg?raw=true "workflow")

Dependencies
============

The script has been devlopped with python 2.7 (tested with 2.7.12)

## Expernal dependencies

* [Mash](https://github.com/marbl/Mash/blob/master/INSTALL.txt) tested with 2.0
* [Dendropy](https://www.dendropy.org/) tested with 4.3.0


Parameters
==========

## Parameters

* -i : tsv file containing paths to reads or/and asssemly or/and sketched files, more than 2 (REQUIRED)
* -o : output tsv matrix name  (default:output)
* -T : maximum number of threads to use (default:1)
* -k : k-mer size for sketching (default:21)
* -s : sketch size = number of k-mer (default:1000)
* -e : output newick tree name (default:output)

## Options

* --S : supress sketch files (defaut:False)
* --NJ : use neighbour joinning algorithm (defaut:UPGMA)

