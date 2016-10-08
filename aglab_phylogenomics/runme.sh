#!/bin/bash

# Find the orthologs based on the OrthoDB mapping results, contained in a *.og file.
# *.og files usually contain the following data.
#
################################################################################
# 1	103372:0012a9	0	3482	4047	7528	22834.2	40.2578	0
# 1	104421:001a20	0	2587	2584	5170	23542.6	40.965	0
# 1	1049336:004382	0	3014	181	3194	3741.4	14.6722	0
# 1	1049336:00465a	0	813	450	1262	6567.3	27.0259	0
# 1	1049336:0046fc	7	1774	3	1776	10107	-1	0
# 1	1049336:0046fd	0	5629	21	5649	8641.6	31.689	0
# 1	112268:001a32	0	106	1272	1377	7322.2	10.7332	0
# 1	1155016:00265a	0	1253	4	1256	1279.2	6.03966	0
# 1	1155016:003741	0	688	21	708	728.4	3.76045	3.2e-212
################################################################################
#
# where the first field contains the name of the cluster
#       the second field contains the gene name. The semicolon ":" splits this field
#          into two parts; the first is the taxon ID of the organism and the second
#          is a gene identifier used in OrthoDB.
#       the other seven fields are irrelevant to the analysis that follows
#
# The second input file, "species.list", is a list of taxon IDs corresponding to the
# species we are interested in. It looks like this:
#
#########
# 6669
# 136037
# 7029
# 121225
# 13037
# OTAUR
# APLAN
########
#
# where each of the above IDs should match the first part of the second field
# of the previous *.og file.
find_orthologs_from_mapping_data.pl orthoDB_mapping.og species.list > Aglab_orthology.txt

# After running the above command, ortholog clusters that contain single-copy genes in all
# of the species in "species.list" will be output in separate files, one for each cluster.
# These files contain only the gene names and have the file extension ".names".

# create a directory and move these files into the directory
mkdir sc_fasta

mv *names sc_fasta

cd sc_fasta

# Generate a fasta file that will contain the sequence for each of the genes that
# appear in these files ("all_species.faa") and run the following script in order
# to extract the sequences
for x in *names; do extract_specific_fasta_seqs_v3.pl $x all_species.faa > `basename $x .names`.fs; done

# create a directory for each fasta file and move it into that directory
for x in *fs; do mkdir $x.raxml; mv $x $x.raxml; mv `basename $x .fs`.names $x.raxml; done

# Run a phylogenetic analysis for each ortholog cluster. This script internally uses:
# "muscle" for multiple sequence alignment,
# "trimal" for trimming the alignments
# "RAxML" for running the phylogenetic analysis.
# See the script for more details.
for x in *raxml; do cd $x; create_tree_from_fasta_file.sh.pl `ls *.fs` 6669; cd ../; done

# Prepare the data for the concatenated phylogenomic analysis.
# More specifically, filter out any ortholog clusters for which the individual phylogenetic analysis
# resulted in a tree that has an average bootstrap value of zero. Concatenate the remaining alignments
# in the file "concat_genes.fs.aln.trimmed", in the "BS_ge_0.01" directory.
filter_trees.pl -p raxml -b 0.01 -o BS_ge_0.01

cd BS_ge_0.01/

# Use the seqret tool (from the EMBOSS suite) to convert this fasta alignment to phylip.
seqret -sformat fasta -sequence concat_genes.fs.aln.trimmed -osformat phylip -outseq concat_genes.fs.aln.trimmed.phy

# Start the RAxML-based phylogenetic analysis on the concatenated alignment
raxmlHPC-PTHREADS-SSE3 -T 12 -f a -x 12345 -p 12345 -N 100 -m PROTGAMMAJTT -s concat_genes.fs.aln.trimmed.phy -w `pwd` -n concat_genes.fs.aln.trimmed.phy.tre -o 6669 > stdout_stderr 2>&1 
