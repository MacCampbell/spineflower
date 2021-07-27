#!/bin/bash -l
#SBATCH -t48:00:00

gunzip 1003_stacks/catalog.fa.gz
samtools faidx 1003_stacks/catalog.fa 
mkdir 1005_paralogs 

ls ${PWD}/1004_alignments/*_sorted_proper_rmdup_ss80k.bam > 1005_paralogs/spineflower_80.bamlist

angsd -bam 1005_paralogs/spineflower_80.bamlist -out 1005_paralogs/spineflower -ref 1003_stacks/catalog.fa -GL 2 -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 10 -minQ 20

gunzip 1005_paralogs/spineflower*.gz

cut -d$'\t' -f1-2  1005_paralogs/spineflower.mafs | sed 1d > 1005_paralogs/spineflower.snp.pos

samtools mpileup -b 1005_paralogs/spineflower_80.bamlist -l 1005_paralogs/spineflower.snp.pos -f 1003_stacks/catalog.fa > 1005_paralogs/spineflower.depth

ngsParalog calcLR -infile 1005_paralogs/spineflower.depth > 1005_paralogs/spineflower.paralogs

#awk '($5 > '100')' 1005_paralogs/spineflower.paralogs | cut -c1-7 | uniq > 1005_paralogs/spineflower.paralogs.list_100
#
#awk '($5 < '0')' 1005_paralogs/spineflower.paralogs | cut -c1-7 | uniq > 1005_paralogs/spineflower.paralogs.list_all
#
#grep '>' 1003_stacks/catalog.fa | cut -c2- | grep -v -f 1005_paralogs/spineflower.paralogs.list_100 | sed 's/$/:/' > 1005_paralogs/spineflower.loci_100
#
#grep '>' 1003_stacks/catalog.fa | cut -c2- | grep -v -f 1005_paralogs/spineflower.paralogs.list_all | sed 's/$/:/' > 1005_paralogs/spineflower.loci_all
#
