#!/bin/bash 

# This script is to generate a PCA and SFS from data with/wihout paralogs

# Name of the job 
#SBATCH -J 0008

#Email myself
#SBATCH --mail-user=maccampbell@ucdavis.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#Job specs

#SBATCH --partition=high

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=2-12:00:00


#SBATCH -o $HOME/spineflower/0008/job-%j.stdout


##Using 0007/nonparalogous-contigs.tsv as a region file let's generate a PCA and a SFS for sanity with and without the flagged regions


## SFS example command
#./angsd -bam bam.filelist -doSaf 1 -out smallFolded -anc  chimpHg19.fa -GL 2 -minMapQ 10 -minQ 20
srun $HOME/angsd/angsd -P 16 -bam 1005_paralogs/spineflower_80.bamlist -doSaf 1 -out 0008/with-paralogs -anc 1003_stacks_catalog.fa -GL 2 -minMapQ 10 -minQ 20
srun $HOME/angsd/angsd -P 16 -bam 1005_paralogs/spineflower_80.bamlist -doSaf 1 -out 0008/with-out-paralogs -anc 1003_stacks_catalog.fa -GL 2 -minMapQ 10 -minQ 20 \
-rf 0007/nonparalogous-contigs.tsv


## Beagle file generation 187 samples
srun $HOME/angsd/angsd -P 16 -bam 1005_paralogs/spineflower_80.bamlist -out 0008/with-paralogs-pca -minInd 168 -GL 1 -doGlf 2  -doMajorMinor 1 -doMaf 2 \
-SNP_pval 1e-6 -minMapQ 10 -minQ 20 

srun $HOME/angsd/angsd -P 16 -bam 1005_paralogs/spineflower_80.bamlist -out 0008/with-out-paralogs-pca -minInd 168 -GL 1 -doGlf 2  -doMajorMinor 1 -doMaf 2 \
-SNP_pval 1e-6 -minMapQ 10 -minQ 20 -rf 0007/nonparalogous-contigs.tsv