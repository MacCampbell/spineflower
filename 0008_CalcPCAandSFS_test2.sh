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

#SBATCH --partition=bigmemm
#SBATCH --mem-per-cpu=10G

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=2-12:00:00

##Setting up test dataset
## head -n 40 1005_paralogs/spineflower_80.bamlist > 0008/40.bamlist

#SBATCH -o /home/maccamp/spineflower/0008/test2.stdout


##Using 0007/nonparalogous-contigs.tsv as a region file let's generate a PCA and a SFS for sanity with and without the flagged regions


## SFS example command
#./angsd -bam bam.filelist -doSaf 1 -out smallFolded -anc  chimpHg19.fa -GL 2 -minMapQ 10 -minQ 20
#srun -p high -t 12:00:00 --mem=16G --nodes=2 --ntasks=6 $HOME/angsd/angsd -P 12 -bam 1005_paralogs/spineflower_80.bamlist -doSaf 1 -out 0008/with-paralogs -anc 1003_stacks/catalog.fa -GL 2 -minMapQ 10 -minQ 20
srun $HOME/angsd/angsd -P 16 -bam 0008/40.bamlist -doSaf 1 -out 0008/test-with-paralogs -anc 1003_stacks/catalog.fa -GL 2 -minMapQ 10 -minQ 20
#srun $HOME/angsd/angsd -P 16 -bam 0008/40.bamlist -doSaf 1 -out 0008/test-with-out-paralogs -anc 1003_stacks/catalog.fa -GL 2 -minMapQ 10 -minQ 20 \
#-rf 0007/nonparalogous-contigs.tsv

srun $HOME/angsd/misc/realSFS 0008/test-with-paralogs.saf.idx -maxIter 100 -P 16 -fold 1 >  0008/test-with-paralogs.sfs
#srun $HOME/angsd/misc/realSFS 0008/test-with-out-paralogs.saf.idx -maxIter 100 -P 16 >  0008/test-with-out-paralogs.sfs

## Beagle file generation 187 samples (for test file, need to change)
srun $HOME/angsd/angsd -P 16 -bam 0008/40.bamlist -out 0008/test-with-paralogs-pca -minInd 168 -GL 1 -doGlf 2  -doMajorMinor 1 -doMaf 2 \
-SNP_pval 1e-6 -minMapQ 10 -minQ 20 

#srun $HOME/angsd/angsd -P 16 -bam 0008/40.bamlist -out 0008/test-with-out-paralogs-pca -minInd 168 -GL 1 -doGlf 2  -doMajorMinor 1 -doMaf 2 \
#-SNP_pval 1e-6 -minMapQ 10 -minQ 20 -rf 0007/nonparalogous-contigs.tsv

## Generating cov matrix
srun python $HOME/pcangsd/pcangsd.py -beagle 0008/test-with-paralogs-pca.beagle.gz -kinship -admix -o 0008/test-with-paralogs-pca -threads 10
#srun python $HOME/pcangsd/pcangsd.py -beagle 0008/test-with-out-paralogs-pca.beagle.gz -kinship -admix -o 0008/test-with-out-paralogs-pca -threads 10

