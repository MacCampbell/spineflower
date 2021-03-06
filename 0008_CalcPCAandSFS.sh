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
#SBATCH --mem=32G


#SBATCH -o /home/maccamp/spineflower/0008/0008.stdout


##Using 0007/nonparalogous-contigs.tsv as a region file let's generate a PCA and a SFS for sanity with and without the flagged regions


## SFS example command
#./angsd -bam bam.filelist -doSaf 1 -out smallFolded -anc  chimpHg19.fa -GL 2 -minMapQ 10 -minQ 20
#srun -p high -t 48:00:00 --mem=16G --nodes=2 --ntasks=6 $HOME/angsd/angsd -P 12 -bam 1005_paralogs/spineflower_80.bamlist -doSaf 1 -out 0008/with-paralogs -anc 1003_stacks/catalog.fa -GL 2 -minMapQ 10 -minQ 20

#srun -p high -t 48:00:00 --mem=16G --nodes=2 $HOME/angsd/angsd -P 12 -bam bamlists/spineflower.bamlist -doSaf 1 -out 0008/with-paralogs -anc 1003_stacks/catalog.fa -GL 2 -minMapQ 10 -minQ 20 > 0008/with-paralogs-saf.stdout 2>0008/with-paralogs-saf.stderr &
#srun -p high -t 48:00:00 --mem=16G --nodes=2 $HOME/angsd/angsd -P 12 -bam bamlists/spineflower.bamlist -doSaf 1 -out 0008/with-out-paralogs -anc 1003_stacks/catalog.fa -GL 2 -minMapQ 10 -minQ 20 -rf 0007/nonparalogous-contigs.tsv > 0008/with-out-paralogs-saf.stdout 2>0008/with-out-paralogs-saf.stderr &


#srun $HOME/angsd/misc/realSFS 0008/with-paralogs.saf.idx -maxIter 100 -P 16 -fold 1 >  0008/with-paralogs.sfs
#srun $HOME/angsd/misc/realSFS 0008/with-out-paralogs.saf.idx -maxIter 100 -P 16 -fold 1 >  0008/with-out-paralogs.sfs

#srun -p bigmemh --mem-per-cpu=10G --mem=32G -t 12:00:00 $HOME/angsd/misc/realSFS 0008/with-out-paralogs.saf.idx -maxIter 100 -P 1 -fold 1 >  0008/with-out-paralogs.sfs
#srun  -p bigmemh --mem-per-cpu=10G --mem=32G -t 12:00:00 $HOME/angsd/misc/realSFS 0008/with-paralogs.saf.idx -maxIter 100 -P 1 -fold 1 >  0008/with-paralogs.sfs

## Beagle file generation 187 samples
#srun -p high -t 24:00:00 $HOME/angsd/angsd -P 16 -bam bamlists/spineflower.bamlist -out 0008/with-paralogs-pca -minInd 168 -GL 1 -doGlf 2  -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 10 -minQ 20 > 0008/with-paralogs-pca.stdout 2> 0008/with-paralogs-pca.stderr &

#srun  -p bigmemh --mem-per-cpu=10G --mem=32G -t 12:00:00 $HOME/angsd/angsd -P 4 -bam bamlists/spineflower.bamlist -out 0008/with-out-paralogs-pca -minInd 168 -GL 1 -doGlf 2  -doMajorMinor 1 -doMaf 2 -SNP_pval 1e-6 -minMapQ 10 -minQ 20 -rf 0007/nonparalogous-contigs.tsv > 0008/with-out-paralogs-pca.stdout 2> 0008/with-out-paralogs-pca.stderr &

## Generating cov matrix
#srun python $HOME/pcangsd/pcangsd.py -beagle 0008/with-paralogs-pca.beagle.gz -kinship -admix -o 0008/with-paralogs-pca -threads 10
#srun -p high -t 01:00:00 python $HOME/pcangsd/pcangsd.py -beagle 0008/with-out-paralogs-pca.beagle.gz -kinship -admix -o 0008/with-out-paralogs-pca -threads 10
#srun -p high -t 01:00:00 python $HOME/pcangsd/pcangsd.py -beagle 0008/with-out-paralogs-pca.beagle.gz -relate 0008/with-out-paralogs-pca.kinship.npy -admix -o 0008/with-out-paralogs-pca-relate -threads 10

