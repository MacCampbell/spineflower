#!/bin/bash
#SBATCH -o 0011/0011_calc_sfs-%j.out
#SBATCH -t 48:00


mkdir 0011
wc=$(wc -l meta/poplist | awk '{print $1}')
x=1
while [ $x -le $wc ]
do
	string="sed -n ${x}p metadata/poplist2"
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1}')
	set -- $var
	pop=$1

echo "#!/bin/bash
#SBATCH --job-name=sfs${x}
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --partition=high
#SBATCH --time=48:00:00
#SBATCH --output=0011/${pop}-%j.slurmout

##############################################

#ls \${PWD}/1004_alignments/${pop}*_ss40k.bam | sed 's/_ss40k//' > 0011/${pop}.bamlist
#This isn't going to work, let's make our own pop.bamlists

nInd=\$(wc -l bamlist/${pop}.bamlist | awk '{print \$1}')
mInd=\$((\${nInd}/2))

#############################################
#Getting sites together (base) maccamp@farm:~/spineflower/0009$ cat selection.sites | perl -pe 's/_/:/g' > sites


angsd -b 0011/${pop}.bamlist -anc 1003_stacks/catalog.fa -ref 1003_stacks/catalog.fa -out 0011/${pop} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -baq 2 -GL 1 -doMajorMinor 1 -doMaf 1 -minInd $mInd -nind $nInd -minMapQ 10 -minQ 20 -doSaf 2 -nThreads 8 -rf 0009/sites
realSFS 0011/${pop}.saf.idx > 0011/${pop}.sfs

Rscript scripts/plotSFS.R 0011/${pop}.sfs
" > sfs_${pop}.sh

sbatch sfs_${pop}.sh
rm sfs_${pop}.sh

x=$(( $x + 1 ))
done

