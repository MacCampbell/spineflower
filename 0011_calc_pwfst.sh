#!/bin/bash
#SBATCH -o 0011/0011_calc_pwst-%j.out
#SBATCH -t 48:00


mkdir 0011

./1000_scripts/list_to_pwcomps.pl meta/poplist > 0011/pw.list

wc=$(wc -l 0011/pw.list | awk '{print $1}')
x=1
while [ $x -le $wc ]
do
	string="sed -n ${x}p 0011/pw.list"
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1,$2}')
	set -- $var
	pop1=$1
	pop2=$2

echo "#!/bin/bash
#SBATCH --job-name=fst${x}
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --partition=high
#SBATCH --time=48:00:00
#SBATCH --output=0011/${pop1}_${pop2}-%j.slurmout
#############################################

realSFS 0011/${pop1}.saf.idx 0011/${pop2}.saf.idx > 0011/${pop1}_${pop2}.2dsfs
#Rscript scripts/plot2DSFS_2D.R 0011/${pop1}_${pop2}.2dsfs \$(wc -l 0011/${pop1}.bamlist | awk '{print \$1}') \$(wc -l 0011/${pop2}.bamlist | awk '{print \$1}') $pop1 $pop2
#Rscript scripts/plot2dSFS.R 0011/${pop1}_${pop2}.2dsfs 0011/${pop1}_${pop2}.2dsfs.pdf ${pop1} ${pop2}
realSFS fst index 0011/${pop1}.saf.idx 0011/${pop2}.saf.idx -sfs 0011/${pop1}_${pop2}.2dsfs -fstout 0011/${pop1}_${pop2} 
realSFS fst stats 0011/${pop1}_${pop2}.fst.idx > 0011/${pop1}_-_${pop2}.fst.stats
" > fst_${pop1}_${pop2}.sh

sbatch fst_${pop1}_${pop2}.sh
rm fst_${pop1}_${pop2}.sh

x=$(( $x + 1 ))
done


#hold til complete
#grep "" 0011/*fst.stats | sed 's/.fst.stats:/      /' | sed 's:0011/::' |  tr "_-_" "\t" | awk '{print $1"_"$2"_"$3"\t"$4"_"$5"_"$6"\t"$8}' > 0011/pwfst.all

#Changed by Mac to accomodate different delimiter
grep "" 0011/*fst.stats | sed 's/.fst.stats:/\t/' | sed 's:0011/::' | perl -pe 's/_-_/\t/g' > 0011/pwfst.all

1000_scripts/pwlist_to_matrix.pl 0011/pwfst.all > 0011/pwfst.fstmatrix
