#!/bin/bash -l
#SBATCH -o slurm_outs/subsample-%j.out
#SBATCH -t 48:00:00 

ls 1004_alignments/*sam | sed 's:1004_alignments/::' | sed 's/.sam//' > indlist

wc=$(wc -l indlist | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p indlist" 
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1}')   
	set -- $var
	c1=$1

echo "#!/bin/bash
	count=\$(samtools view -c 1004_alignments/${c1}_sorted_proper_rmdup.bam)

        if [ 60000 -le \$count ]
        then
                frac=\$(bc -l <<< 60000/\$count)
        samtools view -bs \$frac 1004_alignments/${c1}_sorted_proper_rmdup.bam > 1004_alignments/${c1}_sorted_proper_rmdup_ss80k.bam
                samtools index 1004_alignments/${c1}_sorted_proper_rmdup_ss80k.bam 1004_alignments/${c1}_sorted_proper_rmdup_ss80k.bam.bai
	fi

ss80k_count=\$(samtools view -c 1004_alignments/${c1}_sorted_proper_rmdup_ss80k.bam)

echo \"\${count},\${ss80k_count}\" > 1004_alignments/${c1}_ss.counts" > ss1_${x}.sh

sbatch -p high -t 24:00:00 ss1_${x}.sh
rm ss1_${x}.sh
sleep 2s


	x=$(( $x + 1 ))

done

rm indlist
