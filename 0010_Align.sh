#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -o slurm_outs/align-%j.out

#bwa index 1003_stacks/catalog.fa.gz
mkdir 0010

#ls data/concat/*A.fq | sed 's/_RA.fq//' | sed 's:data/concat/::' > 0010/samplelist

wc=$(wc -l 0010/samplelist | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p 0010/samplelist" 
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1,$2,$3}')   
	set -- $var
	c1=$1
	c2=$2
	c3=$3

echo "#!/bin/bash
#SBATCH -o 0010/align_${c1}-%j.out
bwa mem 1003_stacks/catalog.fa. data/concat/${c1}_RA.fq 1002_samples/${c1}_RB.fq > 0010/${c1}.sam
samtools view -bS 0010/${c1}.sam > 0010/${c1}.bam
samtools sort  0010/${c1}.bam > 0010/${c1}_sorted.bam
samtools view -b -f 0x2 0010/${c1}_sorted.bam > 0010/${c1}_sorted_proper.bam
samtools rmdup 0010/${c1}_sorted_proper.bam 0010/${c1}_sorted_proper_rmdup.bam
samtools index 0010/${c1}_sorted_proper_rmdup.bam 0010/${c1}_sorted_proper_rmdup.bam.bai
reads=\$(samtools view -c 0010/${c1}_sorted.bam)
ppalign=\$(samtools view -c 0010/${c1}_sorted_proper.bam)
rmdup=\$(samtools view -c 0010/${c1}_sorted_proper_rmdup.bam)
echo \"\${reads},\${ppalign},\${rmdup}\" > 0010/${c1}.stats" > aln_${x}.sh
sbatch -p high -t 24:00:00 aln_${x}.sh
#rm aln_${x}.sh
sleep 5s
	x=$(( $x + 1 ))
done

