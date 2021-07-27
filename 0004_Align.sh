#!/bin/bash -l
#SBATCH -t 48:00:00
#SBATCH -o slurm_outs/align-%j.out
bwa index 1003_stacks/catalog.fa.gz
mkdir 1004_alignments/

ls 1002_samples/*1.fq | sed 's/.1.fq//' | sed 's:1002_samples/::' > 1002_samples/samplelist

wc=$(wc -l 1002_samples/samplelist | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p 1002_samples/samplelist" 
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1,$2,$3}')   
	set -- $var
	c1=$1
	c2=$2
	c3=$3

echo "#!/bin/bash
#SBATCH -o slurm_outs/align_${c1}-%j.out
bwa mem 1003_stacks/catalog.fa.gz 1002_samples/${c1}.1.fq 1002_samples/${c1}.2.fq > 1004_alignments/${c1}.sam
samtools view -bS 1004_alignments/${c1}.sam > 1004_alignments/${c1}.bam
samtools sort  1004_alignments/${c1}.bam > 1004_alignments/${c1}_sorted.bam
samtools view -b -f 0x2 1004_alignments/${c1}_sorted.bam > 1004_alignments/${c1}_sorted_proper.bam
samtools rmdup 1004_alignments/${c1}_sorted_proper.bam 1004_alignments/${c1}_sorted_proper_rmdup.bam
samtools index 1004_alignments/${c1}_sorted_proper_rmdup.bam 1004_alignments/${c1}_sorted_proper_rmdup.bam.bai
reads=\$(samtools view -c 1004_alignments/${c1}_sorted.bam)
ppalign=\$(samtools view -c 1004_alignments/${c1}_sorted_proper.bam)
rmdup=\$(samtools view -c 1004_alignments/${c1}_sorted_proper_rmdup.bam)
echo \"\${reads},\${ppalign},\${rmdup}\" > 1004_alignments/${c1}.stats" > aln_${x}.sh
sbatch -p high -t 24:00:00 aln_${x}.sh
rm aln_${x}.sh
sleep 5s
	x=$(( $x + 1 ))
done

