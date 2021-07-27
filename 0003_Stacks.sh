#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -c 8 
#SBATCH -o slurm_outs/stacks-%j.out

module load stacks/2.53
mkdir 1003_stacks/
mkdir 1003_stacks/samples

grep ".1.fq" 1002_samples/wc.txt | head -n 150 | tail -n 100  | awk '{print $2}' | sed 's:.1.fq:	:' | sed 's:1002_samples/::' > 1003_stacks/uStacks.samplelist

grep -f 1003_stacks/uStacks.samplelist reference/pop.metadata > 1003_stacks/uStacks.popmap

input="1003_stacks/uStacks.samplelist"
wc=$(wc -l $input | awk '{print $1}')
x=1
while [ $x -le $wc ]
do

        string="sed -n ${x}p $input"
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1}')
        set -- $var
        sample=$1

	process_radtags -1 1002_samples/${sample}.1.fq -2 1002_samples/${sample}.2.fq -o 1003_stacks/samples/ -q -c --disable-rad-check
        ustacks -f 1003_stacks/samples/${sample}.1.1.fq -o 1003_stacks/ -i $x --name $sample -M 4 -p 8 


        x=$(( $x + 1 ))

done

cstacks -n 6 -P 1003_stacks/ -M 1003_stacks/uStacks.popmap -p 8

sstacks -P 1003_stacks/ -M 1003_stacks/uStacks.popmap -p 8 

tsv2bam -P 1003_stacks/ -M 1003_stacks/uStacks.popmap --pe-reads-dir 1003_stacks/samples/ -t 8 

gstacks -P 1003_stacks/ -M 1003_stacks/uStacks.popmap -t 8 

#populations -P 1003_stacks/ -M 1003_stacks/uStacks.popmap -r 0.65 --vcf --genepop --structure --fstats --hwe -t 8
