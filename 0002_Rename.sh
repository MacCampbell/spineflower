#!/bin/bash
#SBATCH -t 2-00:00:00 
#SBATCH -o slurm_outs/rename-%j.out

mkdir 1002_samples/

input="reference/BMAG049_AGTCAA.metadata"
wc=$(wc -l $input | awk '{print $1}')
x=1
while [ $x -le $wc ]
do

        string="sed -n ${x}p $input"
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1,$2,$3,$4}')
        set -- $var
        c1=$1
        c2=$2
        c3=$3
	c4=$4

	mv 1001_sequence/BMAG050_${c1}_GG${c2}TGCAGG_RA.fastq 1002_samples/${c3}.1.fq
	mv 1001_sequence/BMAG050_${c1}_GG${c2}TGCAGG_RB.fastq 1002_samples/${c3}.2.fq
        x=$(( $x + 1 ))

done


input="reference/BMAG050_GTCCGC.metadata"
wc=$(wc -l $input | awk '{print $1}')
x=1
while [ $x -le $wc ]
do

        string="sed -n ${x}p $input"
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1,$2,$3,$4}')
        set -- $var
        c1=$1
        c2=$2
        c3=$3
        c4=$4

	mv 1001_sequence/BMAG050_${c1}_GG${c2}TGCAGG_RA.fastq 1002_samples/${c3}.1.fq
        mv 1001_sequence/BMAG050_${c1}_GG${c2}TGCAGG_RB.fastq 1002_samples/${c3}.2.fq
        x=$(( $x + 1 ))

done

input="reference/BMAG050_GTGAAA.metadata"
wc=$(wc -l $input | awk '{print $1}')
x=1
while [ $x -le $wc ]
do

        string="sed -n ${x}p $input"
        str=$($string)

        var=$(echo $str | awk -F"\t" '{print $1,$2,$3,$4}')
        set -- $var
        c1=$1
        c2=$2
        c3=$3
        c4=$4

        mv 1001_sequence/BMAG050_${c1}_GG${c2}TGCAGG_RA.fastq 1002_samples/${c3}.1.fq
        mv 1001_sequence/BMAG050_${c1}_GG${c2}TGCAGG_RB.fastq 1002_samples/${c3}.2.fq
        x=$(( $x + 1 ))

done


wc -l 1002_samples/*.fq | sort -nr > 1002_samples/wc.txt
