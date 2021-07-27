#!/bin/bash
#SBATCH -o slurm_outs/0001_Demultiplex-%j.out
#SBATCH -t 48:00 
sbatch -t 48:00:00 1000_scripts/run_BestRadSplit.sh 1001_sequence/plates/BMAG049_R1_AGTCAA.fastq 1001_sequence/plates/BMAG049_R3_AGTCAA.fastq 1001_sequence/BMAG049_AGTCAA_

sbatch -t 48:00:00 1000_scripts/run_BestRadSplit.sh 1001_sequence/plates/BMAG050_R1_GTCCGC.fastq 1001_sequence/plates/BMAG050_R3_GTCCGC.fastq 1001_sequence/BMAG050_GTCCGC_

sbatch -t 48:00:00 1000_scripts/run_BestRadSplit.sh 1001_sequence/plates/BMAG050_R1_GTGAAA.fastq 1001_sequence/plates/BMAG050_R3_GTGAAA.fastq 1001_sequence/BMAG050_GTGAAA_

sbatch -t 48:00:00 1000_scripts/run_BestRadSplit.sh 1001_sequence/plates/BMAG058_R1_TTAGGC.fastq 1001_sequence/plates/BMAG058_R3_TTAGGC.fastq 1001_sequence/BMAG058_TTAGGC_

