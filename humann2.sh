#!/bin/bash
#PBS -m ae
#PBS -M twaddlac@gmail.com
#PBS -N humann2
#PBS -l walltime=72:00:00,nodes=1:ppn=12,mem=64gb
#PBS -joe

module load humann2
module load bowtie2
module load diamond
module load metaphlan

cd "/scratch/at120/virome-pipeline/lhmp/biobakery/humann2"

humann2 --input lhmp.unaligned.fasta --output lhmp.unaligned.humann2 --threads 12