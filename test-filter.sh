#!/bin/bash

# Add MEGAN
# Add bbnorm
# Add Velvet
# Add Spades
# blast after partitioning and after inflating

path="/scratch/at120/virome-pipeline/virome-analysis"
fastq="test"

cd $path

bowtie2=\
$(echo \
	"module load bowtie2 && \
	bowtie2 \
	-p 12 \
	--very-sensitive-local \
	--un-conc $path/$fastq.unconc.fastq \
	-x /scratch/at120/db/kazima-contamination-seqs/kazima-contamination-seqs.fasta \
	-1 $path/$fastq.r1.fastq \
	-2 $path/$fastq.r2.fastq \
	-S $path/$fastq.sam \
	&& rm $path/$fastq.sam"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N bowtie2 -l walltime=72:00:00,nodes=1:ppn=12,mem=12gb)
echo $bowtie2 > ids.txt
echo "Submitted Mapping"

convertBT2=\
$(echo \
	"module load khmer && \
	fastq-to-fasta.py $path/$fastq.unconc.1.fastq > $path/$fastq.unconc.fasta && \
	fastq-to-fasta.py $path/$fastq.unconc.2.fastq >> $path/$fastq.unconc.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N convert -W depend=afterok:$bowtie2 -l walltime=72:00:00,nodes=1:ppn=1,mem=4gb)
echo $convertBT2 >> ids.txt
echo "Submitted Converting to Fasta"

megablastFilter=\
$(echo \
	"module load blast+ && \
	blastn \
	-query $path/$fastq.unconc.fasta \
	-out $path/$fastq.unconc.megablast.tsv \
	-outfmt 6 \
	-evalue 0.00001 \
	-max_target_seqs 1 \
	-culling_limit 2 \
	-num_threads 12 \
	-db /scratch/at120/db/kazima-contamination-seqs/kazima-contamination-seqs.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N megablastFilter -W depend=afterok:$convertBT2 -l walltime=72:00:00,nodes=1:ppn=12,mem=12gb)
echo $megablastFilter >> ids.txt
echo "Submitted Megablast Filter"

unmappedMegablast=\
$(echo \
	"module load biopython && \
	python /scratch/at120/virome-pipeline/get-unmapped.py $path/$fastq.unconc.megablast.tsv $path/$fastq.unconc.fasta $path/$fastq.unconc.megablast.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N blastnFilter -W depend=afterok:$megablastFilter -l walltime=72:00:00,nodes=1:ppn=1,mem=12gb)
echo $unmappedMegablast >> ids.txt
echo "Submitted Unmapped Megablast"

blastnFilter=\
$(echo \
	"module load blast+ && \
	blastn \
	-task blastn \
	-query $path/$fastq.unconc.megablast.fasta \
	-out $path/$fastq.unconc.megablast.blastn.tsv \
	-outfmt 6 \
	-evalue 0.00001 \
	-max_target_seqs 1 \
	-culling_limit 2 \
	-num_threads 12 \
	-db /scratch/at120/db/kazima-contamination-seqs/kazima-contamination-seqs.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N blastnFilter -W depend=afterok:$unmappedMegablast -l walltime=72:00:00,nodes=1:ppn=12,mem=12gb)
echo $blastnFilter >> ids.txt
echo "Submitted blastn Filter"

unmappedBlastn=\
$(echo \
	"module load biopython && \
	python $path/get-unmapped.py $path/$fastq.unconc.megablast.blastn.tsv $path/$fastq.unconc.megablast.fasta $path/$fastq.unconc.megablast.blastn.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N unmappedBlastn -W depend=afterok:$blastnFilter -l walltime=72:00:00,nodes=1:ppn=1,mem=12gb)
echo $unmappedBlastn >> ids.txt
echo "Submitted Unmapped Blastn"

extractPairedReads=\
$(echo \
	"module load khmer && \
	extract-paired-reads.py $path/$fastq.unconc.megablast.blastn.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N extractPairedReads -W depend=afterok:$unmappedBlastn -l walltime=72:00:00,nodes=1:ppn=1,mem=12gb)
echo $extractPairedReads >> ids.txt
echo "Submitted Extract Paired Reads"
