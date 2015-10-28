#!/bin/bash

#PBS -V
#PBS -m ae
#PBS -M twaddlac@gmail.com
#PBS -S /bin/bash
#PBS -j oe -N partition -l walltime=168:00:00,nodes=1:ppn=20,mem=256gb

cd $path

module purge
module load khmer/intel/1.3

cd $path

#########################################
############### Functions ###############
#########################################

run_idba (){

	reads=$1
	job_name=$2
	
	name=$(basename $reads .fasta)

	echo \
		"module load idba && \
		idba_ud \
		--pre_correction \
		--num_threads 20 \
		-r $path/$name.fasta \
		-o $path/$name.idba.d" \
		| qsub -j oe -N $job_name -l walltime=72:00:00,nodes=1:ppn=20,mem=128gb

}

run_blastn (){

	reads=$1
	output=$2
	depend=$3
	job_name=$4

	echo \
		"module load blast+ && \
		blastn \
		-query $path/$reads \
		-out $path/$output \
		-outfmt 5 \
		-evalue 0.00001 \
		-max_target_seqs 20 \
		-culling_limit 2 \
		-num_threads 12 \
		-db /scratch/at120/db/blast/nt/nt" \
		| qsub -j oe -N $job_name -W depend=afterany:$depend -l walltime=168:00:00,nodes=1:ppn=12,mem=24gb

}

run_blastx (){

	reads=$1
	output=$2
	depend=$3
	job_name=$4

	echo \
		"module load blast+ && \
		blastx \
		-query $path/$reads \
		-out $path/$output \
		-outfmt 5 \
		-evalue 0.00001 \
		-max_target_seqs 20 \
		-culling_limit 2 \
		-num_threads 12 \
		-db /scratch/at120/db/blast/nr/nr" \
		| qsub -j oe -N $job_name -W depend=afterany:$depend -l walltime=168:00:00,nodes=1:ppn=12,mem=24gb

}
#########################################
############### Functions ###############
#########################################

echo "Filter below"

filter-below-abund.py $fastq.normC5k20.kh *.fastq.gz 

for i in *below; do
	mv $i $i.fastq;
done

echo "Do Partition" 

do-partition.py -k 32 -x 64e9 --threads 20 kak *.kak.fastq.gz.below.fastq

echo "Extract Partitions" 

extract-partitions.py -X 1000000 kak *.part 

echo "Extract Paired Reads" 

for i in kak*fq
do 
	extract-paired-reads.py $i 
	name=$(basename $i .fq)
	mv $i.se $name.diginorm.se.fastq
	mv $i.pe $name.diginorm.pe.fastq
	fastq-to-fasta.py $name.diginorm.pe.fastq > $name.diginorm.pe.fasta
	
	partitioned_idba_id=$(run_idba $name.diginorm.pe.fasta partitioned_idba_$i)
	echo $partitioned_idba_id >> ids.txt
	echo "Submitted partitioned idba"

	blastn_partitioned_idba_id=$(run_blastn $name.diginorm.pe.idba.d/contig.fa $name.diginorm.pe.idba.d/contig.blastn.nt.xml $partitioned_idba_id blastn_partitioned_idba_$name)
	echo $blastn_partitioned_idba_id >> ids.txt
	echo "Submitted partitioned blastn idba"
:<<'END'
	blastx_partitioned_idba_id=$(run_blastx $name.diginorm.pe.idba.d/contig.fa $name.diginorm.pe.idba.d/contig.blastx.nr.xml $partitioned_idba_id blastx_partitioned_idba_$name)
	echo $blastx_partitioned_idba_id >> ids.txt
	echo "Submitted partitioned blastx idba"
END
	rm $i
	bzip2 -9 $name.diginorm.se.fastq $name.diginorm.pe.fastq
	
done

:<<'END'
#This is to repopulate the reads
echo "sweep reads" 

sweep-reads3.py -x 64e9 kak.group*.fq $fastq.unconc.interleaved.fastq 

echo "extract paired reads"

for i in kak*.sweep3
do 
	extract-paired-reads.py $i 
	name=$(basename $i .fq.sweep3)
	mv $i.se $name.diginorm.inflate.nodn.se.fastq
	mv $i.pe $name.diginorm.inflate.nodn.pe.fastq
	fastq-to-fasta.py $name.diginorm.inflate.nodn.pe.fastq > $name.diginorm.inflate.nodn.pe.fasta
	
	inflated_idba_id=$(run_idba $name.diginorm.inflate.nodn.pe.fasta inflated_idba_$i)
	echo $inflated_idba_id >> ids.txt
	echo "Submitted inflated idba"

	blastn_inflated_idba_id=$(run_blastn $name.diginorm.inflate.nodn.pe.idba.d/contig.fa $name.diginorm.inflate.nodn.pe.idba.d/contig.blastn.nt.xml $inflated_idba_id blastn_inflated_idba_$name)
	echo $blastn_inflated_idba_id >> ids.txt
	echo "Submitted inflated blastn idba"

	blastx_inflated_idba_id=$(run_blastx $name.diginorm.inflate.nodn.pe.idba.d/contig.fa $name.diginorm.inflate.nodn.pe.idba.d/contig.blastx.nr.xml $inflated_idba_id blastx_inflated_idba_$name)
	echo $blastx_inflated_idba_id >> ids.txt
	echo "Submitted inflated blastx idba"
	
done
END
