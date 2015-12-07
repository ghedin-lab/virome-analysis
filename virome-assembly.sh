#!/bin/bash

# Add MEGAN
# Add bbnorm
# Add Velvet
# Add Spades
# blast after partitioning and after inflating

path="/scratch/at120/virome-pipeline/lhmp"
fastq="lhmp"	

cd $path
r
#########################################
############### Functions ###############
#########################################

run_idba (){

	reads=$1
	depend=$2
	job_name=$3
	mem=$4
	
	name=$(basename $reads .fasta)

	echo \
		"module load idba && \
		idba_ud \
		--pre_correction \
		--num_threads 20 \
		-r $path/$name.fasta \
		-o $path/$name.idba.d && \
		bzip2 -9 $path/$name.fasta" \
		| qsub -m ae -M twaddlac@gmail.com -j oe -N $job_name -W depend=afterok:$depend -l walltime=48:00:00,nodes=1:ppn=20,mem=$mem

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
		| qsub -m ae -M twaddlac@gmail.com -j oe -N $job_name -W depend=afterany:$depend -l walltime=168:00:00,nodes=1:ppn=12,mem=24gb

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
		| qsub -m ae -M twaddlac@gmail.com -j oe -N $job_name -W depend=afterany:$depend -l walltime=168:00:00,nodes=1:ppn=12,mem=24gb

}
#########################################
############### Functions ###############
#########################################

####################################
############### Main ###############
####################################

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
	| qsub -m ae -M twaddlac@gmail.com -j oe -N bowtie2 -l walltime=24:00:00,nodes=1:ppn=12,mem=12gb)
echo $bowtie2 > ids.txt
echo "Submitted Mapping"

convertBT2=\
$(echo \
	"module load khmer && \
	fastq-to-fasta.py $path/$fastq.unconc.1.fastq > $path/$fastq.unconc.fasta && \
	fastq-to-fasta.py $path/$fastq.unconc.2.fastq >> $path/$fastq.unconc.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N convert -W depend=afterok:$bowtie2 -l walltime=24:00:00,nodes=1:ppn=1,mem=4gb)
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
	| qsub -m ae -M twaddlac@gmail.com -j oe -N megablastFilter -W depend=afterok:$convertBT2 -l walltime=24:00:00,nodes=1:ppn=12,mem=12gb)
echo $megablastFilter >> ids.txt
echo "Submitted Megablast Filter"

unmappedMegablast=\
$(echo \
	"module load biopython && \
	python /scratch/at120/virome-pipeline/get-unmapped.py $path/$fastq.unconc.megablast.tsv $path/$fastq.unconc.fasta $path/$fastq.unconc.megablast.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N unmappedMegablast -W depend=afterok:$megablastFilter -l walltime=24:00:00,nodes=1:ppn=1,mem=12gb)
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
	| qsub -m ae -M twaddlac@gmail.com -j oe -N blastnFilter -W depend=afterok:$unmappedMegablast -l walltime=24:00:00,nodes=1:ppn=12,mem=12gb)
echo $blastnFilter >> ids.txt
echo "Submitted blastn Filter"

unmappedBlastn=\
$(echo \
	"module load biopython && \
	python get-unmapped.py $path/$fastq.unconc.megablast.blastn.tsv $path/$fastq.unconc.megablast.fasta $path/$fastq.unconc.megablast.blastn.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N unmappedBlastn -W depend=afterok:$blastnFilter -l walltime=24:00:00,nodes=1:ppn=1,mem=12gb)
echo $unmappedBlastn >> ids.txt
echo "Submitted Unmapped Blastn"

extractPairedReads=\
$(echo \
	"module load khmer && \
	extract-paired-reads.py $path/$fastq.unconc.megablast.blastn.fasta"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N extractPairedReads -W depend=afterok:$unmappedBlastn -l walltime=24:00:00,nodes=1:ppn=1,mem=12gb)
echo $extractPairedReads >> ids.txt
echo "Submitted Extract Paired Reads"


# repeatMasker=\
# 	$(echo \
# 		"module load repeat_masker && \
# 		-qq \
# 		-pa 12 \
# 		$path/$fastq.fasta"\
# 		| qsub -m ae -M twaddlac@gmail.com -j oe -N bowtie2 -l walltime=24:00:00,nodes=1:ppn=12,mem=12gb)
# echo $repeatMasker >> ids.txt
# echo "Submitted RepeatMasker"

tagCleaner=\
$(echo \
	"module load tagcleaner && \
	perl tagcleaner.pl \
	-predict \
	-64 \
	-cont \
	-fastq $path/$fastq"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N tagcleaner -W depend=afterok:$blastnFilter -l walltime=24:00:00,nodes=1:ppn=1,mem=12gb)
echo $tagCleaner >> ids.txt
echo "Submitted TagCleaner"

prinseq=\
$(echo \
	"module load prinseq && \
	prinseq-lite.pl \
	-fastq $path/$fastq.unconc.1.fastq \
	-fastq2 $path/$fastq.unconc.2.fastq \
	-out_format 3 \
	-out_good $path/$fastq.prinseq_good \
	-out_bad $path/$fastq.prinseq_bad \
	-lc_threshold 15 \
	-lc_method dust \
	#-graph_data $path/$fastq.prinseq.gd \
	#-graph_stats ld,gc,qd,ns,pt,ts,aq,de,da,sc,dn \
	-no_qual_header \
	-range_len 50,300 \
	-min_qual_mean 25 \
	-ns_max_p 10 \
	-derep 12 \
	-trim_qual_right 20 \
	-trim_qual_type min  \
	-trim_qual_window 1 \
	-trim_qual_step 1 \
	-trim_qual_rule lt \
	-stats_all > $path/$fastq.prinseq.stats.txt && \
	rm $path/$fastq.unconc.?.fastq" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N prinseq -W depend=afterok:$bowtie2 -l walltime=168:00:00,nodes=1:ppn=1,mem=128gb)
echo $prinseq > ids.txt
echo "Submitted Filtering"



:<<'END'
bowtie2=\
$(echo \
	"module load bowtie2 && \
	bowtie2 \
	-p 12 \
	--very-sensitive-local \
	--un-conc $path/$fastq.unconc.fastq \
	--un $path/$fastq.unaligned.fastq \
	-x /scratch/at120/db/kazima-contamination-seqs/kazima-contamination-seqs.fasta \
	-1 $path/$fastq.prinseq_good_1.fastq \
	-2 $path/$fastq.prinseq_good_2.fastq \
	-U $path/$fastq.prinseq_good_1_singletons.fastq,$path/$fastq.prinseq_good_2_singletons.fastq \
	-S $path/$fastq.sam \
	&& rm $path/$fastq.sam"\
	| qsub -m ae -M twaddlac@gmail.com -j oe -N bowtie2 -W depend=afterok:$prinseq -l walltime=24:00:00,nodes=1:ppn=12,mem=12gb)
echo $bowtie2 >> ids.txt
echo "Submitted Mapping"
END

convert=\
$(echo \
	"module load khmer && \
	fastq-to-fasta.py $path/$fastq.prinseq_good_1_singletons.fastq > $path/$fastq.unaligned.fasta && \
	fastq-to-fasta.py $path/$fastq.prinseq_good_2_singletons.fastq >> $path/$fastq.unaligned.fasta && \
	fastq-to-fasta.py $path/$fastq.prinseq_good_1.fastq >> $path/$fastq.unaligned.fasta && \
	fastq-to-fasta.py $path/$fastq.prinseq_good_2.fastq >> $path/$fastq.unaligned.fasta" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N convert -W depend=afterok:$prinseq -l walltime=24:00:00,nodes=1:ppn=1,mem=4gb)
echo $convert >> ids.txt
echo "Submitted All Converting to Fasta"


interleave_filtered=\
$(echo \
	"module load khmer && \
	interleave-reads.py \
	$path/$fastq.prinseq_good_?.fastq \
	> $path/$fastq.unconc.interleaved.fastq" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N interleave_filtered -W depend=afterok:$prinseq -l walltime=12:00:00,nodes=1:ppn=1,mem=4gb)
echo $interleave_filtered >> ids.txt
echo "Submitted interleaving filtered reads"


convert_interleaved=\
$(echo \
	"module load khmer && \
	fastq-to-fasta.py $path/$fastq.unconc.interleaved.fastq > $path/$fastq.unconc.interleaved.fasta" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N convert_interleaved -W depend=afterok:$interleave_filtered -l walltime=12:00:00,nodes=1:ppn=1,mem=4gb)
echo $convert_interleaved >> ids.txt
echo "Submitted convert interleaved"



first_idba=$(run_idba $fastq.unconc.interleaved.fasta $convert_interleaved first_idba 512gb)
echo $first_idba >> ids.txt
echo "Submitted Filtered Reads Assembly"

blastn_first_idba=$(run_blastn $fastq.unconc.interleaved.idba.d/contig.fa $fastq.unconc.interleaved.idba.d/contig.blastn.nt.xml $first_idba blastn_first_idba)
echo $blastn_first_idba >> ids.txt
echo "Submitted Blastn Filtered Reads Assembly"

blastx_first_idba=$(run_blastx $fastq.unconc.interleaved.idba.d/contig.fa $fastq.unconc.interleaved.idba.d/contig.blastx.nr.xml $first_idba blastx_first_idba)
echo $blastx_first_idba >> ids.txt
echo "Submitted Blastx Filtered Reads Assembly"

kraken=\
$(echo \
	"module load kraken && \
	kraken \
    --db /scratch/at120/db/kraken/standard \
    --threads 12 \
    --output $path/$fastq.kraken.out \
    --fasta-input \
    $path/$fastq.unaligned.fasta" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N kraken -W depend=afterok:$convert -l walltime=168:00:00,nodes=1:ppn=12,mem=64gb )
echo $kraken >> ids.txt
echo "Submitted Kraken"

:<<'END'
cdhit=\
$(echo \
	"module load cdhit && \
	cd-hit-est \
	-i $path/$fastq.unaligned.fasta \
	-o $path/$fastq.unaligned.clustered.fasta \
	-d 0 \
	-g 1 \
	-c 0.97 \
	-r 1 \
	-T 12 \
	-M 0" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N cdhit -W depend=afterok:$convert -l walltime=48:00:00,nodes=1:ppn=12,mem=256gb )
echo $cdhit >> ids.txt
echo "Submitted cdhit"

blastn_reads=$(run_blastn $fastq.unaligned.clustered.fasta $fastq.unaligned.clustered.blastn.nt.xml $cdhit blastn_reads)
echo $blastn_reads >> ids.txt
echo "Submitted clustered blastn reads"

blastx_reads=$(run_blastx $fastq.unaligned.clustered.fasta $fastq.unaligned.clustered.blastx.nr.xml $cdhit blastx_reads)
echo $blastx_reads >> ids.txt
echo "Submitted clustered blastx reads"
END

diginorm=$(qsub -v path=$path,fastq=$fastq -W depend=afterok:$interleave_filtered /scratch/at120/virome-pipeline/diginorm.sh)
echo $diginorm >> ids.txt
echo "Submitted Diginorm"

convert_normalized=\
$(echo \
	"module load khmer && \
	fastq-to-fasta.py $path/$fastq.unconc.interleaved.fastq.keep.abundfilt.pe.keep > $path/$fastq.unaligned.pe.kak.fasta" \
	| qsub -m ae -M twaddlac@gmail.com -j oe -N convert_normalized -W depend=afterok:$diginorm -l walltime=12:00:00,nodes=1:ppn=1,mem=4gb)
echo $convert_normalized >> ids.txt
echo "Submitted convert normalized"

#normalized_idba=$(runba $fastq.unaligned.pe.kak.fasta $convert_normalized normalized_idba)
normalized_idba=$(run_idba $fastq.unaligned.pe.kak.fasta $convert_normalized normalized_idba 256gb)
echo $normalized_idba >> ids.txt
echo "Submitted normalized idba"

blastn_normalized_idba=$(run_blastn $fastq.unaligned.pe.kak.idba.d/contig.fa $fastq.unaligned.pe.kak.idba.d/contig.blastn.nt.xml $normalized_idba blastn_normalized_idba)
echo $blastn_normalize_iddba >> ids.txt
echo "Submitted normalized blastn idba"

blastx_normalized_idba=$(run_blastx $fastq.unaligned.pe.kak.keep.idba.d/contig.fa $fastq.unaligned.pe.kak.idba.d/contig.blastx.nr.xml $normalized_idba blastx_normalized_idba)
echo $blastx_normalized_idba >> ids.txt
echo "Submitted normalized blastx idba"

partition=$(qsub -v path=$path,fastq=$fastq -W depend=afterok:$diginorm /scratch/at120/virome-pipeline/partition.sh)
echo $partition >> ids.txt
echo "Sumbitted Partition"

