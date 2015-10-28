#!/bin/bash

#PBS -V
#PBS -m ae
#PBS -M twaddlac@gmail.com	
#PBS -S /bin/bash
#PBS -j oe -N diginorm -l walltime=168:00:00,nodes=1:ppn=1,mem=256gb  

cd $path

module purge
module load khmer/intel/1.3

echo "first norm pe"

normalize-by-median.py -k 20 -C 20 -N 4 -x 64e9 -p --savetable $fastq.normC20k20.kh $fastq.unconc.interleaved.fastq 

echo "second norm se" >> khmer.output

normalize-by-median.py -C 20  --loadtable $fastq.normC20k20.kh --savetable $fastq.normC20k20.kh $fastq.*singletons.fastq 

echo "filter abund">> khmer.output

filter-abund.py -V --threads 12 $fastq.normC20k20.kh *keep 

echo "extract paired" >> khmer.output

extract-paired-reads.py $fastq.unconc.interleaved.fastq.keep.abundfilt 

echo "norm 2 pe" >> khmer.output

normalize-by-median.py -k 20 -C 5 -N 4 -x 64e9 -p  --savetable $fastq.normC5k20.kh $fastq.unconc.interleaved.fastq.keep.abundfilt.pe 

echo "norm 3 se" >> khmer.output

normalize-by-median.py -C 5 --loadtable $fastq.normC5k20.kh --savetable $fastq.normC5k20.kh $fastq.*singletons*.keep.abundfilt $fastq.unconc.interleaved.fastq.keep.abundfilt.se 

echo "zip" >> khmer.output

gzip -9c $fastq.unconc.interleaved.fastq.keep.abundfilt.pe.keep > $fastq.unaligned.pe.kak.fastq.gz
gzip -9c $fastq.*singletons*.keep.abundfilt.keep $fastq.unconc.interleaved.fastq.keep.abundfilt.se.keep > $fastq.unaligned.se.kak.fastq.gz

rm $fastq.normC20k20.kh *.abundfilt *.pe *.se

