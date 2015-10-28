#!/usr/bin/python

from Bio import SeqIO
import sys

seqIDS = []

with open(sys.argv[1]) as blastFile:
		for hit in blastFile.readlines():
				seqIDS.append(hit.split("\t")[0])

with open(sys.argv[2], 'rb') as r1File, open(sys.argv[3], 'rb') as r2File:
		for r1,r2 in SeqIO.parse(r1File, "fastq"),SeqIO.parse(r2File, "fastq"):
			print r2.id
			#if r1.id in seqIDS:
			#	print r1.id