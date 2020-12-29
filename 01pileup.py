#!/usr/bin/python

###################################################
# Summarizes the total number of reads per position
###################################################

import sys
import re
import os
import pysam
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="input file in bam format")
parser.add_argument("-o", "--output", help="output path/to/outputfolder/$name")
parser.add_argument("-rs", "--rstart", help="define starting position", type=int)
parser.add_argument("-re", "--rend", help="define ending position", type=int)
parser.add_argument("-n", "--name", help="samplename")
parser.add_argument("-bq", "--basequality", help="min basequality for filtering", type=float, default=0)
parser.add_argument("-al", "--alignmentquality", help="min basequality for filtering", type=float, default=0)
args = parser.parse_args()




bamfile = args.input
outpre = args.output
base_qual = args.basequality
sample = args.name
alignment_quality = args.alignmentquality
start1 = args.rstart
maxBP = args.rend




start = int(start1) - 1 #if a sam file is provided add 1 due to 1 indexing
n = (int(maxBP) - int(start)) + 1 #region length
n2 = n - 1

# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,n2):
			if(vec[i] > 0):
				V.write(str(i+start+1)+","+sample+","+str(vec[i])+"\n")


def writeSparseMatrix2(mid, vec1, vec2):
	with open(outpre + "."+mid+".txt","w") as V:
		for i in range(0,n2):
			if(vec1[i] > 0):
				V.write(str(i+start+1)+","+sample+","+str(vec1[i])+","+str(vec2[i])+"\n")


# BAQ
# initialize with a pseudo count to avoid dividing by zero
countsAfw = [0.00000001] * n 
countsCfw = [0.00000001] * n 
countsGfw = [0.00000001] * n 
countsTfw = [0.00000001] * n 


countsArv = [0.00000001] * n 
countsCrv = [0.00000001] * n 
countsGrv = [0.00000001] * n 
countsTrv = [0.00000001] * n 




#count for each base on each aligned read in the target region

bam2 = pysam.AlignmentFile(bamfile, "rb")
for read in bam2:
	seq = read.seq
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	for qpos, refpos in read.get_aligned_pairs(True):
		if is_reverse() == False and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
			if(seq[qpos] == "A" and quality[qpos] > base_qual):
				countsAfw[(refpos)-int(start)] += 1
			elif(seq[qpos] == "C" and quality[qpos] > base_qual):
				countsCfw[(refpos)-int(start)] += 1
			elif(seq[qpos] == "G" and quality[qpos] > base_qual):
				countsGfw[(refpos)-int(start)] += 1
			elif(seq[qpos] == "T" and quality[qpos] > base_qual):
				countsTfw[(refpos)-int(start)] += 1
        if is_reverse() == True and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
            if(seq[qpos] == "A" and quality[qpos] > base_qual):
				countsArv[(refpos)-int(start)] += 1
			elif(seq[qpos] == "C" and quality[qpos] > base_qual):
				countsCvr[(refpos)-int(start)] += 1
			elif(seq[qpos] == "G" and quality[qpos] > base_qual):
				countsGvr[(refpos)-int(start)] += 1
			elif(seq[qpos] == "T" and quality[qpos] > base_qual):
				countsTrv[(refpos)-int(start)] += 1


countsAfw = [ int(round(elem)) for elem in countsAfw ]
countsCfw = [ int(round(elem)) for elem in countsCfw ]
countsGfw = [ int(round(elem)) for elem in countsGfw ]
countsTfw = [ int(round(elem)) for elem in countsTfw ]

countsArv = [ int(round(elem)) for elem in countsArv ]
countsCrv = [ int(round(elem)) for elem in countsCrv ]
countsGrv = [ int(round(elem)) for elem in countsGrv ]
countsTrv = [ int(round(elem)) for elem in countsTrv ]


# Allele Counts

writeSparseMatrix2("A", countsAfw, countsArv)
writeSparseMatrix2("C", countsCfw, countsCrv)
writeSparseMatrix2("G", countsGfw, countsGrv)
writeSparseMatrix2("T", countsTfw, countsTrv)

zipped_list = zip(list(countsAfw),list(countsCfw),list(countsGfw),list(countsTrv),countsArv),list(countsCrv),list(countsGrv),list(countsTrv))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)