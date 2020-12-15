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




start = int(start1) - 1
n = (int(maxBP) - int(start)) + 1
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
countsA = [0.00000001] * n 
countsC = [0.00000001] * n 
countsG = [0.00000001] * n 
countsT = [0.00000001] * n 


qualA = [0.0] * n
qualC = [0.0] * n
qualG = [0.0] * n
qualT = [0.0] * n



bam2 = pysam.AlignmentFile(bamfile, "rb")
for read in bam2:
	seq = read.seq
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
			if(seq[qpos] == "A" and quality[qpos] > base_qual):
				qualA[int(refpos)-int(start)] += quality[qpos]
				countsA[(refpos)-int(start)] += 1
			elif(seq[qpos] == "C" and quality[qpos] > base_qual):
				qualC[(refpos)-int(start)] += quality[qpos]
				countsC[(refpos)-int(start)] += 1
			elif(seq[qpos] == "G" and quality[qpos] > base_qual):
				qualG[(refpos)-int(start)] += quality[qpos]
				countsG[(refpos)-int(start)] += 1
			elif(seq[qpos] == "T" and quality[qpos] > base_qual):
				qualT[(refpos)-int(start)] += quality[qpos]
				countsT[(refpos)-int(start)] += 1
			
meanQualA = [round(x/y,1) for x, y in zip(qualA, countsA)]
meanQualC = [round(x/y,1) for x, y in zip(qualC, countsC)]
meanQualG = [round(x/y,1) for x, y in zip(qualG, countsG)]
meanQualT = [round(x/y,1) for x, y in zip(qualT, countsT)]

countsA = [ int(round(elem)) for elem in countsA ]
countsC = [ int(round(elem)) for elem in countsC ]
countsG = [ int(round(elem)) for elem in countsG ]
countsT = [ int(round(elem)) for elem in countsT ]


# Allele Counts

writeSparseMatrix2("A", countsA, meanQualA)
writeSparseMatrix2("C", countsC, meanQualC)
writeSparseMatrix2("G", countsG, meanQualG)
writeSparseMatrix2("T", countsT, meanQualT)

zipped_list = zip(list(countsA),list(countsC),list(countsG),list(countsT))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)