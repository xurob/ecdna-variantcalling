#!/usr/bin/python



import sys
import re
import os
import shutil
import pysam
import argparse
import glob
import pathlib
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileMerger



parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="input directory: please add /")
parser.add_argument("-o", "--output", help="output path/to/outputdir/")
parser.add_argument("-rs", "--rstart", help="define starting position", type=int)
parser.add_argument("-re", "--rend", help="define ending position", type=int)
parser.add_argument("-n", "--name", help="name for main output")
parser.add_argument("-ref", "--reference", help = "please supply a reference file as fasta")
parser.add_argument("-bq", "--basequality", help="min basequality for filtering", type=float, default=0)
parser.add_argument("-ebq", "--excludebasequality", help="exclude basequality from outputfiles, default = true", type=bool, default= True)
parser.add_argument("-al", "--alignmentquality", help="min basequality for filtering", type=float, default=0)
parser.add_argument("-cov", "--coveragefile", help="outputs coverage summary as pdf, default = False", type=bool, default= False)
parser.add_argument("-k", "--keep", help="keep temp files, default = true", type=bool, default= True)

args = parser.parse_args()




file_dir = args.input
outpre = args.output
base_qual = args.basequality
name = args.name
ebq = args.excludebasequality
alignment_quality = args.alignmentquality
start1 = args.rstart
maxBP = args.rend
reffile = args.reference
coverageout = args.coveragefile
keep = args.keep





temppath = outpre + "temp/"
if not os.path.exists(temppath):
		   os.makedirs(temppath)
start = int(start1) - 1 #if a sam file is provided add 1 due to 1 indexing
n = (int(maxBP) - int(start)) + 1 #region length
n2 = n - 1
filelistlength = len(list(glob.iglob(file_dir + '**/*.bam', recursive=True)))
itercounter = 0

# =============================================================================
# piles up the counts for every bam file in dir
# =============================================================================

def writeSparseMatrix(mid, vec, sample1):
		with open(temppath + sample1+"."+mid+".txt","w") as V:
			for i in range(0,n2):
				if(vec[i] > 0):
					V.write(str(i+start+1)+","+sample1+","+str(vec[i])+"\n")
			V.close()
	
	
def writeSparseMatrix2(mid, vec1, vec2, sample1):
	with open(temppath + sample1+ "."+mid+".txt","w") as V:
		for i in range(0,n2):
			if(vec1[i] > 0):
				V.write(str(i+start+1)+","+sample1+","+str(vec1[i])+","+str(vec2[i])+"\n")
		V.close()

def writeSparseMatrix4(mid, vec1, vec2, vec3, vec4, sample1):
	with open(temppath + sample1+ "."+mid+".txt","w") as V:
		for i in range(0,n2):
			if(vec1[i] > 0):
				V.write(str(i+start+1)+","+sample1+","+str(vec1[i])+","+str(vec2[i])+"\n")
		V.close()


def writeSparseMatrix3(mid, vec, sample1):
	with open(temppath + sample1+ "." +mid+".txt","w") as V:
				V.write(sample1+","+str(vec)+"\n")
				V.close()


def outputcountsummary(filename,countsfw,countsrv,countstot,depthstat):
	with PdfPages(temppath + sample+'.pdf') as pdf:
		plt.rcParams['text.usetex'] = False
		fig = plt.figure(figsize=(30, 10))
		plt.plot(countsfw, label='forward counts')
		plt.plot(countsrv, label='reverse counts')
		plt.plot(countstot, label='coverage')
		plt.title(filename + " depth:" + str(depthstat))
		plt.xlabel("position from offset")
		plt.ylabel("counts")
		plt.legend(bbox_to_anchor=(1, 1), loc='upper left', borderaxespad=0.)
		pdf.savefig(fig)  
		
		plt.close()

for filename in glob.iglob(file_dir + '**/*.bam', recursive=True):
	sample = str(os.path.splitext(os.path.basename(filename))[0])
	
	# BAQ
	
	qualAfw = [0.0] * n
	qualCfw = [0.0] * n
	qualGfw = [0.0] * n
	qualTfw = [0.0] * n
	
	qualArv = [0.0] * n
	qualCrv = [0.0] * n
	qualGrv = [0.0] * n
	qualTrv = [0.0] * n
	
	# initialize with a pseudo count to avoid dividing by zero
	countsAfw = [0.00000001] * n 
	countsCfw = [0.00000001] * n 
	countsGfw = [0.00000001] * n 
	countsTfw = [0.00000001] * n 
	
	
	countsArv = [0.00000001] * n 
	countsCrv = [0.00000001] * n 
	countsGrv = [0.00000001] * n 
	countsTrv = [0.00000001] * n 
	
	countsfw = [0.00000001] * n 
	countsrv = [0.00000001] * n 
	countstot = [0.00000001] * n
	
	
	covcount = 0
	
	
	#count for each base on each aligned read in the target region
	
	bam2 = pysam.AlignmentFile(filename, "rb")
	for read in bam2:
		seq = read.seq
		quality = read.query_qualities
		align_qual_read = read.mapping_quality
		readrv = read.is_reverse
		for qpos, refpos in read.get_aligned_pairs(True):
			if readrv == False and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
				if(seq[qpos] == "A" and quality[qpos] > base_qual):
					countsAfw[(refpos)-int(start)] += 1
					covcount +=1
					qualAfw[int(refpos)-int(start)] += quality[qpos]
					countsfw[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
				elif(seq[qpos] == "C" and quality[qpos] > base_qual):
					countsCfw[(refpos)-int(start)] += 1
					covcount +=1
					qualCfw[int(refpos)-int(start)] += quality[qpos]
					countsfw[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
				elif(seq[qpos] == "G" and quality[qpos] > base_qual):
					countsGfw[(refpos)-int(start)] += 1
					covcount +=1
					qualGfw[int(refpos)-int(start)] += quality[qpos]
					countsfw[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
				elif(seq[qpos] == "T" and quality[qpos] > base_qual):
					countsTfw[(refpos)-int(start)] += 1
					covcount +=1
					qualTfw[int(refpos)-int(start)] += quality[qpos]
					countsfw[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
			if readrv == True and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
				if(seq[qpos] == "A" and quality[qpos] > base_qual):
					countsArv[(refpos)-int(start)] += 1
					covcount +=1
					qualArv[int(refpos)-int(start)] += quality[qpos]
					countsrv[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
				elif(seq[qpos] == "C" and quality[qpos] > base_qual):
					countsCrv[(refpos)-int(start)] += 1
					covcount +=1
					qualCrv[int(refpos)-int(start)] += quality[qpos]
					countsrv[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
				elif(seq[qpos] == "G" and quality[qpos] > base_qual):
					countsGrv[(refpos)-int(start)] += 1
					covcount +=1
					qualGrv[int(refpos)-int(start)] += quality[qpos]
					countsrv[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
				elif(seq[qpos] == "T" and quality[qpos] > base_qual):
					countsTrv[(refpos)-int(start)] += 1
					covcount +=1
					qualTrv[int(refpos)-int(start)] += quality[qpos]
					countsrv[int(refpos)-int(start)] +=1
					countstot[int(refpos)-int(start)] +=1
	
	meanQualAfw = [round(x/y,1) for x, y in zip(qualAfw, countsAfw)]
	meanQualCfw = [round(x/y,1) for x, y in zip(qualCfw, countsCfw)]
	meanQualGfw = [round(x/y,1) for x, y in zip(qualGfw, countsGfw)]
	meanQualTfw = [round(x/y,1) for x, y in zip(qualTfw, countsTfw)]
	
	meanQualArv = [round(x/y,1) for x, y in zip(qualArv, countsArv)]
	meanQualCrv = [round(x/y,1) for x, y in zip(qualCrv, countsCrv)]
	meanQualGrv = [round(x/y,1) for x, y in zip(qualGrv, countsGrv)]
	meanQualTrv = [round(x/y,1) for x, y in zip(qualTrv, countsTrv)]


	countsAfw = [ int(round(elem)) for elem in countsAfw ]
	countsCfw = [ int(round(elem)) for elem in countsCfw ]
	countsGfw = [ int(round(elem)) for elem in countsGfw ]
	countsTfw = [ int(round(elem)) for elem in countsTfw ]
	
	countsArv = [ int(round(elem)) for elem in countsArv ]
	countsCrv = [ int(round(elem)) for elem in countsCrv ]
	countsGrv = [ int(round(elem)) for elem in countsGrv ]
	countsTrv = [ int(round(elem)) for elem in countsTrv ]
	
	countsfw = [ int(round(elem)) for elem in countsfw ]
	countsrv = [ (int(round(elem))*-1) for elem in countsrv ]
	countstot = [ int(round(elem)) for elem in countstot ]
	
	# Allele Counts
	if ebq == False:
		writeSparseMatrix4("A", countsAfw, countsArv, meanQualAfw, meanQualArv, sample)
		writeSparseMatrix4("C", countsCfw, countsCrv, meanQualCfw, meanQualCrv, sample)
		writeSparseMatrix4("G", countsGfw, countsGrv, meanQualGfw, meanQualGrv, sample)
		writeSparseMatrix4("T", countsTfw, countsTrv, meanQualTfw, meanQualTrv, sample)
	
	else:
		writeSparseMatrix2("A", countsAfw, countsArv, sample)
		writeSparseMatrix2("C", countsCfw, countsCrv, sample)
		writeSparseMatrix2("G", countsGfw, countsGrv, sample)
		writeSparseMatrix2("T", countsTfw, countsTrv, sample)
	
	zipped_list = zip(list(countsAfw),list(countsCfw),list(countsGfw),list(countsTfw),list(countsArv),list(countsCrv),list(countsGrv),list(countsTrv))
	sums = [sum(item) for item in zipped_list]
	writeSparseMatrix("coverage", sums, sample)
	depth = covcount/n #n = region length
	writeSparseMatrix3("depth", depth, sample)
	if coverageout == True:
			outputcountsummary(sample, countsfw, countsrv, sums, depth)
		
	
	
		
	
# =============================================================================
# make tab delimited reference txt file from fasta
# =============================================================================

shutil.copy2(reffile, temppath)
reffiletemp = temppath + str((os.path.basename(reffile)))
base1 = os.path.splitext(reffiletemp)[0]
os.rename(reffiletemp, base1 + ".txt")
start2 = int(start1)
tempindex = 0
with open(reffile) as fasta:
	with open(outpre + "ref.txt","w") as fasta1:
		for line in fasta:  
			if ">" not in line:
			   for ch in line:
				   if tempindex < maxBP and str(ch) != "\n":  
					   tempindex += 1
					   if(tempindex in range(start2,maxBP+1)):
							   fasta1.write(str(tempindex)+"\t"+str(ch)+"\n")
		if(tempindex > maxBP):
				fasta1.close()
				fasta.close()
						  

# =============================================================================
# combines all pileupfiles
# =============================================================================

dirname = os.path.dirname(os.path.abspath(__file__))
pileupcomb = os.path.join(dirname, '02_merge_pileup_counts.sh')
subprocess.check_call("sh "+pileupcomb + " %s %s %s" % (str(temppath), str(name), str(outpre)), shell=True)

merger = PdfFileMerger()
for pdf in glob.iglob(temppath + '**/*.pdf', recursive=True):
    merger.append(pdf)

merger.write(outpre + "coveragesummary.pdf")
merger.close()
# =============================================================================
# call R Script for creating SE Object
# =============================================================================