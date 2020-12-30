#!/usr/bin/python



import sys
import re
import os
import pysam
import argparse
import glob
import subprocess



parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", help="input directory: please add /")
parser.add_argument("-o", "--output", help="output path/to/outputdir/")
parser.add_argument("-rs", "--rstart", help="define starting position", type=int)
parser.add_argument("-re", "--rend", help="define ending position", type=int)
parser.add_argument("-n", "--name", help="name for main output")
parser.add_argument("-bq", "--basequality", help="min basequality for filtering", type=float, default=0)
parser.add_argument("-ebq", "--excludebasequality", help="exclude basequality from outputfiles, default = true", type=bool, default= True)
parser.add_argument("-al", "--alignmentquality", help="min basequality for filtering", type=float, default=0)
args = parser.parse_args()




file_dir = args.input
outpre = args.output
base_qual = args.basequality
name = args.name
ebq = args.excludebasequality
alignment_quality = args.alignmentquality
start1 = args.rstart
maxBP = args.rend

# =============================================================================
# piles up the counts for every bam file in dir
# =============================================================================

if ebq == True:
	for filename in glob.iglob(file_dir + '**/*.bam', recursive=True):
		start = int(start1) - 1 #if a sam file is provided add 1 due to 1 indexing
		n = (int(maxBP) - int(start)) + 1 #region length
		n2 = n - 1
		sample = string(filename)
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
		
		def writeSparseMatrix3(mid, vec):
			with open(outpre + "."+mid+".txt","w") as V:
				for i in range(0,n2):
					if(vec[i] > 0):
						V.write(sample+","+str(vec[i])+"\n")
		
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
		
		covcount = 0
		
		
		#count for each base on each aligned read in the target region
		
		bam2 = pysam.AlignmentFile(filename, "rb")
		for read in bam2:
			seq = read.seq
			quality = read.query_qualities
			align_qual_read = read.mapping_quality
			for qpos, refpos in read.get_aligned_pairs(True):
				if is_reverse() == False and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
					if(seq[qpos] == "A" and quality[qpos] > base_qual):
						countsAfw[(refpos)-int(start)] += 1
						covcount +=1
					elif(seq[qpos] == "C" and quality[qpos] > base_qual):
						countsCfw[(refpos)-int(start)] += 1
						covcount +=1
					elif(seq[qpos] == "G" and quality[qpos] > base_qual):
						countsGfw[(refpos)-int(start)] += 1
						covcount +=1
					elif(seq[qpos] == "T" and quality[qpos] > base_qual):
						countsTfw[(refpos)-int(start)] += 1
						covcount +=1
		        if is_reverse() == True and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
		            if(seq[qpos] == "A" and quality[qpos] > base_qual):
						countsArv[(refpos)-int(start)] += 1
						covcount +=1
					elif(seq[qpos] == "C" and quality[qpos] > base_qual):
						countsCvr[(refpos)-int(start)] += 1
						covcount +=1
					elif(seq[qpos] == "G" and quality[qpos] > base_qual):
						countsGvr[(refpos)-int(start)] += 1
						covcount +=1
					elif(seq[qpos] == "T" and quality[qpos] > base_qual):
						countsTrv[(refpos)-int(start)] += 1
						covcount +=1
		
		
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
		depth = covcount/n #n = region length
		writeSparseMatrix3("depth", depth)

else:
	for filename in glob.iglob(file_dir + '**/*.bam', recursive=True):
		start = int(start1) - 1 #if a sam file is provided add 1 due to 1 indexing
		n = (int(maxBP) - int(start)) + 1 #region length
		n2 = n - 1
		sample = string(filename)
		# Export Functions
		def writeSparseMatrix(mid, vec):
			with open(outpre + "."+mid+".txt","w") as V:
				for i in range(0,n2):
					if(vec[i] > 0):
						V.write(str(i+start+1)+","+sample+","+str(vec[i])+"\n")
		
		
		def writeSparseMatrix2(mid, vec1, vec2, vec3, vec4):
			with open(outpre + "."+mid+".txt","w") as V:
				for i in range(0,n2):
					if(vec1[i] > 0):
						V.write(str(i+start+1)+","+sample+","+str(vec1[i])+","+str(vec3[i])+","+str(vec2[i])+","+str(vec4[i])+"\n")
		
		def writeSparseMatrix3(mid, vec):
			with open(outpre + "."+mid+".txt","w") as V:
				for i in range(0,n2):
					if(vec[i] > 0):
						V.write(sample+","+str(vec[i])+"\n")
		
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
		
		covcount = 0
		
		
		#count for each base on each aligned read in the target region
		
		bam2 = pysam.AlignmentFile(filename, "rb")
		for read in bam2:
			seq = read.seq
			quality = read.query_qualities
			align_qual_read = read.mapping_quality
			for qpos, refpos in read.get_aligned_pairs(True):
				if is_reverse() == False and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
					if(seq[qpos] == "A" and quality[qpos] > base_qual):
						countsAfw[(refpos)-int(start)] += 1
						covcount +=1
						qualAfw[int(refpos)-int(start)] += quality[qpos]
					elif(seq[qpos] == "C" and quality[qpos] > base_qual):
						countsCfw[(refpos)-int(start)] += 1
						covcount +=1
						qualCfw[int(refpos)-int(start)] += quality[qpos]
					elif(seq[qpos] == "G" and quality[qpos] > base_qual):
						countsGfw[(refpos)-int(start)] += 1
						covcount +=1
						qualGfw[int(refpos)-int(start)] += quality[qpos]
					elif(seq[qpos] == "T" and quality[qpos] > base_qual):
						countsTfw[(refpos)-int(start)] += 1
						covcount +=1
						qualTfw[int(refpos)-int(start)] += quality[qpos]
		        if is_reverse() == True and qpos is not None and refpos is not None and align_qual_read > alignment_quality and int(start) <= int(refpos) <= int(maxBP):
		            if(seq[qpos] == "A" and quality[qpos] > base_qual):
						countsArv[(refpos)-int(start)] += 1
						covcount +=1
						qualArv[int(refpos)-int(start)] += quality[qpos]
					elif(seq[qpos] == "C" and quality[qpos] > base_qual):
						countsCvr[(refpos)-int(start)] += 1
						covcount +=1
						qualCrv[int(refpos)-int(start)] += quality[qpos]
					elif(seq[qpos] == "G" and quality[qpos] > base_qual):
						countsGvr[(refpos)-int(start)] += 1
						covcount +=1
						qualGrv[int(refpos)-int(start)] += quality[qpos]
					elif(seq[qpos] == "T" and quality[qpos] > base_qual):
						countsTrv[(refpos)-int(start)] += 1
						covcount +=1
						qualT[int(refpos)-int(start)] += quality[qpos]
		
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
		
		
		# Allele Counts
		
		writeSparseMatrix2("A", countsAfw, countsArv, meanQualAfw, meanQualArv)
		writeSparseMatrix2("C", countsCfw, countsCrv, meanQualCfw, meanQualCrv)
		writeSparseMatrix2("G", countsGfw, countsGrv, meanQualGfw, meanQualGrv)
		writeSparseMatrix2("T", countsTfw, countsTrv, meanQualTfw, meanQualTrv)
		
		zipped_list = zip(list(countsAfw),list(countsCfw),list(countsGfw),list(countsTrv),countsArv),list(countsCrv),list(countsGrv),list(countsTrv))
		sums = [sum(item) for item in zipped_list]
		writeSparseMatrix("coverage", sums)
		depth = covcount/n #n = region length
		writeSparseMatrix3("depth", depth)


	
# =============================================================================
# combines all pileupfiles
# =============================================================================

dirname = os.path.dirname(os.path.abspath(__file__))
pileupcomb = os.path.join(dirname, '02_merge_pileup_counts.sh')
subprocess.check_call(pileupcomb + " %s %s %s" % (str(file_dir), str(name), str(outpre)), shell=True)


# =============================================================================
# call R Script for creating SE Object
# =============================================================================