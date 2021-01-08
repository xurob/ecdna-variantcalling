
#!/usr/bin/python



import sys
import re
import pysam
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def outputcountsummary(filename,countsfw,countsrv,countstot,depthstat):
	with PdfPages(filename + ".pdf") as pdf:
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


alignment_quality, base_qual = 0
start = 1
maxBP,n  = 16569 #rCRS length


'''
Count matrices for every base on forward and reverse strand
n = length reference genome
'''

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

filename = "PATH/TO/BAMFILE"
#count for each base on each aligned read in the target region

'''
+=1 for every base counted at the respective reference position on the respective strand
'''

bam = pysam.AlignmentFile(filename, "rb")
for read in bam:
	seq = read.seq
	readrv = read.is_reverse
	for qpos, refpos in read.get_aligned_pairs(True):
		if readrv == False:
			if(seq[qpos] == "A"):
				countsAfw[(refpos)-int(start)] += 1
				covcount +=1
				countsfw[int(refpos)-int(start)] +=1
			elif(seq[qpos] == "C"):
				countsCfw[(refpos)-int(start)] += 1
				covcount +=1
				countsfw[int(refpos)-int(start)] +=1
			elif(seq[qpos] == "G"):
				countsGfw[(refpos)-int(start)] += 1
				covcount +=1
				countsfw[int(refpos)-int(start)] +=1
			elif(seq[qpos] == "T"):
				countsTfw[(refpos)-int(start)] += 1
				covcount +=1
				countsfw[int(refpos)-int(start)] +=1
		if readrv == True:
			if(seq[qpos] == "A"):
				countsArv[(refpos)-int(start)] += 1
				covcount +=1
				countsrv[int(refpos)-int(start)] +=1
			elif(seq[qpos] == "C"):
				countsCrv[(refpos)-int(start)] += 1
				covcount +=1
				countsrv[int(refpos)-int(start)] +=1
			elif(seq[qpos] == "G"):
				countsGrv[(refpos)-int(start)] += 1
				covcount +=1
				countsrv[int(refpos)-int(start)] +=1
			elif(seq[qpos] == "T"):
				countsTrv[(refpos)-int(start)] += 1
				covcount +=1
				countsrv[int(refpos)-int(start)] +=1


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



zipped_list = zip(list(countsAfw),list(countsCfw),list(countsGfw),list(countsTfw),list(countsArv),list(countsCrv),list(countsGrv),list(countsTrv))
sums = [sum(item) for item in zipped_list]
depth = covcount/n #n = region length

outputcountsummary(filename, countsfw, countsrv, sums, depth) 
'''

'''
		
		
		
		
					
					
					
					
					
					