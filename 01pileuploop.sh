#!/bin/bash
FILES=/fast/users/xuro_c/scratch/vc_mito_mycn/mito/control/*bam
for f in $FILES; do
	fbname=$(basename "$f" .bam)
	python /fast/users/xuro_c/scratch/vc_mito_mycn/mito/control/01pileup.py -i $f -o /fast/users/xuro_c/scratch/vc_mito_mycn/mito/control/out/$fbname -n $fbname -re 16300 -rs 30
done
