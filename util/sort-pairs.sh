#!/bin/bash
f=$1  # input unsorted pairs file (gzipped)
outprefix=$2

# header
gunzip -c $f | grep "^#" > $outprefix.pairs
    
# sorting 
gunzip -c $f | grep -v '^#' | sort -k2,2 -k4,4 -k3,3g -k5,5g >> $outprefix.pairs
    
# compressing
bgzip -f $outprefix.pairs

# indexing
pairix -f $outprefix.pairs.gz
