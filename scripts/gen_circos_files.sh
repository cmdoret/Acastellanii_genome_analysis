#!/bin/bash

# Generating circos input files.
# cmdoret, 20190603

CIR_DIR="data/input/misc/circos_conf"
REF="$1"
SIGHUNT="$2"
INTERPRO="$3"
CAND="$4"

# Only include scaffolds larger than 100kb
len_cutoff=100000

# Karyotype
awk '/^>/ {if (seqlen){print seqlen}; printf "%s ",$0;seqlen=0;next;}
     {seqlen += length($0)}
     END{print seqlen}' $REF |
  awk -v minlen=$len_cutoff '$2 >= minlen {print "chr","-",$1,$1,"0",$2,"Neff"}' |
  tr -d '>' \
  >${CIR_DIR}/karyotype.txt

# sighunt DIAS
tail -n +2 $SIGHUNT >${CIR_DIR}/sighunt_dias.txt
# Interpro hits (proportion of bacterial hits)
tail -n +2 $INTERPRO |
    awk '{if ($9 == ""){prop=0} else {prop=$9} print $2,$3,$4,prop}' \
    >${CIR_DIR}/interpro_hits.txt 

# Candidate genes
tail -n +2  $CAND |
  awk '{print $1,$2,$3}' \
  >${CIR_DIR}/hgt_candidates.txt

