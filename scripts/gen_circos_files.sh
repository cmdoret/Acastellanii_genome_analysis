#!/bin/bash

# Generating circos input files.
# cmdoret, 20190603

CIR_DIR="data/input/misc/circos_conf"j
REF=$1
SIGHUNT=$2
INTERPRO=$3
CAND=$4
# Karyotype
awk '/^>/ {if (seqlen){print seqlen}; printf "%s ",$0;seqlen=0;next;}
     {seqlen += length($0)}
     END{print seqlen}' $REF |
  awk 'print "chr","-",$1,$1,"0",$2' \
  >${CIR_DIR}/karyotype.txt

# sighunt DIAS
tail -n +2 $SIGHUNT >${CIR_DIR}/sighunt_dias.txt
# Interpro hits (proportion of bacterial hits)
tail -n +2 $INTERPRO |
    awk '{print $2,$3,$4,$9}' \
    >${CIR_DIR}/interpro_hits.txt 
# Candidates

