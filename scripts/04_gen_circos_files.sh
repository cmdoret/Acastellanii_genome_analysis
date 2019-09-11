#!/bin/bash

# Generating circos input files.
# cmdoret, 20190603


REF="$1"
#SIGHUNT="$2"
#INTERPRO="$3"
#CAND="$4"
# MCscanX output prefix
MCSX="$2"
CIR_DIR=$3

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
#tail -n +2 $SIGHUNT >${CIR_DIR}/sighunt_dias.txt
# Interpro hits (proportion of bacterial hits), only keep genes with more than 10 hits
#tail -n +2 $INTERPRO |
#    awk '{print $1,$3,$4,$5}' \
#    >${CIR_DIR}/interpro_hits.txt 

# Candidate genes
#tail -n +2  $CAND |
#  awk '{print $1,$2,$3}' \
#  >${CIR_DIR}/hgt_candidates.txt

# Parsing MCScanX collinearity output into circos links
while read line
do
  # Make links involving HGT candidates red
  g1="${line%% *}"
  g2="${line##* }"
  #attr=$( \
  #  awk -v g1=${g1%%-*} -v g2=${g2%%-*} \
  #      'BEGIN{attr="color=lgrey_a4,z=0"} \
  #      {if($0 ~ g1 || $0 ~ g2){attr="color=dred,z=30";exit}} \
  #      END{print attr}' \
  #      ${CAND} \
  #    )
  attr='lgrey'
  awk -v g1=$g1 -v g2=$g2 -v attr="$attr" \
    'BEGIN{N1=0;N2=0}
     $2 ~ g1 {
       N1+=1;c1=$1;s1=$3;e1=$4}
     $2 ~ g2 {
       N2+=1;c2=$1;s2=$3;e2=$4}
     END{
       if (N1 < 1 || N2 < 1) {exit 1}
       else {print c1,s1,e1,c2,s2,e2,attr}}' "$MCSX.gff"
  # Convert gene IDs to coordinates
done < <(grep -v "^#" "$MCSX.collinearity" | tr -d ' ' | tr '\t' ' ' | cut -d$' ' -f2,3)  \
     > "${CIR_DIR}/mcsx.txt"
