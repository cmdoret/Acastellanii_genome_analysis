#!/bin/bash

# Generating circos input files.
# cmdoret, 20190603


REF="$1"
# MCscanX output prefix
MCSX="$2"
CIR_DIR=$3
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

# Candidate genes
#tail -n +2  $CAND |
#  awk '{print $1,$2,$3}' \
#  >${CIR_DIR}/hgt_candidates.txt

# Parsing MCScanX collinearity output into circos links
n_genes=$(grep -v "^#" "$MCSX.collinearity" | wc -l)
n_done=0
while read line
do
  # Make links involving HGT candidates red
  g1="${line%% *}"
  g2="${line##* }"
#  attr=$( \
#    awk -v g1=${g1%%-*} -v g2=${g2%%-*} \
#        'BEGIN{attr="color=lgrey_a4,z=0"} \
#        {if($0 ~ g1 || $0 ~ g2){attr="color=dred,z=30";exit}} \
#        END{print attr}' \
#        ${CAND} \
#  )
  attr='lgrey'
  # Get coordinate of gene ID from GFF file
  awk -v g1=$g1 -v g2=$g2 -v attr="$attr" \
    'BEGIN{N1=0;N2=0}
     $2 ~ g1 {
       N1+=1;c1=$1;s1=$3;e1=$4}
     $2 ~ g2 {
       N2+=1;c2=$1;s2=$3;e2=$4}
     END{
       if (N1 < 1 || N2 < 1) {exit 1}
       else {print c1,s1,e1,c2,s2,e2,attr}}' "$MCSX.gff"
  echo "$n_done / $n_genes         \r" >2
  (( n_done++ ))
  # Convert gene IDs to coordinates
done < <(grep -v "^#" "$MCSX.collinearity" | tr -d ' ' | tr '\t' ' ' | cut -d$' ' -f2,3)  \
     > "${CIR_DIR}/mcsx.txt"
