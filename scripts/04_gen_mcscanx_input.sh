#!/user/bin/env bash

# This script allows to quickly setup input files for MCScanX. It requires a set
# of gene coordinates in GFF format and a genome assembly in fasta format.
# 12.11.2017
# Cyril Matthey-Doret


# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -g genes -o output_folder -r ref [-h]
   -g   gtf file with genes coordinates
   -o   output folder for MCScanX files
   -f   reference genome
   -c   number of cpus to use for BLAST
   -h   displays this help
EOF
   exit 0
}


N_CPU=1
# Parsing CL arguments
while getopts ":g:o:f:c:h" opt; do
   case $opt in
   g )  GFF=${OPTARG} ;;
   o )  OUT_F=${OPTARG} ;;
   f )  FASTA=${OPTARG};;
   c )  N_CPU=${OPTARG} ;;
   h )  usage ;;
   \?)  usage ;;
   esac
done
# Testing if mandatory arguments have been provided
if [ "x" == "x$FASTA" ] || [ "x" == "x$GFF" ] || \
   [ "x" == "x$OUT_F" ];
then
  echo "Error: A reference genome, input GFF file and output folder \
  are required."
  usage
  exit 0
fi


#1: format GFF file for MCScanX, only keeping CDS
echo "Formatting GFF..."
OUT_GFF="$OUT_F/MCScanX_genes.gff"
awk 'BEGIN{OFS="\t"}
     $3 ~ "CDS" {id_match=match($9,/ID=[^;]*/)
      id=substr($9,RSTART+3,RLENGTH-3)
      print $1,id,$4,$5}' $GFF > $OUT_GFF
echo "GFF formatted !"

#2: build blast database from sequences and all vs all blast
BASE_OUT=$(basename $FASTA)
BASE_OUT=${BASE_OUT%%.*}
makeblastdb -in $FASTA -dbtype prot -out "${OUT_F}/${BASE_OUT}"

echo "Blasting transcriptome against itself."
blastp -query $FASTA \
       -db "${OUT_F}/${BASE_OUT}" \
       -outfmt 6 \
       -max_target_seqs 5 \
       -num_threads $N_CPU \
       -out "${OUT_F}/${BASE_OUT}.blast"

#3 Renaming chromosomes according to MCScanX conventions
MC_IN="$OUT_F/MCScanX_in"
# Use two first characters of fasta file as organism abbreviation
#SP_CODE=${BASE_OUT:0:2}
#sed "s/^scaffold_\([0-9]*\)/$SP_CODE\1/g" | 
sed 's/\.cds//g' "$OUT_GFF" > "$MC_IN.gff"
#sed "s/scaffold_\([0-9]*\)/$SP_CODE\1/g" | 
sed 's/::[^	]*//g' "${OUT_F}/${BASE_OUT}.blast" > "$MC_IN.blast"


echo "All input files are ready:"
echo "  - BLAST output: $MC_IN.blast"
echo "  - GFF file: $MC_IN.gff"

