#!/user/bin/env bash

# This script allows to quickly setup input files for MCScanX. It requires a set
# of gene coordinates in GFF format and a genome assembly in fasta format.
# 12.11.2017
# Cyril Matthey-Doret


# Help message
function usage () {
   cat <<EOF
Usage: `basename $0` -g genes -o output_folder -r ref [-p prev_ref] [-c corresp] [-l] [-h]
   -g   gtf file with genes coordinates
   -o   output folder for MCScanX files
   -r   reference genome
   -h   displays this help
EOF
   exit 0
}



# Parsing CL arguments
while getopts ":g:o:r:p:c:lh" opt; do
   case $opt in
   o )  OUT_F=${OPTARG} ;;
   r )  REF=${OPTARG};;
   h )  usage ;;
   \?)  usage ;;
   esac
done

# Testing if mandatory arguments have been provided
if [ "x" == "x$REF" ] || [ "x" == "x$GFF" ] || \
   [ "x" == "x$OUT_F" ];
then
  echo "Error: A reference genome, input GFF file and output folder \
  are required."
  usage
  exit 0
fi


#1: Extract records for features of interest from gff file
echo -n "Extracting genes coordinates from the gff file..."
OUT_GFF="$OUT_F/MCScanX_genes.gff"
awk '$3 ~ "gene" {print $0}' "$GFF" > "$OUT_GFF"
echo "genes extracted !"


#2: convert GFF to BED format and extract sequence
# bed format: chrom start end id 0 strand
echo "Converting GFF to BED."
OUT_BED="${OUT_GFF%.*}.bed"
awk 'BEGIN{OFS="\t"}
     {id_match=match($9,/ID=[^;]*/)
      id=substr($9,RSTART+3,RLENGTH-3)
      print $1,$4,$5,id,0,$7}' $OUT_GFF > $OUT_BED

awk 'BEGIN{OFS="\t"}{print $1,$4,$2,$3}' "$OUT_BED" > "$OUT_GFF"
OUT_SEQ="$OUT_BED.fasta"


echo "Extracting gene sequences from reference genome"
bedtools getfasta -fi $REF -bed $OUT_BED -fo $OUT_SEQ -name

#4: build blast database from sequences and all vs all blast
makeblastdb -in $OUT_SEQ -dbtype nucl

echo "Blasting transcriptome against itself."
blastn -query $OUT_SEQ \
       -db $OUT_SEQ \
       -outfmt 6 \
       -max_target_seqs 5 \
       -out $OUT_SEQ.blast

# Renaming chromosomes according to MCScanX conventions
MC_IN="$OUT_F/MCScanX_in"
sed 's/^chr\([0-9]*\)/lf\1/' "$OUT_GFF" > "$MC_IN.gff"
sed 's/chr\([0-9]*\)/lf\1/g' "$OUT_SEQ.blast" \
    | sed 's/::[^	]*//g' > "$MC_IN.blast"


echo "All input files are ready:"
echo "  - BLAST output: $MC_IN.blast"
echo "  - GFF file: $MC_IN.gff"

echo "Running MCScanX"
MCScanX -s 3 $MC_IN

echo "Generating graphics control file for circle plotter"
# 800 pixels, displaying all chromosomes in the input GFF file
echo "1920" > $OUT_F/circle.ctl
cut -f1 "$MC_IN.gff" | uniq | paste -s -d, - >> $OUT_F/circle.ctl