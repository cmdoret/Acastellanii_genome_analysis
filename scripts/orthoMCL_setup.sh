#!bin/bash
# Using grankek's docker container to setup orthoMCL
# Configuration script adapted from kadeng/OrthoMCL-docker
# cmdoret, 20190515

FASTA=$1
N_CPU=$2
RUN_DIR=${3:-$PWD}
MNT=/host_dir
# Start MySQL
echo "Starting MySQL container..."
docker run --name orthomcl-mysql -e MYSQL_ROOT_PASSWORD=asdf1234 -e MYSQL_USER=orthomcl_user -e MYSQL_PASSWORD=shhh_this_is_secret -e MYSQL_DATABASE=orthomcl -d mysql | tee .ortho_mcl_db_container_id | head -n1

# Only keep container id in temporary file
echo `tail -n 1 .ortho_mcl_db_container_id` > .ortho_mcl_db_container_id

sleep 10

# Start OrthoMCL container
echo "Starting OrthoMCL container..."
echo "Mounting local directory ${RUN_DIR} to $MNT within the OrthoMCL container"
docker run -it --name orthomcl-run --link orthomcl-mysql:mysql -d -v ${RUN_DIR}:${MNT} granek/orthomcl 

# Generate orthomcl config (env variables exist in container, output is written to file in host path)
docker exec orthomcl-run echo """dbVendor=mysql
 dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1:$MYSQL_PORT_3306_TCP_ADDR:$MYSQL_PORT_3306_TCP_PORT
 dbLogin=orthomcl_user
 dbPassword=shhh_this_is_secret
 similarSequencesTable=SimilarSequences
 orthologTable=Ortholog
 inParalogTable=InParalog
 coOrthologTable=CoOrtholog
 interTaxonMatchView=InterTaxonMatch
 percentMatchCutoff=50
 evalueExponentCutoff=-5
 oracleIndexTblSpc=NONE""" > "$RUN_DIR/orthoMCL.config"

# Run orthoMCL setup inside the container
docker exec orthomcl-run orthomclInstallSchema $MNT/orthoMCL.config
# Note I did not filter proteins, could be added
# Make blastdb for input fasta file
#makeblastdb -in $FASTA -out $RUN_DIR/orthoMCL.blastdb -dbtype prot
echo "## Creating diamond index of all proteins"
diamond makedb --in $FASTA --db $RUN_DIR/orthoMCL.diamond
# Run all-vs-all blastp on all proteins
#docker exec orthomcl-run blastall -p blastp -a $N_CPU -i $MNT/fasta/all_proteins_orthoMCL.fa -m 8 -d $MNT/orthoMCL.blastdb -o $MNT/orthoMCL.blastp
echo "## Using diamond blastp on all vs all proteins"
diamond blastp --threads $N_CPU \
               --query $FASTA \
               --db $RUN_DIR/orthoMCL.diamond \
               -f 6 -o $RUN_DIR/orthoMCL.blastp
echo "====== Starting BlastParser"
# Parse BLAST tabular output  using orthomclBlastParser
docker exec orthomcl-run orthomclBlastParser $MNT/orthoMCL.blastp $MNT/fasta/ >> $RUN_DIR/orthoMCL.parsed.txt
echo "====== Starting LoadBlast"
# Load output into the database using orthomclLoadBlast
docker exec orthomcl-run orthomclLoadBlast $MNT/orthoMCL.config $MNT/orthoMCL.parsed.txt
echo "====== Starting mclPairs"
# Compute pairwise relationships
docker exec orthomcl-run orthomclPairs $MNT/orthoMCL.config $MNT/pairs.log cleanup=all
# Dump files of pairs from the database. Generates mclInput file and pairs folder
docker exec orthomcl-run orthomclDumpPairsFiles $MNT/orthoMCL.config
#(12) run the mcl program on the mcl_input.txt file created in Step 11.
echo "## Running Markov clustering on blast output"
docker exec orthomcl-run mcl 
#(13) run orthomclMclToGroups to convert mcl output to groups.txt
docker exec orthomcl-run orthomclMclToGroups

docker stop orthomcl-run orthomcl-mysql
docker rm orthomcl-run orthomcl-mysql