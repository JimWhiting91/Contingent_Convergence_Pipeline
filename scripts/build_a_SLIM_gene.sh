#!/bin/bash

#---------------------------------------------
# THIS SCRIPT BUILDS A REALISTIC GUPPY 'GENE' APPROXIMATION FOR USE IN SLIM SIMULATIONS
# SCRIPT USES INFO FROM A GFF TO BUILD THE GENE

usage: bash build_a_SLIM_gene.sh path/to/gff gene_ID
gene_ID is usually the job array index when running over an HPC
#---------------------------------------------

# Path the GFF is defined in the command line
GFF=$1
ID=$2

# Create a temp GFF with only CDS and only INTRONs
awk '$3 ~ /CDS/' ${GFF} > gff_codon.tmp
awk '$3 ~ /intron/' ${GFF} > gff_intron.tmp

# How many exons will the gene have, now capped at 10 max for speed purposes?
EXON_N=$(sed -n 's/.*gene_id=//p' gff_codon.tmp | cut -c 1-19 | uniq -c | awk '{print $1}' | awk '$1 < 11' | shuf -n1)

INTRON_N=$(($EXON_N - 1))

# Extract EXON_N exon lengths
awk '{print $5-$4}' gff_codon.tmp | shuf -n${EXON_N} > exon_lengths.tmp

# Extract 1-EXON_N intron lengths
awk '{print $5-$4}' gff_intron.tmp | shuf -n${INTRON_N} > intron_lengths.tmp

# Initialise the genome elements and write them in SLIM format
END_BP=1000
for i in $(seq 1 $INTRON_N)
do

EXON_SIZE=$(sed "${i}q;d" exon_lengths.tmp)
INTRON_SIZE=$(sed "${i}q;d" intron_lengths.tmp)

# Add the exon
echo "initializeGenomicElement(g2, $END_BP, $(($END_BP+EXON_SIZE)));" >> 2_Chunk_Slim_mod_$ID
echo $EXON_SIZE >> exons.txt

# Add the intron
echo "initializeGenomicElement(g1, $(($END_BP+EXON_SIZE+1)), $(($END_BP+EXON_SIZE+$INTRON_SIZE)));" >> 2_Chunk_Slim_mod_$ID

END_BP=$(($END_BP+EXON_SIZE+$INTRON_SIZE+1))
done

# Add the final exon and the final flank
echo "initializeGenomicElement(g2, $END_BP, $(($END_BP+$(sed "${EXON_N}q;d" exon_lengths.tmp))));" >> 2_Chunk_Slim_mod_$ID
echo $(sed "${EXON_N}q;d" exon_lengths.tmp) >> exons.txt
echo "initializeGenomicElement(g1, $(($END_BP+$(sed "${EXON_N}q;d" exon_lengths.tmp)+1)), $(($END_BP+$(sed "${EXON_N}q;d" exon_lengths.tmp)+1000)));" >> 2_Chunk_Slim_mod_$ID
echo "}" >> 2_Chunk_Slim_mod_$ID

rm -f ./*.tmp

