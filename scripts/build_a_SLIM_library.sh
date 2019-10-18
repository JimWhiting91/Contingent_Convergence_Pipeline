#---------------------------------------------
# THIS SCRIPT BUILDS A REALISTIC GUPPY 'GENE' APPROXIMATION FOR USE IN SLIM SIMULATIONS
# SCRIPT USES INFO FROM A GFF TO BUILD THE GENE

# usage e.g. bash scripts/build_a_SLIM_library.sh some.gff3

#---------------------------------------------

# Path the GFF is defined in the command line
GFF=$1

# ID should correspond to HPC array variable. Otherwise can be run in a forloop or similar
ID=$MOAB_JOBARRAYINDEX
OUTPUT_DIR=outputs/slim_genes
WINDOW_SIZE=25000

mkdir $OUTPUT_DIR
cd $OUTPUT_DIR

# Create a temp GFF with only CDS and only INTRONs
awk '$3 ~ /CDS/' ${GFF} > gff_codon_${ID}.tmp
awk '$3 ~ /intron/' ${GFF} > gff_intron_${ID}.tmp

# We use an until loop to make sure our gene doesn't exceed 20kb
GENE_SIZE=25001
until (( $GENE_SIZE < 25000 ))
do

# How many exons will the gene have, now capped at 10 max for speed purposes?
EXON_N=$(sed -n 's/.*gene_id=//p' gff_codon_${ID}.tmp | cut -c 1-19 | uniq -c | awk '{print $1}' | shuf -n1)

INTRON_N=$(($EXON_N - 1))

# Extract EXON_N exon lengths
awk '{print $5-$4}' gff_codon_${ID}.tmp | shuf -n${EXON_N} > exon_lengths_${ID}.tmp

# Extract 1-EXON_N intron lengths
awk '{print $5-$4}' gff_intron_${ID}.tmp | shuf -n${INTRON_N} > intron_lengths_${ID}.tmp

# How big is the 'gene'?
GENE_SIZE=$(cat exon_lengths_${ID}.tmp intron_lengths_${ID}.tmp | awk '{s+=$1} END {print s}')
echo $GENE_SIZE
# End until loop
done

# Initialise the genome elements and write them in SLIM format

# Calculate the difference between start/end of window. This puts our gene in the centre of the window.
FLANK=$((12500-$GENE_SIZE/2-1))
END_BP=$(($FLANK+1))

echo "initializeGenomicElement(g1, 0, $FLANK);" >> 2_Chunk_Slim_mod_$ID
for i in $(seq 1 $INTRON_N)
do

EXON_SIZE=$(sed "${i}q;d" exon_lengths_${ID}.tmp)
INTRON_SIZE=$(sed "${i}q;d" intron_lengths_${ID}.tmp)

# Add the exon
echo "initializeGenomicElement(g2, $END_BP, $(($END_BP+EXON_SIZE)));" >> 2_Chunk_Slim_mod_$ID
echo $EXON_SIZE >> exons_$ID.txt

# Add the intron
echo "initializeGenomicElement(g1, $(($END_BP+EXON_SIZE+1)), $(($END_BP+EXON_SIZE+$INTRON_SIZE)));" >> 2_Chunk_Slim_mod_$ID

END_BP=$(($END_BP+EXON_SIZE+$INTRON_SIZE+1))
done

# Add the final exon and the final flank
echo "initializeGenomicElement(g2, $END_BP, $(($END_BP+$(sed "${EXON_N}q;d" exon_lengths_${ID}.tmp))));" >> 2_Chunk_Slim_mod_$ID
echo $(sed "${EXON_N}q;d" exon_lengths_${ID}.tmp) >> exons_$ID.txt
echo "initializeGenomicElement(g1, $(($END_BP+$(sed "${EXON_N}q;d" exon_lengths_${ID}.tmp)+1)), 25000);" >> 2_Chunk_Slim_mod_$ID
echo "}" >> 2_Chunk_Slim_mod_$ID

# Finally pass all info onto the information file
echo -e "${ID}\t${GENE_SIZE}\t${EXON_N}\t$(cat exon_lengths_${ID}.tmp | awk '{s+=$1} END {print s}')" >> slim_gene_info.txt

rm -f ./*_${ID}.tmp exons_$ID.txt

