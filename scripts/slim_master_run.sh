#---------------------------
# SLIM MASTER SCRIPT

# Parameters are entered at the top of scripts
# Designed to be run over an HPC, with $MOAB_JOBARRAYINDEX corresponding to array variable. Can be replaced with forloop or similar although this would take substantially longer.

# This script performs simulations with both diverged (PhenoDiv) and common (PhenoNull) phenotypic optima.

#---------------------------

# Set up working directory
MASTER=~/slim_dir

# INPUT BASIC PARAMETERS
# Size of Pop and mutation rate
POP_N=1000
MUTATION="4.89e-6"
REC_RATE="1e-6"
MUT_EFFECT=1
OUT_DIR="BASIC_OUTPUT"
SLIM_LIBRARY=$MASTER/outputs/slim_genes

# Make output folder if it does not already exist
mkdir $MASTER/demographics_model/$OUT_DIR

# Make directory to work in
mkdir $MASTER/demographics_model/$OUT_DIR/GENE_${MOAB_JOBARRAYINDEX}
cd $MASTER/demographics_model/$OUT_DIR/GENE_${MOAB_JOBARRAYINDEX}

# Get Exon size and selection coef
EXON_SIZE=$(awk -v GENE_N="${MOAB_JOBARRAYINDEX}" '$1==GENE_N {print $4}' $SLIM_LIBRARY/slim_gene_info.txt)
SEL_COEF=$(awk -v GENE_N="${MOAB_JOBARRAYINDEX}" '$1==GENE_N {print $5}' $SLIM_LIBRARY/slim_gene_info.txt)

# Modify 1_Chunk to include new selection coefficient and pop_size and mutation rate
sed "s/selection_coef/${SEL_COEF}/g" $MASTER/scripts/1_Chunk_Spatial_Slim.txt > 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/pop_coef/${POP_N}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/ID_coef/${MOAB_JOBARRAYINDEX}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/MUT_RATE/${MUTATION}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/REC_RATE/${REC_RATE}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/EXON_L/${EXON_SIZE}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}

# Modify 1.5_Chunk to include new selection coefficient and pop_size and mutation rate
sed "s/selection_coef/${SEL_COEF}/g" $MASTER/scripts/1.5_Chunk_Spatial_Slim.txt > 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/pop_coef/${POP_N}/g" 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/ID_coef/${MOAB_JOBARRAYINDEX}/g" 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/MUT_RATE/${MUTATION}/g" 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/REC_RATE/${REC_RATE}/g" 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/EXON_L/${EXON_SIZE}/g" 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}

# Modify 1_Chunk to include new selection coefficient and pop_size and mutation rate
sed "s/selection_coef/${SEL_COEF}/g" $MASTER/scripts/1_Chunk_Spatial_Slim_NullPheno.txt > 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}
sed -i "s/pop_coef/${POP_N}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}
sed -i "s/ID_coef/${MOAB_JOBARRAYINDEX}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}
sed -i "s/MUT_RATE/${MUTATION}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}
sed -i "s/REC_RATE/${REC_RATE}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}
sed -i "s/EXON_L/${EXON_SIZE}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}

# We also need to seed the dXY function with the number of sites
SITE_N=25000
sed -i "s/siteN/${SITE_N}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}
sed -i "s/siteN/${SITE_N}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}

# Modify 3_Chunks so burn-in time = 5 x 2 N
BURN_GEN=$((5*2*POP_N))
cp $MASTER/scripts/3_Chunk_Spatial_Slim_FIRSTHALF.txt 3_Chunk_Slim_mod_SIM_FIRST_HALF_${MOAB_JOBARRAYINDEX}
cp $MASTER/scripts/3.5_Chunk_Spatial_Slim_SECONDHALF.txt 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${MOAB_JOBARRAYINDEX}
sed -i "s/siteN/${SITE_N}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${MOAB_JOBARRAYINDEX}

##############
# Run burn_in
###############

# Build script
cat 1.5_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX} $SLIM_LIBRARY/2_Chunk_Slim_mod_${MOAB_JOBARRAYINDEX} 3_Chunk_Slim_mod_SIM_FIRST_HALF_${MOAB_JOBARRAYINDEX} > burn_in_SLIM_${MOAB_JOBARRAYINDEX}.slim

# Run
slim -d MUT_EFFECT=${MUT_EFFECT} burn_in_SLIM_${MOAB_JOBARRAYINDEX}.slim

###############
# Script now cycles through treatments
###############

# Treatment arrays 32 treatment total
bottleneck_array=(100 $POP_N 100 $POP_N 100 $POP_N 100 $POP_N 100 $POP_N 100 $POP_N 100 $POP_N 100 $POP_N)
finalsize_array=(0.01 0.01 0.1 0.1 0.5 0.5 1.0 1.0 0.01 0.01 0.1 0.1 0.5 0.5 1.0 1.0)
migration_array=(0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002)

# We also now iterate over genes, which starts here
for iter in {1..5}
do

#----------------------
# For loop to write scripts
for treat in $(seq 1 16)
do

# Define the treatment variables
BOTTLE_COEF=${bottleneck_array[${treat}-1]}
SIZE_COEF=${finalsize_array[${treat}-1]}
MIGRATION_COEF=${migration_array[${treat}-1]}

# Modify 1_Chunk so that coefficients are correct
sed "s/migration_coef/$MIGRATION_COEF/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX} > 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/bottleneck_coef/$BOTTLE_COEF/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/finalsize_coef/$SIZE_COEF/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/run_coef/${treat}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/ITERATION_N/${iter}/g" 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}_${treat}

# Modify Null Chunk
sed "s/migration_coef/$MIGRATION_COEF/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX} > 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/bottleneck_coef/$BOTTLE_COEF/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/finalsize_coef/$SIZE_COEF/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/run_coef/${treat}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}_${treat}
sed -i "s/ITERATION_N/${iter}/g" 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}_${treat}

# Build script for treatment
cat 1_Chunk_Slim_mod_demographics_${MOAB_JOBARRAYINDEX}_${treat} $SLIM_LIBRARY/2_Chunk_Slim_mod_${MOAB_JOBARRAYINDEX} 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${MOAB_JOBARRAYINDEX} > treatment_${treat}.slim
rm -f 1_Chunk_Slim_mod_demographics_1_*

# Build null scripts
# Slight edit for end segment to alter output
sed 's/.txt/NULL.txt/g' 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${MOAB_JOBARRAYINDEX} > 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_NULL_${MOAB_JOBARRAYINDEX}
sed -i 's/burned_geneNULL.txt/burned_gene.txt/g' 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_NULL_${MOAB_JOBARRAYINDEX} 

cat 1_Chunk_Slim_mod_NULL_${MOAB_JOBARRAYINDEX}_${treat} $SLIM_LIBRARY/2_Chunk_Slim_mod_${MOAB_JOBARRAYINDEX} 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_NULL_${MOAB_JOBARRAYINDEX} > treatment_${treat}_NULL.slim
rm -f 1_Chunk_Slim_mod_NULL_1_*

done
#----------------------

# Can now use parallel to run each treatment over a CPU and output result
parallel "slim -d MUT_EFFECT=${MUT_EFFECT}" ::: $(ls treatment*)

# End the iteration loop
done

# Concatenate all the outputs together
cat $(ls GENE_${MOAB_JOBARRAYINDEX}_*NULL.txt) > GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM_NULL.out
rm -f GENE_${MOAB_JOBARRAYINDEX}_*NULL.txt
cat $(ls GENE_${MOAB_JOBARRAYINDEX}_*.txt) > GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM.out

# Add a header
echo -e "Gene_ID\tIteration\tRun_ID\tSelection_Coef\tBottleneck\tPop2_Size\tMigration\tGene_Size\tGeneration\tFST\tdXY\tPop1_PI\tPop2_PI\tPop1_Mean\tPop1_sd\tPop2_Mean\tPop2_sd\tEvolv_Gen\tExon_Length_Total" | cat - GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM.out > tmp && mv tmp GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM.out
mv GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM.out $MASTER/demographics_model/$OUT_DIR/
sed -i '/^$/d' $MASTER/demographics_model/$OUT_DIR/GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM.out

# Add a header
echo -e "Gene_ID\tIteration\tRun_ID\tSelection_Coef\tBottleneck\tPop2_Size\tMigration\tGene_Size\tGeneration\tFST\tdXY\tPop1_PI\tPop2_PI\tPop1_Mean\tPop1_sd\tPop2_Mean\tPop2_sd\tEvolv_Gen\tExon_Length_Total" | cat - GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM_NULL.out > tmp && mv tmp GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM_NULL.out
mv GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM_NULL.out $MASTER/demographics_model/$OUT_DIR/
sed -i '/^$/d' $MASTER/demographics_model/$OUT_DIR/GENE_${MOAB_JOBARRAYINDEX}_allTreatments_SLIM_NULL.out


# Kill the folder
rm -rf $MASTER/demographics_model/sim_genes/${MOAB_JOBARRAYINDEX}
rm -f $MASTER/scripts/logs/slim_master_run_isca.txt.*.${MOAB_JOBARRAYINDEX}


