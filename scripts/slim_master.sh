#!/bin/bash

#---------------------------
# NOTES:
# This pipeline has been written to run as a batch array on an HPC.
# It therefore applies an ID automatically to each simulated gene, running all demographic treatments per gene ID.
# Depending on the system, the identifier for this variable will differ. Here it is set through MOAB_JOBARRAYINDEX, but this can be set in the settings section below.
# If not using an HPC, this script can be modified to a forloop, in which SLIM_RUN_ID is set for each iteration of the loop

# Script requires GNU Parallel [Tange, O. (2018). Gnu parallel 2018. Lulu. com.]

# Example of PBS flags for HPC is given below:
# Example spawns array with 1000 jobs (throttled to 10), giving each 3 hrs walltime, 16 processors per job
#---------------------------

##PBS -d . 
##PBS -l walltime=3:00:00
##PBS -l nodes=1:ppn=16
##PBS -N slim_run 
##PBS -e logs/slim_run.err.txt 
##PBS -o logs/slim_run.out.txt 
##PBS -V
##PBS -t 1-1000%10

#---------------------------
# SET THE FOLLOWING:
# Working directory
# Path to GFF
# Population size (Ancestral) eg. 1000
# Mutation rate eg. "4.89e-6"
# Name of output directory
#---------------------------

## WORKING DIRECTORY
MASTER=

## BASIC PARAMETERS
GFF=
POP_N=1000
MUTATION="4.89e-6"
OUT_DIR=

## JOB ID IDENTIFIER - System-specific (see note above)
SLIM_RUN_ID=${MOAB_JOBARRAYINDEX}

#---------------------------

# Make output folder if it does not already exist
mkdir $MASTER/demographics_model/$OUT_DIR

# Run the analysis for each gene in a contained directory
mkdir $MASTER/demographics_model/sim_genes/${SLIM_RUN_ID}
cd $MASTER/demographics_model/sim_genes/${SLIM_RUN_ID}

# Run build_a_SLIM_gene.sh to construct the gene to run through SLIM
bash $MASTER/scripts/build_a_SLIM_gene.sh $GFF ${SLIM_RUN_ID}

# Get Exon size
EXON_SIZE=$(awk '{ sum += $1 } END { print sum }' exons.txt)

# Assign a selection coefficient at random, selection coefficient corresponds to the sd of the dnorm function for fitness effects of distance from the optimum
selection_array=(0.1 0.5 1.0 5.0 10)
SEL_COEF=${selection_array[$RANDOM % ${#selection_array[@]} ]}

# Modify 1_Chunk to include new selection coefficient and pop_size and mutation rate
sed "s/selection_coef/${SEL_COEF}/g" $MASTER/scripts/1_Chunk_Spatial_Slim.txt > 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/pop_coef/${POP_N}/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/ID_coef/${SLIM_RUN_ID}/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/MUT_RATE/${MUTATION}/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/EXON_L/${EXON_SIZE}/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}

# Modify 1.5_Chunk to include new selection coefficient and pop_size and mutation rate
sed "s/selection_coef/${SEL_COEF}/g" $MASTER/scripts/1.5_Chunk_Spatial_Slim.txt > 1.5_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/pop_coef/${POP_N}/g" 1.5_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/ID_coef/${SLIM_RUN_ID}/g" 1.5_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/MUT_RATE/${MUTATION}/g" 1.5_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}
sed -i "s/EXON_L/${EXON_SIZE}/g" 1.5_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}

# We also need to seed the dXY function with the number of sites
SITE_N=$(tail -2 2_Chunk_Slim_mod_${SLIM_RUN_ID} | head -n1 | awk '{print $3}' | sed "s/);//g")
sed -i "s/siteN/${SITE_N}/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}

# Modify 3_Chunks so burn-in time = 5 x 2 N
BURN_GEN=$((5*2*POP_N))
sed "s/BURN_GEN/${BURN_GEN}/g" $MASTER/scripts/3_Chunk_Spatial_Slim_FIRSTHALF.txt > 3_Chunk_Slim_mod_SIM_FIRST_HALF_${SLIM_RUN_ID}

# Add Generation times to 3_5 file. Generation times are defined
GEN_1_1=$(echo "(((0.05 * 2 * $POP_N) - 50))/1" | bc)
GEN_1_2=$(echo "((0.05*2*$POP_N))/1" | bc)
GEN_2_1=$(echo "(((0.158*2*$POP_N)-50))/1" | bc)
GEN_2_2=$(echo "((0.158*2*$POP_N))/1" | bc)
GEN_3_1=$(echo "(((1.581*2*$POP_N)-50))/1" | bc)
GEN_3_2=$(echo "((1.581*2*$POP_N))/1" | bc)
GEN_4_1=$(echo "(((5*2*$POP_N)-50))/1" | bc)
GEN_4_2=$(echo "((5*2*$POP_N))/1" | bc)

if [ $GEN_1_1 -eq 0 ]
then 
GEN_1_1=2
GEN_1_2=52
fi

sed "s/GEN_1_1/${GEN_1_1}/g" $MASTER/scripts/3.5_Chunk_Spatial_Slim_SECONDHALF.txt > 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_1_2/${GEN_1_2}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_2_1/${GEN_2_1}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_2_2/${GEN_2_2}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_3_1/${GEN_3_1}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_3_2/${GEN_3_2}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_4_1/${GEN_4_1}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/GEN_4_2/${GEN_4_2}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}
sed -i "s/siteN/${SITE_N}/g" 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID}

##############
# Run burn_in
###############

# Build script
cat 1.5_Chunk_Slim_mod_demographics_${SLIM_RUN_ID} 2_Chunk_Slim_mod_${SLIM_RUN_ID} 3_Chunk_Slim_mod_SIM_FIRST_HALF_${SLIM_RUN_ID} > burn_in_SLIM_${SLIM_RUN_ID}.slim

# Run
slim burn_in_SLIM_${SLIM_RUN_ID}.slim

###############
# Script now cycles through treatments
###############

# Treatment arrays 32 treatment total
bottleneck_array=(2 20 100 $POP_N 2 20 100 $POP_N 2 20 100 $POP_N 2 20 100 $POP_N 2 20 100 $POP_N 2 20 100 $POP_N 2 20 100 $POP_N 2 20 100 $POP_N)
finalsize_array=(0.01 0.01 0.01 0.01 0.1 0.1 0.1 0.1 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 0.01 0.01 0.01 0.01 0.1 0.1 0.1 0.1 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0)
migration_array=(0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002 0.002)

#----------------------
# For loop to write scripts
for treat in $(seq 1 32)
do

# Define the treatment variables
BOTTLE_COEF=${bottleneck_array[${treat}-1]}
SIZE_COEF=${finalsize_array[${treat}-1]}
MIGRATION_COEF=${migration_array[${treat}-1]}

# Modify 1_Chunk so that coefficients are correct
sed "s/migration_coef/$MIGRATION_COEF/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID} > 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}_${treat}
sed -i "s/bottleneck_coef/$BOTTLE_COEF/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}_${treat}
sed -i "s/finalsize_coef/$SIZE_COEF/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}_${treat}
sed -i "s/run_coef/${treat}/g" 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}_${treat}

# Build script for treatment
cat 1_Chunk_Slim_mod_demographics_${SLIM_RUN_ID}_${treat} 2_Chunk_Slim_mod_${SLIM_RUN_ID} 3.5_Chunk_Slim_mod_SIM_SECOND_HALF_${SLIM_RUN_ID} > treatment_${treat}.slim
rm -f 1_Chunk_Slim_mod_demographics_1_*
done
#----------------------

# Can now use parallel to run each treatment over a CPU and output result
parallel "slim" ::: $(ls treatment*)

# Concatenate all the outputs together
cat $(ls GENE_${SLIM_RUN_ID}_*) > GENE_${SLIM_RUN_ID}_allTreatments_SLIM.out

# Add a header
echo -e "Gene_ID\tRun_ID\tSelection_Coef\tBottleneck\tPop2_Size\tMigration\tGene_Size\tGeneration\tFST\tdXY\tPop1_PI\tPop2_PI\tPop1_Mean\tPop1_sd\tPop2_Mean\tPop2_sd\tEvolv_Gen\tExon_Length_Total" | cat - GENE_${SLIM_RUN_ID}_allTreatments_SLIM.out > tmp && mv tmp GENE_${SLIM_RUN_ID}_allTreatments_SLIM.out
mv GENE_${SLIM_RUN_ID}_allTreatments_SLIM.out $MASTER/demographics_model/$OUT_DIR/
sed -i '/^$/d' $MASTER/demographics_model/$OUT_DIR/GENE_${SLIM_RUN_ID}_allTreatments_SLIM.out

# Retain the gene chunk
mv 2_Chunk_Slim_mod_${SLIM_RUN_ID} $MASTER/demographics_model/$OUT_DIR/

# Kill the folder
rm -rf $MASTER/demographics_model/sim_genes/${SLIM_RUN_ID}
