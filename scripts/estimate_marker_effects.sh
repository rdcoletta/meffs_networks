#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=40gb
#SBATCH -J estimate_marker_effects
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

echo "file with markers: ${MARKERS}"
echo "file with phenotypes: ${BLUES}"
echo "output folder: ${OUTFOLDER}"
echo "marker effect model: ${MEFFMODEL}"
echo ""

# estimate effects - all data
Rscript scripts/estimate_marker_effects.R ${MARKERS} ${BLUES} ${OUTFOLDER} \
                                          --marker-eff-model=${MEFFMODEL}

echo ""

# estimate effects - only non-missing genotypes
Rscript scripts/estimate_marker_effects.R ${MARKERS} ${BLUES} ${OUTFOLDER} \
                                          --marker-eff-model=${MEFFMODEL} \
                                          --no-missing-genotypes
