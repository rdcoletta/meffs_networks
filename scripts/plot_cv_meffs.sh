#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=20gb
#SBATCH -J plot_cv_meffs
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# create output folder
OUTFOLDER=analysis/networks/${TRAIT}/meff_${MEFF_MODEL}/norm_${NORM_METHOD}
mkdir -p ${OUTFOLDER}

# plot cv distribution
Rscript scripts/plot_cv_meffs.R ${MEFF_FILE} ${OUTFOLDER} --norm-method=${NORM_METHOD}
