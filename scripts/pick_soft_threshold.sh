#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J pick_soft_threshold
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# pick soft threshold
Rscript scripts/pick_soft_threshold.R ${MEFF_FILE} ${MARKERS} ${OUTFOLDER} --norm-method=${NORM_METHOD} --cv-threshold=${CV_THRESHOLD}
