#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=80gb
#SBATCH -J summarize_ld_per_network
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# summarize ld
Rscript scripts/summarize_ld_per_network.R ${FOLDER} --meff-model=${MEFF_MODEL} --norm-method=${NORM_METHOD} --minsize=${MINSIZE} --pamStage=${PAM}
