#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=100gb
#SBATCH -J build_meff_network
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

echo "file with R variables from pick_soft_threshold.R script: ${RDATA}"
echo "output folder: ${OUTFOLDER}"
echo "lowest power for which the scale-free topology fit index curve: ${SFT}"
echo ""

# estimate effects - all data
Rscript scripts/build_meff_network.R ${RDATA} ${OUTFOLDER} --soft-threshold=${SFT}
