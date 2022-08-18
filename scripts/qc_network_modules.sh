#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J qc_network_modules
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

echo "file with R variables from define_network_modules.R script: ${RDATA}"
echo "output folder: ${OUTFOLDER}"
echo "lowest power for which the scale-free topology fit index curve: ${SFT}"
echo "number of most important markers to print for each module: ${NHUBS}"

# qc network modules
Rscript scripts/qc_network_modules.R ${RDATA} ${OUTFOLDER} ${SFT} --n-hub-markers=${NHUBS}
