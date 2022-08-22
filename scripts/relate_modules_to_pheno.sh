#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J relate_modules_to_pheno
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

echo "file with R variables from define_network_modules.R script: ${RDATA}"
echo "file containing trait values per environment: ${TRAIT_FILE}"
echo "output folder: ${OUTFOLDER}"
echo ""

# define network modules
Rscript scripts/relate_modules_to_pheno.R ${RDATA} ${TRAIT_FILE} ${OUTFOLDER}
