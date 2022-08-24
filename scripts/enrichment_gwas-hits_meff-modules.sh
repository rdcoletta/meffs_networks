#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J enrichment_gwas-hits_meff-modules
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# check gwas hits in modules
Rscript scripts/enrichment_gwas-hits_meff-modules.R ${MEFF_FILE} ${MODSUMMARY} ${MODPVALS} ${MODFOLDER} ${OUTFOLDER} --soft-threshold=${SFT} --edge-threshold=${EDGE}
