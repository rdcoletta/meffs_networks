#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=180gb
#SBATCH -J qc_hybrid_hmp
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue


module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# create folder to save results
mkdir -p ${FOLDER}
# qc hapmap file
Rscript scripts/qc_hybrid_hmp.R ${HMP} ${FOLDER}
