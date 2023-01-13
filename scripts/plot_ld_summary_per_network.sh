#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=120gb
#SBATCH -J plot_ld_summary_per_network
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# summarize ld
Rscript scripts/plot_ld_summary_per_network.R ${LDFILE} ${OUTFOLDER}
