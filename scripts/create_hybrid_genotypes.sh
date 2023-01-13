#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=110gb
#SBATCH -J create_hybrid_genotypes
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

Rscript scripts/create_hybrid_genotypes.R ${HYBINFO} ${RILGENO} ${OUT} --het-parents=${HET} --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P
