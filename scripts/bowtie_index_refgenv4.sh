#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=60gb
#SBATCH -J bowtie_index_refgenv4
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load bowtie/1.1.2

# go to refgen v4 folder
cd ~/projects/marker-effects_networks/data/refgen/B73v4

# build index
bowtie-build chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa index
