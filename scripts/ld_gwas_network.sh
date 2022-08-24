#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J ld_gwas_network
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# define plk filename
PLK=${LDFOLDER}/markers_gwas_network
# define ld output name
OUT=${LDFOLDER}/ld_markers_gwas_network
# hmp2plk
run_pipeline.pl -Xmx40g -importGuess ${GENODATA} \
                -includeSiteNamesInFile ${MARKERSLIST} \
                -export ${PLK} \
                -exportType Plink > /dev/null
# calculate LD within chromosomes -- but only display markers in high LD
plink --file ${PLK}.plk \
      --make-founders \
      --r2 gz dprime with-freqs \
      --ld-snp-list ${GWASHITS} \
      --ld-window-r2 0.9 \
      --ld-window 400000000 \
      --ld-window-kb 400000 \
      --out ${OUT}

# remove unnecessary files
rm ${LDFOLDER}/*.ped
rm ${LDFOLDER}/*.map
rm ${LDFOLDER}/*.nosex
