#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J ld_markers_modules
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/marker-effects_networks

# get module names
modules=$(awk '$7 == "adjacency"' ${MODSUMMARY} | cut -f 6 | sort | uniq)
# loop through modules
for mod in ${modules[@]}; do
  # create list of markers to calculate LD
  MODMARKERS=${MODFOLDER}/list_markers_${mod}.txt
  awk '$7 == "adjacency"' ${MODSUMMARY} | awk -v mod="${mod}" '$6 == mod' | cut -f 1 > ${MODMARKERS}
  # define plk filename
  PLK=${MODFOLDER}/markers_geno_${mod}
  # define ld output name
  OUT=${MODFOLDER}/ld_markers_${mod}
  # hmp2plk
  run_pipeline.pl -Xmx40g -importGuess ${MARKERS} \
                  -includeSiteNamesInFile ${MODMARKERS} \
                  -export ${PLK} \
                  -exportType Plink > /dev/null
  # calculate LD across chromosomes
  plink --file ${PLK}.plk \
        --make-founders \
        --r2 gz dprime with-freqs inter-chr \
        --ld-window-r2 0 \
        --out ${OUT}
done

# remove unnecessary files
rm ${MODFOLDER}/*.ped
rm ${MODFOLDER}/*.map
rm ${MODFOLDER}/*.nosex
