#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J define_network_modules
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# check whether --pamStage option was requested or not -- if not, set it to blank
[[ ${PAM} == "on" ]] && PAM="--pamStage" || PAM=""

# go to project folder
cd ~/projects/marker-effects_networks

echo "file with R variables from build_meff_network.R script: ${RDATA}"
echo "output folder: ${OUTFOLDER}"
echo "minimum number of markers to be in a module: ${MINSIZE}"
echo "threshold to merge similar modules based on their eigengenes: ${MEDISS}"
if [[ ${PAM} == "--pamStage" ]]; then
  echo "pamStage: ON"
else
  echo "pamStage: OFF"
fi
echo ""

# define network modules
Rscript scripts/define_network_modules.R ${RDATA} ${OUTFOLDER} --min-mod-size=${MINSIZE} --ME-diss-threshold=${MEDISS} --soft-threshold=${SFT} ${PAM}
