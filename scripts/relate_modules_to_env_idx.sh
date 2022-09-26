#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=50gb
#SBATCH -J relate_modules_to_env_idx
#SBATCH -o /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/marker-effects_networks/analysis/MSI_dump/%x_%j.err
#SBATCH --mail-type=NONE
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# check whether --per-intervals option was requested or not -- if not, set it to blank
[[ ${IDX_TYPE} == "intervals" ]] && INT="--per-intervals" || INT=""

# go to project folder
cd ~/projects/marker-effects_networks

echo "file with R variables from define_network_modules.R script: ${RDATA}"
echo "file containing environmental indices per environment: ${EC_FILE}"
echo "output folder: ${OUTFOLDER}"
if [[ ${INT} == "--per-intervals" ]]; then
  echo "per intervals: TRUE"
else
  echo "per intervals: FALSE"
fi
echo ""

# define network modules
Rscript scripts/relate_modules_to_env_idx.R ${RDATA} ${EC_FILE} ${OUTFOLDER} ${INT}
