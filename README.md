# Marker effects networks
by Rafael Della Coletta and Candice Hirsch

> The objective of this project is to determine whether we can create networks of marker effects from a maize population across multiple environments. If this is possible, then the question becomes whether or not we can extract useful information about which markers (or modules of markers) are more associated with environmental adaptability.




## Project folder

All data, scripts, analyses, notes and other things related to this project is located on my MSI account:

```bash
cd /home/hirschc1/della028/projects/marker-effects_networks
mkdir -p {analysis,data,scripts,tests}
mkdir -p analysis/MSI_dump
```




## Genotypic data

We will estimate marker effects (SNPs and SVs) of markers of hybrid maize lines used on my genomic prediction project. The data was already LD pruned with PLINK with different window sizes (`--indep-pairwise [10kb|100kb|1000kb] 1 0.9 --geno 0.25`).

```bash
# transfer data
cd ~/projects
cp genomic_prediction/hybrids/data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-*kb.geno-miss-0.25.hmp.txt marker-effects_networks/data/
cp genomic_prediction/hybrids/data/SNPs_IDs_poly_hybrids.txt marker-effects_networks/data/
cp genomic_prediction/hybrids/data/SVs_IDs_poly.txt marker-effects_networks/data/

# go back to project folder
cd ~/projects/marker-effects_networks
```

Number of markers remaining after pruning the genotypic dataset with different window sizes:

| window size | SVs remaining | SNPs remaining |
| ----------- | ------------- | -------------- |
| 10 kb       | 3,110         | 359,990        |
| 100 kb      | 2,300         | 103,295        |
| 1000 kb     | 1,067         | 18,606         |



## Phenotypic data

We also need phenotypic data to estimate the effects markers. For this purpose, I will use phenotypic data from my genomic prediction project where five traits (EHT, Moisture, PHT, TWT, YLD) of ~400 hybrids were evaluated in 4-11 environments.

```bash
# transfer data
cd ~/projects
cp genomic_prediction/hybrids/data/1stStage_BLUEs.*-per-env.txt marker-effects_networks/data/

# go back to project folder
cd ~/projects/marker-effects_networks
```




## Estimate marker effects

We will estimate marker effects for the five traits using four different models to see how much the network changes based on the type of model selected. The rationale for using more than one model is that different models may be better for a certain genetic architecture, which is unknown from empirical data. In addition, we chose to use marker dataset pruned in 100kb windows to balance the number of markers vs computation time.

```bash
# estimate effects
MARKERS=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt
for trait in EHT Moisture PHT TWT YLD; do
  BLUES=data/1stStage_BLUEs.${trait}-per-env.txt
  OUTFOLDER=analysis/marker_effects/${trait}
  for marker_eff_model in rrblup bayescpi mrr gwas; do
    sbatch --export=MARKERS=${MARKERS},BLUES=${BLUES},OUTFOLDER=${OUTFOLDER},MEFFMODEL=${marker_eff_model} scripts/estimate_marker_effects.sh
  done
done

# plot results for qc
module load R/3.6.0
Rscript scripts/plot_marker_effects.R analysis/marker_effects \
                                      --traits=EHT,Moisture,PHT,TWT,YLD \
                                      --models=rrblup,bayescpi,mrr,gwas
Rscript scripts/plot_marker_effects.R analysis/marker_effects \
                                      --traits=EHT,Moisture,PHT,TWT,YLD \
                                      --models=rrblup,bayescpi,mrr,gwas \
                                      --no-missing-genotypes
```
