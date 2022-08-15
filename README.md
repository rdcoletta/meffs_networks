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

> BayesCpi model didn't perform as well as the other models. It will be interesting to see the networks that come out of the marker effects from this model. It may be a control for a "bad" input to build networks.




## Build networks

Due to the exploratory nature of this project, we'll test different parameters when building networks with [WGCNA](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html). This R package was developed and optimized to build co-expression networks of genes. While networks can be built with a different type of data, it is still unknown if marker effects networks would behave similarly to co-expression networks (i.e. have a scale-free topology). Another concern is how to transform/normalize the marker effect data (if needed in the first place), since the distribution of marker effects are not skewed as gene expression data.


### Step 1. Choose a coefficient of variation filter

Markers with very low coefficient of variation (CV) across samples are not very informative when building networks and should be removed. The `scripts/plot_cv_meffs.R` plots the distribution of CV scores for all the markers so I can decide how many markers to remove. The script will also filters out duplicated markers (i.e. only one marker from a group of markers with exact same effects will be kept) to reduce redundancy, reduce total amount of memory required to build the network and also speed up network construction.

```bash
# trait to build network
TRAIT=YLD
# model used to calculate marker effects
for meff_model in rrblup gwas; do
  # marker effects file
  MEFF_FILE=analysis/marker_effects/${TRAIT}/marker_effects.${meff_model}.txt
  # type of normalization method to use
  for norm_method in minmax zscore none; do
    sbatch --export=TRAIT=${TRAIT},MEFF_FILE=${MEFF_FILE},MEFF_MODEL=${meff_model},NORM_METHOD=${norm_method} scripts/plot_cv_meffs.sh
  done
done
```

#### Step 2. Pick a soft threshold

In order to build networks with WGCNA, I need to specify a power to model the scale-free topology of the network. The `scripts/pick_soft_threshold.R` plots the scale-free topology fit index for different powers so I can decide which one is more adequate for my data.

```bash
# marker file
MARKERS=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt
# trait to build network
TRAIT=YLD

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # marker effects file
  MEFF_FILE=analysis/marker_effects/${TRAIT}/marker_effects.${meff_model}.txt
  # type of normalization method to use
  for norm_method in minmax zscore none; do
    # adjust cv threshold to type of normalization
    [[ ${norm_method} == 'minmax' ]] && CV_THRESHOLD=0.12
    [[ ${norm_method} == 'zscore' ]] && CV_THRESHOLD=0.75
    [[ ${norm_method} == 'none' ]] && CV_THRESHOLD=0.75
    # create output folder
    OUTFOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}
    mkdir -p ${OUTFOLDER}
    # submit job
    sbatch --export=MARKERS=${MARKERS},MEFF_FILE=${MEFF_FILE},OUTFOLDER=${OUTFOLDER},NORM_METHOD=${norm_method},CV_THRESHOLD=${CV_THRESHOLD} scripts/pick_soft_threshold.sh
  done
done
```

> Based on preliminary testing, using a more stringent CV cutoff resulted in better scale-free topology fit. Also, Z-score and no normalization yield the same results here and in my preliminary tests, so I'll just run `minmax`- and `zscore`-normalized data.
