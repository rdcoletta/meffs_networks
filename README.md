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


### Step 2. Pick a soft threshold

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



### Step 3. Build the network

Depending on how the marker effects were calculated, I will need to use different powers as the soft threshold for WGCNA (`20` for rrblup model and `24` for gwas model). Despite being different, both of them are much higher powers that what's tipically used in gene co-expression networks, suggesting that building marker effects networks will require tuning additional parameters downstream the pipeline. The `scripts/build_meff_network.R` builds the network for a given soft threshold.

```bash
# trait to build network
TRAIT=YLD

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # folder with results
    FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}
    # file with saved R variables from step 2
    RDATA=${FOLDER}/pick_soft_threshold.RData
    # lowest power for which the scale-free topology fit index curve
    [[ ${meff_model} == 'rrblup' ]] && SFT=20
    [[ ${meff_model} == 'gwas' ]] && SFT=24
    # submit job
    sbatch --export=RDATA=${RDATA},OUTFOLDER=${FOLDER},SFT=${SFT} scripts/build_meff_network.sh
  done
done
```


### Step 4. Define network modules

Once networks are built, I need to define the modules (or clusters) based by cutting the dendrogram at a certain height. The `cutreeDynamic()` function of WGCNA implemented in `scripts/define_network_modules.R` There are two main arguments in this function that can have a big impact on how modules are identified: `DeepSplit` (higher number results in more smaller clusters) and `pamStage` (if `FALSE` there will be more markers not assigned to any module - they will be part of the grey module). The first one may be harder to tune so I'll use WGCNA's default, but I will test how turning `pamStage = TRUE` and `pamStage = FALSE` affects module definition. Additionally, I need to determine what is the minimum number of markers to use in each module and how that affects the quality of the network.

```bash
# trait to build network
TRAIT=YLD

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # folder with results
    FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}
    # file with saved R variables from step 2
    RDATA=${FOLDER}/build_meff_network.RData
    # lowest power for which the scale-free topology fit index curve
    [[ ${meff_model} == 'rrblup' ]] && SFT=20
    [[ ${meff_model} == 'gwas' ]] && SFT=24
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # threshold to merge similar modules based on their eigengenes
        MEDISS=0.25
        # create new output folder
        OUTFOLDER=${FOLDER}/min_mod_size_${minsize}/pamStage_${pam}
        mkdir -p ${OUTFOLDER}
        # submit job
        sbatch --export=RDATA=${RDATA},OUTFOLDER=${OUTFOLDER},MINSIZE=${minsize},MEDISS=${MEDISS},SFT=${SFT},PAM=${pam} scripts/define_network_modules.sh
      done
    done
  done
done
```

> MDS plots from GWAS effects networks seem to be separating clusters much better than those from rrblup effects networks. Also, modules from the two networks (`pamStage` on and off) of `analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_100` couldn't be defined. Actually only one module was identified, which is not useful.



### Step 5. QC network modules

We checked how "noisy" the modules defined in the previous step were by looking at the distributions of `kDiff` (i.e. number of connections of a marker with other makers in the module - number of connections of same marker with markers outside its module) and the `clustering coefficient` (i.e. tendency of a marker to associate with its own module). We also calculated some network quality metrics for each module (mean connectivity, density, degree centralization, and heterogeneity).

```bash
# trait to build network
TRAIT=YLD

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # file with saved R variables from step 2
        RDATA=${FOLDER}/define_network_modules.RData
        # lowest power for which the scale-free topology fit index curve
        [[ ${meff_model} == 'rrblup' ]] && SFT=20
        [[ ${meff_model} == 'gwas' ]] && SFT=24
        # number of most important markers to print for each module
        NHUBS=10
        # submit job
        sbatch --export=RDATA=${RDATA},OUTFOLDER=${FOLDER},SFT=${SFT},NHUBS=${NHUBS} scripts/qc_network_modules.sh
      done
    done
  done
done
```

> In general, kDiff plots from GWAS effects networks seem to be have better quality modules (i.e. markers with more connections within their own module than with other modules) than those from rrblup effects networks. Turning off the `pamStage` when assigning modules also seem to result in slightly better kDiff, cluster coefficient (i.e. tendency to associate with only a select group) and TOM plots than when this option is on (more visible when GWAS effects were used).

In addition, I wrote `scripts/compare_two_networks.R` to check which modules in one network corresponds to another module in a different network. This will help us understand how different parameters chosen during network construction affects the clustering of markers, and also help identify correponding modules when correlating networks with traits.

**PAM stage on vs off**

```bash
module load R/3.6.0

# trait to build network
TRAIT=YLD
# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # folder with results
      FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}
      echo -e "\n${FOLDER}\n"
      # output folder
      OUT=analysis/networks/${TRAIT}/network_comparisons/pamStage_on-vs-off/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}
      mkdir -p ${OUT}
      # compare networks
      Rscript scripts/compare_two_networks.R \
              ${FOLDER}/pamStage_off/define_network_modules.RData \
              ${FOLDER}/pamStage_on/define_network_modules.RData \
              ${FOLDER}/pamStage_off/module_membership.txt \
              ${FOLDER}/pamStage_on/module_membership.txt \
              ${FOLDER}/pamStage_off/kDiff_per_module.txt \
              ${FOLDER}/pamStage_on/kDiff_per_module.txt \
              ${OUT} \
              --name-net1="pamStage_FALSE" \
              --name-net2="pamStage_TRUE" 2> /dev/null
    done
  done
done
```

> As expected, using `pamStage = FALSE` results in more markers assigned to the grey module compared to `pamStage = TRUE`, but the remaining modules have a pretty high module membership correlation between networks. Usually, there's a 1:1 to relationship (i.e. one module in one network can be identified in just one module from another network), but sometimes a module from `pamStage = FALSE` can be found in two or more modules from `pamStage = TRUE` (especially for networks with high total number of modules).

**Minimum number of markers per module**

```bash
module load R/3.6.0

# trait to build network
TRAIT=YLD
# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # define modules with and without the PAM stage
    for pam in on off; do
      # folder with results
      FOLDER25=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_25/pamStage_${pam}
      FOLDER50=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_50/pamStage_${pam}
      FOLDER100=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_100/pamStage_${pam}
      # output folders
      OUT25v50=analysis/networks/${TRAIT}/network_comparisons/minsize_25-vs-50/meff_${meff_model}/norm_${norm_method}/pamStage_${pam}
      mkdir -p ${OUT25v50}
      OUT50v100=analysis/networks/${TRAIT}/network_comparisons/minsize_50-vs-100/meff_${meff_model}/norm_${norm_method}/pamStage_${pam}
      mkdir -p ${OUT50v100}
      OUT25v100=analysis/networks/${TRAIT}/network_comparisons/minsize_25-vs-100/meff_${meff_model}/norm_${norm_method}/pamStage_${pam}
      mkdir -p ${OUT25v100}
      # compare 25 vs 50
      echo -e "\n${FOLDER25}\n${FOLDER50}\n"
      Rscript scripts/compare_two_networks.R \
              ${FOLDER25}/define_network_modules.RData \
              ${FOLDER50}/define_network_modules.RData \
              ${FOLDER25}/module_membership.txt \
              ${FOLDER50}/module_membership.txt \
              ${FOLDER25}/kDiff_per_module.txt \
              ${FOLDER50}/kDiff_per_module.txt \
              ${OUT25v50} \
              --name-net1="min_mod_size_25" \
              --name-net2="min_mod_size_50" 2> /dev/null
      # compare 50 vs 100
      echo -e "\n${FOLDER50}\n${FOLDER100}\n"
      Rscript scripts/compare_two_networks.R \
              ${FOLDER50}/define_network_modules.RData \
              ${FOLDER100}/define_network_modules.RData \
              ${FOLDER50}/module_membership.txt \
              ${FOLDER100}/module_membership.txt \
              ${FOLDER50}/kDiff_per_module.txt \
              ${FOLDER100}/kDiff_per_module.txt \
              ${OUT50v100} \
              --name-net1="min_mod_size_50" \
              --name-net2="min_mod_size_100" 2> /dev/null
      # compare 25 vs 100
      echo -e "\n${FOLDER25}\n${FOLDER100}\n"
      Rscript scripts/compare_two_networks.R \
              ${FOLDER25}/define_network_modules.RData \
              ${FOLDER100}/define_network_modules.RData \
              ${FOLDER25}/module_membership.txt \
              ${FOLDER100}/module_membership.txt \
              ${FOLDER25}/kDiff_per_module.txt \
              ${FOLDER100}/kDiff_per_module.txt \
              ${OUT25v100} \
              --name-net1="min_mod_size_25" \
              --name-net2="min_mod_size_100" 2> /dev/null
    done
  done
done
```

> These comparisons are harder to interpret, but overall it looks like there are more modules in one network that are split in the other (and vice-versa), especially when `pamStage = TRUE` (in this case, you can see that some of the split modules don't have a very high correlation with its counterpart(s) in the other network). When `pamStage = FALSE`, you generally see more modules represented in the grey module of the other network.

**Minmax vs zscore normalization**

```bash
module load R/3.6.0

# trait to build network
TRAIT=YLD
# marker effects from rrblup
for meff_model in rrblup gwas; do
  # minimum number of markers per module
  for minsize in 25 50 100; do
    # define modules with and without the PAM stage
    for pam in on off; do
      # folder with results
      FOLDER_M=analysis/networks/${TRAIT}/meff_${meff_model}/norm_minmax/min_mod_size_${minsize}/pamStage_${pam}
      FOLDER_Z=analysis/networks/${TRAIT}/meff_${meff_model}/norm_zscore/min_mod_size_${minsize}/pamStage_${pam}
      echo -e "\n${FOLDER_M}\n${FOLDER_Z}\n"
      # output folder
      OUT=analysis/networks/${TRAIT}/network_comparisons/minmax-vs-zscore/meff_${meff_model}/min_mod_size_${minsize}/pamStage_${pam}
      mkdir -p ${OUT}
      # compare networks
      Rscript scripts/compare_two_networks.R \
              ${FOLDER_M}/define_network_modules.RData \
              ${FOLDER_Z}/define_network_modules.RData \
              ${FOLDER_M}/module_membership.txt \
              ${FOLDER_Z}/module_membership.txt \
              ${FOLDER_M}/kDiff_per_module.txt \
              ${FOLDER_Z}/kDiff_per_module.txt \
              ${OUT} \
              --name-net1="norm_minmax" \
              --name-net2="norm_zscore" \
              --diff-n-markers 2> /dev/null
    done
  done
done
```

> Comparisons are even harder here because different normalizations had different CV cutoffs and not many markers overlap between the networks (~20-40%). Module relationship rarely follow a 1:1 ratio, so other metrics (e.g. kDiff) may be more informative about which network has better quality.

**rrBLUP vs GWAS effects**

```bash
module load R/3.6.0

# trait to build network
TRAIT=YLD
# type of normalization method to use
for norm_method in minmax zscore; do
  # minimum number of markers per module
  for minsize in 25 50 100; do
    # define modules with and without the PAM stage
    for pam in on off; do
      # folder with results
      FOLDER_R=analysis/networks/${TRAIT}/meff_rrblup/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
      FOLDER_G=analysis/networks/${TRAIT}/meff_gwas/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
      echo -e "\n${FOLDER_R}\n${FOLDER_G}\n"
      # output folder
      OUT=analysis/networks/${TRAIT}/network_comparisons/rrblup-vs-gwas/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
      mkdir -p ${OUT}
      # compare networks
      Rscript scripts/compare_two_networks.R \
              ${FOLDER_R}/define_network_modules.RData \
              ${FOLDER_G}/define_network_modules.RData \
              ${FOLDER_R}/module_membership.txt \
              ${FOLDER_G}/module_membership.txt \
              ${FOLDER_R}/kDiff_per_module.txt \
              ${FOLDER_G}/kDiff_per_module.txt \
              ${OUT} \
              --name-net1="meff_rrblup" \
              --name-net2="meff_gwas" \
              --diff-n-markers 2> /dev/null
    done
  done
done
```

> There's more overlap of markers between different marker effect type (~40-60%) than between normalization methods (~20-40%), but there's still a lot of modules that are split between marker effect types (i.e. not a 1:1 ratio).




## Relate modules to phenotypic and environmental data

To validate whether the modules I found have any biological meaning, we need to test if there are any modules that are stastically associated with the trait of interest. The `scripts/relate_modules_to_pheno.sh` does just that and outputs the p-values of correlation tests (corrected for multiple testing) and of permutations.

```bash
# trait to build network
TRAIT=YLD
# phenotypic data to correlate to modules
TRAIT_FILE=data/1stStage_BLUEs.YLD-per-env.txt

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # file with saved R variables from step 2
        RDATA=${FOLDER}/define_network_modules.RData
        # submit job
        sbatch --export=RDATA=${RDATA},TRAIT_FILE=${TRAIT_FILE},OUTFOLDER=${FOLDER} scripts/relate_modules_to_pheno.sh
      done
    done
  done
done
```

> Overall, correlation tests were more stringent than permutations (meaning there were more modules significantly associated with trait when doing permutations).


### Calculate LD for markers within each module

Another important thing to check how many markers in a module are in LD with each other. If all of them are in perfect LD, then the marker effect networks would be just picking up LD without any other biological reason for grouping those markers. The `scripts/ld_markers_modules.sh` calculates LD all markers within a module (even if they are in different chromosomes).

```bash
# trait to build network
TRAIT=YLD
# genotypic data of markers
MARKERS=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # create folder to keep list of markers per module
        MODFOLDER=${FOLDER}/modules_ld
        mkdir -p ${MODFOLDER}
        # get module summary
        MODSUMMARY=${FOLDER}/kDiff_per_module.txt
        # submit job
        sbatch --export=MODSUMMARY=${MODSUMMARY},MODFOLDER=${MODFOLDER},MARKERS=${MARKERS} scripts/ld_markers_modules.sh
      done
    done
  done
done
```

In addition, I want to make sure that all (or most) GWAS hits and top 100 non-significant hits are represented in the network (i.e. if the actual GWAS hit is not present in the network due to filtering, I want to see if at least one marker in LD with that hit is present).

```bash
# get non-redudant set of gwas hits and top non-significant markers for each environment
GWASHITS=data/list_gwas_top-peaks_YLD-per-env.txt
cut -f 2 ~/projects/genomic_prediction/hybrids/analysis/gwas/gwas_top-peaks_YLD-per-env.txt | sed 1d | sort | uniq > ${GWASHITS}

# trait to build network
TRAIT=YLD
# genotypic data of markers
GENODATA=data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # create folder to keep list of markers per module
        LDFOLDER=${FOLDER}/gwas_network_ld
        mkdir -p ${LDFOLDER}
        # get list of network modules
        NETMARKERS=$(cut -f 1 ${FOLDER}/kDiff_per_module.txt | sed 1d | sort | uniq)
        # merge top hits with all markers in network
        MARKERSLIST=${LDFOLDER}/markers_gwas_network_merged.txt
        echo $(cat ${GWASHITS}) ${NETMARKERS} | tr " " "\n" | sort | uniq > ${MARKERSLIST}
        # submit job
        sbatch --export=LDFOLDER=${LDFOLDER},GENODATA=${GENODATA},GWASHITS=${GWASHITS},MARKERSLIST=${MARKERSLIST} scripts/ld_gwas_network.sh
      done
    done
  done
done
```



### Check in which modules GWAS hits for trait are located

Finally, we can test whether or not modules associated with a trait are also enriched for GWAS hits. If that's the case, then it would be another evidence that the modules identified in the networks have biological meaning and could potentially be used to identify markers associated with trait variability across environments.

First, I'll run `scripts/check_ld_network_gwas-hits.sh` to check which network markers are in LD with gwas hits or top non-significant hits.

```bash
# trait to build network
TRAIT=YLD
# non-redudant set of gwas hits and top non-significant markers for each environment
GWASHITS=~/projects/genomic_prediction/hybrids/analysis/gwas/gwas_top-peaks_YLD-per-env.txt

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # get module summary
        MODSUMMARY=${FOLDER}/kDiff_per_module.txt
        # get ld file per module
        LDMOD=${FOLDER}/gwas_network_ld/ld_markers_gwas_network.ld.gz
        # name of output file
        OUT=${FOLDER}/kDiff_per_module.gwas-status.txt
        if [[ -e ${MODSUMMARY} ]]; then
          # submit job
          sbatch --export=MODSUMMARY=${MODSUMMARY},LDMOD=${LDMOD},GWASHITS=${GWASHITS},OUT=${OUT} scripts/check_ld_network_gwas-hits.sh
        fi
      done
    done
  done
done
```

Then, the `scripts/enrichment_gwas-hits_meff-modules.sh` perform enrichment tests for GWAS hits and plot both the connections of each marker within a module (highlighting the GWAS hits and hub markers or which markers are in LD with each other). For plotting, only connections of markers with adjacency values above 0.1 will be shown (adjacency values vary between 0 and 1).

```bash
# trait to build network
TRAIT=YLD

# marker effects from rrblup
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # marker effects file
        MEFF_FILE=analysis/marker_effects/${TRAIT}/marker_effects.${meff_model}.txt
        # lowest power for which the scale-free topology fit index curve
        [[ ${meff_model} == 'rrblup' ]] && SFT=20
        [[ ${meff_model} == 'gwas' ]] && SFT=24
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # folder with ld per module
        MODFOLDER=${FOLDER}/modules_ld
        # get module summary
        MODSUMMARY=${FOLDER}/kDiff_per_module.gwas-status.txt
        # get p-values module-trait association
        MODPVALS=${FOLDER}/module-pheno_pvals.txt
        # adjacency threshold for including edges in the network plots
        EDGE=0.1
        if [[ -e ${MODSUMMARY} ]]; then
          # submit job
          sbatch --export=MEFF_FILE=${MEFF_FILE},MODSUMMARY=${MODSUMMARY},MODPVALS=${MODPVALS},MODFOLDER=${MODFOLDER},OUTFOLDER=${FOLDER},SFT=${SFT},EDGE=${EDGE} scripts/enrichment_gwas-hits_meff-modules.sh
        fi
      done
    done
  done
done
```
