# Marker effects networks
by Rafael Della Coletta and Candice Hirsch

> The objective of this project was to determine whether or not we could create networks of marker effects from a maize population across multiple environments, and identify which markers (or modules of markers) were associated with environmental adaptability. For more details about this analysis, please read the manuscript: https://www.biorxiv.org/content/10.1101/2023.01.19.524532v1

<!-- TOC START min:1 max:3 link:true asterisk:false update:true -->
- [Marker effects networks](#marker-effects-networks)
  - [Important notes!](#important-notes)
  - [Genotypic data](#genotypic-data)
    - [Hapmap format](#hapmap-format)
    - [Data QC](#data-qc)
    - [Convert SNP coordinates from B73v2 to B73v4](#convert-snp-coordinates-from-b73v2-to-b73v4)
    - [Correct miscalls with sliding window approach](#correct-miscalls-with-sliding-window-approach)
    - [Generate hybrid genotypes](#generate-hybrid-genotypes)
    - [Prune markers by LD](#prune-markers-by-ld)
  - [Phenotypic data](#phenotypic-data)
  - [Estimate marker effects](#estimate-marker-effects)
  - [Build networks](#build-networks)
    - [Step 1. Choose a coefficient of variation filter](#step-1-choose-a-coefficient-of-variation-filter)
    - [Step 2. Pick a soft threshold](#step-2-pick-a-soft-threshold)
    - [Step 3. Build the network](#step-3-build-the-network)
    - [Step 4. Define network modules](#step-4-define-network-modules)
    - [Step 5. QC network modules](#step-5-qc-network-modules)
    - [Step 6. Calculate LD for markers within each module](#step-6-calculate-ld-for-markers-within-each-module)
  - [Relate modules to environmental covariables](#relate-modules-to-environmental-covariables)
<!-- TOC END -->



## Important notes!

All data necessary to replicate this analysis are available in the Data Repository for the University of Minnesota (DRUM): **ADD DOI HERE!**. When submitting the files for publication, I named them Files S1-S6 for easier reference within the manuscript. However, please rename them to their original filenames before running any of the scripts according to this table:

| File | Name                                                                 |
| ---- | -------------------------------------------------------------------- |
| S1   | NIFA_CompleteDataset.csv                                             |
| S2   | 1stStage_BLUEs.YLD-per-env.txt                                       |
| S3   | usda_22kSNPs_parents.hmp.txt                                         |
| S4   | usda_22kSNPs_rils.hmp.txt                                            |
| S5   | usda_hybrids_SNP-chip.maf-filter.pruned-100kb.geno-miss-0.25.hmp.txt |
| S6   | env_covariables_means_per_intervals.txt                              |

All scripts necessary for data analysis are available at the `scripts` folder in this GitHub repository, and I recommend that you create a folder called `data` to keep the files above and another called `analysis` to keep any other files generated throughout this analysis.

If you want to skip the genotypic and phenotypic data quality control steps, and start building the marker effect networks, then jump to the [Estimate marker effects](#estimate-marker-effects) section. You will need files S2, S5 and S6 (properly renamed according to the table above).

All the analysis, unless indicated otherwise, were run using the Minnesota Supercomputing Institute (MSI) clusters, which uses the SLURM scheduler. You may need to adapt some scripts if you use a different scheduler.

Finally, here is a list of software necessary to run the analysis and their respective versions:

| Software | Version |
| -------- | ------- |
| R        | 3.6.0   |
| Python   | 3.6.6   |
| GNU bash | 4.2.46  |
| Tassel   | 5.2.56  |
| Plink    | 1.9     |
| Bowtie   | 1.1.2   |



## Genotypic data

We will estimate marker effects of SNPs from the hybrid maize lines, but we only have genotypic data (20k Infinium chip data) for the 333 RILs used to generate the hybrids. Thus, we have to create hybrid genotypes based on genotypes of parental (RIL) data *in silico*.


### Hapmap format

Before doing any analysis, we need to transform the genotypic data that comes from the SNP chip into the [HapMap format](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load) for more info about the format). I used the `scripts/usda_geno2hmp.R` script to transform all parental and RIL data to hapmap format:

> For the manuscript, I provided the hapmap version of the raw genotypic datasets, so you can skip the first two steps (i.e., running `scripts/usda_geno2hmp.R` and `scripts/remove_extra_geno-data.R`) and start at the Tassel filtering step.

```bash
module load R/3.6.0

Rscript scripts/usda_geno2hmp.R data/SNP_chip/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv \
                                data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318_corrected.csv \
                                data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318_corrected.csv \
                                data
```

> RILs have their genotype ID instead of name because some RILs had multiple IDs. That's why I generated a list relating names to IDs. Also, I converted `--` to `NN` to represent missing data for compatibility with TASSEL 5 software.

Since there are genotypic data for RILs than used in this project, I need to filter out some RILs. The file `data/usda_RILs.txt` contains all the RILs used, and the `scripts/remove_extra_geno-data.R` script removes extra genotypic data.

```bash
module load R/3.6.0

Rscript scripts/remove_extra_geno-data.R data/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt \
                                         data/usda_22kSNPs_parents.hmp.txt \
                                         data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318_corrected.hmp.txt \
                                         data/usda_22kSNPs_rils.hmp.txt \
                                         data/usda_RILs.txt
```

I will use the software [Tassel 5](https://www.maizegenetics.net/tassel) for some basic QC, because this software can generate summaries (including allele frequency) pretty quickly. Tassel requires that the genotypic data is sorted by position and chromosome number. Thus, I used `SortGenotypeFilePlugin` from Tassel to quickly sort the genotypic data and create another hapmap file called `data/usda_22kSNPs_rils.sorted.hmp.txt`. Then I ran Tassel again to transform the sorted data into diploid hapmap format. The file `data/usda_22kSNPs_rils.sorted.diploid.hmp.txt` file will be used for the rest of QC analysis.

```bash
# sort ril data
run_pipeline.pl -Xmx2g -SortGenotypeFilePlugin \
                       -inputFile data/usda_22kSNPs_rils.hmp.txt \
                       -outputFile data/usda_22kSNPs_rils.sorted.hmp.txt \
                       -fileType Hapmap
run_pipeline.pl -Xmx2g -importGuess data/usda_22kSNPs_rils.sorted.hmp.txt \
                       -export data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                       -exportType HapmapDiploid

# sort parents
run_pipeline.pl -Xmx2g -SortGenotypeFilePlugin \
                       -inputFile data/usda_22kSNPs_parents.hmp.txt \
                       -outputFile data/usda_22kSNPs_parents.sorted.hmp.txt \
                       -fileType Hapmap
run_pipeline.pl -Xmx2g -importGuess data/usda_22kSNPs_parents.sorted.hmp.txt \
                       -export data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                       -exportType HapmapDiploid
```


### Data QC

The QC of parental and RIL data will be done on a family basis. Thus, I have to find out which RILs belong to the same biparental cross. I wrote the `scripts/create_list_biparental-crosses.py` script to generate a table with a biparental cross in a column and all the RIL IDs in another, and ran this in the command line to create the table:

```bash
python3 scripts/create_list_biparental-crosses.py data/id_table_Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318_corrected.txt \
                                                 data/usda_biparental-crosses.txt
```

I wrote `scripts/markers_summary.R` to summarize the number of missing and heterozygous markers for the parents and RILs and save them at the `analysis/qc` folder.

```bash
module load R/3.6.0

# create directory to store qc results
mkdir -p analysis/qc
# run qc
Rscript scripts/markers_summary.R data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_parents.txt \
                                  data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_rils.txt \
                                  data/usda_biparental-crosses.txt
```

The PHG35 parent has much more heterozygotes than expected for a fully inbred line. This might very likely be due to polen contamination in the seeds used to do the genotyping. I can use the marker genotypes in the RIL data from all PHG35 progeny and figure out the PHG35 haplotype based on the genotype of the other parent used to develop that progeny. To do that, I wrote `scripts/reconstruct_PHG35_from_RIL_data.R`.

```bash
module load R/3.6.0

Rscript scripts/reconstruct_PHG35_from_RIL_data.R data/usda_biparental-crosses.txt \
                                                  data/usda_22kSNPs_parents.sorted.diploid.hmp.txt \
                                                  data/usda_22kSNPs_rils.sorted.diploid.hmp.txt
```

To make sure the reconstruction worked, I ran again `scripts/markers_summary.R`. Plots at `analysis/qc` folder show that, despite higher mising data relatively to other parents, the number of heterozygotes is now zero.

```bash
Rscript scripts/markers_summary.R data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt \
                                  analysis/qc/summary_markers_parents_PHG35-reconstructed.txt \
                                  data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                                  analysis/qc/summary_markers_rils.txt \
                                  data/usda_biparental-crosses.txt
```

By looking at the boxplot of missing data of RILs, what stands out is that cross `B73xPHG47` had the highest variation and amount of missing data than other populations, and few populations showed very few RILs with more extreme missing data. Although the maximum amount of missing data for a RIL was about 9%, the median amount of missing data per cross was below 2.5%. For heterozygous markers, they were also mostly between 0 and 3% (which would be the expected for an F6 population), however population `B73xPHG47` also displayed a lot of variation and higher values of heterozygous markers. Very few RILs displayed high levels of heterozygosity, but those values were much higher than the expected (between 10 to 30%). Finally, I cannot see any parent that contributes with more heterozygosity or missing data in the RILs.

With RILs from a biparental cross, we have the expectation that each segregating locus will have only two alleles with frequency ~50% each. Large deviations from this expectation may indicate some errors during genotyping and the locus should be removed from analysis.

```bash
# save current field delimiter
curr_IFS=${IFS}
{
  # skip header of "usda_biparental-crosses.txt" file
  read
  # set delimeter to tab
  IFS="\t"
  # read file line by line
  while read -r line; do
    # split line and assign each column to a variable
    # also change "*" by "x" in the cross name
    cross=$(echo ${line} | cut -f1 | tr "*" "x")
    ril_list=$(echo ${line} | cut -f2)
    # check if directory exists; if it doesn't, create one to store results
    [[ -d analysis/qc/${cross} ]] || mkdir -p analysis/qc/${cross}
    # run tassel (added "\" at the end of the line just to improve readability)
    run_pipeline.pl -Xmx6g -importGuess data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                    -FilterTaxaBuilderPlugin -taxaList ${ril_list} -endPlugin \
                    -GenotypeSummaryPlugin -endPlugin \
                    -export analysis/qc/${cross}/${cross}_OverallSummary,analysis/qc/${cross}/${cross}_AlleleSummary,analysis/qc/${cross}/${cross}_SiteSummary,analysis/qc/${cross}/${cross}_TaxaSummary
  done
} < "data/usda_biparental-crosses.txt" > "analysis/qc/tassel_log.txt"
# set delimiter back to what it was
IFS=${curr_IFS}

# remove empty folders (RILs without genotypic data)
# note that although i'm using *, the folder won't be deleted if it's not empty
rmdir analysis/qc/*
```

In addition, I will use the sorted diploid hapmap files from parents and RILs to find markers with segregation distortion (i.e. allele frequencies are too different from what is expected in RILs) and estimate recombination frequency for each biparental population. For this purpose, I wrote the `scripts/recomb-freq_biparental-crosses.R` script.

```bash
# find segregation distortion
Rscript scripts/recomb-freq_biparental-crosses.R data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt \
                                                 data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                                                 data/usda_biparental-crosses.txt \
                                                 analysis/qc \
                                                 data/biparental-crosses

# create file with problematic snps from all crosses
cat analysis/qc/*/markers-seg* analysis/qc/*/markers-missing* | sort | uniq > analysis/qc/rqtl_snps-to-remove_all-crosses.txt

# # number of markers to be removed
# wc -l analysis/qc_geno/rqtl_snps-to-remove_all-crosses.txt
# # 1174

# filter parental file
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt \
                -excludeSiteNamesInFile analysis/qc/rqtl_snps-to-remove_all-crosses.txt \
                -export data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt \
                -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_rils.sorted.diploid.hmp.txt \
                -excludeSiteNamesInFile analysis/qc/rqtl_snps-to-remove_all-crosses.txt \
                -export data/usda_22kSNPs_rils.sorted.diploid.filtered.hmp.txt \
                -exportType HapmapDiploid
```


### Convert SNP coordinates from B73v2 to B73v4

The SNP chip data was designed from the B73 v2 reference assembly, so I need to convert the marker coordinates to a more up-to-date version: B73 v4.

```bash
# make directory to store refgen sequence, and go to that directory
mkdir -p data/refgen/B73v2
cd data/refgen/B73v2

# download v2 sequences for all chromosomes
for chr in {1..10}; do
  wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/Zea_mays.AGPv2.dna.chromosome.$chr.fa
done
# download README
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-7/fasta/zea_mays/dna/README

# return to project's home directory
cd ~/projects/marker-effects_networks
# extract probes sequences
python3 scripts/extract_probe_seqs.py data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt \
                                     data/refgen/B73v2/ \
                                     100 \
                                     data/probes-100bp_22kSNPchip_B73v2.fa
```

Once I have a fasta file with all probe sequences, I can align them to the refgen v4 assembly using `bowtie`. But first, I need to download the v4 assembly and build an index (by running `scripts/bowtie_index_refgenv4.sh`).

```bash
# make directory to store refgen sequence, and go to that directory
mkdir -p data/refgen/B73v4
cd data/refgen/B73v4

# download sequences for all chromosomes
for chr in {1..10}; do
  wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr$chr.fna.gz
done
# decompress them
gunzip *.gz
# change extension name .fna to .fa for bowtie compatibility
for file in *.fna; do
    mv -- "$file" "${file%.fna}.fa"
done

# return to project's home directory
cd ~/projects/marker-effects_networks

# build index for refgen v4
sbatch scripts/bowtie_index_refgenv4.sh
```

Once the index is built, mapping reads with bowtie is pretty straightforward:

```bash
# load bowtie
module load bowtie/1.1.2

# map probes to refgen v4
# bowtie [options] [index_basename] [unpaired_reads] [output_file]
bowtie -f -v 0 -m 1 --un data/probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa -S \
       data/refgen/B73v4/index \
       data/probes-100bp_22kSNPchip_B73v2.fa \
       data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam 2> data/probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt
```

> Bowtie options:
> * `-f`: input is a fasta file
> * `-v 0`: map end-to-end with 0 mismatches allowed
> * `-m 1`: don't report alignments if read map in more than 1 place
> * `--un data/probes-100bp_22kSNPchip_not-aligned-to-B73v4.fa`: unmapped reads are written to this file
> * `-S`: output is SAM
> * `data/refgen/B73v4/index`: basename of the ref gen v4 index
> * `data/probes-100bp_22kSNPchip_B73v2.fa`: reads to map
> * `data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam`: output name
> * `2> data/probes-100bp_22kSNPchip_aligned-to-B73v4.stats.txt`: this will print alignment statistics from std err to this file

Here are the stats of the alignment:

| Reads                   | Stats  |
| ----------------------- | ------ |
| Processed               | 18,965 |
| Uniquely mapped         | 18,694 |
| Multimapped (supressed) | 64     |
| Unmapped                | 207    |

Now, I need to correct the positions of the SNPs in the hapmap files based on new coordinates from SAM file. I wrote a python script called `scripts/convert_hmp_v2-to-v4.py` to do that, and another file with SNP coordinates in both v2 and v4. This script also removes reads that uniquely mapped with no mismatch to a different chromosome in v4, since it's very unlikely that there would be such major differences between v2 and v4 assemblies (i.e., these are probably mismapped).

```bash
python3 scripts/convert_hmp_v2-to-v4.py data/probes-100bp_22kSNPchip_aligned-to-B73v4.sam \
                                       data/SNP_positions_v2-to-v4_probes-100bp.txt \
                                       data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.hmp.txt,data/usda_22kSNPs_rils.sorted.diploid.filtered.hmp.txt
# 22 reads discarded due to mapping into different chromosome in v4
```

Inspecting the file, I realized that some markers that didn't match between v2 and v4 and that these mismatches were due to mapping to complementary strand. Thus, I corrected the strand of those markers in the hapmap files by running `scripts/correct_SNP_strands.R`. Finally, I corrected the alleles' column of these hapmap files using Tassel.

```bash
# correct markers not phased
Rscript scripts/correct_SNP_strands.R data/SNP_positions_v2-to-v4_probes-100bp.txt \
                                      data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                                      data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt

# correct alleles column
run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                        -export data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                        -exportType HapmapDiploid

run_pipeline.pl -Xmx10g -importGuess data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt \
                        -export data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt \
                        -exportType HapmapDiploid
```


### Correct miscalls with sliding window approach

It's also important to minimize the number of miscalls in the SNP chip data for RILs because they can affect downstream analysis. We will do that by using the sliding window approach described by [Huang et al. (2009)](https://genome.cshlp.org/content/19/6/1068.abstract) for each RIL family separately.

```bash
# divide hmp by cross
# parents
Rscript scripts/divide_hmp_by_cross.R data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      data/hapmap_by_cross \
                                      --parents
# rils
Rscript scripts/divide_hmp_by_cross.R data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt \
                                      data/usda_biparental-crosses.txt \
                                      data/hapmap_by_cross \
                                      --rils

# run sliding window approach
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in ${crosses}; do
  Rscript scripts/sliding_window_approach.R ${cross} \
                                            data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.hmp.txt \
                                            data/hapmap_by_cross/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.${cross}.hmp.txt \
                                            --window_size=15 --window_step=1 --min_snps_per_window=5
done

# merge files
crosses=$(ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq)

for cross in ${crosses}; do
  echo $cross
  if [[ ${cross} == "B73xLH82" ]]; then
    HMP1=data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.B73xLH82.sliding-window.hmp.txt
    HMP2=data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.B73xPH207.sliding-window.hmp.txt
  elif [[ ${cross} == "B73xPH207" ]]; then
    continue
  else
    HMP1=data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.sliding-window.hmp.txt
    HMP2=data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.${cross}.sliding-window.hmp.txt
  fi
  echo ${HMP1}
  echo ${HMP2}
  echo ${OUT}
  OUT=data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.sliding-window.hmp.txt
  run_pipeline.pl -fork1 -h ${HMP1} \
                  -fork2 -h ${HMP2} \
                  -combine3 -input1 -input2 \
                  -mergeGenotypeTables \
                  -export ${OUT} \
                  -exportType HapmapDiploid \
                  -runfork1 -runfork2 > /dev/null
done
```


### Generate hybrid genotypes

Use information available on `data/NIFA_CompleteDataset.csv` to get pedigree of hybrids, and then create hybrid genotypes based on genotypes of parental (RIL) data generated above.

```bash
# create hybrid genotypes
HYBINFO=data/NIFA_CompleteDataset.csv
RILGENO=data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.sliding-window.hmp.txt
OUT=data/usda_hybrids_SNP-chip.hmp.txt
HET=missing
# create hybrids
sbatch --export=HYBINFO=${HYBINFO},RILGENO=${RILGENO},OUT=${OUT},HET=${HET} scripts/create_hybrid_genotypes.sh

# remove monomorphic markers (maf < 0.05)
run_pipeline.pl -Xmx40g \
                -importGuess data/usda_hybrids_SNP-chip.hmp.txt \
                -FilterSiteBuilderPlugin \
                -siteMinAlleleFreq 0.05 \
                -endPlugin \
                -export data/usda_hybrids_SNP-chip.maf-filter,data/usda_hybrids_SNP-chip.maf-filter2 -exportType HapmapDiploid
# remove unnecessary file
rm data/usda_hybrids_SNP-chip.maf-filter2.json.gz
```

Check how much having missing data the genotypic dataset currently has, and also plot marker distribution along the chromosomes.

```bash
# qc
sbatch --export=HMP=data/usda_hybrids_SNP-chip.maf-filter.hmp.txt,FOLDER=analysis/qc/hybrid_hmp scripts/qc_hybrid_hmp.sh
```

| missing filter | marker # |
| -------------- | -------- |
| no filter      | 14,676   |
| 0.2            | 13,757   |
| 0.3            | 14,210   |
| 0.4            | 14,379   |
| 0.5            | 14,466   |


### Prune markers by LD

We need to prune the full dataset (`data/usda_hybrids_SNP-chip.maf-filter.hmp.txt`) based on LD to decrease the number of perfectly correlated markers. I tested three different window sizes and three different missing data thresholds to calculate LD to see how many markers would be removed.

```bash
# hapmap file to filter
HMP=data/usda_hybrids_SNP-chip.maf-filter.hmp.txt
# prunning parameters
VARCOUNT=1
R2=0.9

for genomiss in 0.1 0.25 0.5; do
  for winsize in 10kb 100kb 1000kb; do
    # folder to save intermediate files
    FOLDER=analysis/prune_markers_ld/winsize_${winsize}
    mkdir -p ${FOLDER}
    # output name
    OUT=data/usda_hybrids_SNP-chip.maf-filter.pruned-${winsize}.geno-miss-${genomiss}.hmp.txt
    # set intermediate filenames
    PLK=$(basename ${HMP} .hmp.txt)
    PRUNNED=$(echo ${PLK}.marker-ids)
    # transform hmp to plk
    run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                    -export ${FOLDER}/${PLK} \
                    -exportType Plink
    # prune plink file by ld
    plink --file ${FOLDER}/${PLK}.plk \
          --indep-pairwise ${winsize} ${VARCOUNT} ${R2} \
          --geno ${genomiss} \
          --out ${FOLDER}/${PRUNNED} \
          --allow-extra-chr \
          --make-founders
    # filter hapmap
    run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                    -includeSiteNamesInFile ${FOLDER}/${PRUNNED}.prune.in \
                    -export ${OUT} \
                    -exportType HapmapDiploid
  done
done

# count how many markers were left after prunning
wc -l data/usda_hybrids_SNP-chip.maf-filter.pruned-*.geno-miss-*.hmp.txt
```

Number of markers remaining after pruning the genotypic dataset with different window sizes:

| missing data filter | window size | markers remaining |
| ------------------- | ----------- | ----------------- |
| 0.1                 | 10 kb       | 9,482             |
| 0.1                 | 100 kb      | 8,377             |
| 0.1                 | 1000 kb     | 5,670             |
| 0.25                | 10 kb       | 11,777            |
| 0.25                | 100 kb      | 10,334            |
| 0.25                | 1000 kb     | 6,834             |
| 0.5                 | 10 kb       | 12,120            |
| 0.5                 | 100 kb      | 10,615            |
| 0.5                 | 1000 kb     | 6,978             |

We decided to use `data/usda_hybrids_SNP-chip.maf-filter.pruned-100kb.geno-miss-0.25.hmp.txt` to balance missing data, window size and number of markers to build networks.

> This is the File S5 I mention in the manuscript.



## Phenotypic data

We also need phenotypic data to estimate the effects markers. For this purpose, I will use grain yield data from my genomic prediction project where ~400 hybrids were evaluated in 9 environments from file `data/NIFA_CompleteDataset.csv`. First we run ANOVA to see if these population has enough variability for the trait across these environments, and then we obtain BLUEs for each hybrid in each environment.

```bash
module load R/3.6.0

# run anova
Rscript scripts/anova_pheno_hybrids.R \
        data/NIFA_CompleteDataset.csv \
        BEC-BL19,BEC-BL20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
        UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P \
        YLD \
        data

# get BLUes
Rscript scripts/get_BLUEs_from_empirical-data.R \
        data/NIFA_CompleteDataset.csv \
        YLD \
        data/1stStage_BLUEs.YLD-per-env.txt \
        analysis/pheno_BLUES_qc \
        --envs=BEC-BL19,BEC-BL20,COR19,COR20,MIN19,MIN20,SYN19,SYN20,URB19 \
        --checks=UIUC-CHECK1,UIUC-CHECK2,UIUC-CHECK3,UIUC-CHECK4,6049V2P

# count number of hybrids per environment
grep -v NA data/1stStage_BLUEs.YLD-per-env.txt | cut -f 2 | sed 1d | sort | uniq -c
```

Total number of hybrids evaluated per environment:

| env      | hybrids |
| -------- | ------- |
| BEC-BL19 | 367     |
| BEC-BL20 | 398     |
| COR19    | 370     |
| COR20    | 366     |
| MIN19    | 356     |
| MIN20    | 383     |
| SYN19    | 361     |
| SYN20    | 397     |
| URB19    | 372     |

> Originally it was 10 environments, but after preliminary analysis, we were not able to get good marker effects estimation for one environment and decided to drop it.



## Estimate marker effects

We will estimate marker effects for the trait using two different models to see how much the network changes based on the type of model selected. The rationale for using more than one model is that different models may be better for a certain genetic architecture, which is unknown from empirical data. In addition, we chose to use marker dataset pruned in 100kb windows to balance the number of markers vs computation time.

```bash
# estimate effects
MARKERS=data/usda_hybrids_SNP-chip.maf-filter.pruned-100kb.geno-miss-0.25.hmp.txt
for trait in YLD; do
  BLUES=data/1stStage_BLUEs.${trait}-per-env.txt
  OUTFOLDER=analysis/marker_effects/${trait}
  for marker_eff_model in rrblup gwas; do
    sbatch --export=MARKERS=${MARKERS},BLUES=${BLUES},OUTFOLDER=${OUTFOLDER},MEFFMODEL=${marker_eff_model} scripts/estimate_marker_effects.sh
  done
done

# plot results for qc
module load R/3.6.0
Rscript scripts/plot_marker_effects.R analysis/marker_effects \
                                      analysis/marker_effects/qc \
                                      --traits=YLD \
                                      --models=rrblup,gwas
```

Correlation between real phenotypes and predicted GEBVs:

| env      | rrBLUP | GWAS |
| -------- | ------ | ---- |
| BEC-BL19 | 0.86   | 0.86 |
| BEC-BL20 | 0.72   | 0.72 |
| COR19    | 0.77   | 0.76 |
| COR20    | 0.82   | 0.82 |
| MIN19    | 0.78   | 0.79 |
| MIN20    | 0.72   | 0.7  |
| SYN19    | 0.82   | 0.83 |
| SYN20    | 0.75   | 0.74 |
| URB19    | 0.77   | 0.77 |



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
    # create output folder
    OUTFOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}
    mkdir -p ${OUTFOLDER}
    # submit job
    sbatch --export=TRAIT=${TRAIT},MEFF_FILE=${MEFF_FILE},MEFF_MODEL=${meff_model},NORM_METHOD=${norm_method},OUTFOLDER=${OUTFOLDER} scripts/plot_cv_meffs.sh
  done
done
```


### Step 2. Pick a soft threshold

In order to build networks with WGCNA, I need to specify a power to model the scale-free topology of the network. The `scripts/pick_soft_threshold.R` plots the scale-free topology fit index for different powers so I can decide which one is more adequate for my data.

```bash
# marker file
MARKERS=data/usda_hybrids_SNP-chip.maf-filter.pruned-100kb.geno-miss-0.25.hmp.txt
# trait to build network
TRAIT=YLD

# test different cv thresholds depending on normalization method

# marker effects model
for meff_model in rrblup gwas; do
  # marker effects file
  MEFF_FILE=analysis/marker_effects/${TRAIT}/marker_effects.${meff_model}.txt
  # type of normalization method to use
  for norm_method in minmax; do
    # adjust cv threshold to type of normalization
    for cv_threshold in 0.05 0.075 0.1 0.125 0.15 0.175 0.2 0.25; do
      # create output folder
      OUTFOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/cv_thresholds
      mkdir -p ${OUTFOLDER}
      # submit job
      sbatch --export=MARKERS=${MARKERS},MEFF_FILE=${MEFF_FILE},OUTFOLDER=${OUTFOLDER},NORM_METHOD=${norm_method},CV_THRESHOLD=${cv_threshold} scripts/pick_soft_threshold.sh
    done
  done
done

# marker effects model
for meff_model in rrblup gwas; do
  # marker effects file
  MEFF_FILE=analysis/marker_effects/${TRAIT}/marker_effects.${meff_model}.txt
  # type of normalization method to use
  for norm_method in zscore none; do
    # adjust cv threshold to type of normalization
    for cv_threshold in 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 1; do
      # create output folder
      OUTFOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/cv_thresholds
      mkdir -p ${OUTFOLDER}
      # submit job
      sbatch --export=MARKERS=${MARKERS},MEFF_FILE=${MEFF_FILE},OUTFOLDER=${OUTFOLDER},NORM_METHOD=${norm_method},CV_THRESHOLD=${cv_threshold} scripts/pick_soft_threshold.sh
    done
  done
done
```

The CV threshold that maximized scale-free topology model fit and mantained about the same number of markers was 0.05 for minmax normalization and 0.5 for Z-score normalization. Also, Z-score and no normalization yield the same results, so I'll just use `minmax`- and `zscore`-normalized data for downstream analysis.

```bash
# marker effects model
for meff_model in rrblup gwas; do
  # marker effects file
  MEFF_FILE=analysis/marker_effects/${TRAIT}/marker_effects.${meff_model}.txt
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # adjust cv threshold to type of normalization
    [[ ${norm_method} == 'minmax' ]] && CV_THRESHOLD=0.05
    [[ ${norm_method} == 'zscore' ]] && CV_THRESHOLD=0.5
    # create output folder
    OUTFOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}
    mkdir -p ${OUTFOLDER}
    # submit job
    sbatch --export=MARKERS=${MARKERS},MEFF_FILE=${MEFF_FILE},OUTFOLDER=${OUTFOLDER},NORM_METHOD=${norm_method},CV_THRESHOLD=${CV_THRESHOLD} scripts/pick_soft_threshold.sh
  done
done
```



### Step 3. Build the network

I will use power of `24` as the soft threshold for WGCNA as it seems to have a good fit for all normalization and model effect scenarios. It's worth noting that this is a much higher power than what's tipically used in gene co-expression networks, suggesting that building marker effects networks will require tuning additional parameters downstream the pipeline. The `scripts/build_meff_network.R` builds the network for a given soft threshold.

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
    # power for which the scale-free topology fit index curve
    SFT=24
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
    # power for which the scale-free topology fit index curve
    SFT=24
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

> Modules from `gwas` networks with `100` minimum number of markers couldn't be defined.



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
        # power for which the scale-free topology fit index curve
        SFT=24
        # number of most important markers to print for each module
        NHUBS=10
        # submit job
        sbatch --export=RDATA=${RDATA},OUTFOLDER=${FOLDER},SFT=${SFT},NHUBS=${NHUBS} scripts/qc_network_modules.sh
      done
    done
  done
done

# summarize network quality results
Rscript scripts/summarize_qc_networks.R analysis/networks/YLD
```

> In general, kDiff plots from rrBLUP effects networks seem to be have better quality modules (i.e. markers with more connections within their own module than with other modules) than those from GWAS effects networks. It's hard to tell if turning off the `pamStage` when assigning modules had any meaningful impact on kDiff, cluster coefficient and TOM plots than when this option is on.

In addition, I wrote `scripts/compare_two_networks.R` to check which modules in one network corresponds to another module in a different network. This will help us understand how different parameters chosen during network construction affects the clustering of markers, and also help identify correponding modules when correlating networks with traits. To do that, I calculate the correlation between the module memberships (the correlation of the module eigengene and the marker effect profile) of two different networks.

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




### Step 6. Calculate LD for markers within each module

An important thing to check is how many markers in a module are in LD with each other. If all of them are in perfect LD, then the marker effect networks would be just picking up LD without any other biological reason for grouping those markers. The `scripts/ld_markers_modules.sh` calculates LD all markers within a module (even if they are in different chromosomes).

```bash
# trait to build network
TRAIT=YLD
# genotypic data of markers
MARKERS=data/usda_hybrids_SNP-chip.maf-filter.pruned-100kb.geno-miss-0.25.hmp.txt

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

# summarize network quality results
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}/modules_ld
        # submit job
        sbatch --export=FOLDER=${FOLDER},MEFF_MODEL=${meff_model},NORM_METHOD=${norm_method},MINSIZE=${minsize},PAM=${pam} scripts/summarize_ld_per_network.sh
      done
    done
  done
done

# save header before merging summary of all networks
head -n 1 analysis/networks/${TRAIT}/meff_rrblup/norm_minmax/min_mod_size_25/pamStage_on/modules_ld/summary_ld.txt > analysis/networks/${TRAIT}/summary_ld_per_network.txt
for meff_model in rrblup gwas; do
  # type of normalization method to use
  for norm_method in minmax zscore; do
    # minimum number of markers per module
    for minsize in 25 50 100; do
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}/modules_ld
        # combine results from different networks and plot summary
        sed 1d ${FOLDER}/summary_ld.txt >> analysis/networks/${TRAIT}/summary_ld_per_network.txt
      done
    done
  done
done
# compress file
gzip analysis/networks/YLD/summary_ld_per_network.txt
# plot summary
LDFILE=analysis/networks/YLD/summary_ld_per_network.txt.gz
OUTFOLDER=analysis/networks/YLD/qc_networks
sbatch --export=LDFILE=${LDFILE},OUTFOLDER=${OUTFOLDER} scripts/plot_ld_summary_per_network.sh
```



## Relate modules to environmental covariables

Since we are interested in understanding phenotypic plasticity across environments, we also want to see if the network modules correlate with environmental covariables. To do that, we first extract these covariables running `scripts/get_env_types.R`. I created `data/env_sites_coord.csv` with environmental coordinates and planting/harvest dates.

> The information for file `data/env_sites_coord.csv` is available as Table S2 in the manuscript. Running `scripts/get_env_types.R` will create the `data/env_covariables/env_covariables_means_per_intervals.txt` which is the File S6, so you can skip this step if you want.

```bash
module load R/3.6.0

# get env covariables
Rscript scripts/get_env_types.R data/env_sites_coord.csv \
                                data/env_covariables \
                                --country=USA \
                                --interval-window=3
```

> Note: had to run this script locally due to some sort of incompatibility with MSI server. I transfered the files to MSI afterwards using FileZilla.

Due to high number of time points and environmental covariables, we'll also do a PCA to reduce dimensionality.

```bash
Rscript scripts/pca_env_idx.R \
        data/env_covariables/env_covariables_means_per_intervals.txt \
        data/env_covariables \
        --per-intervals
```

In addition, we were also interested in correlating the module eigenvalues with environmental indices given by Finlay-Wilkson regression.

```bash
module load R/3.6.0

# get env indices from FW
Rscript scripts/get_FW_env-idx.R data/1stStage_BLUEs.YLD-per-env.txt \
                                 analysis/FW/YLD
```

Then, we correlate these covariables with the eigenvalues of each module with `scripts/relate_modules_to_env_idx.sh` and summarize their correlation with `scripts/summarize_mod-env-idx_relationship.R`.

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
      # define modules with and without the PAM stage
      for pam in on off; do
        # folder with results
        FOLDER=analysis/networks/${TRAIT}/meff_${meff_model}/norm_${norm_method}/min_mod_size_${minsize}/pamStage_${pam}
        # file with saved R variables from step 2
        RDATA=${FOLDER}/define_network_modules.RData
        # env indices to correlate to modules
        for idx_type in means intervals pca fw; do
          # means over entire season
          [[ ${idx_type} == "means" ]] && EC_FILE=data/env_covariables/env_covariables_means.txt
          # idx per 3-day intervals
          [[ ${idx_type} == "intervals" ]] && EC_FILE=data/env_covariables/env_covariables_means_per_intervals.txt
          # principal components
          [[ ${idx_type} == "pca" ]] && EC_FILE=data/env_covariables/pca_env_idx.txt
          # FW indices
          [[ ${idx_type} == "fw" ]] && EC_FILE=analysis/FW/YLD/FW_env_idx.txt
          # define output folder
          OUTFOLDER=${FOLDER}/mod_env-idx_${idx_type}
          # submit job
          sbatch --export=RDATA=${RDATA},EC_FILE=${EC_FILE},OUTFOLDER=${OUTFOLDER},IDX_TYPE=${idx_type} scripts/relate_modules_to_env_idx.sh
        done
      done
    done
  done
done

# summarize correlations
Rscript scripts/summarize_mod-env-idx_relationship.R analysis/networks/YLD
```

Finally, plot the relationship between environmental indices and module eigeinvalues, the markers from different networks correlated with same environmental index, and the top covariables that contributes the most to each principal component of the PCA analysis.

```bash
module load R/3.6.0

for idx_type in means intervals pca fw; do

  # means over entire season
  [[ ${idx_type} == "means" ]] && EC_FILE=data/env_covariables/env_covariables_means.txt
  # idx per 3-day intervals
  [[ ${idx_type} == "intervals" ]] && EC_FILE=data/env_covariables/env_covariables_means_per_intervals.txt
  # principal components
  [[ ${idx_type} == "pca" ]] && EC_FILE=data/env_covariables/pca_env_idx.txt
  # FW indices
  [[ ${idx_type} == "fw" ]] && EC_FILE=analysis/FW/YLD/FW_env_idx.txt

  # define p-values
  if [[ ${idx_type} == "means" ]]; then
    PVAL=0.1
  elif [[ ${idx_type} == "intervals" ]]; then
    PVAL=0.001
  else
    PVAL=0.1
  fi

  # plot env_idx-module relationships
  Rscript scripts/plot_mod-env-idx.R \
          analysis/networks/YLD \
          analysis/networks/YLD/module-env-idx_per_network.${idx_type}.txt \
          ${EC_FILE} \
          analysis/networks/YLD/plots_mod-env-idx \
          --p-value=${PVAL}

  # check markers in modules that are correlated to the same env idx
  Rscript scripts/markers_correlated_same_env_idx.R \
          analysis/networks/YLD \
          analysis/networks/YLD/module-env-idx_per_network.${idx_type}.txt \
          data/usda_hybrids_SNP-chip.maf-filter.pruned-100kb.geno-miss-0.25.hmp.txt \
          analysis/networks/YLD/overlap_markers \
          --p-value=${PVAL}

done

# plot top covariables contributing to a PC
Rscript scripts/plot_pc_contributions.R \
        data/env_covariables/pca_contributions.txt \
        analysis/networks/YLD/pc_contributions \
        --prop-covariables=0.1
```
