library(data.table)
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

usage <- function() {
  cat("
description: run PCA on genotypic data and generate diagnostic plots to help
             decide how many principal components should be used in GWAS.

usage: Rscript pca_prior_gwas.R [geno_file] [output_folder] [...]

positional arguments:
  geno_file           file with numeric genotypic data
  gwas_model          type of gwas that will be run: 'additive' (default) or 'dominant'
  output_folder       folder to save files

optional argument:
  --help              show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 3) stop(usage(), "missing positional argument(s)")

# get positional arguments
geno_file <- args[1]
gwas_model <- args[2]
output_folder <- args[3]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)


#### run PCA ----

# load data without header for compatibility with GAPIT
geno <- fread(geno_file, head = FALSE, data.table = FALSE)

# adjust numericalization parameters according to model
if (gwas_model == "additive") {
  marker_effect = "Add"
  marker_impute = "Major"
}

if (gwas_model == "dominant") {
  marker_effect = "Dom"
  marker_impute = "Middle"
}

# transform to numeric
geno_num <- GAPIT.HapMap(G = geno, SNP.effect = marker_effect, SNP.impute = marker_impute)
rm(geno)
# format numeric marker data to be compatible to GAPIT
markers <- data.frame(geno_num$GD)
colnames(markers) <- geno_num$GI$SNP
markers <- cbind(taxa = geno_num$GT, markers)
# generate genetic map with marker names, chr and position
map <- geno_num$GI
rm(geno_num)

# working directory
init_folder <- getwd()
# change directories for plotting
setwd(paste0(init_folder, "/", output_folder))
# plot pca files
pca <- GAPIT.PCA(markers[, -1], taxa = markers[, 1], PCA.total = 10)
# return to initial directory
setwd(init_folder)



#### debug ----

# geno_file <- "data/g2f_hybrids_geno_filtered.pruned.hmp.txt"
# output_folder <- "analysis/pca_prior_gwas"
# gwas_model <- "additive"
