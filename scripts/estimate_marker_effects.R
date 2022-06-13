library(data.table)
library(rrBLUP)
library(bWGR)
library(tidyr)
library(tibble)
library(ggplot2)
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

usage <- function() {
  cat("
description: estimate effects of markers for a trait across multiple environments.

usage: Rscript estimate_marker_effects.R [markers_file] [blues_file] [output_folder] [...]

positional arguments:
  markers_file                hapmap file containing marker data
  blues_file                  file with BLUEs in long format (3 columns: genotype, environment, value)
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --marker-eff-model          the type of model to calculate effects of each marker across environments.
                              Available models are 'rrblup' (univariate rrBLUP; default), 'bayescpi' (univariate
                              BayesCpi), 'mrr' (multivariate rigde-regression), 'gwas' (univariate Q+K GWAS)
  --impute-effect=VALUE       the marker effect ('Add' or 'Dom') when imputing missing data at the hapmap
                              numericalization step (default: 'Add')
  --impute-type=VALUE         the marker type ('Major', 'Middle' or 'Minor') when imputing missing data at
                              the hapmap numericalization step (default: 'Middle')
  --no-missing-genotypes      exclude genotypes with missing data in any environment


"
  )
}

getArgValue <- function(arg) {

  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])

}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# get positional arguments
markers_file <- args[1]
blues_file <- args[2]
output_folder <- args[3]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
marker_eff_model <- "rrblup"
impute_effect <- "Add"
impute_type <- "Middle"
no_missing_genotypes <- FALSE

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--marker-eff-model", "--impute-effect", "--impute-type", "--no-missing-genotypes")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")

  # change default based on the argument provided
  for (argument in opt_args_allowed) {
    if (any(grepl(argument, opt_args_requested))) {
      arg_name <- gsub("-", "_", gsub("--", "", argument))
      arg_value <- getArgValue(opt_args[grep(argument, opt_args)])
      assign(arg_name, arg_value)
    }
  }

}

# make sure optional arguments are valid
if (!marker_eff_model %in% c("rrblup", "bayescpi", "mrr", "gwas") ) {
  stop("Optional argument '--marker-eff-model' should be 'rrblup', 'bayescpi', 'mrr', or 'gwas'")
}

if (!impute_effect %in% c("Add", "Dom", "Minor")) {
  stop("Optional argument '--impute-effect' should be 'Add' or 'Dom'")
}

if (!impute_type %in% c("Major", "Middle", "Minor")) {
  stop("Optional argument '--impute-type' should be 'Major', 'Middle' or 'Minor'")
}




#### load data ----

# "header = FALSE" for compatibility with GAPIT
hapmap <- fread(markers_file, header = FALSE, data.table = FALSE)
# numericalize hapmap
hapmap <- GAPIT.HapMap(G = hapmap, SNP.effect = impute_effect, SNP.impute = impute_type)
#   hm$GT: vector with line names
#   hm$GD: matrix where each row is a line and each column is a marker (the numbers of the cells are the numeric genotypes)
#   hm$GI: data frame with marker information (name, chromosome and position)

# length(hapmap$GI$SNP) == NCOL(markers)
# length(hapmap$GT) == NROW(markers)
markers <- data.frame(hapmap$GD - 1)
colnames(markers) <- hapmap$GI$SNP
rownames(markers) <- hapmap$GT
if (marker_eff_model == "gwas") map <- hapmap$GI
rm(hapmap)

# remove duplicated genotypes to avoid singularities on relationship matrix
markers <- markers[!duplicated(markers), ]

# load and format phenotypic data
means <- fread(blues_file, header = TRUE, data.table = FALSE)

# transform data to wide format
colnames(means) <- c("genotype", "environment", "pheno_value")
means <- pivot_wider(means, names_from = environment, values_from = pheno_value)
means <- data.frame(means)
means <- column_to_rownames(means, var = "genotype")
if (no_missing_genotypes) means <- data.frame(na.omit(means))

# using only phenotypes with snp information
means <- means[intersect(rownames(markers), rownames(means)), ]
markers <- markers[intersect(rownames(markers), rownames(means)), ]
if (!all(rownames(markers) == rownames(means))) stop("IDs don't match between geno and pheno data")



#### estimate marker effects ----

if (marker_eff_model == "rrblup") {

  cat("estimating marker effects by univariate rrBLUP\n")

  # create empty df to store results
  marker_effects <- data.frame(matrix(ncol = 0, nrow = ncol(markers)))
  rownames(marker_effects) <- colnames(markers)
  gebvs <- data.frame(matrix(ncol = 0, nrow = nrow(means)))
  rownames(gebvs) <- rownames(means)

  # estimate effects for each environment separately
  for (env in colnames(means)) {

    # run rrblup
    rr_model <- mixed.solve(y = means[, env], Z = markers)
    # get effects
    marker_effects_env <- data.frame(rr_model$u)
    colnames(marker_effects_env) <- env
    # calculate gebvs
    gebvs_env <- data.frame(as.numeric(rr_model$beta) + as.matrix(markers) %*% as.matrix(marker_effects_env))
    # calculate cor obs vs pred
    cor_obs_pred <- cor(as.numeric(means[, env]), as.numeric(gebvs_env[, env]), use = "complete.obs")
    cat("  ", env, " - correlation obs vs pred:", round(cor_obs_pred, digits = 2), "\n")
    # append effects to main df
    marker_effects <- cbind(marker_effects, marker_effects_env)
    gebvs <- cbind(gebvs, gebvs_env)

  }
  rm(env, rr_model, gebvs_env, cor_obs_pred, marker_effects_env)

}

if (marker_eff_model == "bayescpi") {

  cat("estimating marker effects by univariate BayesCpi\n")

  # create empty df to store results
  marker_effects <- data.frame(matrix(ncol = 0, nrow = ncol(markers)))
  rownames(marker_effects) <- colnames(markers)
  gebvs <- data.frame(matrix(ncol = 0, nrow = nrow(means)))
  rownames(gebvs) <- rownames(means)

  # estimate effects for each environment separately
  for (env in colnames(means)) {

    # run BayesCpi
    bcpi_model <- wgr(y = means[, env], X = as.matrix(markers), pi = 0.95, iv = FALSE)
    # get effects
    marker_effects_env <- data.frame(bcpi_model$b)
    colnames(marker_effects_env) <- env
    # calculate gebvs
    gebvs_env <- data.frame(bcpi_model$hat)
    colnames(gebvs_env) <- env
    # calculate cor obs vs pred
    cor_obs_pred <- cor(as.numeric(means[, env]), as.numeric(gebvs_env[, env]), use = "complete.obs")
    cat("  ", env, " - correlation obs vs pred:", round(cor_obs_pred, digits = 2), "\n")
    # append effects to main df
    marker_effects <- cbind(marker_effects, marker_effects_env)
    gebvs <- cbind(gebvs, gebvs_env)

  }
  rm(env, bcpi_model, gebvs_env, cor_obs_pred, marker_effects_env)

}

if (marker_eff_model == "mrr") {

  cat("estimating marker effects by multivariate ridge-regression\n")

  # estimate effects
  mrr_model <- mrr(Y = as.matrix(means), X = as.matrix(markers))
  # get effects
  marker_effects <- data.frame(mrr_model$b)
  colnames(marker_effects) <- colnames(means)
  rownames(marker_effects) <- colnames(markers)
  # calculate gebvs
  gebvs <- data.frame(mrr_model$hat)
  colnames(gebvs) <- colnames(means)
  rownames(gebvs) <- rownames(means)
  # calculate cor obs vs pred
  cor_obs_pred <- diag(cor(as.matrix(means), as.matrix(gebvs), use = "complete.obs"))
  for (env in names(cor_obs_pred)) {
    cat("  ", env, "- correlation obs vs pred:", round(as.numeric(cor_obs_pred[env]), digits = 2), "\n")
  }

}

if (marker_eff_model == "gwas") {

  cat("estimating marker effects by GWAS Q+K model\n")

  # create empty df to store results
  marker_effects <- data.frame(matrix(ncol = 0, nrow = ncol(markers)))
  rownames(marker_effects) <- colnames(markers)
  gebvs <- data.frame(matrix(ncol = 0, nrow = nrow(means)))
  rownames(gebvs) <- rownames(means)

  # estimate effects for each environment separately
  for (env in colnames(means)) {

    # run gwas q+k
    pheno <- rownames_to_column(means[, env, drop = FALSE], var = "genotype")
    geno <- rownames_to_column(markers + 1, var = "genotype")
    sink(paste0(output_folder, "/gwas.log"))
    gwas_model <- GAPIT(Y = pheno,
                        GD = geno,
                        GM = map,
                        kinship.algorithm = "VanRaden",
                        PCA.total = 5,
                        model = "MLM",
                        file.output = FALSE)
    sink()
    # get effects
    if (!all(gwas_model$GWAS$SNP == map$SNP)) stop("markers don't match")
    marker_effects_env <- data.frame(gwas_model$GWAS$effect)
    colnames(marker_effects_env) <- env
    # calculate gebvs
    gebvs_env <- gwas_model$Pred
    gebvs_env <- gebvs_env[match(pheno$genotype, gebvs_env$Taxa), "Prediction", drop = FALSE]
    colnames(gebvs_env) <- env
    # calculate cor obs vs pred
    cor_obs_pred <- cor(as.numeric(means[, env]), as.numeric(gebvs_env[, env]), use = "complete.obs")
    cat("  ", env, " - correlation obs vs pred:", round(cor_obs_pred, digits = 2), "\n")
    # append effects to main df
    marker_effects <- cbind(marker_effects, marker_effects_env)
    gebvs <- cbind(gebvs, gebvs_env)

  }
  rm(env, pheno, geno, gwas_model, gebvs_env, cor_obs_pred, marker_effects_env)

}

# write effects
marker_effects <- rownames_to_column(marker_effects, var = "marker")
if (no_missing_genotypes) {
  outfile_effects <- paste0(output_folder, "/marker_effects.", marker_eff_model, ".no-missing-genos.txt")
} else {
  outfile_effects <- paste0(output_folder, "/marker_effects.", marker_eff_model, ".txt")
}
fwrite(marker_effects, file = outfile_effects, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# write gebvs
gebvs <- rownames_to_column(gebvs, var = "genotype")
if (no_missing_genotypes) {
  outfile_gebvs <- paste0(output_folder, "/GEBVs.", marker_eff_model, ".no-missing-genos.txt")
} else {
  outfile_gebvs <- paste0(output_folder, "/GEBVs.", marker_eff_model, ".txt")
}
fwrite(gebvs, file = outfile_gebvs, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# markers_file <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt"
# # blues_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
# # output_folder <- "analysis/marker_effects/YLD"
# blues_file <- "data/1stStage_BLUEs.EHT-per-env.txt"
# output_folder <- "analysis/marker_effects/EHT"
# marker_eff_model <- "rrblup"
# # marker_eff_model <- "mrr"
# impute_effect <- "Add"
# impute_type <- "Middle"
# no_missing_genotypes <- FALSE
# # no_missing_genotypes <- TRUE
