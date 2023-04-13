library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

usage <- function() {
  cat("
description: summarize LD results for each module of each network.

usage: Rscript qc_permutation_results.R [mod_env_idx_file] [mod_env_idx_perm_file] [...]

positional arguments:
  mod_env_idx_file            file with correlations between modules and env indices for real network
  mod_env_idx_perm_file       file with correlations between modules and env indices for permuted networks

optional argument:
  --help                      show this helpful message
  --meff-model=[VALUE]        name of marker effect models
  --norm-method=[VALUE]       name of normalization methods
  --minsize=[VALUE]           name of minimum module size
  --pamStage=[VALUE]          name of pamStage values
  --pval-threshold=[VALUE]    p-value threshold to filter results (default: 0.05)

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
mod_env_idx_file <- args[1]
mod_env_idx_perm_file <- args[2]

# set default
meff_model <- "rrblup,gwas"
norm_method <- "minmax,zscore"
minsize <- "25,50,100"
pamStage <- "on,off"
pval_threshold <- "0.05"

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {
  
  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--meff-model", "--norm-method", "--minsize", "--pamStage", "--pval-threshold")
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

# adjust format from optional arguments
model <- unlist(strsplit(meff_model, split = ","))
norm <- unlist(strsplit(norm_method, split = ","))
size <- as.numeric(unlist(strsplit(minsize, split = ",")))
pam <- unlist(strsplit(pamStage, split = ","))
pval_threshold <- as.numeric(pval_threshold)



#### qc permutations ----

# load correlation file from observed network
mod_env_idx <- fread(mod_env_idx_file, header = TRUE, data.table = FALSE)
# summarize correlation values
obs_cor_vals <- mod_env_idx %>%
  filter(meff_model == model & norm_method == norm & minsize == size & pamStage == pam) %>%
  mutate(n_tests = n(),
         n_mods = length(unique(module)),
         n_idx = length(unique(env_idx))) %>%
  filter(pval < pval_threshold) %>%
  summarize(net_type = "observed",
            n_tests = unique(n_tests),
            n_mods = unique(n_mods),
            n_idx = unique(n_idx),
            n_sig_mods = length(unique(module)),
            n_sig_idx = length(unique(env_idx)),
            sig_cor = round(min(abs(cor)), digits = 2),
            sig_pval = round(max(pval), digits = 4))

# get minimum correlation value below 0.05
cor_threshold <- obs_cor_vals$sig_cor

# load correlation file from permutated networks
mod_env_idx_perm <- fread(mod_env_idx_perm_file, header = TRUE, data.table = FALSE)
# get summary of permutated networks
perm_summary <- mod_env_idx_perm %>%
  filter(minsize == size & pamStage == pam) %>%
  mutate(iter = factor(iter)) %>%
  group_by(iter) %>%
  mutate(n_tests = n(),
         n_mods = length(unique(module)),
         n_idx = length(unique(env_idx))) %>%
  summarize(net_type = "permutation",
            n_tests = unique(n_tests),
            n_mods = unique(n_mods),
            n_idx = unique(n_idx))
# get number of modules above correlation threshold
perm_cor_vals <- mod_env_idx_perm %>%
  filter(minsize == size & pamStage == pam) %>%
  mutate(iter = factor(iter)) %>%
  group_by(iter) %>%
  filter(cor >= cor_threshold | cor <= -cor_threshold) %>%
  summarize(n_sig_mods = length(unique(module)),
            n_sig_idx = length(unique(env_idx)),
            sig_cor = round(min(abs(cor)), digits = 2),
            sig_pval = round(max(pval), digits = 4))
# merge summary
perm_cor_vals <- merge(x = perm_summary, y = perm_cor_vals, by = "iter", all.x = TRUE)
perm_cor_vals$n_sig_mods[is.na(perm_cor_vals$n_sig_mods)] <- 0
perm_cor_vals$n_sig_idx[is.na(perm_cor_vals$n_sig_idx)] <- 0
rm(perm_summary)

# check how many modules are associated with env PCs in each permutation
mods_pcs_perm <- data.frame(table(perm_cor_vals$n_sig_mods))
colnames(mods_pcs_perm) <- c("n_sig_mods", "count")

# write file
outfile <- unlist(strsplit(mod_env_idx_perm_file, split = "/"))
outfile <- paste0(outfile[-length(outfile)], collapse = "/")
outfile <- paste0(outfile, "/modules_sig-cor_with_pcs.txt")
fwrite(mods_pcs_perm, file = outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# sum(perm_cor_vals$n_sig_mods >= obs_cor_vals$n_sig_mods)
# sum(!is.na(perm_cor_vals$sig_cor))



#### get distribution of module memberships ----

# get eigenmarker and correlate with each marker effect within module
# then plot distribution of correlation values within modules (and obtain mean and std dev)

# get significant modules from observed network
mod_env_idx_sig <- mod_env_idx %>%
  filter(meff_model == model & norm_method == norm & minsize == size & pamStage == pam & pval < pval_threshold)

# load data
load(paste0("analysis/networks/YLD/meff_", model, "/norm_", norm,
            "/min_mod_size_", size, "/pamStage_", pam, "/define_network_modules.RData"))
rm(geneTree, marker_info)
# calculate module membership (i.e. correlation marker effect and module eigenmarker)
MMs <- signedKME(marker_effects, MEs, outputColumnName = "MM")
colnames(MMs) <- gsub("^MM", "", colnames(MMs), perl = TRUE)

# create empty df
obs_MM_vals <- data.frame(stringsAsFactors = FALSE)
# get mean and std err for each module
for (mod in unique(mod_env_idx_sig$module)) {
  
  # adjust module name to match names in vector and MM data frame
  mod <- gsub("^ME", "", mod, perl = TRUE)
  # get which markers are in each module
  markers_mod <- moduleColors == mod
  # filter module membership
  sig_mod_MM <- MMs[markers_mod, mod]
  # get MM values
  obs_MM_vals <- rbind(obs_MM_vals,
                       data.frame(network = "observed",
                                  module = mod,
                                  n_markers = length(sig_mod_MM),
                                  MM_mean = mean(abs(sig_mod_MM)),
                                  MM_se = sd(abs(sig_mod_MM)) / sqrt(length(sig_mod_MM))))
  
}


# get significant modules from permutated networks
mod_env_idx_perm_sig <- mod_env_idx_perm %>%
  filter(minsize == size & pamStage == pam & (cor >= cor_threshold | cor <= -cor_threshold))

# create empty df
perm_MM_vals <- data.frame(stringsAsFactors = FALSE)
# get mean and std err for each module in each iteration
for (iter in unique(mod_env_idx_perm_sig$iter)) {
  
  # get file name
  meff_mod_Rdata_iter <- paste0("analysis/networks/YLD/meff_", model, "/norm_", norm, "/permutation/iter", iter,
                                "/min_mod_size_", size, "/pamStage_", pam, "/define_network_modules.RData")
  
  if (file.exists(meff_mod_Rdata_iter)) {
    
    # load data
    load(meff_mod_Rdata_iter)
    rm(geneTree, marker_info)
    # calculate module membership (i.e. correlation marker effect and module eigenmarker)
    MMs <- signedKME(marker_effects, MEs, outputColumnName = "MM")
    colnames(MMs) <- gsub("^MM", "", colnames(MMs), perl = TRUE)
    # get modules in iteration
    mods_iter <- unique(mod_env_idx_perm_sig[mod_env_idx_perm_sig$iter == iter, "module"])
    
    for (mod in mods_iter) {
      
      # adjust module name to match names in vector and MM data frame
      mod <- gsub("^ME", "", mod, perl = TRUE)
      # get which markers are in each module
      markers_mod <- moduleColors == mod
      # filter module membership
      sig_mod_MM <- MMs[markers_mod, mod]
      # get MM values
      perm_MM_vals <- rbind(perm_MM_vals,
                            data.frame(network = paste0("iter", iter),
                                       module = mod,
                                       n_markers = length(sig_mod_MM),
                                       MM_mean = mean(abs(sig_mod_MM)),
                                       MM_se = sd(abs(sig_mod_MM)) / sqrt(length(sig_mod_MM))))
      
    }
  }
}

# merge files
MM_vals <- rbind(obs_MM_vals, perm_MM_vals)

# write file
outfile <- unlist(strsplit(mod_env_idx_perm_file, split = "/"))
outfile <- paste0(outfile[-length(outfile)], collapse = "/")
outfile <- paste0(outfile, "/mod-membership_sig-cor_with_pcs.txt")
fwrite(MM_vals, file = outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)




#### debug ----

# mod_env_idx_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
# mod_env_idx_perm_file <- "analysis/networks/YLD/meff_gwas/norm_minmax/permutation/permutation_module-env-idx_per_network.pca.txt"
# model <- "gwas"
# norm <- "minmax"
# size <- 25
# pam <- "off"
# pval_threshold <- 0.05
# # mod_env_idx_perm_file <- "analysis/networks/YLD/meff_gwas/norm_zscore/permutation/permutation_module-env-idx_per_network.pca.txt"
# # model <- "gwas"
# # norm <- "zscore"
# # size <- 10
# # pam <- "on"
# # pval_threshold <- 0.1