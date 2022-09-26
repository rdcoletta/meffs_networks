library(data.table)
library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

usage <- function() {
  cat("
description: get correlations of module eigenvalues with mean trait values across environments.

usage: Rscript relate_modules_to_env_idx.R [meff_mod_Rdata] [env_idx_file] [output_folder] [...]

positional arguments:
  meff_mod_Rdata              .RData file containing R variables from define_network_modules.R script
  env_idx_file                file containing environmental covariables (3 columns: env, covariable, value)
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --per-intervals             add this option if environment index was calculated per intervals instead of
                              the mean of the season. If this option is on, then 'env_idx_file' argument should
                              have 4 columns: env, covariable, intervals, value.


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
meff_mod_Rdata <- args[1]
env_idx_file <- args[2]
output_folder <- args[3]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
per_intervals <- FALSE

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--per-intervals")
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



#### correlate trait data ----

# load module data
load(meff_mod_Rdata)

# load trait data
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)

if (per_intervals & ncol(env_idx) != 4) {
  stop("expecting input file with 4 columns: env, covariable, intervals, value")
} 

# format env idx data
if (!per_intervals) {
  env_idx <- pivot_wider(env_idx, names_from = "covariable", values_from = "value")
} else {
  env_idx <- pivot_wider(env_idx, names_from = c("covariable", "intervals"), names_sep = "-", values_from = "value")
}

# format module data
rownames(MEs) <- gsub("trait_", "", rownames(MEs))
rownames(MEs) <- gsub(".", "-", rownames(MEs), fixed = TRUE)
# make sure samples match
env_idx <- env_idx[match(rownames(MEs), env_idx$env), ]
env_idx <- column_to_rownames(env_idx, var = "env")

# calculate correlation between module and trait
moduleTraitCor <- cor(MEs, env_idx, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(env_idx))
# correct multiple testing
moduleTraitPvalue <- apply(moduleTraitPvalue, MARGIN = 2, function(x) p.adjust(x, method = "fdr"))

# plot all env idx at once if not using intervals
if (!per_intervals) {
  
  # display correlations and their p-values
  pdf(file = paste0(output_folder, "/env_idx_module_correlations.pdf"), width = 15, height = 12)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(env_idx),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-env_idx relationships\n(FDR-corrected p-values)"))
  dev.off()
  
} else {
  
  # if using intervals, plot heatmap for each idx separately
  
  # get env idx names (without the intervals)
  idx_names <- sapply(names(env_idx), function(x) {
    x <- unlist(strsplit(x, split = "-"))
    x <- x[-length(x)]
    x <- paste0(x, collapse = "-")
    return(x)
  })
  idx_names <- unique(idx_names)
  
  for (idx in idx_names) {
    
    # keep only idx of interest
    moduleTraitCor_idx <- as.data.frame(moduleTraitCor) %>% 
      select(matches(paste0("^", idx, "-"), perl = TRUE)) %>% 
      as.matrix()
    moduleTraitPvalue_idx <- as.data.frame(moduleTraitPvalue) %>% 
      select(matches(paste0("^", idx, "-"), perl = TRUE)) %>% 
      as.matrix()
    
    # display correlations and their p-values
    pdf(file = paste0(output_folder, "/env_idx_module_correlations.", idx, "_intervals.pdf"),
        width = 15, height = 12)
    textMatrix = paste(signif(moduleTraitCor_idx, 2), "\n(",
                       signif(moduleTraitPvalue_idx, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor_idx)
    par(mar = c(6, 8.5, 3, 3));
    # display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor_idx,
                   xLabels = colnames(moduleTraitCor_idx),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-env_idx relationships\n(FDR-corrected p-values)"))
    dev.off()
    
  }
}

# reformat cor table
moduleTraitCor <- data.frame(moduleTraitCor, type = "cor", check.names = FALSE)
moduleTraitCor <- rownames_to_column(moduleTraitCor, var = "module")
# reformat pval table
moduleTraitPvalue <- data.frame(moduleTraitPvalue, type = "pval", check.names = FALSE)
moduleTraitPvalue <- rownames_to_column(moduleTraitPvalue, var = "module")
# merge cor and pval tables
pval_table <- rbind(moduleTraitCor, moduleTraitPvalue)

if (!per_intervals) {
  
  # reformat table again
  pval_table <- pivot_longer(pval_table,
                             cols = -c(module, type),
                             names_to = "env_idx",
                             values_to = "value")
  pval_table <- pivot_wider(pval_table,
                            names_from = "type",
                            values_from = "value")
  
} else {
  
  # reformat table again
  pval_table <- pivot_longer(pval_table,
                             cols = -c(module, type),
                             names_to = c("env_idx", "interval"),
                             names_sep = "-",
                             values_to = "value")
  pval_table <- pivot_wider(pval_table,
                            names_from = "type",
                            values_from = "value")
  
}

# reorder by pval
pval_table <- pval_table[order(pval_table$pval), ]
# write table
fwrite(pval_table, file = paste0(output_folder, "/module-env-idx_pvals.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# meff_mod_Rdata <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_off/define_network_modules.RData"
# # env_idx_file <- "data/env_covariables/env_covariables_means.txt"
# env_idx_file <- "data/env_covariables/env_covariables_means_per_intervals.txt"
# output_folder <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_off"
# per_intervals <- TRUE
# # per_intervals <- FALSE
