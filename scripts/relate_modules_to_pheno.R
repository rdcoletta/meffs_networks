library(data.table)
library(WGCNA)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

usage <- function() {
  cat("
description: get correlations of module eigenvalues with mean trait values across environments.

usage: Rscript relate_modules_to_pheno.R [meff_mod_Rdata] [trait_file] [output_folder] [...]

positional arguments:
  meff_mod_Rdata              .RData file containing R variables from define_network_modules.R script
  trait_file                  file containing trait values per environment (3 columns: genotype, env, trait_value)
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --norm-method=VALUE         method to normalize means across multiple environments. Available methods are
                              'minmax' (default), 'zscore', 'none' (i.e. no normalization)
  --n-perm=VALUE              number of permutations to do (default: 1000)
  --seed=VALUE                seed number for shuffling trait data


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
trait_file <- args[2]
output_folder <- args[3]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
norm_method <- "minmax"
seed <- 8871
n_perm <- 1000

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--norm-method", "--seed")
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
if (!norm_method %in% c("minmax", "zscore", "none") ) {
  stop("Optional argument '--norm-method' should be 'minmax', 'zscore' or 'none'")
}



#### correlate trait data ----

# load module data
load(meff_mod_Rdata)

# load trait data
pheno <- fread(trait_file, header = TRUE, data.table = FALSE)

# normalize trait values, if requested
if (norm_method == "minmax") {
  # get max and min values of the data
  max <- max(pheno$real_pheno, na.rm = TRUE)
  min <- min(pheno$real_pheno, na.rm = TRUE)
  # normalize data
  pheno$real_pheno <- sapply(pheno$real_pheno, function(x) (x - min)/(max - min))
  rm(max, min) 
}
if (norm_method == "zscore") {
  # get mean and sd values of the data
  mean <- mean(pheno$real_pheno, na.rm = TRUE)
  sd <- sd(pheno$real_pheno, na.rm = TRUE)
  # normalize data
  pheno$real_pheno <- sapply(pheno$real_pheno, function(x) (x - mean)/sd)
  rm(mean, sd)
}

# shuffle data for permutation test
set.seed(seed)
for (i in 1:n_perm) {
  pheno[, paste0("ctrl_", i)] <- sample(pheno$real_pheno, replace = FALSE)
}

# get mean for each environment
means <- pheno %>% 
  group_by(env) %>% 
  summarize(across(-hybrid, ~ mean(.x, na.rm = TRUE))) %>% 
  ungroup() %>% 
  as.data.frame()

# format module data
rownames(MEs) <- gsub("trait_", "", rownames(MEs))
rownames(MEs) <- gsub(".", "-", rownames(MEs), fixed = TRUE)
# make sure samples match
means <- means[match(rownames(MEs), means$env), ]
means <- column_to_rownames(means, var = "env")
colnames(means)[1] <- "trait"

# calculate correlation between module and trait
moduleTraitCor <- cor(MEs, means$trait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(means))
# correct multiple testing
moduleTraitPvalue <- apply(moduleTraitPvalue, MARGIN = 2, function(x) p.adjust(x, method = "fdr"))

# display correlations and their p-values
pdf(file = paste0(output_folder, "/trait_module_correlations.pdf"), width = 10, height = 6)
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "trait",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships\n(FDR-corrected p-values)"))
dev.off()

# plot module vs trait pattern across envs
cor_plot <- rownames_to_column(MEs, "env") %>% 
  pivot_longer(-env, names_to = "mod_name", values_to = "mod_eigen") %>% 
  group_by(mod_name) %>% 
  mutate(trait_value = means$trait) %>% 
  ungroup() %>% 
  pivot_longer(mod_eigen:trait_value, names_to = "trait", values_to = "value") %>% 
  mutate(mod_name = factor(mod_name,
                           levels = names(moduleTraitCor[, 1]),
                           labels = paste0(names(moduleTraitCor[, 1]),
                                           " (", round(moduleTraitCor[, 1], digits = 2), ")"))) %>% 
  ggplot(aes(x = env, y = value, color = trait, group = trait)) +
  geom_line() +
  facet_wrap(~ mod_name) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = paste0(output_folder, "/trait_module_patterns.pdf"),
       plot = cor_plot, device = "pdf")



#### run permutation test with pheno data ----

# credits: https://dgarcia-eu.github.io/SocialDataScience/5_SocialNetworkPhenomena/056_PermutationTests/PermutationTests

# get correlations for permutated data
cor_perm <- cor(MEs, means[, -1], use = "p")

# prepare matrix to keep pvalues of permutation test
pval_perm <- matrix(nrow = nrow(moduleTraitCor), ncol = 1)
rownames(pval_perm) <- rownames(moduleTraitCor)

# calculate pvalue of correlations for each module
for (mod in rownames(moduleTraitCor)) {
  
  # # plot histogram
  # hist(cor_perm[mod, ], xlim = c(-1, 1))
  # abline(v = moduleTraitCor[mod, ], col="red")
  
  # get pvalue
  pval_perm[mod, ] <- (sum(abs(cor_perm[mod, ]) >= abs(moduleTraitCor[mod, ])) + 1) / length(cor_perm[mod, ])
  
}

# # ???
# # correct for multiple testing
# pval_perm <- apply(pval_perm, MARGIN = 2, function(x) p.adjust(x, method = "fdr"))
# # ???

# write pvalues from different methods
pval_table <- data.frame(cor_test = moduleTraitPvalue, perm_pheno = pval_perm)
pval_table <- rownames_to_column(pval_table, var = "module")
pval_table <- pval_table[order(pval_table$cor_test, pval_table$perm_pheno), ]
fwrite(pval_table, file = paste0(output_folder, "/module-pheno_pvals.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### run permutation test with module data ----

# set.seed(seed)
# # create an empty matrix to store correlation results
# shuffledMEsCor <- matrix(nrow = nrow(moduleTraitCor), ncol = n_perm)
# rownames(shuffledMEsCor) <- rownames(moduleTraitCor)
# 
# for (i in 1:n_perm) {
#   
#   # break marker-module association
#   shuffledColors <- sample(moduleColors, replace = FALSE)
#   # calculate eigengenes
#   shuffledMEList <- moduleEigengenes(marker_effects, colors = shuffledColors)
#   shuffledMEs <- shuffledMEList$eigengenes
#   
#   # format module data
#   rownames(shuffledMEs) <- gsub("trait_", "", rownames(shuffledMEs))
#   rownames(shuffledMEs) <- gsub(".", "-", rownames(shuffledMEs), fixed = TRUE)
#   # make sure samples match
#   shuffledMEs <- shuffledMEs[match(rownames(means), rownames(shuffledMEs)),
#                              match(rownames(shuffledMEsCor), colnames(shuffledMEs))]
#   
#   # calculate correlation between module and trait
#   shuffledMEsCor[, i] <- cor(means$trait, shuffledMEs)
#   
# }
# 
# # prepare matrix to keep pvalues of permutation test
# shuffledMEsPvalue <- matrix(nrow = nrow(moduleTraitCor), ncol = 1)
# rownames(shuffledMEsPvalue) <- rownames(moduleTraitCor)
# 
# # calculate pvalue of correlations for each module
# for (mod in rownames(shuffledMEsCor)) {
#   
#   
#   # #### CAN I MAKE THIS HISTOGRAM AS ABSOLUTE VALUES TO IGNORE SIGN OF CORRELATION?
#   # # plot histogram
#   # hist(abs(shuffledMEsCor[mod, ]), xlim = c(-1, 1))
#   # abline(v = abs(moduleTraitCor[mod, ]), col="red")
#   # ####
#   
#   # get pvalue
#   shuffledMEsPvalue[mod, ] <- (sum(abs(shuffledMEsCor[mod, ]) >= abs(moduleTraitCor[mod, ])) + 1) / length(shuffledMEsCor[mod, ])
#   
# }
# 
# # write pvalues from different methods
# pval_table <- data.frame(cor_test = moduleTraitPvalue, perm_pheno = pval_perm, perm_markers = shuffledMEsPvalue)
# pval_table <- rownames_to_column(pval_table, var = "module")
# fwrite(pval_table, file = paste0(output_folder, "/module-pheno_pvals.txt"),
#        quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# meff_mod_Rdata <- "tests/networks/YLD/define_network_modules.RData"
# trait_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
# output_folder <- "tests/networks/YLD"
# norm_method <- "minmax"
# seed <- 8871
# n_perm <- 1000