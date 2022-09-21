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

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--norm-method")
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
colnames(pheno) <- c("genotype", "env", "trait_value")

# normalize trait values, if requested
if (norm_method == "minmax") {
  # get max and min values of the data
  max <- max(pheno$trait_value, na.rm = TRUE)
  min <- min(pheno$trait_value, na.rm = TRUE)
  # normalize data
  pheno$trait_value <- sapply(pheno$trait_value, function(x) (x - min)/(max - min))
  rm(max, min)
}
if (norm_method == "zscore") {
  # get mean and sd values of the data
  mean <- mean(pheno$trait_value, na.rm = TRUE)
  sd <- sd(pheno$trait_value, na.rm = TRUE)
  # normalize data
  pheno$trait_value <- sapply(pheno$trait_value, function(x) (x - mean)/sd)
  rm(mean, sd)
}

# get mean for each environment
means <- pheno %>%
  group_by(env) %>%
  summarize(across(-genotype, ~ mean(.x, na.rm = TRUE))) %>%
  ungroup() %>%
  as.data.frame()

# format module data
rownames(MEs) <- gsub("trait_", "", rownames(MEs))
rownames(MEs) <- gsub(".", "-", rownames(MEs), fixed = TRUE)
# make sure samples match
means <- means[match(rownames(MEs), means$env), ]
rownames(means) <- NULL
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


# reformat cor table
moduleTraitCor <- data.frame(trait_avg = moduleTraitCor, type = "cor", check.names = FALSE)
moduleTraitCor <- rownames_to_column(moduleTraitCor, var = "module")
# reformat pval table
moduleTraitPvalue <- data.frame(trait_avg = moduleTraitPvalue, type = "pval", check.names = FALSE)
moduleTraitPvalue <- rownames_to_column(moduleTraitPvalue, var = "module")
# merge cor and pval tables
pval_table <- rbind(moduleTraitCor, moduleTraitPvalue)


# reformat table again
pval_table <- pivot_longer(pval_table, cols = -c(module, type), names_to = "env_idx", values_to = "value")
pval_table <- pivot_wider(pval_table, names_from = "type", values_from = "value")
# reorder by pval
pval_table <- pval_table[order(pval_table$pval), ]
# write table
fwrite(pval_table, file = paste0(output_folder, "/module-pheno_pvals.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# plot module vs trait pattern across envs
cor_plot <- rownames_to_column(MEs, "env") %>%
  pivot_longer(-env, names_to = "mod_name", values_to = "mod_eigen") %>% 
  full_join(y = rownames_to_column(means, "env"), by = "env") %>% 
  rename(module = mod_name, trait_value = trait) %>%
  full_join(y = pval_table[, colnames(pval_table) != "env_idx"], by = "module") %>%
  pivot_longer(mod_eigen:trait_value, names_to = "type", values_to = "value") %>% 
  mutate(module = factor(module,
                         levels = pval_table[, "module", drop = TRUE],
                         labels = apply(pval_table, MARGIN = 1, function(mod) {
                           paste0(mod["module"], "\n(cor = ",
                                  round(as.numeric(mod["cor"]), 2), ", pval = ",
                                  round(as.numeric(mod["pval"]), 2), ")")
                         }))) %>% 
  ggplot(aes(x = env, y = value, color = type, group = type)) +
  geom_line() +
  facet_wrap(~ module, nrow = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.minor = element_blank())

ggsave(filename = paste0(output_folder, "/trait_module_patterns.pdf"), plot = cor_plot,
       device = "pdf", height = 12, width = 3 * ceiling(length(MEs) / 4))



#### debug ----

# meff_mod_Rdata <- "tests/networks/YLD/define_network_modules.RData"
# trait_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
# output_folder <- "tests/networks/YLD"
# norm_method <- "minmax"
