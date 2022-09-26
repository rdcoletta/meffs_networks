library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(factoextra)

usage <- function() {
  cat("
description: calculate PCA of environmental indices and respective contributions to the PCs.

usage: Rscript pca_env_idx.R [env_idx_file] [output_folder] [...]

positional arguments:
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
env_idx_file <- args[1]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
per_intervals <- FALSE

# assert to have the correct optional arguments
pos_args <- 2
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



#### PCA ----

# load data
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)

# format matrix
env_idx <- unite(env_idx, covariable:intervals, col = "covariable", sep = "_") %>% 
  pivot_wider(names_from = "covariable", values_from = "value") %>% 
  column_to_rownames(var = "env") %>% 
  as.matrix()

# perform PCA
pca <- prcomp(env_idx, scale = TRUE)

# write PCA results
pca_results <- pca$x[match(rownames(env_idx), rownames(pca$x)), ]
pca_results <- rownames_to_column(data.frame(pca_results), var = "env")
pca_results <- pivot_longer(pca_results, -env, names_to = "covariable", values_to = "value")
fwrite(pca_results, file = paste0(output_folder, "/pca_env_idx.txt"),
       quote = FALSE, sep = "\t", row.names = FALSE)

# calculate percent variance explained for each PC
pca_pve <- fviz_eig(pca)
ggsave(pca_pve, filename = paste0(output_folder, "/pca_scree_plot.pdf"),
       device = "pdf", height = 8, width = 10)

# plot PCA of environments
pca_envs <- fviz_pca_ind(pca, col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
ggsave(pca_envs, filename = paste0(output_folder, "/pca_envs.pdf"),
       device = "pdf", height = 12, width = 12)

# plot pca biplot with top 25 env idx contributing to PCs
pca_biplot <- fviz_pca_biplot(pca, col.ind = "contrib", col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                              select.var = list(contrib = 25), repel = TRUE, max.overlap = 100)
ggsave(pca_biplot, filename = paste0(output_folder, "/pca_biplot.pdf"),
       device = "pdf", height = 12, width = 12)


# calculate contributions (percent variance) of env covariables to each PC
pca_contrib <- data.frame(stringsAsFactors = FALSE)
for (pc in 1:length(unique(pca_results$covariable))) {
  
  # calculate contributions to PC
  pca_contrib_pc <- fviz_contrib(pca, choice = "var", axes = pc)
  pca_contrib_pc <- pca_contrib_pc$data
  # adjust row and col names
  dimnames(pca_contrib_pc) <- list(1:nrow(pca_contrib_pc), c("env_idx", paste0("PC", pc)))
  
  # add results to main df
  if (pc == 1) {
    pca_contrib <- pca_contrib_pc
  } else {
    pca_contrib <- merge(x = pca_contrib, y = pca_contrib_pc, by = "env_idx")
  }
  
}
# reorder df based on env idx
pca_contrib <- pca_contrib[order(pca_contrib$env_idx), ]
# write contributions
fwrite(pca_contrib, file = paste0(output_folder, "/pca_contributions.txt"),
       quote = FALSE, sep = "\t", row.names = FALSE)



#### debug ----

# env_idx_file <- "data/env_covariables/env_covariables_means_per_intervals.txt"
# output_folder <- "data/env_covariables"
# per_intervals <- TRUE