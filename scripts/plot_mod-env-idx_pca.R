library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)


usage <- function() {
  cat("
description: build marker effect networks.

usage: Rscript plot_mod-env-idx_pca.R [folder_base] [mod_env_idx_cor_file] [pca_env_idx_file]
                                      [pca_contrib_file] [output_folder] [...]

positional arguments:
  folder_base                 path to folder with results of module-env idx relationship
  mod_env_idx_cor_file        file with correlations between module and principal components
  pca_env_idx_file            PCA results across environments
  pca_contrib_file            PCA loadings of env covariables
  output_folder               folder to output plots

optional argument:
  --help                      show this helpful message


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
folder_base <- args[1]
mod_env_idx_cor_file <- args[2]
pca_env_idx_file <- args[3]
pca_contrib_file <- args[4]
output_folder <- args[5]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# assert to have the correct optional arguments
pos_args <- 5
if (length(args) != pos_args) stop(usage(), "missing positional argument(s)")




#### plot module relationship with PCs ----

# load module-env index relationship
mod_env_idx_cor <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)
# select only significant associations
mod_env_idx_cor <- subset(mod_env_idx_cor, pval < 0.05)
rownames(mod_env_idx_cor) <- 1:nrow(mod_env_idx_cor)
mod_env_idx_cor$cor <- sapply(mod_env_idx_cor$cor, function(x) round(x, digits = 2))
mod_env_idx_cor$pval <- sapply(mod_env_idx_cor$pval, function(x) round(x, digits = 2))

# load pca results
pca_env_idx <- fread(pca_env_idx_file, header = TRUE, data.table = FALSE)
pca_contrib <- fread(pca_contrib_file, header = TRUE, data.table = FALSE)

for (row in 1:nrow(mod_env_idx_cor)) {
  
  # load module data
  meff_mod_Rdata <- paste0(folder_base,
                           "/meff_", mod_env_idx_cor[row, "meff_model"],
                           "/norm_", mod_env_idx_cor[row, "norm_method"],
                           "/min_mod_size_", mod_env_idx_cor[row, "minsize"],
                           "/pamStage_", mod_env_idx_cor[row, "pamStage"],
                           "/define_network_modules.RData")
  load(meff_mod_Rdata)
  
  # get module eigenvalues
  MEs <- rownames_to_column(MEs, var = "env") %>% 
    select(env, mod_env_idx_cor[row, "module"]) %>% 
    pivot_longer(-env, names_to = "type", values_to = "value") %>% 
    mutate(env = gsub(".", "-", env, fixed = TRUE),
           type = gsub("^ME", "", type, perl = TRUE))
 
  # get environments to keep
  envs_to_keep <- unique(MEs$env)
  
  # get PCA results 
  pca_env_idx_net <- subset(pca_env_idx, env %in% envs_to_keep & covariable == mod_env_idx_cor[row, "env_idx"]) %>% 
    rename(type = "covariable") %>%
    select(env, type, value)
  
  # plot ME and PCA in their own scales
  plot_mod_pc <- rbind(MEs, pca_env_idx_net) %>%  
    ggplot(aes(x = env, y = value, color = type)) +
    geom_line(aes(group = type), show.legend = FALSE) +
    facet_wrap(~ type, nrow = 2, scales = "free_y") +
    labs(title = bquote("Module-PC correlation =" ~ bold(.(as.character(mod_env_idx_cor[row, "cor"])))
                        ~ .(paste0("(pval: ", mod_env_idx_cor[row, "pval"], ")"))),
         subtitle = bquote("meff model:" ~ bold(.(mod_env_idx_cor[row, "meff_model"])) ~ "/"
                           ~ "normalization:" ~ bold(.(mod_env_idx_cor[row, "norm_method"])) ~ "/"
                           ~ "min module size:" ~ bold(.(as.character(mod_env_idx_cor[row, "minsize"]))) ~ "/"
                           ~ "pamStage:" ~ bold(.(mod_env_idx_cor[row, "pamStage"])))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(plot_mod_pc, filename = paste0(output_folder, "/mod_pc.", paste0(mod_env_idx_cor[row, 1:6], collapse = "_"), ".pdf"),
         device = "pdf", height = 10, width = 12)
  
  # normalize to same scale
  MEs$value <- sapply(MEs$value, function(x) (x - min(MEs$value))/(max(MEs$value) - min(MEs$value)))
  pca_env_idx_net$value <- sapply(pca_env_idx_net$value, function(x) (x - min(pca_env_idx_net$value))/(max(pca_env_idx_net$value) - min(pca_env_idx_net$value)))
  
  # plot them together
  plot_mod_pc_norm <- rbind(MEs, pca_env_idx_net) %>%  
    ggplot(aes(x = env, y = value, color = type)) +
    geom_line(aes(group = type)) +
    labs(title = bquote("Module-PC correlation =" ~ bold(.(as.character(mod_env_idx_cor[row, "cor"])))
                        ~ .(paste0("(pval: ", mod_env_idx_cor[row, "pval"], ")"))),
         subtitle = bquote("meff model:" ~ bold(.(mod_env_idx_cor[row, "meff_model"])) ~ "/"
                           ~ "normalization:" ~ bold(.(mod_env_idx_cor[row, "norm_method"])) ~ "/"
                           ~ "min module size:" ~ bold(.(as.character(mod_env_idx_cor[row, "minsize"]))) ~ "/"
                           ~ "pamStage:" ~ bold(.(mod_env_idx_cor[row, "pamStage"]))),
         y = "normalized value") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(plot_mod_pc_norm, filename = paste0(output_folder, "/mod_pc_norm.", paste0(mod_env_idx_cor[row, 1:6], collapse = "_"), ".pdf"),
         device = "pdf", height = 10, width = 12)
  
  # get contributions to pc
  pca_contrib_net <- pca_contrib[, c("env_idx", mod_env_idx_cor[row, "env_idx"])]
  colnames(pca_contrib_net)[2] <- "contrib"
  # keep only top 20 covariables
  pca_contrib_net <- pca_contrib_net[order(pca_contrib_net$contrib, decreasing = TRUE)[1:20], ]
  pca_contrib_net$env_idx <- factor(pca_contrib_net$env_idx,
                                    levels = pca_contrib_net$env_idx[order(pca_contrib_net$contrib, decreasing = TRUE)])
  
  # plot contributions
  plot_pc_contrib <- ggplot(pca_contrib_net) +
    geom_col(aes(x = env_idx, y = contrib)) +
    labs(title = bquote("Contributions to" ~ bold(.(mod_env_idx_cor[row, "env_idx"]))),
         x = "Environmental covariables",
         y = "Contributions (%)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(plot_pc_contrib, filename = paste0(output_folder, "/pc_contrib.", mod_env_idx_cor[row, "env_idx"], ".pdf"),
         device = "pdf", height = 10, width = 16)
  
}



#### debug ----

# folder_base <- "analysis/networks/YLD"
# mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
# pca_env_idx_file <- "data/env_covariables/pca_env_idx.txt"
# pca_contrib_file <- "data/env_covariables/pca_contributions.txt"
# output_folder <- "analysis/networks/YLD/mod_pc_contrib"
