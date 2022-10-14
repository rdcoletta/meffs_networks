library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)


usage <- function() {
  cat("
description: plot the relationship betweenmodule eigeinvalues and environmental indices.

usage: Rscript plot_mod-env-idx.R [folder_base] [mod_env_idx_cor_file] [env_idx_file]
                                  [output_folder] [...]

positional arguments:
  folder_base                 path to folder with results of module-env idx relationship
  mod_env_idx_cor_file        file with correlations between module and env indices
  env_idx_file                environmental indices across environments
  output_folder               folder to output plots

optional argument:
  --help                      show this helpful message
  --p-value=VALUE             p-value threshold to filter correlations (default: 0.05)


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
env_idx_file <- args[3]
output_folder <- args[4]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default
p_value <- "0.05"

# assert to have the correct optional arguments
pos_args <- 4
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--p-value")
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
if (suppressWarnings(!is.na(as.numeric(p_value)))) {
  p_value <- as.numeric(p_value)
} else {
  stop("Optional argument '--p-value' should be a number")
}




#### plot module relationship with env idx ----

# load module-env index relationship
mod_env_idx_cor <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)
# select only significant associations
mod_env_idx_cor <- subset(mod_env_idx_cor, pval < p_value)
rownames(mod_env_idx_cor) <- 1:nrow(mod_env_idx_cor)
# round numbers
mod_env_idx_cor$cor <- sapply(mod_env_idx_cor$cor, function(x) round(x, digits = 2))
mod_env_idx_cor$pval <- sapply(mod_env_idx_cor$pval, function(x) round(x, digits = 2))

# load env idx results
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)

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
  
  # get env idx results 
  env_idx_net <- subset(env_idx, env %in% envs_to_keep & covariable == mod_env_idx_cor[row, "env_idx"]) %>% 
    rename(type = "covariable") %>%
    select(env, type, value)
  
  # plot ME and env idx in their own scales
  plot_mod_env_idx <- rbind(MEs, env_idx_net) %>%  
    ggplot(aes(x = env, y = value, color = type)) +
    geom_line(aes(group = type), show.legend = FALSE) +
    facet_wrap(~ type, nrow = 2, scales = "free_y") +
    labs(title = bquote("Module-Environmental index correlation =" ~ bold(.(as.character(mod_env_idx_cor[row, "cor"])))
                        ~ .(paste0("(pval: ", mod_env_idx_cor[row, "pval"], ")"))),
         subtitle = bquote("meff model:" ~ bold(.(mod_env_idx_cor[row, "meff_model"])) ~ "/"
                           ~ "normalization:" ~ bold(.(mod_env_idx_cor[row, "norm_method"])) ~ "/"
                           ~ "min module size:" ~ bold(.(as.character(mod_env_idx_cor[row, "minsize"]))) ~ "/"
                           ~ "pamStage:" ~ bold(.(mod_env_idx_cor[row, "pamStage"])))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(plot_mod_env_idx, filename = paste0(output_folder, "/", paste0(mod_env_idx_cor[row, 1:6], collapse = "_"), ".pdf"),
         device = "pdf", height = 10, width = 12)
  
  # normalize to same scale
  MEs$value <- sapply(MEs$value, function(x) (x - min(MEs$value))/(max(MEs$value) - min(MEs$value)))
  env_idx_net$value <- sapply(env_idx_net$value, function(x) (x - min(env_idx_net$value))/(max(env_idx_net$value) - min(env_idx_net$value)))
  
  # plot them together
  plot_mod_env_idx_norm <- rbind(MEs, env_idx_net) %>%  
    ggplot(aes(x = env, y = value, color = type)) +
    geom_line(aes(group = type)) +
    labs(title = bquote("Module-Environmental index correlation =" ~ bold(.(as.character(mod_env_idx_cor[row, "cor"])))
                        ~ .(paste0("(pval: ", mod_env_idx_cor[row, "pval"], ")"))),
         subtitle = bquote("meff model:" ~ bold(.(mod_env_idx_cor[row, "meff_model"])) ~ "/"
                           ~ "normalization:" ~ bold(.(mod_env_idx_cor[row, "norm_method"])) ~ "/"
                           ~ "min module size:" ~ bold(.(as.character(mod_env_idx_cor[row, "minsize"]))) ~ "/"
                           ~ "pamStage:" ~ bold(.(mod_env_idx_cor[row, "pamStage"]))),
         y = "normalized value") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(plot_mod_env_idx_norm, filename = paste0(output_folder, "/", paste0(mod_env_idx_cor[row, 1:6], collapse = "_"), ".norm.pdf"),
         device = "pdf", height = 10, width = 12)
  
}



#### debug ----

# folder_base <- "analysis/networks/YLD"
# # mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
# mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.means.txt"
# # env_idx_file <- "data/env_covariables/env_idx.txt"
# env_idx_file <- "data/env_covariables/env_covariables_means.txt"
# output_folder <- "analysis/networks/YLD/plots_mod-env-idx"
# p_value <- 0.05