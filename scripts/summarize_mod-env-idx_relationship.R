library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(ggrepel)

usage <- function() {
  cat("
description: summarize gwas enrichment results of marker effect networks.

usage: Rscript summarize_mod-env-idx_relationship.R [folder_base] [...]

positional arguments:
  folder_base                 path to folder with results of gwas enrichment

optional argument:
  --help                      show this helpful message
  --meff-model=[LIST]         comma-separated list of marker effect models
  --norm-method=[LIST]        comma-separated list of normalization methods
  --minsize=[LIST]            comma-separated list of minimum module size
  --pamStage=[LIST]           comma-separated list of pamStage values
  --idx-type=[LIST]           comma-separated list of types of env indices

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

# set default
meff_model <- "rrblup,gwas"
norm_method <- "minmax,zscore"
minsize <- "25,50,100"
pamStage <- "on,off"
idx_type <- "means,intervals,pca,fw"

# assert to have the correct optional arguments
if (length(args) < 1) stop(usage(), "missing positional argument(s)")

if (length(args) > 1) {
  
  opt_args <- args[-1]
  opt_args_allowed <- c("--meff-model", "--norm-method", "--minsize", "--pamStage", "--idx-type")
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
meff_model <- unlist(strsplit(meff_model, split = ","))
norm_method <- unlist(strsplit(norm_method, split = ","))
minsize <- as.numeric(unlist(strsplit(minsize, split = ",")))
pamStage <- unlist(strsplit(pamStage, split = ","))
idx_type <- unlist(strsplit(idx_type, split = ","))



#### summarize module-env idx relationships ----

for (type in idx_type) {
  
  # create empty df to store results
  mod_env_idx_results <- data.frame(stringsAsFactors = FALSE)
  
  for (model in meff_model) {
    for (norm in norm_method) {
      for (size in minsize) {
        for (pam in pamStage) {
          
          # get folder with enrichment results for a network setting
          network_mod_env <- paste0(folder_base, "/meff_", model, "/norm_", norm,
                                    "/min_mod_size_", size, "/pamStage_", pam, 
                                    "/mod_env-idx_", type, "/module-env-idx_pvals.txt")
          
          # get results
          network_mod_env <- try(fread(network_mod_env, header = TRUE, data.table = FALSE))
          
          if (class(network_mod_env) != "try-error") {
            # add network settings
            network_mod_env <- data.frame(meff_model = model, norm_method = norm,
                                          minsize = size, pamStage = pam,
                                          network_mod_env, stringsAsFactors = FALSE)
            # get accuracy results
            mod_env_idx_results <- rbind(mod_env_idx_results, network_mod_env)
          }
          
        }
      }
    }
  }
  
  # write summary
  fwrite(mod_env_idx_results, file = paste0(folder_base, "/module-env-idx_per_network.", type, ".txt"))
  
  # adjust factor levels
  mod_env_idx_results$module <- gsub("^ME", "", mod_env_idx_results$module, perl = TRUE)
  mod_env_idx_results$module <- factor(mod_env_idx_results$module,
                                       levels = unique(mod_env_idx_results$module))
  
  # create base plot based on env idx type
  if (type == "intervals") {
    plot_mod_env_idx <- ggplot(data = mod_env_idx_results,
                               aes(x = pval, y = cor, fill = module,
                                   label = paste0(env_idx, "_", interval, "\n", module)))
  } else {
    plot_mod_env_idx <- ggplot(data = mod_env_idx_results,
                               aes(x = pval, y = cor, fill = module,
                                   label = paste0(env_idx, "\n", module)))
  }
  # add remaining layers to plot
  plot_mod_env_idx <- plot_mod_env_idx +
    facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_vline(xintercept = 0.5, color = "gray50") +
    geom_point(shape = 21, color = "black", show.legend = FALSE) +
    geom_vline(xintercept = 0.05, color = "firebrick") +
    geom_vline(xintercept = 0.1, color = "firebrick", linetype = "dashed") +
    geom_label_repel(data = subset(mod_env_idx_results, pval < 0.1), max.overlaps = 50,
                     fill = "white", size = 2, show.legend = FALSE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) +
    scale_fill_manual(values = levels(mod_env_idx_results$module)) +
    labs(title = "Correlation between module eigenvalues and env index",
         x = "p-values", y = "correlation") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # save plot
  ggsave(filename = paste0(folder_base, "/mod-env-idx_correlation.", type, ".pdf"),
         plot = plot_mod_env_idx, device = "pdf", width = 18, height = 14)
  
  # save only most significant correlations
  mod_env_idx_results_sig <- filter(mod_env_idx_results, pval < 0.2)
  fwrite(mod_env_idx_results_sig, file = paste0(folder_base, "/mod-env-idx_most-sig-corr.", type, ".txt"))
  
}




#### debug ----

# folder_base <- "analysis/networks/YLD"
# meff_model <- c("rrblup", "gwas")
# norm_method <- c("minmax", "zscore")
# minsize <- c(25, 50, 100)
# pamStage <- c("on", "off")
# idx_type <- c("means", "intervals", "pca", "fw")