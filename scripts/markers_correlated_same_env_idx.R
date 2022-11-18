library(data.table)
library(gtools)
library(ggplot2)
library(UpSetR)

usage <- function() {
  cat("
description: check markers in modules that are correlated to the same env idx..

usage: Rscript markers_correlated_same_env_idx.R [folder_base] [mod_env_idx_cor_file] [output_folder] [...]

positional arguments:
  folder_base                 path to folder with results of module-env idx relationship
  mod_env_idx_cor_file        file with correlations between module and principal components
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
output_folder <- args[3]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default
p_value <- "0.05"

# assert to have the correct optional arguments
pos_args <- 3
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




#### plot markers correlated to same env idx ----

# load module-env index relationship
mod_env_idx_cor <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)
# select only significant associations
mod_env_idx_cor <- subset(mod_env_idx_cor, pval < p_value)
rownames(mod_env_idx_cor) <- 1:nrow(mod_env_idx_cor)
# round numbers
mod_env_idx_cor$cor <- sapply(mod_env_idx_cor$cor, function(x) round(x, digits = 2))
mod_env_idx_cor$pval <- sapply(mod_env_idx_cor$pval, function(x) round(x, digits = 2))

# split data to loop through groups
if ("interval" %in% colnames(mod_env_idx_cor)) {
  mod_env_idx_cor_split <- split(mod_env_idx_cor, mod_env_idx_cor[, c("env_idx", "interval")])
} else {
  mod_env_idx_cor_split <- split(mod_env_idx_cor, mod_env_idx_cor[, "env_idx"])
}

# plot relationship per env_idx
for (idx in names(mod_env_idx_cor_split)) {

  # get only modules associated with that env_idx
  sig_env_idx <- mod_env_idx_cor_split[[idx]]

  if (nrow(sig_env_idx) > 0) {

    # create empty dataframe to store markers from different networks
    markers_mod_idx <- data.frame(stringsAsFactors = FALSE)

    for (row in 1:nrow(sig_env_idx)) {

      # get details about network
      model <- sig_env_idx[row, "meff_model"]
      norm <- sig_env_idx[row, "norm_method"]
      size <- sig_env_idx[row, "minsize"]
      pam <- sig_env_idx[row, "pamStage"]
      module <- gsub("^ME", "", sig_env_idx[row, "module"], perl = TRUE)

      # get markers in module for a specific network
      network_mod_markers <- paste0(folder_base, "/meff_", model, "/norm_", norm,
                                    "/min_mod_size_", size, "/pamStage_", pam,
                                    "/kDiff_per_module.txt")
      network_mod_markers <- fread(network_mod_markers, header = TRUE, data.table = FALSE)
      network_mod_markers <- network_mod_markers[network_mod_markers$source == "TOM", c("module", "marker")]
      network_mod_markers <- network_mod_markers[network_mod_markers$module == module, ]

      # add network settings
      network_mod_markers <- data.frame(network = paste(model, norm, size, pam, sep = "_"),
                                        network_mod_markers, stringsAsFactors = FALSE)
      network_mod_markers <- tidyr::unite(network_mod_markers, network:module, col = "network_module", sep = "-")
      # get accuracy results
      markers_mod_idx <- rbind(markers_mod_idx, network_mod_markers)

    }
    rm(network_mod_markers)

    # transform df into list
    list_markers_mod_idx <- list()
    for (net in unique(markers_mod_idx$network_module)) {
      module_name <- unique(markers_mod_idx[markers_mod_idx$network == net, "network_module"])
      list_markers_mod_idx[[module_name]] <- markers_mod_idx[markers_mod_idx$network == net, "marker"]
    }
    # names(list_markers_mod_idx) <- sapply(names(list_markers_mod_idx), function(x) unlist(strsplit(x, split = "-"))[2])
    names(list_markers_mod_idx) <- sapply(names(list_markers_mod_idx), function(x) gsub("-", "\n", x))

    if (length(list_markers_mod_idx) > 1) {

      # plot intersections of markers
      upset_plot <- upset(fromList(list_markers_mod_idx), order.by = "freq", mb.ratio = c(0.55, 0.45), nsets = 100)
      pdf(file = paste0(output_folder, "/markers_cor_", idx, ".pdf"), onefile = FALSE, width = 12, height = 10)
      print(upset_plot)
      dev.off()

    }
  }
}



#### debug ----

# folder_base <- "analysis/networks/YLD"
# mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
# # mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.per-intervals.txt"
# output_folder <- "analysis/networks/YLD/plots_mod-env-idx"
# p_value <- 0.1
