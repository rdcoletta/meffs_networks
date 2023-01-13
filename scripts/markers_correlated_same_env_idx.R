library(data.table)
library(gtools)
library(ggplot2)
library(UpSetR)

usage <- function() {
  cat("
description: check markers in modules that are correlated to the same env idx while considering LD among them.

usage: Rscript markers_correlated_same_env_idx.R [folder_base] [mod_env_idx_cor_file] [geno_data] [output_folder] [...]

positional arguments:
  folder_base                 path to folder with results of module-env idx relationship
  mod_env_idx_cor_file        file with correlations between module and principal components
  geno_data                   hapmap file to calculate LD among markers
  output_folder               folder to output plots

optional argument:
  --help                      show this helpful message
  --p-value=VALUE             p-value threshold to filter correlations (default: 0.05)
  --tassel-path=PATH          absolute path to TASSEL 5
                              (default: '/home/hirschc1/della028/software/tassel-5-standalone')
  --plink-path=PATH           absolute path to PLINK 1.9
                              (default: '/home/hirschc1/della028/software/plink_linux_x86_64_20200219')


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
geno_data <- args[3]
output_folder <- args[4]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default
p_value <- "0.05"
tassel_path <- "/home/hirschc1/della028/software/tassel-5-standalone"
plink_path <- "/home/hirschc1/della028/software/plink_linux_x86_64_20200219"

# assert to have the correct optional arguments
pos_args <- 4
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--p-value", "--tassel-path", "--plink-path")
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

    cat ("calculating LD for markers correlated with ", idx, "...\n", sep = "")

    # get list of marker names to calculate LD
    list_markers_for_ld <- paste0(output_folder, "/markers_correlated_", idx, ".txt")
    fwrite(x = data.frame(unique(unlist(list_markers_mod_idx))), file = list_markers_for_ld,
           quote = FALSE, na = NA, row.names = FALSE, col.names = FALSE)
    # create plink file to calculate LD
    commands_hmp2plk <- paste0(tassel_path, "/run_pipeline.pl",
                               " -Xmx40g -importGuess ", geno_data,
                               " -includeSiteNamesInFile ", list_markers_for_ld,
                               " -export ", gsub(".txt", "", list_markers_for_ld),
                               " -exportType Plink > /dev/null")
    system(commands_hmp2plk)
    # calculate LD across chromosomes
    commands_plink <- paste0(plink_path, "/plink",
                             " --file ", gsub(".txt", ".plk", list_markers_for_ld),
                             " --make-founders",
                             " --r2 gz dprime with-freqs inter-chr",
                             " --ld-window-r2 0.9",
                             " --out ", gsub(".txt", "", list_markers_for_ld))
    system(commands_plink)

    cat ("...done!\n", sep = "")

    # read LD file
    ld_file <- gsub(".txt", ".ld.gz", list_markers_for_ld)
    
    if (file.exists(ld_file)) {
      
      # load ld file
      ld_file <- fread(ld_file, header = TRUE, data.table = FALSE)
      
      # substitute the name of markers in LD to each other with a representative marker
      for (mod in names(list_markers_mod_idx)) {
        
        for (marker in 1:length(list_markers_mod_idx[[mod]])) {
          # get marker name
          marker_name <- list_markers_mod_idx[[mod]][marker]
          # check which markers are in LD with it
          marker_ld <- ld_file[ld_file$SNP_A == marker_name | ld_file$SNP_B == marker_name, c("SNP_A", "SNP_B")]
          # if there's at least one marker in LD
          if (nrow(marker_ld) > 0) {
            # get the name of the first marker in the group
            # (it will always be the same because file is ordered by position)
            marker_ld <- unique(unlist(marker_ld))[1]
            # change the marker name to the representative marker
            list_markers_mod_idx[[mod]][marker] <- marker_ld
          }
        }
        
        # keep redundant marker names
        list_markers_mod_idx[[mod]] <- unique(list_markers_mod_idx[[mod]])
        
      }
      
      if (length(list_markers_mod_idx) > 1) {
        
        # plot intersections of markers
        upset_plot <- upset(fromList(list_markers_mod_idx), order.by = "freq", mb.ratio = c(0.55, 0.45), nsets = 100)
        pdf(file = paste0(output_folder, "/markers_cor_", idx, ".pdf"), onefile = FALSE, width = 12, height = 10)
        print(upset_plot)
        dev.off()
        
      }
      
    } else {
      cat("\n--- The file ", ld_file, "does not exist because there were no markers in LD ---\n\n")
    }
    
  }
}



#### debug ----

# folder_base <- "analysis/networks/YLD"
# mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
# geno_data <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt"
# output_folder <- "analysis/networks/YLD/overlap_markers"
# p_value <- 0.1
# tassel_path <- "/home/hirschc1/della028/software/tassel-5-standalone"
# plink_path <- "/home/hirschc1/della028/software/plink_linux_x86_64_20200219"
