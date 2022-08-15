library(WGCNA)
library(data.table)
library(tibble)
library(tidyr)
library(ggplot2)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

usage <- function() {
  cat("
description: build marker effect networks.

usage: Rscript build_meff_network.R [sft_Rdata] [output_folder] [...]

positional arguments:
  sft_Rdata                   .RData file containing R variables from pick_soft_threshold.R script
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --soft-threshold=VALUE      the lowest power for which the scale-free topology fit index curve (default: 10)


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
sft_Rdata <- args[1]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
soft_threshold <- "10"

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--soft-threshold")
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
if (suppressWarnings(!is.na(as.integer(soft_threshold)))) {
  soft_threshold <- as.integer(soft_threshold)
} else {
  stop("Optional argument '--soft-threshold' should be an integer")
}



#### build weighted network ----

# load data from previous step
load(sft_Rdata)

# calculate the adjacencies, using the soft thresholding power 10
adjacency <- adjacency(marker_effects, power = soft_threshold)

# plot scale free topology
pdf(file = paste0(output_folder, "/scale-free_topology_plot.pdf"), width = 12, height = 9)
par(mfrow = c(1, 2))
hist(rowSums(adjacency), main = paste0("Connectivity distribution\n(power = ", soft_threshold, ")"), xlab = "Connectivity")
scaleFreePlot(adjacency, nBreaks = 10, truncated = TRUE, removeFirst = FALSE)
dev.off()

# to minimize effects of noise and spurious associations, we transform the adjacency into
# topological Overlap Matrix, and calculate the corresponding dissimilarity

# trn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
rm(adjacency)
dissTOM <- 1 - TOM
rm(TOM)

# use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes

# call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# plot the resulting clustering tree (dendrogram)
pdf(file = paste0(output_folder, "/marker_clustering_TOM-diss.pdf"), width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Marker clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# save data for next step
save(marker_effects, marker_info, geneTree, dissTOM, file = paste0(output_folder, "/build_meff_network.RData"))



#### debug ----

# sft_Rdata <- "tests/networks/YLD/pick_soft_threshold.RData"
# output_folder <- "tests/networks/YLD"
# soft_threshold <- 12
