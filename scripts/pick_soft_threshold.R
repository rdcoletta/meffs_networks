library(WGCNA)
library(data.table)
library(tibble)
library(tidyr)
library(ggplot2)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

usage <- function() {
  cat("
description: determine soft-threshold to use when building marker effect networks.

usage: Rscript pick_soft_threshold.R [marker_effects_file] [marker_hmp_file] [output_folder] [...]

positional arguments:
  marker_effects_file         file containing marker effects for a trait (3 columns: environment, marker, effect)
  marker_hmp_file             hapmap file containing marker data
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --norm-method               method to normalize effects across multiple environments. Available methods are
                              'minmax' (default), 'zscore', 'none' (i.e. no normalization)
  --cv-threshold              markers with coefficient of variation (CV) below this threshold will be removed.
                              If number is provided, then an absolute cut-off will be applied. Alternatively,
                              you can choose between 'Q1', 'Q2' or 'Q3' (default) to remove markers with CV below
                              the first, second or third quantiles.
  --debug                     add this option to randomly select 10,000 markers and speed up computation time


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
marker_effects_file <- args[1]
marker_hmp_file <- args[2]
output_folder <- args[3]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
norm_method <- "minmax"
cv_threshold <- "Q3"
debug <- FALSE

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--norm-method", "--cv-threshold", "--debug")
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

if (suppressWarnings(!is.na(as.numeric(cv_threshold)))) {
  cv_threshold <- as.numeric(cv_threshold)
} else {
  if (!cv_threshold %in% c("Q1", "Q2", "Q3")) {
    stop("Optional argument '--cv-threshold' should be a number or one of these quantiles: 'Q1', 'Q2' or 'Q3'")
  } 
}



#### filter marker data ----

# load marker effects
marker_effects <- fread(marker_effects_file, header = TRUE, data.table = FALSE)
marker_effects <- pivot_wider(marker_effects, names_from = "env", values_from = "effect")
colnames(marker_effects) <- gsub(".", "-", colnames(marker_effects), fixed = TRUE)

# load marker info
marker_info <- fread(marker_hmp_file, header = TRUE, data.table = FALSE)
marker_info <- marker_info[, 1:4]

# transform to rownames
marker_effects <- column_to_rownames(marker_effects, "marker")
marker_info <- column_to_rownames(marker_info, "rs#")

if (debug) {
  # randomly sample markers for testing
  set.seed(2726)
  marker_effects <- marker_effects[sort(sample(1:nrow(marker_effects), size = 10000, replace = FALSE)), ]
  marker_info <- marker_info[match(rownames(marker_effects), rownames(marker_info)), ]
}

# normalize marker effects 
if (norm_method == "minmax") {
  # get max and min values of the data
  max <- max(as.matrix(marker_effects))
  min <- min(as.matrix(marker_effects))
  # normalize data
  marker_effects <- data.frame(apply(marker_effects, MARGIN = c(1,2), function(x) (x - min)/(max - min)))
  rm(max, min) 
}
if (norm_method == "zscore") {
  # get mean and sd values of the data
  mean <- mean(as.matrix(marker_effects))
  sd <- sd(as.matrix(marker_effects))
  # normalize data
  marker_effects <- data.frame(apply(marker_effects, MARGIN = c(1,2), function(x) (x - mean)/sd))
  rm(mean, sd)
}

# calculate coefficient of variation
markers_cv <- apply(marker_effects,  MARGIN = 1, function(marker) abs(sd(marker) / mean(marker)))
# set CV treshold if a quantile (and not absolute) value was provided
if (cv_threshold %in% c("Q1", "Q2", "Q3")) {
  cv_threshold <- as.numeric(gsub("Q", "", cv_threshold))
  cv_threshold <- as.numeric(quantile(markers_cv)[cv_threshold + 1])
}
# plot cv distribution
plot_cv_cutoff <- ggplot(data.frame(cv = markers_cv)) +
  geom_histogram(aes(x = cv)) +
  geom_vline(xintercept = cv_threshold, color = "firebrick") #+ coord_cartesian(ylim = c(0,100))
ggsave(filename = paste0(output_folder, "/cv_cutoff.pdf"), plot = plot_cv_cutoff,
       device = "pdf", height = 10, width = 8)

# get markers that pass CV threshold
high_cv_markers <- names(which(markers_cv >= as.numeric(cv_threshold)))
n_markers_removed <- nrow(marker_effects) - length(high_cv_markers)
cat("removed ", n_markers_removed, " markers (", round((n_markers_removed / NROW(marker_effects) * 100), digits = 2),
    "%) with CV < ", cv_threshold, "\n", sep = "")
# remove markers with low CV
marker_effects <- marker_effects[high_cv_markers, ]
rm(markers_cv, cv_threshold, high_cv_markers, n_markers_removed)

# match gene/sample orders
marker_order <- intersect(rownames(marker_info), rownames(marker_effects))
marker_info <- marker_info[marker_order, ]
marker_effects <- marker_effects[marker_order, ]
if (!all(rownames(marker_info) == rownames(marker_effects))) stop("marker names don't match")

# make samples as rows and markers as columns
marker_effects <- data.frame(t(marker_effects))
marker_info <- data.frame(t(marker_info))



#### qc ----

# check for markers and samples with too many missing values
gsg <- goodSamplesGenes(marker_effects, verbose = 3)
# if 'gsg$allOK' returns TRUE, all genes have passed the cuts
# if not, we remove the offending genes and samples from the data:
if (!gsg$allOK) {
  # optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(marker_effects)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(marker_effects)[!gsg$goodSamples], collapse = ", ")));
  # remove the offending genes and samples from the data:
  marker_effects <- marker_effects[gsg$goodSamples, gsg$goodGenes]
}
# match marker/sample orders
marker_order <- intersect(colnames(marker_info), colnames(marker_effects))
marker_info <- marker_info[, marker_order]
marker_effects <- marker_effects[, marker_order]

# cluster the samples (in contrast to clustering genes that will come later) to see if there
# are any obvious outliers
sample_tree <- hclust(dist(marker_effects), method = "average")
sizeGrWindow(12,9)
pdf(file = paste0(output_folder, "/sample_clustering.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
# remove unnecessary variables
rm(gsg, sample_tree, marker_order)



#### pick soft threshold  ----

# choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# call the network topology analysis function
sft <- pickSoftThreshold(marker_effects, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
pdf(file = paste0(output_folder, "/scale-free_topology_fit_index.pdf"), width = 9, height = 6)
par(mfrow = c(1, 2))
cex1 <- 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# save data for next step
save(marker_effects, marker_info, sft, file = paste0(output_folder, "/pick_soft_threshold.RData"))



#### debug ----

# marker_effects_file <- "analysis/marker_effects/YLD/marker_effects.txt"
# marker_hmp_file <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt"
# output_folder <- "tests/networks/YLD"
# # norm_method <- "none"
# norm_method <- "minmax"
# # norm_method <- "zscore"
# cv_threshold <- "Q3"
# # cv_threshold <- 0.2
# # debug <- FALSE