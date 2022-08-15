library(data.table)
library(tibble)
library(tidyr)
library(ggplot2)

usage <- function() {
  cat("
description: plot distribution of coefficient of variation of marker effects.

usage: Rscript plot_cv_meffs.R [marker_effects_file] [output_folder] [...]

positional arguments:
  marker_effects_file         file containing marker effects for a trait (first column: marker names,
                              remaining columns: environments)
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --norm-method               method to normalize effects across multiple environments. Available methods are
                              'minmax' (default), 'zscore', 'none' (i.e. no normalization)
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
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
norm_method <- "minmax"
debug <- FALSE

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--norm-method", "--debug")
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




#### plot cv ----

# load marker effects
marker_effects <- fread(marker_effects_file, header = TRUE, data.table = FALSE)
# marker_effects <- pivot_wider(marker_effects, names_from = "env", values_from = "effect")
colnames(marker_effects) <- gsub(".", "-", colnames(marker_effects), fixed = TRUE)

# transform to rownames
marker_effects <- column_to_rownames(marker_effects, "marker")

# remove markers with exact same effects across envs
marker_effects <- marker_effects[!duplicated(marker_effects), ]

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

# using zscore normalization (or no normalization) can have negative marker effect values,
# which will end up inflating the coefficient of variation. Thus, I'm shiffiting the scale
# of each marker based on the most negative value
if (norm_method == "zscore" | norm_method == "none") {
  marker_effects <- apply(marker_effects,  MARGIN = 1, function(marker) marker - min(marker))
  marker_effects <- data.frame(t(marker_effects), stringsAsFactors = FALSE)
} 

# calculate coefficient of variation
markers_cv <- apply(marker_effects,  MARGIN = 1, function(marker) abs(sd(marker) / mean(marker)))
cv_quantiles <- as.numeric(quantile(markers_cv))

# plot cv distribution
plot_cv_cutoff <- ggplot(data.frame(cv = markers_cv)) +
  geom_histogram(aes(x = cv)) +
  geom_vline(xintercept = cv_quantiles[2], color = "firebrick") +
  geom_vline(xintercept = cv_quantiles[3], color = "firebrick") +
  geom_vline(xintercept = cv_quantiles[4], color = "firebrick") +
  annotate("text", x = cv_quantiles[2], y = -Inf, label = "Q1",
           vjust = -0.5, hjust = 1, color = "firebrick", fontface = "bold") +
  annotate("text", x = cv_quantiles[3], y = -Inf, label = "Q2",
           vjust = -0.5, hjust = 1, color = "firebrick", fontface = "bold") +
  annotate("text", x = cv_quantiles[4], y = -Inf, label = "Q3",
           vjust = -0.5, hjust = 1, color = "firebrick", fontface = "bold") +
  labs(caption = paste0("Q1: ", sum(markers_cv > cv_quantiles[2]), " markers left (cv = ",
                        round(cv_quantiles[2], digits = 2), ")\n",
                        "Q2: ", sum(markers_cv > cv_quantiles[3]), " markers left (cv = ",
                        round(cv_quantiles[3], digits = 2), ")\n",
                        "Q3: ", sum(markers_cv > cv_quantiles[4])," markers left (cv = ",
                        round(cv_quantiles[4], digits = 2), ")\n"))
ggsave(filename = paste0(output_folder, "/cv_cutoff.pdf"), plot = plot_cv_cutoff,
       device = "pdf", height = 10, width = 8)



#### debug ----

# marker_effects_file <- "analysis/marker_effects/YLD/marker_effects.rrblup.txt"
# # marker_effects_file <- "analysis/marker_effects/YLD/marker_effects.rrblup.no-missing-genos.txt"
# output_folder <- "tests/networks/YLD"
# # norm_method <- "none"
# norm_method <- "minmax"
# # norm_method <- "zscore"
# debug <- FALSE
