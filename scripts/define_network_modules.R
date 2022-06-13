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

usage: Rscript define_network_modules.R [meff_net_Rdata] [output_folder] [...]

positional arguments:
  meff_net_Rdata              .RData file containing R variables from build_meff_network.R script
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --min-mod-size              minimum number of markers to be in a module (default: 50)
  --ME-diss-threshold         threshold to merge similar modules based on their eigengenes
                              (default: 0.25, i.e. correlation of 0.75)


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
meff_net_Rdata <- args[1]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
min_mod_size <- "50"
ME_diss_threshold <- 0.25

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) != pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--min-mod-size", "--ME-diss-threshold")
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
if (suppressWarnings(!is.na(as.integer(min_mod_size)))) {
  min_mod_size <- as.integer(min_mod_size)
} else {
  stop("Optional argument '--min-mod-size' should be an integer")
}

if (suppressWarnings(!is.na(as.numeric(ME_diss_threshold)))) {
  ME_diss_threshold <- as.numeric(ME_diss_threshold)
} else {
  if (ME_diss_threshold < 0 | ME_diss_threshold > 1) {
    stop("Optional argument '--ME-diss-threshold' should be a number between 0 and 1") 
  }
}



#### define network modules ----

# load data from previous step
load(meff_net_Rdata)

# in the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a marker
# branches of the dendrogram group together densely interconnected markers.
# module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”).
# there are several methods for branch cutting - WGCNA standard method is the Dynamic Tree Cut from the package dynamicTreeCut.

# module identification using dynamic tree cut
# the function returns modules labeled from the largest to smallest (label 0 is reserved for unassigned genes)
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = min_mod_size)

# # lists the sizes of the modules
# table(dynamicMods)

# convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
# table(dynamicColors)

# plot the module assignment under the gene dendrogram
sizeGrWindow(8,6)
pdf(file = paste0(output_folder, "/marker_dendrogram.pdf"), width = 12, height = 9)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Marker dendrogram and module colors")
dev.off()

# the Dynamic Tree Cut may identify modules whose marker effect profiles are very similar. It may be
# prudent to merge such modules since their markers are highly interconnexted. To quantify similarity
# of entire modules, we calculate their eigengenes and cluster them on their correlation

# calculate eigengenes
MEList <- moduleEigengenes(marker_effects, colors = dynamicColors)
MEs <- MEList$eigengenes
# calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

# plot the result
sizeGrWindow(7, 6)
pdf(file = paste0(output_folder, "/eigengenes_clustering.pdf"), width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# plot the cut line into the dendrogram
abline(h = ME_diss_threshold, col = "red")
dev.off()

# call an automatic merging function
merge <- mergeCloseModules(marker_effects, dynamicColors, cutHeight = ME_diss_threshold, verbose = 3)
# the merged module colors
mergedColors <- merge$colors
# eigengenes of the new merged modules:
mergedMEs <- merge$newMEs

# to see what the merging did to our module colors, we plot the gene dendrogram again,
# with the original and merged module colors underneath
sizeGrWindow(12, 9)
pdf(file = paste0(output_folder, "/marker_dendrogram_merged.pdf"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# in the subsequent analysis, we will use the merged module colors in mergedColors

# rename to moduleColors
moduleColors <- mergedColors
# construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors())
moduleLabels <- match(moduleColors, colorOrder)
MEs <- mergedMEs

# save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = paste0(outfolder, "/define_network_modules.RData"))



#### debug ----

# meff_net_Rdata <- "tests/networks/YLD/build_meff_network.RData"
# output_folder <- "tests/networks/YLD"
# min_mod_size <- 50
# ME_diss_threshold <- 0.25