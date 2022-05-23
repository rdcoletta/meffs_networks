library(WGCNA)
library(data.table)
library(tibble)
library(tidyr)
library(ggplot2)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)


marker_effects_file <- "analysis/marker_effects/YLD/marker_effects.txt"
marker_info_file <- "data/usda_hybrids_projected-SVs-SNPs.poly.low-missing.pruned-100kb.geno-miss-0.25.hmp.txt"
output_folder <- "analysis/networks/YLD"
# norm_method <- "none"
norm_method <- "minmax"
# norm_method <- "zscore"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)



#### input file ----

# load marker effects
marker_effects <- fread(marker_effects_file, header = TRUE, data.table = FALSE)
marker_effects <- pivot_wider(marker_effects, names_from = "env", values_from = "effect")

# # randomly sample markers to speed up
# set.seed(8173)
# marker_effects <- marker_effects[sample(1:nrow(marker_effects), size = 50000, replace = FALSE), ]

# load marker info
marker_info <- fread(marker_info_file, header = TRUE, data.table = FALSE)
marker_info <- marker_info[, 1:4]

# transform to rownames
marker_effects <- column_to_rownames(marker_effects, "marker")
marker_info <- column_to_rownames(marker_info, "rs#")

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
# set a CV treshold
cv_threshold <- as.numeric(quantile(markers_cv)[4])
# cv_threshold <- 0.2
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
# all(rownames(marker_info) == rownames(marker_effects))

# make samples as rows and markers as columns
marker_effects <- data.frame(t(marker_effects))
marker_info <- data.frame(t(marker_info))


# check for genes and samples with too many missing values
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
# are any obvious outliers.
sample_tree <- hclust(dist(marker_effects), method = "average")
sizeGrWindow(12,9)
pdf(file = paste0(output_folder, "/sample_clustering.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()
# don't seem to have a sample much different than all the others
# remove unnecessary variables
rm(gsg, sample_tree, marker_order)

# # write files
# marker_effects_filtered <- rownames_to_column(marker_effects)
# marker_info_filtered <- rownames_to_column(marker_info)
# fwrite(marker_effects_filtered, "tgca-hnsc_expr.tumor.filtered.txt", sep = "\t", quote = FALSE, na = NA, row.names = FALSE)
# fwrite(marker_info_filtered, "genes_marker_info.tumor.filtered.txt", sep = "\t", quote = FALSE, na = NA, row.names = FALSE)
# rm(marker_effects_filtered, marker_info_filtered)



#### pick soft threshold  ----

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
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
save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna_tumor_network.RData")




####################################
# ADD CODE BELOW TO ANOTHER SCRIPT #
####################################


#### build weighted network ----

# sft$powerEstimate
# # choose the power 10, which is the lowest power for which the scale-free topology fit index curve
# # flattens out upon reaching a high value (in this case, ~0.8).
sft$powerEstimate <- 10

# calculate the adjacencies, using the soft thresholding power 10
adjacency <- adjacency(marker_effects, power = sft$powerEstimate)

sizeGrWindow(12,9)
pdf(file = paste0(output_folder, "/scale-free_topology_plot.pdf"), width = 12, height = 9)
scaleFreePlot(adjacency, nBreaks = 10, truncated = TRUE, removeFirst = FALSE)
dev.off()

# To minimize effects of noise and spurious associations, we transform the adjacency into
# Topological Overlap Matrix, and calculate the corresponding dissimilarity:

# Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM
rm(adjacency)


# use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf(file = paste0(output_folder, "/marker_clustering_TOM-diss.pdf"), width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Marker clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# TOMplot(diss = dissTOM, dendro = geneTree)

# save data for next step
save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna_tumor_network.RData")




####################################
# ADD CODE BELOW TO ANOTHER SCRIPT #
####################################

#### define network modules ----

# In the clustering tree (dendrogram), each leaf, that is a short vertical line, corresponds to a gene.
# Branches of the dendrogram group together densely interconnected, highly co-expressed genes.
# Module identification amounts to the identification of individual branches (”cutting the branches off the dendrogram”).
# There are several methods for branch cutting; our standard method is the Dynamic Tree Cut from the package dynamicTreeCut.
# The next snippet of code illustrates its use

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

# The function returned 44 modules labeled 1–44 largest to smallest. Label 0 is reserved for unassigned genes.
# The above command lists the sizes of the modules

# We now plot the module assignment under the gene dendrogram:

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf(file = paste0(output_folder, "/marker_dendrogram.pdf"), width = 12, height = 9)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Marker dendrogram and module colors")
dev.off()

# The Dynamic Tree Cut may identify modules whose expression profiles are very similar. It may be
# prudent to merge such modules since their genes are highly co-expressed. To quantify co-expression
# similarity of entire modules, we calculate their eigengenes and cluster them on their correlation:

# Calculate eigengenes
MEList = moduleEigengenes(marker_effects, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge:
MEDissThres = 0.25

# Plot the result
sizeGrWindow(7, 6)
pdf(file = paste0(output_folder, "/eigenmarkers_clustering.pdf"), width = 7, height = 6)
plot(METree, main = "Clustering of module eigenmarkers",
     xlab = "", sub = "")
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(marker_effects, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


# To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
sizeGrWindow(12, 9)
pdf(file = paste0(output_folder, "/marker_dendrogram_merged.pdf"), width = 12, height = 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# In the subsequent analysis, we will use the merged module colors in mergedColors. We save the
# relevant variables for use in subsequent parts of the tutorial:

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "wgcna_tumor_network.RData")

# resulted in 37 modules
