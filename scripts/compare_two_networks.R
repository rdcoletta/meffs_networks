library(WGCNA)
library(data.table)
library(tidyr)
library(ggplot2)

usage <- function() {
  cat("
description: check if GWAS hits are present in significant marker effect modules.

usage: Rscript compare_two_networks.R [meff_mod_Rdata1] [meff_mod_Rdata2] [MM_net1_file] [MM_net2_file]
                                      [modules_net1_file] [modules_net2_file] [output_folder] [...]

positional arguments:
  meff_mod_Rdata1             .RData file from network 1 containing R variables from build_meff_network.R script
  meff_mod_Rdata2             .RData file from network 2 containing R variables from build_meff_network.R script
  MM_net1_file                file from network 1 with module membership matrix
  MM_net2_file                file from network 2 with module membership matrix
  modules_net1_file           file from network 1 with markers assigned to modules
  modules_net2_file           file from network 2 with markers assigned to modules
  output_folder               name of folder to save results

optional argument:
  --help                      show this helpful message
  --name-net1=[VALUE]         name of network 1 (default: network1)
  --name-net2=[VALUE]         name of network 2 (default: network2)
  --diff-n-markers            add this option if comparing networks with different number of markers


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
meff_mod_Rdata1 <- args[1]
meff_mod_Rdata2 <- args[2]
MM_net1_file <- args[3]
MM_net2_file <- args[4]
modules_net1_file <- args[5]
modules_net2_file <- args[6]
output_folder <- args[7]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
name_net1 <- "network1"
name_net2 <- "network2"
diff_n_markers <- FALSE

# assert to have the correct optional arguments
pos_args <- 7
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--name-net1", "--name-net2", "--diff-n-markers")
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



#### compare dendros ----

# only compare dendrograms if networks have the same markers
if (!diff_n_markers) {
  
  # load module definition from network 1
  load(meff_mod_Rdata1)
  # adjust names
  net1_MEs <- MEs
  net1_moduleLabels <- moduleLabels
  net1_moduleColors <- moduleColors
  net1_geneTree <- geneTree
  rm(MEs, moduleLabels, moduleColors, geneTree)
  
  # load module definition from network 2
  load(meff_mod_Rdata2)
  # adjust names
  net2_MEs <- MEs
  net2_moduleLabels <- moduleLabels
  net2_moduleColors <- moduleColors
  net2_geneTree <- geneTree
  rm(MEs, moduleLabels, moduleColors, geneTree)
  
  # plot dendrogram with the two networks
  pdf(file = paste0(output_folder, "/marker_dendrogram_with_modules.pdf"), width = 12, height = 9)
  plotDendroAndColors(dendro = net1_geneTree,
                      colors = cbind(net1_moduleColors, net2_moduleColors),
                      groupLabels = c(name_net1, name_net2),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
}



#### compare module membership ----

# load module membership matrices for each network
MM_net1 <- fread(MM_net1_file, header = TRUE, data.table = FALSE)
MM_net2 <- fread(MM_net2_file, header = TRUE, data.table = FALSE)

# if networks have different number of markers...
if (diff_n_markers) {
  # ... keep only common markers
  total_markers_net1 <- nrow(MM_net1)
  total_markers_net2 <- nrow(MM_net2)
  common_markers <- intersect(MM_net1$marker, MM_net2$marker)
  MM_net1 <- subset(MM_net1, marker %in% common_markers)
  MM_net2 <- subset(MM_net2, marker %in% common_markers)
}

# calculate correlation between module memberships of the two networks
MM_cor <- cor(MM_net1[, -1], MM_net2[, -1])
rownames(MM_cor) <- gsub("^MM", "", rownames(MM_cor), perl = TRUE)
colnames(MM_cor) <- gsub("^MM", "", colnames(MM_cor), perl = TRUE)

# load data with markers assigned to modules for each network
modules_net1 <- fread(modules_net1_file, header = TRUE, data.table = FALSE)
modules_net2 <- fread(modules_net2_file, header = TRUE, data.table = FALSE)

# if networks have different number of markers...
if (diff_n_markers) {
  # ... keep only common markers
  modules_net1 <- subset(modules_net1, marker %in% common_markers)
  modules_net2 <- subset(modules_net2, marker %in% common_markers)
}

# merge module assignment from both networks
merged_modules <- merge(x = modules_net1[modules_net1$source == "TOM", c("marker", "module")],
                        y = modules_net2[modules_net2$source == "TOM", c("marker", "module")],
                        by = "marker")
colnames(merged_modules) <- c("marker", "network1", "network2")

# create empty df to store results
mods_relationship <- data.frame(stringsAsFactors = FALSE)
for (mod_net1 in unique(modules_net1$module)) {
  
  # identify module relationship between networks
  mods_net2_in_net1 <- unique(merged_modules[merged_modules$network1 == mod_net1, "network2"])
  mods_net2_in_net1 <- t(MM_cor[mod_net1, mods_net2_in_net1, drop = FALSE])
  mods_net2_in_net1 <- data.frame(network1 = mod_net1,
                                  network2 = rownames(mods_net2_in_net1),
                                  MM_cor = as.numeric(mods_net2_in_net1))
  # add results to main data frame
  mods_relationship <- rbind(mods_relationship, mods_net2_in_net1)
  
}

# reorder factor levels
mods_relationship$network1 <- factor(mods_relationship$network1,
                                     levels = c("grey", sort(unique(modules_net1[modules_net1$module != "grey", "module"]))))
mods_relationship$network2 <- factor(mods_relationship$network2,
                                     levels = c("grey", sort(unique(modules_net2[modules_net2$module != "grey", "module"]))))

# plot module relationships
plot_mod_relations <- ggplot(mods_relationship) +
  geom_tile(aes(x = network1, y = network2, fill = MM_cor), color = "black") +
  # scale_fill_viridis_c(option = "B", limits = c(-1, 1), begin = 0, end = 0.9, name = "MM cor") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-1, 1), name = "MM cor") +
  theme_light() +
  theme(panel.grid.major.y = element_line(linetype = "dotted"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Module relationships between networks", x = name_net1, y = name_net2)

# if networks have different number of markers...
if (diff_n_markers) {
  # ...add a subtitle with percentage of markers kept from each network
  markers_net1_left <- round(length(common_markers) / total_markers_net1 * 100, digits = 1)
  markers_net2_left <- round(length(common_markers) / total_markers_net2 * 100, digits = 1)
  plot_mod_relations <- plot_mod_relations +
    labs(subtitle = paste0("- Markers remaining from '", name_net1, "' network: ", markers_net1_left, "%\n",
                           "- Markers remaining from '", name_net2, "' network: ", markers_net2_left, "%"))
}

ggsave(filename = paste0(output_folder, "/mod_relationships_between_nets.pdf"),
       plot = plot_mod_relations, device = "pdf", width = 12, height = 12)



#### debug ----

# meff_mod_Rdata1 <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_off/define_network_modules.RData"
# # meff_mod_Rdata2 <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_on/define_network_modules.RData"
# meff_mod_Rdata2 <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/define_network_modules.RData"
# MM_net1_file <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_off/module_membership.txt"
# # MM_net2_file <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_on/module_membership.txt"
# MM_net2_file <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/module_membership.txt"
# modules_net1_file <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_off/kDiff_per_module.txt"
# # modules_net2_file <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/pamStage_on/kDiff_per_module.txt"
# modules_net2_file <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/kDiff_per_module.txt"
# output_folder <- "analysis/networks/YLD/meff_gwas/norm_zscore/min_mod_size_50/compare_pamStage_on-vs-off"
# name_net1 <- "stagePAM_FALSE"
# name_net2 <- "stagePAM_TRUE"
# # diff_n_markers <- FALSE
# diff_n_markers <- TRUE
