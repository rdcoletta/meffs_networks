library(WGCNA)
library(data.table)
library(tibble)
library(tidyr)
library(ggplot2)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

usage <- function() {
  cat("
description: QC modules of a marker effect network.

usage: Rscript qc_network_modules.R [meff_mod_Rdata] [output_folder] [sft] [...]

positional arguments:
  meff_mod_Rdata              .RData file containing R variables from build_meff_network.R script
  output_folder               name of folder to save results
  sft                         soft threshold used to build network

optional argument:
  --help                      show this helpful message
  --n-hub-markers=VALUE       number of most important markers to print for each module (default: 10)


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
meff_mod_Rdata <- args[1]
output_folder <- args[2]
sft <- as.numeric(args[3])
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
n_hub_markers <- "10"

# assert to have the correct optional arguments
pos_args <- 3
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--n-hub-markers")
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
if (suppressWarnings(!is.na(as.integer(n_hub_markers)))) {
  n_hub_markers <- as.integer(n_hub_markers)
} else {
  stop("Optional argument '--n-hub-markers' should be an integer")
}




#### qc modules ----

# load data
load(meff_mod_Rdata)

# calculate adjacency and TOM again to define connection strength among markers
cat(paste0(Sys.time(), " - calculating adjacency and TOM matrices...\n"))
adjacency <- adjacency(marker_effects, power = sft)
TOM <- TOMsimilarity(adjacency)
rownames(TOM) <- rownames(adjacency)
colnames(TOM) <- rownames(adjacency)
cat(paste0(Sys.time(), " - done!\n\n"))

# get connectivity among markers
degrees_adj <- cbind(intramodularConnectivity(adjacency, moduleColors),
                     module = moduleColors, source = "adjacency")
degrees_TOM <- cbind(intramodularConnectivity(TOM, moduleColors),
                     module = moduleColors, source = "TOM")

# get clustering coefficients
cat(paste0(Sys.time(), " - getting clustering coefficients...\n"))
cat(paste0(Sys.time(), " -   adjacency matrix\n"))
degrees_adj$clusterCoeff  <- clusterCoef(adjacency)
cat(paste0(Sys.time(), " -   TOM matrix\n"))
degrees_TOM$clusterCoeff <- clusterCoef(TOM)
cat(paste0(Sys.time(), " - done!\n\n"))

# merge types of connection strength
degrees_adj <- rownames_to_column(degrees_adj, var = "marker")
degrees_TOM <- rownames_to_column(degrees_TOM, var = "marker")
degrees <- rbind(degrees_adj, degrees_TOM)
rm(degrees_adj, degrees_TOM)

# plot kDiff among markers
# (i.e. difference between number of connections within vs outside modules)
degrees$module <- factor(degrees$module)
plot_kDiff <- ggplot(degrees) + 
  facet_wrap(~ source, nrow = 2) +
  geom_boxplot(aes(x = module, y = kDiff, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(degrees$module)) +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0(output_folder, "/kDiff_per_module.pdf"),
       plot = plot_kDiff, device = "pdf")

plot_clusterCoeff <- ggplot(degrees) + 
  facet_wrap(~ source, nrow = 2) +
  geom_boxplot(aes(x = module, y = clusterCoeff, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(degrees$module)) +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0(output_folder, "/clusterCoeff_per_module.pdf"),
       plot = plot_clusterCoeff, device = "pdf")

# also write table with results
fwrite(degrees, file = paste0(output_folder, "/kDiff_per_module.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# get module membership values
# (i.e. the correlation of the module eigengene and the marker effects profile)
mod_membership <- signedKME(marker_effects, MEs, outputColumnName = "MM")
mod_membership <- rownames_to_column(mod_membership, var = "marker")
# write file
fwrite(mod_membership, file = paste0(output_folder, "/module_membership.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# find hub markers and calculate some network stats
for (strength_type in c("adjacency", "TOM")) {
  
  # subset by type of connection strength used
  degrees_source <- subset(degrees, source == strength_type)
  degrees_source <- degrees_source[, -ncol(degrees_source)]
  
  # get network stats
  mean_k <- mean(degrees_source$kTotal)
  density <- mean_k / (nrow(degrees_source) - 1)
  degree_centralization <- (max(degrees_source$kTotal) / nrow(degrees_source)) - density
  heterogeneity <- sd(degrees_source$kTotal) / mean_k
  # write file
  stats_net <- data.frame(mean_k = mean_k,
                          density = density,
                          degree_centralization = degree_centralization,
                          heterogeneity = heterogeneity)
  fwrite(stats_net, file = paste0(output_folder, "/network_stats.", strength_type, ".txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
  
  # create empty data frames to store
  stats_mod <- data.frame(stringsAsFactors = FALSE)
  hub_markers <- data.frame(stringsAsFactors = FALSE)
  
  for (mod in unique(moduleColors)) {
    
    # subset by module
    degrees_source_mod <- subset(degrees_source, module == mod)
    
    # get module stats
    mean_k <- mean(degrees_source_mod$kTotal)
    density <- mean_k / (nrow(degrees_source_mod) - 1)
    degree_centralization <- (max(degrees_source_mod$kTotal) / nrow(degrees_source_mod)) - density
    heterogeneity <- sd(degrees_source_mod$kTotal) / mean_k
    stats_mod <- rbind(stats_mod,
                       data.frame(module = mod,
                                  mean_k = mean_k,
                                  density = density,
                                  degree_centralization = degree_centralization,
                                  heterogeneity = heterogeneity))
    
    # reorder based on highest connections within module
    degrees_source_mod <- degrees_source_mod[order(degrees_source_mod$kWithin, decreasing = TRUE), ]
    # get marker hubs for that module
    mod_hub_markers <- degrees_source_mod[1:n_hub_markers, -ncol(degrees_source_mod)]
    # append to data frame
    hub_markers <- rbind(hub_markers, mod_hub_markers)
    
  }
  
  # write results
  fwrite(stats_mod, file = paste0(output_folder, "/module_stats.", strength_type, ".txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
  fwrite(hub_markers, file = paste0(output_folder, "/hub_markers_per_module.", strength_type, ".txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
  
  # plot results
  stats_mod <- pivot_longer(stats_mod, -module, names_to = "metric", values_to = "value")
  stats_mod$module <- factor(stats_mod$module)
  plot_stats_mod <- ggplot(stats_mod, aes(x = module, y = value, fill = module)) +
    facet_wrap(~ metric, scales = "free_y") +
    geom_col(color = "black", show.legend = FALSE) +
    scale_fill_manual(values = levels(stats_mod$module)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "Summary of network metrics",
         subtitle = paste0("(degree_centralization: ", round(stats_net$degree_centralization, 2), " / ",
                           "density: ", round(stats_net$density, 2), " / ",
                           "heterogeneity: ", round(stats_net$heterogeneity, 2), " / ",
                           "mean_k: ", round(stats_net$mean_k, 2), ")"))
  
  ggsave(filename = paste0(output_folder, "/summary_stats.", strength_type, ".pdf"),
         plot = plot_stats_mod, device = "pdf", width = 12, height = 8)
  
}



#### debug ----

# meff_mod_Rdata <- "tests/networks/YLD/define_network_modules.RData"
# output_folder <- "tests/networks/YLD"
# sft <- 12
# n_hub_markers <- 10

