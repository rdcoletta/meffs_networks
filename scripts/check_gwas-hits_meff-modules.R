library(WGCNA)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(network)
library(sna)
library(ggnetwork)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

usage <- function() {
  cat("
description: check if GWAS hits are present in significant marker effect modules.

usage: Rscript check_gwas-hits_meff-modules.R [gwas_filename] [modules_filename] [mod_pheno_filename]
                                              [mod_ld_folder][output_folder] [...]

positional arguments:
  gwas_filename                 file with gwas hits
  meff_filename                 file containing marker effects for a trait (first column: marker names,
                                remaining columns: environments)
  modules_filename              file with markers assigned to modules and with networks stats (e.g. kDiff)
  mod_pheno_filename            file with p-values of module-phenotype associations
  mod_ld_folder                 folder with ld files of each module
  output_folder                 name of folder to save results

optional argument:
  --help                        show this helpful message
  --type-connect-str=[VALUE]    choose the type of connection strenght between nodes: 'TOM' (default) or
                                'adjacency'
  --soft-threshold=[VALUE]      the lowest power for which the scale-free topology fit index curve (default: 10)
  --edge-threshold=[VALUE]      adjacency threshold for including edges in the network plots




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
gwas_filename <- args[1]
meff_filename <- args[2]
modules_filename <- args[3]
mod_pheno_filename <- args[4]
mod_ld_folder <- args[5]
output_folder <- args[6]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
type_connect_str <- "TOM"
soft_threshold <- "12"
edge_threshold <- "0.25"

# assert to have the correct optional arguments
pos_args <- 6
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--type-connect-str", "--soft-threshold", "--edge-threshold")
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
if (!type_connect_str %in% c("TOM", "adjacency")) {
  stop("Optional argument '--type-connect-str' should be either 'TOM' or 'adjacency'")
}

if (suppressWarnings(!is.na(as.integer(soft_threshold)))) {
  soft_threshold <- as.integer(soft_threshold)
} else {
  stop("Optional argument '--soft-threshold' should be an integer")
}

if (suppressWarnings(!is.na(as.numeric(edge_threshold)))) {
  edge_threshold <- as.numeric(edge_threshold)
} else {
  stop("Optional argument '--edge-threshold' should be a number")
}



#### check for gwas hits ----

# load data
gwas <- fread(gwas_filename, header = TRUE, data.table = FALSE)
meff <- fread(meff_filename, header = TRUE, data.table = FALSE)
modules <- fread(modules_filename, header = TRUE, data.table = FALSE)
mod_pheno <- fread(mod_pheno_filename, header = TRUE, data.table = FALSE)

# keep either adjacency or TOM values
modules <- subset(modules, source == type_connect_str)

# get gwas hits -- this includes hits in ld with each other
gwas_hits <- unique(gwas[!is.na(gwas[, 3]), 3])
# keep only hits available in this network
gwas_hits <- gwas_hits[gwas_hits %in% modules$marker]

# get order of hub markers in each module
modules <- group_by(modules, module) %>%
  mutate(hub_pos = rank(desc(kDiff)),
         n_mod = n()) %>%
  arrange(hub_pos, .by_group = TRUE)

# add p-values to modules
mod_pheno$module <- gsub("^ME", "", mod_pheno$module, perl = TRUE)
colnames(mod_pheno)[2:3] <- c("pval_cor_test", "pval_perm")
modules <- merge(x = modules, y = mod_pheno, by = "module")
rm(mod_pheno)

# find which modules have these gwas hits
modules$is_gwas_hit <- modules$marker %in% gwas_hits

# how many of all gwas hits are in a significant module
cat("gwas hit in a significant module: ",
    round(sum(modules$is_gwas_hit & modules$pval_cor_test < 0.05) * 100 / length(gwas_hits), 2), "% (cor test) or ",
    round(sum(modules$is_gwas_hit & modules$pval_perm < 0.05) * 100 / length(gwas_hits), 2), "% (perm)\n",
    sep = "")

# get name of modules associated with trait
modules_cor_trait <- unique(modules[modules$pval_cor_test < 0.05 | modules$pval_perm < 0.05, "module"])

cat("analyzing signficant modules...\n")
# test for enrichment of gwas hits in modules associated with the trait
for (mod in modules_cor_trait) {
  
  cat(paste0("  ", mod, " (", which(modules_cor_trait == mod), "/", length(modules_cor_trait), ")\n"))
  
  # subset data
  mod_subset <- subset(modules, module == mod)

  # run hypergeometric test
  pval_gwas_hyper <- phyper(q = sum(mod_subset$is_gwas_hit) - 1,
                            m = length(gwas_hits),
                            n = nrow(modules) - length(gwas_hits),
                            k = nrow(mod_subset),
                            lower.tail = FALSE, log.p = FALSE)

  # check effects of gwas hits across envs
  meff_gwas_plot <- filter(meff, marker %in% gwas_hits) %>%
    mutate(in_module = factor(if_else(marker %in% mod_subset$marker, true = TRUE, false = FALSE),
                              levels = c(FALSE, TRUE))) %>%
    pivot_longer(-c(marker, in_module), names_to = "env", values_to = "effects") %>%
    ggplot() +
    geom_line(aes(x = env, y = effects, group = marker, color = in_module)) +
    scale_color_manual(values = c("#8080804D", "firebrick")) +
    labs(title = paste0("Effects of GWAS hits in ", mod, " module across envs"),
         subtitle = paste0("(markers: ", nrow(mod_subset), "/", nrow(modules), " - ",
                           "gwas hits: ", sum(mod_subset$is_gwas_hit), "/", length(gwas_hits), " - ",
                           "pval hyper: ", round(pval_gwas_hyper, 3), ")")) +
    theme_bw() +
    theme(panel.grid = element_blank())

  # print(meff_gwas_plot)
  ggsave(filename = paste0(output_folder, "/meff_gwas.", mod, ".pdf"),
         plot = meff_gwas_plot, device = "pdf", width = 10, height = 6)


  # load file
  ld_module <- paste0(mod_ld_folder, "/ld_markers_", mod, ".ld.gz")
  ld_module <- fread(ld_module, header = TRUE, data.table = FALSE)

  # plot distribution r2 for all markers
  ld_dist_plot <- ggplot(ld_module, aes(x = R2)) +
    geom_histogram(color = "black", fill = mod) +
    labs(title = paste0("LD distribution of ", mod, " module")) +
    theme_bw()

  # print(ld_dist_plot)
  ggsave(filename = paste0(output_folder, "/ld_dist.", mod, ".pdf"),
         plot = ld_dist_plot, device = "pdf", width = 6, height = 6)

  # select markers in that module
  mod_markers <- modules[modules$module == mod & modules$source == type_connect_str, "marker"]

  # get effects for markers in the module
  mod_meff <- meff[meff$marker %in% mod_markers, ]
  rownames(mod_meff) <- NULL
  mod_meff <- column_to_rownames(mod_meff, var = "marker")
  mod_meff <- t(mod_meff)

  # calculate TOM for that module
  mod_TOM <- TOMsimilarityFromExpr(mod_meff, power = soft_threshold, TOMType = "unsigned")
  dimnames(mod_TOM) <- list(colnames(mod_meff), colnames(mod_meff))
  # reformat TOM matrix to create a network
  mod_TOM <- exportNetworkToCytoscape(mod_TOM,
                                      weighted = TRUE,
                                      threshold = edge_threshold,
                                      nodeNames = colnames(mod_meff))

  # get edges of that module
  mod_edges <- data.frame(lapply(mod_TOM$edgeData, as.character), stringsAsFactors = FALSE)
  mod_edges$weight <- as.numeric(mod_edges$weight)
  # get ld between two markers in network
  mod_edges$R2 <- apply(mod_edges[, c("fromNode", "toNode")], MARGIN = 1, function(marker_pair, ld_module) {

    r2_value <- ld_module[(ld_module$SNP_A == marker_pair["fromNode"] & ld_module$SNP_B == marker_pair["toNode"]) |
                            (ld_module$SNP_A == marker_pair["toNode"] & ld_module$SNP_B == marker_pair["fromNode"]), "R2"]
    if (length(r2_value) == 0) r2_value <- NA

    return(r2_value)

  }, ld_module = ld_module)
  # add categorical variable relating two markers --> are they in LD or not?
  mod_edges$LD <- ifelse(mod_edges$R2 > 0.9, yes = TRUE, no = FALSE)

  # create network data
  mod_network <- network(mod_edges, directed = FALSE, matrix.type = "edgelist")

  # get node names
  node_names <- get.vertex.attribute(mod_network, attrname = "vertex.names")
  # add kWithin to node attributes
  node_kWithin <- modules[modules$module == mod & modules$source == type_connect_str, c("marker", "kWithin")]
  node_kWithin <- node_kWithin[match(node_names, node_kWithin$marker), ]
  set.vertex.attribute(mod_network,
                       attrname = "kWithin",
                       value = node_kWithin$kWithin)
  # add whether marker is gwas hit or not
  node_gwas_hit <- node_names %in% gwas_hits
  set.vertex.attribute(mod_network,
                       attrname = "is_gwas_hit",
                       value = node_gwas_hit)

  # add vertex layout to network
  mod_network <- ggnetwork(mod_network, layout = "fruchtermanreingold", weights = "weight", cell.jitter = 0.5)
  # View(mod_network)

  # visualize network by hub markers
  plot_net_hub <- ggplot(mod_network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "#8080800D") +
    geom_nodes(aes(size = kWithin, color = is_gwas_hit)) +
    scale_color_manual(values = c("black", "firebrick")) +
    theme_blank() +
    labs(title = paste0(mod, " module"),
         subtitle = paste0("(markers: ", nrow(mod_subset), "/", nrow(modules), " - ",
                           "gwas hits: ", sum(mod_subset$is_gwas_hit), "/", length(gwas_hits), " - ",
                           "pval hyper: ", round(pval_gwas_hyper, 3), ")"))

  ggsave(filename = paste0(output_folder, "/mod_viz.", mod, ".by-hubs.pdf"),
         plot = plot_net_hub, device = "pdf", width = 10, height = 12)

  # visualize network by LD status
  plot_net_ld <- ggplot(mod_network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(color = LD)) +
    geom_nodes(color = "black") +
    scale_color_manual(values = c("#8080800D", "firebrick")) +
    theme_blank() +
    labs(title = paste0(mod, " module"),
         subtitle = paste0("(markers: ", nrow(mod_subset), "/", nrow(modules), " - ",
                           "gwas hits: ", sum(mod_subset$is_gwas_hit), "/", length(gwas_hits), " - ",
                           "pval hyper: ", round(pval_gwas_hyper, 3), ")"))

  ggsave(filename = paste0(output_folder, "/mod_viz.", mod, ".by-ld.pdf"),
         plot = plot_net_ld, device = "pdf", width = 10, height = 12)

}
cat("done!\n\n")


# phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
#
# where:
# q = (number of markers from gwas hits in module) - 1;
# m = number of markers in gwas hits
# n = total number of markers on network - number of markers in gwas hits;
# k = number of markers in module.

# Why "overlap - 1" in phyper()?
#   By default, the R function phyper computes the inclusive lower tail of the distribution: P(X ≤ x).
#   With the option “lower.tail=FALSE”, phyper() returns the exclusive upper tail P(X&gtx).
#   We want the *inclusive* upper tail : P-value = P(X≥x). For this, we can compute the exclusive upper tail of the value just below x.
#   Indeed, since the distribution is discrete, P(X >x-1) = P(X ≥x).
#   Source: <http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html>

# Networks are generally plotted using a force-directed algorithm, which is a class of algorithms
# for drawing graphs. The default function in igraph is called the Fruchterman-Reingold (1991) algorithm.
# Typically, force-directed algorithms use a physical simulation where some kind of attractive force
# (imagine a spring) are used to attract nodes connected by edges together. So ‘tightly’ connected clusters
# of nodes will show up close to each other, and those that are ‘loosely’ connected will be repulsed towards
# the outside. However, the algorithm does not specify where any node has to be other than these constraints.
# (source: https://dshizuka.github.io/networkanalysis/03_plots.html)
# Note: ggnetwork use same default as igraph.


# Hexadecimal color code for transparency --> https://gist.github.com/lopspower/03fb1cc0ac9f32ef38f4
# ggnetwork tutorial --> https://cran.r-project.org/web/packages/ggnetwork/vignettes/ggnetwork.html


#### debug ----

# gwas_filename <- "../genomic_prediction/hybrids/analysis/gwas/summary_gwas.txt"
# meff_filename <- "analysis/marker_effects/YLD/marker_effects.rrblup.txt"
# modules_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/kDiff_per_module.txt"
# mod_pheno_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/module-pheno_pvals.txt"
# mod_ld_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/modules_ld"
# output_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50"
# type_connect_str <- "adjacency"
# soft_threshold <- 12
# edge_threshold <- 0.25