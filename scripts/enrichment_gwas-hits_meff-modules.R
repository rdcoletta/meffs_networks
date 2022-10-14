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
description: check if GWAS hits are present in marker effect modules.

usage: Rscript enrichment_gwas-hits_meff-modules.R [meff_filename] [modules_filename] [mod_pheno_filename]
                                                   [mod_ld_folder] [output_folder] [...]

positional arguments:
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
  --edge-threshold=[VALUE]      adjacency threshold for including edges in the network plots (default: 0.1)


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
meff_filename <- args[1]
modules_filename <- args[2]
mod_pheno_filename <- args[3]
mod_ld_folder <- args[4]
output_folder <- args[5]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default of optional args
type_connect_str <- "TOM"
soft_threshold <- "10"
edge_threshold <- "0.1"

# assert to have the correct optional arguments
pos_args <- 5
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
meff <- fread(meff_filename, header = TRUE, data.table = FALSE)
modules <- fread(modules_filename, header = TRUE, data.table = FALSE)
mod_pheno <- fread(mod_pheno_filename, header = TRUE, data.table = FALSE)

# keep either adjacency or TOM values
modules <- subset(modules, source == type_connect_str)

# get order of hub markers in each module
modules <- group_by(modules, module) %>%
  mutate(hub_pos = rank(desc(kDiff)),
         n_mod = n()) %>%
  arrange(hub_pos, .by_group = TRUE)

# add p-values to modules
mod_pheno$module <- gsub("^ME", "", mod_pheno$module, perl = TRUE)
modules <- merge(x = modules, y = mod_pheno, by = "module")
rm(mod_pheno)

# how many of all gwas hits are in a significant module
gwas_hits <- unique(modules[modules$gwas_status == "gwas_hit", "marker"])
cat("gwas hit in a significant module: ",
    round(sum(modules$gwas_status == "gwas_hit" & modules$pval < 0.05) * 100 / length(gwas_hits), 2), "%",
    sep = "")

# how many of all gwas hits are in a significant module
gwas_top_non_sigs <- unique(modules[modules$gwas_status == "gwas_top_non_sig", "marker"])
cat("gwas top non-significant hit in a significant module: ",
    round(sum(modules$gwas_status == "gwas_top_non_sig" & modules$pval < 0.05) * 100 / length(gwas_top_non_sigs), 2), "%",
    sep = "")

# create empty df to store stats for each module of network
net_stats <- data.frame(stringsAsFactors = FALSE)

cat("analyzing signficant modules...\n")
# test for enrichment of gwas hits in modules associated with the trait
for (mod in unique(modules$module)) {

  cat(paste0("  ", mod, " (", which(unique(modules$module) == mod), "/", length(unique(modules$module)), ")\n"))

  # subset data
  mod_subset <- subset(modules, module == mod)

  # network numbers
  n_markers_net <- nrow(modules)
  n_gwas_hits_net <- length(gwas_hits)
  n_gwas_top_ns_net <- length(gwas_top_non_sigs)

  # module numbers
  n_markers_mod <- nrow(mod_subset)
  n_gwas_hits_mod <- sum(mod_subset$gwas_status == "gwas_hit")
  n_gwas_top_ns_mod <- sum(mod_subset$gwas_status == "gwas_top_non_sig")
  pval_trait_mod_cor <- unique(mod_subset$pval)
  sig_module <- pval_trait_mod_cor < 0.05

  # run hypergeometric test on gwas hits
  pval_gwas_hyper <- phyper(q = n_gwas_hits_mod - 1,
                            m = n_gwas_hits_net,
                            n = n_markers_net - n_gwas_hits_net,
                            k = n_markers_mod,
                            lower.tail = FALSE, log.p = FALSE)
  # run hypergeometric test on gwas top non-significant hits
  pval_gwas_top_ns_hyper <- phyper(q = n_gwas_top_ns_mod - 1,
                                   m = n_gwas_top_ns_net,
                                   n = n_markers_net - n_gwas_top_ns_net,
                                   k = n_markers_mod,
                                   lower.tail = FALSE, log.p = FALSE)

  # check effects of gwas hits across envs
  meff_gwas_plot <- filter(meff, marker %in% gwas_hits) %>%
    mutate(in_module = factor(if_else(marker %in% mod_subset$marker, true = TRUE, false = FALSE),
                              levels = c(FALSE, TRUE))) %>%
    pivot_longer(-c(marker, in_module), names_to = "env", values_to = "effects") %>%
    ggplot() +
    geom_line(aes(x = env, y = effects, group = marker, color = in_module)) +
    scale_color_manual(values = c("#8080804D", "firebrick")) +
    labs(title = bquote("Effects of GWAS hits in" ~ bold(.(mod)) ~ "module across envs"),
         subtitle = paste0("markers - ", n_markers_mod, "/", n_markers_net, "\n",
                           "gwas hits - ", n_gwas_hits_mod, "/", n_gwas_hits_net, " (",
                           "pval hyper: ", round(pval_gwas_hyper, 3), ")\n",
                           "gwas top non-sig hits - ", n_gwas_top_ns_mod, "/", n_gwas_top_ns_net, " (",
                           "pval hyper: ", round(pval_gwas_top_ns_hyper, 3), ")\n")) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.subtitle = element_text(color = "gray40"))

  # print(meff_gwas_plot)
  ggsave(filename = paste0(output_folder, "/meff_gwas.", mod, ".pdf"),
         plot = meff_gwas_plot, device = "pdf", width = 10, height = 6)


  # load file
  ld_module <- paste0(mod_ld_folder, "/ld_markers_", mod, ".ld.gz")
  # use try() function to not break code if module doesn't have an LD file
  # (grey modules may not have enough markers -- not subjected to minimum module size)
  ld_module <- try(fread(ld_module, header = TRUE, data.table = FALSE))

  if (class(ld_module) != "try-error") {

    # get ld stats
    ld_stats <- data.frame(t(as.numeric(summary(ld_module$R2)[2:5])), stringsAsFactors = FALSE)
    colnames(ld_stats) <- c("ld_Q1", "ld_median", "ld_mean", "ld_Q3")

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
    # add whether marker is gwas hit, top non-sig hit, or not a gwas hit
    node_gwas_status <- mod_subset[mod_subset$marker %in% node_names, c("marker", "gwas_status")]
    node_gwas_status <- node_gwas_status[match(node_names, node_gwas_status$marker), ]
    set.vertex.attribute(mod_network,
                         attrname = "gwas_status",
                         value = node_gwas_status$gwas_status)

    # add vertex layout to network
    mod_network <- ggnetwork(mod_network, layout = "fruchtermanreingold", weights = "weight", cell.jitter = 0.5)
    # View(mod_network)

    # visualize network by hub markers
    mod_network$gwas_status <- factor(mod_network$gwas_status,
                                      levels = c("gwas_hit", "gwas_top_non_sig", "not_gwas_hit"))
    plot_net_hub <- ggplot(mod_network, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "#8080800D") +
      geom_nodes(aes(size = kWithin, color = gwas_status)) +
      scale_color_manual(values = c("firebrick", "darkblue", "black"), drop = FALSE) +
      theme_blank() +
      labs(title = paste0(mod, " module"),
           subtitle = paste0("markers - ", nrow(mod_subset), "/", nrow(modules), "\n",
                             "gwas hits - ", sum(mod_subset$gwas_status == "gwas_hit"), "/", length(gwas_hits), " (",
                             "pval hyper: ", round(pval_gwas_hyper, 3), ")\n",
                             "gwas top non-sig hits - ", sum(mod_subset$gwas_status == "gwas_top_non_sig"), "/", length(gwas_top_non_sigs), " (",
                             "pval hyper: ", round(pval_gwas_top_ns_hyper, 3), ")\n")) +
      theme(plot.subtitle = element_text(color = "gray40"))

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
                             "gwas hits: ", sum(mod_subset$gwas_status == "gwas_hit"), "/", length(gwas_hits), " - ",
                             "pval hyper: ", round(pval_gwas_hyper, 3), ")"))

    ggsave(filename = paste0(output_folder, "/mod_viz.", mod, ".by-ld.pdf"),
           plot = plot_net_ld, device = "pdf", width = 10, height = 12)

  } else {

    # set ld stats to NA
    ld_stats <- data.frame(ld_Q1 = NA, ld_median = NA, ld_mean = NA, ld_Q3 = NA,
                           stringsAsFactors = FALSE)

  }

  # add module numbers to final df
  mod_stats <- data.frame(mod = mod,
                          n_markers_mod = n_markers_mod,
                          pval_trait_mod_cor = pval_trait_mod_cor,
                          sig_module = sig_module,
                          n_gwas_hits_mod = n_gwas_hits_mod,
                          n_gwas_top_ns_mod = n_gwas_top_ns_mod,
                          pval_gwas_hyper = pval_gwas_hyper,
                          pval_gwas_top_ns_hyper = pval_gwas_top_ns_hyper,
                          ld_stats,
                          stringsAsFactors = FALSE)

  net_stats <- rbind(net_stats, mod_stats)

}

# write file
fwrite(net_stats, file = paste0(output_folder, "/gwas_enrichment_per_module.txt"))
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

# meff_filename <- "analysis/marker_effects/YLD/marker_effects.rrblup.txt"
# modules_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/kDiff_per_module.gwas-status.txt"
# mod_pheno_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/module-pheno_pvals.txt"
# mod_ld_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/modules_ld"
# output_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off"
# type_connect_str <- "TOM"
# soft_threshold <- 20
# edge_threshold <- 0.1
