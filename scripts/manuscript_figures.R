library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(WGCNA)
library(network)
library(sna)
library(ggnetwork)
library(factoextra)
library(ggh4x)
library(ggrepel)
library(UpSetR)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)
dir.create("figures")



#### figure 1 ----

##### a #####

input_folder <- "analysis/marker_effects/YLD"
blues_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
output_folder <- "figures"

# create empty data frames to store results across models
gebvs_all <- data.frame(stringsAsFactors = FALSE)

for (model in c("rrblup", "gwas")) {
  
  # get filename
  file_gebvs <- paste0(input_folder, "/GEBVs.", model, ".txt")
  # load file
  gebvs <- fread(file_gebvs, header = TRUE, data.table = FALSE)
  # transform to long format
  gebvs <- pivot_longer(gebvs, -genotype, names_to = "env", values_to = "gebv")
  # append results to main df
  gebvs_all <- rbind(gebvs_all, data.frame(gebvs, model = model))
  
}

# read file with genotype means for each environment
means <- fread(blues_file, header = TRUE, data.table = FALSE)
means <- data.frame(means, model = "real_pheno")
colnames(means) <- colnames(gebvs_all)
# correct env names for compatibility
means$env <- gsub("-", ".", means$env)
# append means to gebv df
gebvs_all <- rbind(gebvs_all, means)

# compare distribution of gebvs with real data
plot_boxplot_gebvs <- ggplot(gebvs_all) +
  geom_boxplot(aes(x = env, y = gebv, fill = model)) +
  labs(x = "Environments", y = "Yield", fill = "Data type") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    labels = c("RR-BLUP", "GWAS", "Observed")) +
  theme_bw()
ggsave(filename = paste0(output_folder, "/fig_1a.png"), plot = plot_boxplot_gebvs,
       device = "png", dpi = 300, units = "in", width = 8, height = 4)

# clean R environment
rm(list = ls())

##### b #####

# output from 'scripts/pick_soft_threshold.R'
# network selected:
# - meff model: rrblup
# - norm method: zscore

##### c #####

# output from 'scripts/build_meff_network.R'
# network selected:
# - meff model: rrblup
# - norm method: zscore



#### figure 2 ----

##### a #####

# output from 'scripts/define_network_modules.R'
# network selected:
# - meff model: rrblup
# - norm method: zscore
# - min mod size: 50
# - panStage: off

##### b #####

input_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off"
module_stats_file <- paste0(input_folder, "/kDiff_per_module.txt")
output_folder <- "figures"
type_connections <- "TOM"

# read module stats
module_stats <- fread(module_stats_file, header = TRUE, data.table = FALSE)
module_stats$module <- factor(module_stats$module)
# format table
module_stats <- module_stats %>% 
  filter(source == type_connections) %>% 
  select(marker, module, kDiff, kRatio, clusterCoeff) %>% 
  pivot_longer(-c(marker, module), names_to = "stat", values_to = "value") %>% 
  mutate(stat = factor(stat, levels = c("kDiff", "kRatio", "clusterCoeff")))

# plot kDiff among markers
# (i.e. difference between number of connections within vs outside modules)
plot_mod_stats <- ggplot(module_stats) +
  facet_wrap(~ stat, nrow = 3, scales = "free_y") +
  geom_boxplot(aes(x = module, y = value, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(module_stats$module)) +
  labs(x = "Modules", y = "Values") +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0(output_folder, "/fig_2b.png"), plot = plot_mod_stats,
       device = "png", dpi = 300, units = "in", width = 8, height = 6)

# clean R environment
rm(list = ls())



#### figure 3 ----

input_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off"
output_folder <- "figures"
meff_mod <- paste0(input_folder, "/define_network_modules.RData")

# load data
load(meff_mod)

# plot heatmap
png(paste0(output_folder, "/fig_3.png"), units = "in", res = 300, width = 8, height = 8)
pheatmap(mat = t(MEs), angle_col = 0, main = "Modules eigenvalues per environment")
dev.off()

# clean R environment
rm(list = ls())



#### figure 4 ----

##### a #####

input_folder <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off"
meff_mod_file <- paste0(input_folder, "/define_network_modules.RData")
modules_filename <- paste0(input_folder, "/kDiff_per_module.txt")
mod_ld_folder <- paste0(input_folder, "/modules_ld")
output_folder <- "figures"
type_connect_str <- "TOM"
soft_threshold <- 24
edge_threshold <- 0.01
mod <- "darkolivegreen"
ld_module <- paste0(mod_ld_folder, "/ld_markers_", mod, ".ld.gz")

# load data
load(meff_mod_file)
rm(geneTree, marker_info, moduleColors, moduleLabels)
modules <- fread(modules_filename, header = TRUE, data.table = FALSE)
# keep either adjacency or TOM values
modules <- subset(modules, source == type_connect_str)
# get order of hub markers in each module
modules <- group_by(modules, module) %>%
  mutate(hub_pos = rank(desc(kDiff)),
         n_mod = n()) %>%
  arrange(hub_pos, .by_group = TRUE)

# subset data
mod_subset <- subset(modules, module == mod)
# select markers in that module
mod_markers <- pull(modules[modules$module == mod & modules$source == type_connect_str, "marker"])

# get effects for markers in the module
mod_meff <- marker_effects[, mod_markers]
# get eigenvalues for that module
mod_ME <- MEs[, paste0("ME", mod), drop = FALSE]
mod_ME <- rownames_to_column(mod_ME, var = "env")
colnames(mod_ME)[2] <- "ME"

# plot module meffs and eigenvalues across environments
meff_plot <- data.frame(t(mod_meff)) %>% 
  rownames_to_column(var = "marker") %>% 
  pivot_longer(-marker, names_to = "env", values_to = "effects") %>% 
  ggplot() +
  geom_line(aes(x = env, y = effects, group = marker), alpha = 0.2, color = "grey50") +
  geom_line(data = mod_ME, aes(x = env, y = ME, group = 1), size = 2, color = mod) +
  labs(title = bquote("Effects of markers in the" ~ bold(.(mod)) ~ "module across envs"),
       x = "Environments", y = "Marker effects") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = paste0(output_folder, "/fig_4a.png"), plot = meff_plot,
       device = "png", dpi = 300, units = "in", width = 7, height = 4)

##### b #####

# load ld file
ld_module <- fread(ld_module, header = TRUE, data.table = FALSE)

# plot distribution r2 for all markers
ld_dist_plot <- ggplot(ld_module, aes(x = R2)) +
  geom_histogram(data = subset(ld_module, R2 >= 0.9), fill = "firebrick") +
  geom_histogram(data = subset(ld_module, R2 < 0.9), fill = "grey50") +
  labs(title = paste0("LD distribution of ", mod, " module"),
       x = bquote(R^2), y = "Count") +
  theme_bw()
ggsave(filename = paste0(output_folder, "/fig_4b.png"), plot = ld_dist_plot,
       device = "png", dpi = 300, units = "in", width = 4, height = 4)

##### c #####

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
# remove edges with no LD info
mod_edges <- mod_edges[!is.na(mod_edges$R2), ]

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

# add vertex layout to network
mod_network <- ggnetwork(mod_network, layout = "fruchtermanreingold", weights = "weight", cell.jitter = 0.5)
# View(mod_network)

# visualize network by LD status
plot_net_ld <- ggplot() +
  geom_edges(data = subset(mod_network, LD == FALSE),
             aes(x = x, y = y, xend = xend, yend = yend), color = "grey50", alpha = 0.2) +
  geom_edges(data = subset(mod_network, LD == TRUE),
             aes(x = x, y = y, xend = xend, yend = yend), color = "firebrick") +
  geom_nodes(data = mod_network,
             aes(x = x, y = y), color = "black") +
  theme_blank()
ggsave(filename = paste0(output_folder, "/fig_4c.png"), plot = plot_net_ld,
       device = "png", dpi = 300, units = "in", width = 6, height = 6)

# clean R environment
rm(list = ls())



#### figure 5 ----

env_idx_file <- "data/env_covariables/env_covariables_means_per_intervals.txt"
output_folder <- "figures"
example_idx <- "T2M"

##### a #####

# load and format data
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)
env_idx <- env_idx %>% 
  filter(covariable == example_idx) %>% 
  select(-covariable)

# plot env idx across intervals
env_idx_plot <- ggplot(env_idx, aes(x = intervals, y = value, color = env)) +
  geom_line(aes(group = env), show.legend = FALSE) +
  facet_wrap(~ env, nrow = 3, ncol = 3) +
  labs(x = "Time intervals", y = "Temperature (centered and scaled)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave(filename = paste0(output_folder, "/fig_5a.png"), plot = env_idx_plot,
       device = "png", dpi = 300, units = "in", width = 6, height = 6)

##### b #####

# load and format data
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)
env_idx <- data.frame(env = env_idx$env, idx = paste(env_idx$covariable, env_idx$intervals, sep = "_"), value = env_idx$value)
env_idx <- pivot_wider(env_idx, names_from = "idx", values_from = "value")
env_idx <- column_to_rownames(env_idx, var = "env")

# correlate variables
cor_env_idx <- cor(as.matrix(env_idx))

# plot heatmap
png(paste0(output_folder, "/fig_5b.png"), units = "in", res = 300, width = 6, height = 4)
pheatmap(mat = cor_env_idx, angle_col = 0, main = "Correlation among environmental parameters",
         show_colnames = FALSE, show_rownames = FALSE)
dev.off()

##### c #####

# load data
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)
env_idx <- env_idx %>% 
  unite(covariable:intervals, col = "covariable", sep = "_") %>% 
  pivot_wider(names_from = "covariable", values_from = "value") %>% 
  column_to_rownames(var = "env") %>% 
  as.matrix()

# perform PCA
pca <- prcomp(env_idx, scale = TRUE)
pca_results <- pca$x[match(rownames(env_idx), rownames(pca$x)), ]
pca_results <- rownames_to_column(data.frame(pca_results), var = "env")
pca_results <- pivot_longer(pca_results, -env, names_to = "covariable", values_to = "value")

# plot PCA
png(paste0(output_folder, "/fig_5c.png"), units = "in", res = 300, width = 5, height = 3)
fviz_pca_ind(pca, col.ind = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  scale_x_continuous(breaks = seq(-30, 30, 15)) +
  scale_y_continuous(breaks = seq(-20, 20, 10)) +
  scale_color_viridis_c(option = "C", breaks = seq(0, 30, 5), end = 0.8) +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-20, 20)) +
  labs(title = "PCA of environmental parameters", x = "PC1", y = "PC2", color = "Contributions") +
  theme_bw() +
  theme(legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(size = 8), #change legend title font size
        legend.text = element_text(size = 8)) #change legend text font size
dev.off()

# clean R environment
rm(list = ls())



#### figure 6 ----

model <- "rrblup"
norm <- "zscore"
size <- "50"
pam <- "off"
meff_mod_Rdata <- paste0("analysis/networks/YLD/meff_", model, "/norm_", norm, "/min_mod_size_", size,
                         "/pamStage_", pam, "/define_network_modules.RData")
mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
env_idx_file <- "data/env_covariables/pca_env_idx.txt"
output_folder <- "figures"

##### a ######

# load data
mod_env_idx_results <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)

# select network of interest
mod_env_idx_results <- subset(mod_env_idx_results, meff_model == model & norm_method == norm &
                              minsize == size & pamStage == pam)
rownames(mod_env_idx_results) <- 1:nrow(mod_env_idx_results)

# adjust factor levels
mod_env_idx_results$module <- gsub("^ME", "", mod_env_idx_results$module, perl = TRUE)
mod_env_idx_results$module <- factor(mod_env_idx_results$module,
                                     levels = unique(mod_env_idx_results$module))

# plot association
plot_cor_mod_pc <- ggplot(mod_env_idx_results, aes(x = pval, y = cor, fill = module, label = paste0(env_idx, "\n", module))) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_vline(xintercept = 0.5, color = "gray50") +
  geom_point(data = head(mod_env_idx_results, n = 3), shape = 21, fill = levels(mod_env_idx_results$module)[1:3], show.legend = FALSE) +
  geom_point(data = tail(mod_env_idx_results, n = nrow(mod_env_idx_results) - 3), color = "grey60", show.legend = FALSE) +
  geom_label_repel(data = head(mod_env_idx_results, n = 3), max.overlaps = 50,
                   fill = "white", size = 2, show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) +
  labs(x = "p-values", y = "Correlation") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(filename = paste0(output_folder, "/fig_6a.png"), plot = plot_cor_mod_pc,
       device = "png", dpi = 300, units = "in", width = 4, height = 4)


##### b #####

# load module-env index relationship
mod_env_idx_cor <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)
# select network of interest
mod_env_idx_cor <- subset(mod_env_idx_cor, meff_model == model & norm_method == norm &
                          minsize == size & pamStage == pam)
rownames(mod_env_idx_cor) <- 1:nrow(mod_env_idx_cor)
# select top 3 most significant modules
mod_env_idx_cor <- head(mod_env_idx_cor[order(mod_env_idx_cor$pval), ], n = 3)
# round numbers
mod_env_idx_cor$cor <- sapply(mod_env_idx_cor$cor, function(x) round(x, digits = 2))
mod_env_idx_cor$pval <- sapply(mod_env_idx_cor$pval, function(x) round(x, digits = 2))

# load env idx results
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)

for (row in 1:nrow(mod_env_idx_cor)) {
  
  # load module data
  load(meff_mod_Rdata)
  
  # get module eigenvalues
  MEs <- rownames_to_column(MEs, var = "env") %>% 
    select(env, mod_env_idx_cor[row, "module"]) %>% 
    pivot_longer(-env, names_to = "type", values_to = "value") %>% 
    mutate(env = gsub(".", "-", env, fixed = TRUE),
           type = gsub("^ME", "", type, perl = TRUE))
  
  # get environments to keep
  envs_to_keep <- unique(MEs$env)
  
  # get env idx results
  env_idx_net <- subset(env_idx, env %in% envs_to_keep & covariable == mod_env_idx_cor[row, "env_idx"]) %>% 
    rename(type = "covariable")
  
  # # plot ME and env idx in their own scales
  # plot_mod_env_idx <- rbind(MEs, env_idx_net) %>%
  #   ggplot(aes(x = env, y = value, color = type)) +
  #   geom_line(aes(group = type), show.legend = FALSE) +
  #   facet_wrap(~ type, nrow = 2, scales = "free_y") +
  #   labs(title = bquote("Module-Environmental index correlation =" ~ bold(.(as.character(mod_env_idx_cor[row, "cor"])))
  #                       ~ .(paste0("(pval: ", mod_env_idx_cor[row, "pval"], ")"))),
  #        subtitle = bquote("meff model:" ~ bold(.(mod_env_idx_cor[row, "meff_model"])) ~ "/"
  #                          ~ "normalization:" ~ bold(.(mod_env_idx_cor[row, "norm_method"])) ~ "/"
  #                          ~ "min module size:" ~ bold(.(as.character(mod_env_idx_cor[row, "minsize"]))) ~ "/"
  #                          ~ "pamStage:" ~ bold(.(mod_env_idx_cor[row, "pamStage"])))) +
  #   theme_bw() +
  #   theme(panel.grid.minor = element_blank())
  
  # normalize to same scale
  MEs$value <- sapply(MEs$value, function(x) (x - min(MEs$value))/(max(MEs$value) - min(MEs$value)))
  env_idx_net$value <- sapply(env_idx_net$value, function(x) (x - min(env_idx_net$value))/(max(env_idx_net$value) - min(env_idx_net$value)))
  
  # plot them together
  plot_mod_env_idx_norm <- rbind(MEs, env_idx_net) %>%  
    mutate(type = factor(type, levels = unique(type))) %>% 
    ggplot(aes(x = env, y = value, color = type)) +
    geom_line(aes(group = type), size = 1.5) +
    scale_color_manual(values = c(gsub("^ME", "",  mod_env_idx_cor[row, "module"], perl = TRUE), "grey70")) + 
    labs(title = bquote("Module-Environmental index correlation =" ~ bold(.(as.character(mod_env_idx_cor[row, "cor"])))
                        ~ .(paste0("(pval: ", mod_env_idx_cor[row, "pval"], ")"))),
         subtitle = bquote("meff model:" ~ bold(.(mod_env_idx_cor[row, "meff_model"])) ~ "/"
                           ~ "normalization:" ~ bold(.(mod_env_idx_cor[row, "norm_method"])) ~ "/"
                           ~ "min module size:" ~ bold(.(as.character(mod_env_idx_cor[row, "minsize"]))) ~ "/"
                           ~ "pamStage:" ~ bold(.(mod_env_idx_cor[row, "pamStage"]))),
         x = "Environments", y = "Normalized value", color = "Eigenvalues") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  ggsave(filename = paste0(output_folder, "/fig_6b", row, ".png"), plot = plot_mod_env_idx_norm,
         device = "png", dpi = 300, units = "in", width = 8, height = 3)
  
}


##### c #####

# output from 'plot_pc_contributions.R'

# clean R environment
rm(list = ls())



#### figure 7 ----

mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
env_idx_file <- "data/env_covariables/pca_env_idx.txt"
output_folder <- "figures"
p_value <- 0.1
idx <- "PC5"
list_markers_for_ld <- paste0("analysis/networks/YLD/overlap_markers/markers_correlated_", idx, ".txt")

# load module-env index relationship
mod_env_idx_cor <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)
# select only significant associations
mod_env_idx_cor <- subset(mod_env_idx_cor, pval < p_value)
rownames(mod_env_idx_cor) <- 1:nrow(mod_env_idx_cor)
# round numbers
mod_env_idx_cor$cor <- sapply(mod_env_idx_cor$cor, function(x) round(x, digits = 2))
mod_env_idx_cor$pval <- sapply(mod_env_idx_cor$pval, function(x) round(x, digits = 2))
# correct module name
mod_env_idx_cor$module <- gsub("^ME", "", mod_env_idx_cor$module, perl = TRUE)

# split data to loop through groups
mod_env_idx_cor_split <- split(mod_env_idx_cor, mod_env_idx_cor[, "env_idx"])

# get only modules associated with that env_idx
sig_env_idx <- mod_env_idx_cor_split[[idx]]

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
  network_mod_markers <- paste0("analysis/networks/YLD/meff_", model, "/norm_", norm,
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
names(list_markers_mod_idx) <- sapply(names(list_markers_mod_idx), function(x) gsub("-", "_", x))

# read LD file
ld_file <- fread(gsub(".txt", ".ld.gz", list_markers_for_ld), header = TRUE, data.table = FALSE)

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

# rename marker sets
sig_env_idx$name <- pull(unite(sig_env_idx[, c("meff_model", "norm_method", "minsize", "pamStage", "module")], col = "name", sep = "_"))
sig_env_idx <- sig_env_idx[match(names(list_markers_mod_idx), sig_env_idx$name), ]
sig_env_idx$alias <- NA
for (row in 1:nrow(sig_env_idx)) {
  sig_env_idx[row, "alias"] <- paste0("network ", row, "\n(cor: ", sig_env_idx[row, "cor"],
                                      ", p-value: ", sig_env_idx[row, "pval"], ")")
}
names(list_markers_mod_idx) <- sig_env_idx[, "alias"]

# plot intersections of markers
upset_plot <- upset(fromList(list_markers_mod_idx), sets = rev(sig_env_idx[, "alias"]), keep.order = TRUE,
                    order.by = "freq", mb.ratio = c(0.55, 0.45), nsets = 100)
png(paste0(output_folder, "/fig_7.png"), units = "in", res = 300, width = 5, height = 3)
print(upset_plot)
dev.off()

# clean R environment
rm(list = ls())



#### supplemental figure 1 ----

# output from 'scripts/pick_soft_threshold.R'
