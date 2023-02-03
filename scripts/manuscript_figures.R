library(data.table)
library(readxl)
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
library(ggridges)
library(viridis)
library(UpSetR)
enableWGCNAThreads()
allowWGCNAThreads()
options(stringsAsFactors = FALSE)
dir.create("figures")



#### figure 1 ----

##### a #####

input_folder <- "analysis/marker_effects/YLD"
blues_file <- "data/1stStage_BLUEs.YLD-per-env.txt"
env_idx_file <- "data/env_covariables/env_covariables_means_per_intervals.txt"
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
gebvs_all$env <- gsub(".", "-", gebvs_all$env, fixed = TRUE)
# append means to gebv df
gebvs_all <- rbind(gebvs_all, means)
# adjust factor levels for plotting
gebvs_all$model <- factor(gebvs_all$model, levels = c("rrblup", "gwas", "real_pheno"))

# compare distribution of gebvs with real data
plot_boxplot_gebvs <- ggplot(gebvs_all) +
  geom_boxplot(aes(x = env, y = gebv, fill = model)) +
  labs(x = "Environments", y = "Grain yield", fill = "Data type") +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    labels = c("RR-BLUP", "GWAS", "Observed")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/fig_1a.png"), plot = plot_boxplot_gebvs,
       device = "png", dpi = 350, units = "in", width = 8, height = 2)


##### b #####

# load and format data
env_idx <- fread(env_idx_file, header = TRUE, data.table = FALSE)
env_idx <- data.frame(env = env_idx$env, idx = paste(env_idx$covariable, env_idx$intervals, sep = "_"), value = env_idx$value)
env_idx <- pivot_wider(env_idx, names_from = "idx", values_from = "value")
env_idx <- column_to_rownames(env_idx, var = "env")

# correlate variables
cor_env_idx <- cor(as.matrix(env_idx))

# plot heatmap
png(paste0(output_folder, "/fig_1b.png"), units = "in", res = 350, width = 4, height = 3)
pheatmap(mat = cor_env_idx, angle_col = 0, show_colnames = FALSE, show_rownames = FALSE, legend_breaks = seq(-1, 1, 0.5), legend_labels = c(-1, -0.5, 0, 0.5, 1))
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
eigs <- pca$sdev^2

# plot PCA
png(paste0(output_folder, "/fig_1c.png"), units = "in", res = 350, width = 4, height = 3)
fviz_pca_ind(pca, col.ind = "grey30", repel = TRUE, labelsize = 3) +
  scale_x_continuous(breaks = seq(-30, 30, 15)) +
  scale_y_continuous(breaks = seq(-20, 20, 10)) +
  scale_color_viridis_c(option = "C", breaks = seq(0, 30, 5), end = 0.8) +
  coord_cartesian(xlim = c(-30, 30), ylim = c(-20, 20)) +
  labs(x = paste0("PC1 (", round(eigs[1] / sum(eigs), digits = 3) * 100, "%)"),
       y = paste0("PC2 (", round(eigs[2] / sum(eigs), digits = 3) * 100, "%)"),
       color = "Contributions") +
  theme_bw() +
  theme(title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.key.size = unit(0.2, 'in'),
        legend.title = element_text(size = 8), #change legend title font size
        legend.text = element_text(size = 8)) #change legend text font size
dev.off()

# clean R environment
rm(list = ls())



#### figure 2 ----

results_folder <- "analysis/networks/YLD"
qc_networks_file <- paste0(results_folder, "/summary_qc_networks.txt")
net_alias_file <- "supp_materials/supp_table6.xlsx"
output_folder <- "figures"

##### a #####

# create empty df to store results
sft_all <- data.frame(stringsAsFactors = FALSE)

for (model in c("rrblup", "gwas")) {
  for (norm in c("minmax", "zscore")) {

    if (norm == "minmax") {
      for (cv in c(0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25)) {

        # load file
        sft_file <- paste0(results_folder, "/meff_", model, "/norm_", norm,
                           "/cv_thresholds/scale-free_topology_fit_index.cv-", cv, ".txt")
        sft_results <- fread(sft_file, header = TRUE, data.table = FALSE)
        # get results
        sft_all <- rbind(sft_all, data.frame(model = model, norm = norm,  cv = cv, sft_results))

      }
    }

    if (any(norm == "zscore" | norm == "none")) {
      for (cv in c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1)) {

        # load file
        sft_file <- paste0(results_folder, "/meff_", model, "/norm_", norm,
                           "/cv_thresholds/scale-free_topology_fit_index.cv-", cv, ".txt")
        sft_results <- fread(sft_file, header = TRUE, data.table = FALSE)
        # get results
        sft_all <- rbind(sft_all, data.frame(model = model, norm = norm,  cv = cv, sft_results))

      }
    }

  }
}
rm(sft_results)

# add signed R2 values
sft_all$signed_R2 <- -sign(sft_all$slope) * sft_all$SFT.R.sq

# plot power distribution
sft_plot <- sft_all %>%
  mutate(norm = factor(norm, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         model = factor(model, levels = c("gwas", "rrblup"), labels = c("GWAS", "RR-BLUP"))) %>%
  ggplot(aes(x = Power, y = signed_R2, group = as.factor(cv), color = as.factor(cv))) +
  facet_grid(norm ~ model) +
  geom_line() +
  geom_point(size = 1.5) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_color_viridis_d(end = 0.9, option = "B", direction = -1) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  scale_y_continuous(breaks = seq(-1, 1, 1)) +
  labs(x = "Power", y = bquote("Scale free topology model fit, signed" ~ R^2), color = "CV cut-off") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.direction = "vertical") +
guides(color = guide_legend(ncol = 2))
ggsave(filename = paste0(output_folder, "/fig_2a.png"), plot = sft_plot,
       device = "png", dpi = 350, units = "in", width = 8, height = 3)


##### b #####

plot_total_markers <- sft_all %>%
  mutate(norm = factor(norm, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         model = factor(model, levels = c("gwas", "rrblup"), labels = c("GWAS", "RR-BLUP"))) %>%
  group_by(model, norm, cv) %>%
  summarize(total.markers = mean(total.markers, na.rm = TRUE)) %>%
  ggplot(aes(x = as.factor(cv), y = total.markers)) +
  facet_grid2(norm ~ model, scales = "free_x", independent = "x") +
  geom_col() +
  # scale_y_continuous(breaks = seq(0, 10000, 1000)) +
  labs(x = "CV cut-off", y = "Total markers in network") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))
ggsave(filename = paste0(output_folder, "/fig_2b.png"), plot = plot_total_markers,
       device = "png", dpi = 350, units = "in", width = 4, height = 3)

##### c #####

# get network aliases
net_alias <- read_excel(net_alias_file)
net_alias$order <- paste(net_alias$Model, net_alias$Normalization, net_alias$`Minimum module size`, net_alias$pamStage, sep = "_")

# load qc file
qc_results <- fread(qc_networks_file, header = TRUE, data.table = FALSE)
# adjust data to match aliases order
qc_results <- qc_results %>%
  filter(source == "adjacency") %>%
  mutate(norm_method = factor(norm_method, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         meff_model = factor(meff_model, levels = c("gwas", "rrblup"), labels = c("GWAS", "RR-BLUP"))) %>%
  unite(col = "order", meff_model:pamStage, remove = FALSE) %>%
  mutate(order = match(order, net_alias$order)) %>%
  arrange(order)

# plot markers per module
plot_markers_per_mod <- qc_results %>%
  group_by(order, meff_model, norm_method, minsize, pamStage, module) %>%
  summarize(n_markers = n()) %>%
  ungroup() %>%
  group_by(order, meff_model, norm_method, minsize, pamStage) %>%
  mutate(n_modules = n()) %>%
  ungroup() %>%
  mutate(n_markers_no_grey = NA,
         n_markers_no_grey = case_when(module != "grey" ~ n_markers)) %>%
  pivot_longer(c(n_markers, n_markers_no_grey, n_modules), names_to = "metric", values_to = "values") %>%
  mutate(order = factor(order),
         metric = factor(metric, levels = c("n_modules", "n_markers", "n_markers_no_grey"),
                         labels = c("Total modules", "Unassigned markers", "Markers per module"))) %>%
  ggplot() +
  facet_wrap(~ metric, nrow = 3, scales = "free_y") +
  geom_col(data = function(x) subset(x, metric == "Total modules"),
           aes(x = order, y = values), fill = "grey40", show.legend = FALSE) +
  geom_col(data = function(x) subset(x, metric == "Unassigned markers"  & module == "grey"),
           aes(x = order, y = values), fill = "grey", show.legend = FALSE) +
  geom_boxplot(data = function(x) subset(x, metric == "Markers per module"),
               aes(x = order, y = values), fill = "grey60", outlier.size = 0.25, show.legend = FALSE) +
  labs(x = "Networks", y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 8))
ggsave(filename = paste0(output_folder, "/fig_2c.png"), plot = plot_markers_per_mod,
       device = "png", dpi = 350, units = "in", width = 4, height = 3)

# clean R environment
rm(list = ls())



#### figure 3 ----

results_folder <- "analysis/networks/YLD"
qc_networks_file <- paste0(results_folder, "/summary_qc_networks.txt")
net_alias_file <- "supp_materials/supp_table6.xlsx"
output_folder <- "figures"

# get network aliases
net_alias <- read_excel(net_alias_file)
net_alias$order <- paste(net_alias$Model, net_alias$Normalization, net_alias$`Minimum module size`, net_alias$pamStage, sep = "_")

# load qc file
qc_results <- fread(qc_networks_file, header = TRUE, data.table = FALSE)
# adjust data to match aliases order
qc_results <- qc_results %>%
  filter(source == "adjacency") %>%
  mutate(norm_method = factor(norm_method, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         meff_model = factor(meff_model, levels = c("gwas", "rrblup"), labels = c("GWAS", "RR-BLUP"))) %>%
  unite(col = "order", meff_model:pamStage, remove = FALSE) %>%
  mutate(order = match(order, net_alias$order)) %>%
  arrange(order)

##### a #####

plot_k_metrics <- qc_results %>%
  filter(module != "grey") %>%
  group_by(order, meff_model, norm_method, minsize, pamStage, module) %>%
  mutate(n_markers = n()) %>%
  ungroup() %>%
  pivot_longer(c(kDiff, clusterCoeff, kRatio), names_to = "metric", values_to = "values") %>%
  mutate(module = factor(module),
         order = factor(order)) %>%
  ggplot() +
  facet_wrap(~ metric, nrow = 3, scales = "free_y") +
  geom_boxplot(aes(x = order, y = values), fill = "grey80", outlier.size = 0.25, show.legend = FALSE) +
  facetted_pos_scales(y = list(metric == "clusterCoeff" ~ scale_y_continuous(breaks = seq(0, 0.4, 0.2)),
                               metric == "kDiff" ~ scale_y_continuous(breaks = seq(-20, 20, 20)),
                               metric == "kRatio" ~ scale_y_continuous(breaks = seq(-1, 1, 1)))) +
  labs(x = "Networks", y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/fig_3a.png"), plot = plot_k_metrics,
       device = "png", dpi = 350, units = "in", width = 4, height = 3)

##### b #####

plot_pos_kDiff <- qc_results %>%
  filter(module != "grey") %>%
  group_by(order, meff_model, norm_method, minsize, pamStage, module) %>%
  summarize(n_markers = n(),
            median_kDiff = median(kDiff),
            mean_kDiff = mean(kDiff),
            se_kDiff = mean_kDiff / sqrt(n_markers)) %>%
  ungroup() %>%
  group_by(order, meff_model, norm_method, minsize, pamStage) %>%
  summarize(n_modules = n(),
            positive_mean_kDiff = round(sum(mean_kDiff > 0) / n_modules, digits = 2)) %>%
  ungroup() %>%
  pivot_longer(positive_mean_kDiff, names_to = "metric", values_to = "values") %>%
  mutate(order = factor(order)) %>%
  ggplot() +
  geom_col(aes(x = order, y = values), fill = "grey40", show.legend = FALSE) +
  labs(x = "Networks", y = "Proportion of modules\nwith mean kDiff > 0") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/fig_3b.png"), plot = plot_pos_kDiff,
       device = "png", dpi = 350, units = "in", width = 4, height = 2)

# clean R environment
rm(list = ls())



#### figure 4 -----

input_folder <- "analysis/networks/YLD/meff_gwas/norm_minmax/min_mod_size_25/pamStage_off"
meff_mod_file <- paste0(input_folder, "/define_network_modules.RData")
modules_filename <- paste0(input_folder, "/kDiff_per_module.txt")
mod_ld_folder <- paste0(input_folder, "/modules_ld")
output_folder <- "figures"
type_connect_str <- "TOM"
soft_threshold <- 24
edge_threshold <- 0.01
mod <- "firebrick2"
ld_module <- paste0(mod_ld_folder, "/ld_markers_", mod, ".ld.gz")

##### a #####

# load ld file
ld_module <- fread(ld_module, header = TRUE, data.table = FALSE)

# plot distribution r2 for all markers
ld_dist_plot <- ggplot(ld_module, aes(x = R2)) +
  geom_histogram(data = subset(ld_module, R2 >= 0.9), color = "black", fill = "gold2", binwidth = 0.05) +
  geom_histogram(data = subset(ld_module, R2 < 0.9), color = "black", fill = "grey50", binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = bquote("Linkage disequilibrium" ~ (r^2)), y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/fig_4a.png"), plot = ld_dist_plot,
       device = "png", dpi = 350, units = "in", width = 4, height = 2)

##### b #####

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
             aes(x = x, y = y, xend = xend, yend = yend), color = "grey50", size = 1, alpha = 0.2) +
  geom_edges(data = subset(mod_network, LD == TRUE),
             aes(x = x, y = y, xend = xend, yend = yend), color = "gold3", size = 1) +
  geom_nodes(data = mod_network,
             aes(x = x, y = y), color = "black", size = 3) +
  theme_blank()
ggsave(filename = paste0(output_folder, "/fig_4b.png"), plot = plot_net_ld,
       device = "png", dpi = 350, units = "in", width = 4, height = 3)

##### c (MSI) #####
# ran this on MSI

# get network aliases
net_alias_file <- "supp_materials/supp_table6.xlsx"
net_alias <- read_excel(net_alias_file)
net_alias$order <- paste(net_alias$Model, net_alias$Normalization, net_alias$`Minimum module size`, net_alias$pamStage, sep = "_")

# file name
ld_results_file <- "analysis/networks/YLD/summary_ld_per_network.txt.gz"
# ld_results_file <- "analysis/networks/YLD/summary_ld_per_network.TEST.txt"

# load file
ld_results <- fread(ld_results_file, header = TRUE, data.table = FALSE)
# adjust files to match network aliases
ld_results <- ld_results %>%
  mutate(norm_method = factor(norm_method, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         meff_model = factor(meff_model, levels = c("gwas", "rrblup"), labels = c("GWAS", "RR-BLUP"))) %>%
  unite(col = "order", meff_model:pamStage, remove = FALSE) %>%
  mutate(order = match(order, net_alias$order)) %>%
  arrange(order)

# create plot
ld_dist <- ld_results %>%
  mutate(order = factor(order)) %>%
  ggplot(aes(x = R2, y = order, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 8, rel_min_height = 0.001, size = 0.1) +
  scale_fill_binned(breaks = seq(0, 1, 0.1), guide = guide_bins(show.limits = TRUE)) +
  labs(x = bquote("Linkage disequilibrium" ~ (r^2)), y = "Networks") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  binned_scale(aesthetics = "fill", scale_name = "custom",
               palette = ggplot2:::binned_pal(scales::manual_pal(values = c(rep("#80808080", 8), rep("gold2", 2)))),
               guide = "bins",
               breaks = seq(0, 1, 0.1), limits = c(0, 1), show.limits = TRUE)
ggsave(filename = paste0(output_folder, "/fig_4c.pdf"), plot = ld_dist,
       device = "pdf", units = "in", width = 4, height = 6)

# converted to "png" on my Mac (300dpi)

# clean R environment
rm(list = ls())



#### figure 5 ----

input_folder <- "analysis/networks/YLD/meff_gwas/norm_minmax/min_mod_size_25/pamStage_off"
output_folder <- "figures"
meff_mod <- paste0(input_folder, "/define_network_modules.RData")

# load data
load(meff_mod)
rownames(MEs) <- gsub(".", "-", rownames(MEs), fixed = TRUE)
colnames(MEs) <- gsub("^ME", "", colnames(MEs), perl = TRUE)

# plot heatmap
png(paste0(output_folder, "/fig_5.png"), units = "in", res = 350, width = 8, height = 4)
pheatmap(mat = MEs, angle_col = 90, fontsize_col = 6)
dev.off()

png(paste0(output_folder, "/figure5.png"), units = "in", res = 350, width = 8, height = 4)
pheatmap(mat = MEs, angle_col = 90, fontsize_col = 6)
dev.off()

# clean R environment
rm(list = ls())



#### figure 6 ----

model <- "gwas"
norm <- "minmax"
size <- "25"
pam <- "off"
meff_mod_Rdata <- paste0("analysis/networks/YLD/meff_", model, "/norm_", norm, "/min_mod_size_", size,
                         "/pamStage_", pam, "/define_network_modules.RData")
mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
env_idx_file <- "data/env_covariables/pca_env_idx.txt"
output_folder <- "figures"
pca_contrib_file <- "data/env_covariables/pca_contributions.txt"

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
  geom_point(data = head(mod_env_idx_results, n = 3), shape = 21, size = 3,
             fill = mod_env_idx_results$module[1:3], show.legend = FALSE) +
  geom_point(data = tail(mod_env_idx_results, n = nrow(mod_env_idx_results) - 3),
             size = 3, color = "grey60", alpha = 0.4, show.legend = FALSE) +
  geom_label_repel(data = head(mod_env_idx_results, n = 3), max.overlaps = 10, nudge_y = c(-0.3, 0.5, 0.3), nudge_x = c(0, 0, 0.2),
                   fill = "white", size = 3, show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1)) +
  labs(x = "p-values", y = "Correlation") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/fig_6a.png"), plot = plot_cor_mod_pc,
       device = "png", dpi = 350, units = "in", width = 4, height = 5)

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
    geom_line(aes(group = type), size = 1) +
    scale_color_manual(values = c(gsub("^ME", "",  mod_env_idx_cor[row, "module"], perl = TRUE), "grey70")) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = "Environments", y = "Normalized value", color = "Eigenvalues") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.title.y = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
          legend.direction = "horizontal",
          legend.position = "top",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.margin = margin(0,0,0,0))
  ggsave(filename = paste0(output_folder, "/fig_6b", row, ".png"), plot = plot_mod_env_idx_norm,
         device = "png", dpi = 350, units = "in", width = 4, height = 2)

}


##### c #####

# load pc contributions
pca_contrib <- fread(pca_contrib_file, header = TRUE, data.table = FALSE)

pca_contrib[, c("idx", "interval")] <- do.call(rbind, lapply(pca_contrib$env_idx, function(idx) {
  idx <- unlist(strsplit(idx, split = "_"))
  if (length(idx) > 2) idx <- c(paste0(idx[1:(length(idx) - 1)], collapse = "_"), idx[length(idx)])
  return(idx)
}))

pca_contrib <- pca_contrib %>%
  select(-env_idx) %>%
  pivot_longer(-c(idx, interval), names_to = "pc", values_to = "contrib") %>%
  mutate(interval = as.numeric(interval))

for (row in 1:nrow(mod_env_idx_cor)) {

  # get pc
  PC <- mod_env_idx_cor[row, "env_idx"]

  # plot summary contributions
  plot_summary_pc <- pca_contrib %>%
    filter(pc == PC) %>%
    arrange(desc(contrib)) %>%
    ggplot(aes(x = interval, y = contrib)) +
    facet_wrap(~ idx, ncol = 6) +
    geom_smooth(se = TRUE, color = "grey30") +
    scale_x_continuous(breaks = c(25, 125)) +
    labs(x = "Time interval", y = "Total contributions (%)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 6, margin = margin(1, 0, 1, 0)))
  plot_summary_pc <- plot_summary_pc +
  scale_y_continuous(breaks = seq(0, round(max(ggplot_build(plot_summary_pc)$data[[1]]["y"]), digits = 1), length.out = 2))
  ggsave(filename = paste0(output_folder, "/fig_6c", row, "_", PC, ".png"), plot = plot_summary_pc,
         device = "png", dpi = 350, units = "in", width = 4, height = 3)

}

# clean R environment
rm(list = ls())



#### figure 7 ----

mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
env_idx_file <- "data/env_covariables/pca_env_idx.txt"
net_alias_file <- "supp_materials/supp_table6.xlsx"
output_folder <- "figures"
p_value <- 0.11
indices <- c("PC4", "PC5", "PC7", "PC9")

# load file with network aliases
net_alias <- read_excel(net_alias_file)
net_alias$Model <- factor(net_alias$Model, levels = c("GWAS", "RR-BLUP"), labels = c("gwas", "rrblup"))
net_alias$Normalization <- factor(net_alias$Normalization, levels = c("Minmax", "Z-score"), labels = c("minmax", "zscore"))

for (idx in indices) {

  list_markers_for_ld <- paste0("analysis/networks/YLD/overlap_markers/markers_correlated_", idx, ".txt")

  # load module-env index relationship
  mod_env_idx_cor <- fread(mod_env_idx_cor_file, header = TRUE, data.table = FALSE)
  # select only significant associations
  mod_env_idx_cor <- subset(mod_env_idx_cor, pval <= p_value)
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
    net_name <- subset(net_alias, Model == model & Normalization == norm & `Minimum module size` == size & pamStage == pam)[, 1]
    network_mod_markers <- data.frame(network = net_name,
                                      network_mod_markers, stringsAsFactors = FALSE)
    network_mod_markers <- tidyr::unite(network_mod_markers, Network:module, col = "network_module", sep = "-")
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
  # reorder list
  list_markers_mod_idx <- list_markers_mod_idx[order(names(list_markers_mod_idx))]

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

  # rename networks for plotting
  names(list_markers_mod_idx) <- sapply(names(list_markers_mod_idx), function(x) {
    paste0("network ", gsub("_", "\n(", x), ")")
  })
  # # uncomment if want to use cor and pval values
  # sig_env_idx$name <- apply(sig_env_idx, MARGIN = 1, function(sig_cor, net_alias) {
  #   name <- subset(net_alias, Model == sig_cor["meff_model"]
  #                  & Normalization == sig_cor["norm_method"]
  #                  & `Minimum module size` == sig_cor["minsize"]
  #                  & pamStage == sig_cor["pamStage"])[, 1]
  #   name <- paste0(name, "_", sig_cor["module"])
  #   return(name)
  # }, net_alias = net_alias)
  # sig_env_idx <- sig_env_idx[match(names(list_markers_mod_idx), sig_env_idx$name), ]
  # sig_env_idx$alias <- apply(sig_env_idx, MARGIN = 1, function(x) {
  #   paste0("network ", gsub("_", " - ", x["name"]), "\n(cor: ", x["cor"], ", pval: ", x["pval"], ")")
  # })
  # names(list_markers_mod_idx) <- sig_env_idx[, "alias"]

  # plot intersections of markers
  upset_plot <- upset(fromList(list_markers_mod_idx), sets = rev(names(list_markers_mod_idx)), keep.order = TRUE,
                      mainbar.y.label = paste0("Shared\nmarkers\nassociated\nwith ", idx), sets.x.label = "Markers per module",
                      order.by = "freq", mb.ratio = c(0.4, 0.6), nsets = 100, text.scale = c(0.8, 1, 1, 0.8, 1, 1))
  # text.scale = c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  png(paste0(output_folder, "/fig_7", letters[which(idx == indices)], ".png"), units = "in", res = 350, width = 4, height = 2.5)
  print(upset_plot)
  dev.off()

}

# clean R environment
rm(list = ls())


#### supp figure 1 ----

# output from 'scripts/get_BLUEs_from_empirical-data.R'



#### supp figure 2 ----

# output from 'scripts/estimate_marker_effects.R'



#### supp figure 3 ----

# output from 'scripts/estimate_marker_effects.R'



#### supp figure 4 ----

# output from 'scripts/pca_env_idx.R'


#### supp figure 5 ----

results_folder <- "analysis/networks/YLD"
output_folder <- "figures"

# create empty df to store results
sft_all <- data.frame(stringsAsFactors = FALSE)

for (model in c("rrblup", "gwas")) {
  for (norm in c("minmax", "zscore")) {

    if (norm == "minmax") {
      for (cv in c(0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25)) {

        # load file
        sft_file <- paste0(results_folder, "/meff_", model, "/norm_", norm,
                           "/cv_thresholds/scale-free_topology_fit_index.cv-", cv, ".txt")
        sft_results <- fread(sft_file, header = TRUE, data.table = FALSE)
        # get results
        sft_all <- rbind(sft_all, data.frame(model = model, norm = norm,  cv = cv, sft_results))

      }
    }

    if (any(norm == "zscore" | norm == "none")) {
      for (cv in c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1)) {

        # load file
        sft_file <- paste0(results_folder, "/meff_", model, "/norm_", norm,
                           "/cv_thresholds/scale-free_topology_fit_index.cv-", cv, ".txt")
        sft_results <- fread(sft_file, header = TRUE, data.table = FALSE)
        # get results
        sft_all <- rbind(sft_all, data.frame(model = model, norm = norm,  cv = cv, sft_results))

      }
    }

  }
}
rm(sft_results)

# add signed R2 values
sft_all$signed_R2 <- -sign(sft_all$slope) * sft_all$SFT.R.sq

##### a -----

mean_k_plot <- sft_all %>%
  mutate(norm = factor(norm, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         model = factor(model, levels = c("gwas", "rrblup"), labels = c("GWAS", "rrBLUP"))) %>%
  ggplot(aes(x = Power, y = mean.k., group = as.factor(cv), color = as.factor(cv))) +
  facet_grid(model ~ norm) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(end = 0.9, option = "B", direction = -1) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Power", y = "Mean connectivity", color = "CV cut-off") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/supp_fig_4a.png"), plot = mean_k_plot,
       device = "png", dpi = 350, units = "in", width = 10, height = 6)

##### b -----

median_k_plot <- sft_all %>%
  mutate(norm = factor(norm, levels = c("minmax", "zscore"), labels = c("Minmax", "Z-score")),
         model = factor(model, levels = c("gwas", "rrblup"), labels = c("GWAS", "rrBLUP"))) %>%
  ggplot(aes(x = Power, y = median.k., group = as.factor(cv), color = as.factor(cv))) +
  facet_grid(model ~ norm) +
  geom_line() +
  geom_point() +
  scale_color_viridis_d(end = 0.9, option = "B", direction = -1) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = "Power", y = "Median connectivity", color = "CV cut-off") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/supp_fig_4b.png"), plot = median_k_plot,
       device = "png", dpi = 350, units = "in", width = 10, height = 6)

# clean R environment
rm(list = ls())




#### supp figure 6 (MSI) ----
# note: ran this on MSI

input_folder <- "analysis/networks/YLD"
output_folder <- "figures"
i <- 1

for (model in c("gwas", "rrblup")) {
  for (norm in c("minmax", "zscore")) {
    for (size in c(25, 50, 100)) {
      for (pam in c("off", "on")) {

        # get filename
        meff_mod <- paste0(input_folder, "/meff_", model, "/norm_", norm,
                           "/min_mod_size_", size, "/pamStage_", pam,
                           "/define_network_modules.RData")

        # exclude networks with gwas and 100 min mod size
        if (!all(model == "gwas" & size == 100)) {

          # load data
          load(meff_mod)
          # plot heatmap
          pdf(paste0(output_folder, "/supp-fig_5-net_", i, ".pdf"), width = 8, height = 10)
          pheatmap(mat = t(MEs), angle_col = 0, main = paste0("Network ", i))
          dev.off()
          # add count
          i <- i + 1

        }

      }
    }
  }
}

# clean R environment
rm(list = ls())



#### supp figure 7 (MSI) ----
# note: ran this on MSI

mod_env_idx_cor_file <- "analysis/networks/YLD/module-env-idx_per_network.pca.txt"
env_idx_file <- "data/env_covariables/pca_env_idx.txt"
output_folder <- "figures"
i <- 1

for (model in c("gwas", "rrblup")) {
  for (norm in c("minmax", "zscore")) {
    for (size in c(25, 50, 100)) {
      for (pam in c("off", "on")) {

        # get filename
        meff_mod_Rdata <- paste0("analysis/networks/YLD/meff_", model, "/norm_", norm, "/min_mod_size_", size,
                                 "/pamStage_", pam, "/define_network_modules.RData")

        # exclude networks with gwas and 100 min mod size
        if (!all(model == "gwas" & size == 100)) {

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
            labs(title = paste0("Network ", i), x = "p-values", y = "Correlation") +
            theme_bw() +
            theme(panel.grid = element_blank())
          ggsave(filename = paste0(output_folder, "/supp-fig_7-net_", i, ".pdf"), plot = plot_cor_mod_pc,
                 device = "pdf", units = "in", width = 4, height = 4)

          # add count
          i <- i + 1

        }

      }
    }
  }
}

# clean R environment
rm(list = ls())
