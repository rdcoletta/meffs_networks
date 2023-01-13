library(data.table)
library(dplyr)
library(ggplot2)
library(ggridges)
library(viridis)
library(ggh4x)

usage <- function() {
  cat("
description: plot distribution of LD for all networks.

usage: Rscript plot_ld_summary_per_network.R [ld_results_file] [output_folder] [...]

positional arguments:
  marker_effects_file         file containing LD results from all networks
  output_folder               name of file to save plot

optional argument:
  --help                      show this helpful message


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
ld_results_file <- args[1]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) != pos_args) stop(usage(), "missing positional argument(s)")



#### plot ld distribution ----

# load file
ld_results <- fread(ld_results_file, header = TRUE, data.table = FALSE)

# create plot
ld_dist <- ld_results %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>%
  mutate(network = cur_group_id()) %>%
  ungroup() %>%
  mutate(network = factor(network)) %>%
  ggplot(aes(x = R2, y = network, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 8, rel_min_height = 0.001) +
  scale_fill_binned(breaks = seq(0, 1, 0.1), guide = guide_bins(show.limits = T)) +
  labs(title = "LD distribution across all modules of each network",
       x = bquote(R^2), y = "Networks") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        axis.text.y = element_blank()) +
  binned_scale(aesthetics = "fill", scale_name = "custom", 
               palette = ggplot2:::binned_pal(scales::manual_pal(values = c(rep("#80808080", 8), rep("firebrick", 2)))),
               guide = "bins",
               breaks = seq(0, 1, 0.1), limits = c(0, 1), show.limits = T)
ggsave(plot = ld_dist, filename = paste0(output_folder, "/ld_dist_all_networks.pdf"),
       device = "pdf", width = 8, height = 8)

# calculate proportion markers in ld per module
prop_markers_ld_per_mod <- ld_results %>% 
  group_by(meff_model, norm_method, minsize, pamStage, module) %>% 
  summarize(prop_markers_ld = sum(R2 >= 0.9) / n()) %>% 
  ungroup()
fwrite(prop_markers_ld_per_mod, file = paste0(output_folder, "/n_markers_ld_networks.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# plot proportion markers in ld per module
prop_markers_ld_plot <- prop_markers_ld_per_mod %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>% 
  mutate(network = cur_group_id()) %>%
  ungroup() %>% 
  mutate(network = factor(network)) %>%
  ggplot() +
  facet_nested(~ meff_model + norm_method + minsize + pamStage, scales = "free",
               nest_line = element_line(linetype = 1, color = "gray80", lineend = "butt")) +
  geom_boxplot(aes(x = network, y = prop_markers_ld), outlier.size = 0.25, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Networks", y = "Proportion of marker pairs in LD within a module") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0, "lines"))
ggsave(plot = prop_markers_ld_plot, filename = paste0(output_folder, "/n_markers_ld_networks.pdf"),
       device = "pdf", width = 10, height = 6)



#### debug ----

# ld_results_file <- "analysis/networks/YLD/summary_ld_per_network.TEST.txt"
# output_folder <- "analysis/networks/YLD/qc_networks"


# ld_summary <- fread(paste0(output_folder, "/n_markers_ld_networks.txt"), header = TRUE, data.table = FALSE)
# 
# ld_summary %>%
#   filter(meff_model == "gwas"
#          & norm_method == "minmax"
#          & minsize == 25
#          & pamStage == "off"
#          & module == "firebrick2")
# 
# ld_summary %>%
#   group_by(meff_model, norm_method, minsize, pamStage) %>%
#   summarize(avg_prop_ld = mean(prop_markers_ld),
#             se_prop_ld = avg_prop_ld / sqrt(n())) %>%
#   ungroup()
# 
# ld_summary %>%
#   summarize(avg_prop_ld = mean(prop_markers_ld),
#             se_prop_ld = avg_prop_ld / sqrt(n()))
