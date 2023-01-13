library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)

usage <- function() {
  cat("
description: summarize QC results of marker effect networks.

usage: Rscript summarize_qc_networks.R [folder_base] [...]

positional arguments:
  folder_base                 path to folder with QC results

optional argument:
  --help                      show this helpful message
  --meff-model=[LIST]         comma-separated list of marker effect models
  --norm-method=[LIST]        comma-separated list of normalization methods
  --minsize=[LIST]            comma-separated list of minimum module size
  --pamStage=[LIST]           comma-separated list of pamStage values

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
folder_base <- args[1]

# set default
meff_model <- "rrblup,gwas"
norm_method <- "minmax,zscore"
minsize <- "25,50,100"
pamStage <- "on,off"

# assert to have the correct optional arguments
if (length(args) < 1) stop(usage(), "missing positional argument(s)")

if (length(args) > 1) {

  opt_args <- args[-1]
  opt_args_allowed <- c("--meff-model", "--norm-method", "--minsize", "--pamStage")
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

# adjust format from optional arguments
meff_model <- unlist(strsplit(meff_model, split = ","))
norm_method <- unlist(strsplit(norm_method, split = ","))
minsize <- as.numeric(unlist(strsplit(minsize, split = ",")))
pamStage <- unlist(strsplit(pamStage, split = ","))



#### summarize gwas qc results ----

# create empty df to store results
qc_results <- data.frame(stringsAsFactors = FALSE)

for (model in meff_model) {
  for (norm in norm_method) {
    for (size in minsize) {
      for (pam in pamStage) {

        # get folder with qc results for a network setting
        folder_qc <- paste0(folder_base, "/meff_", model, "/norm_", norm,
                                "/min_mod_size_", size, "/pamStage_", pam)

        # get results
        network_qc <- paste0(folder_qc, "/kDiff_per_module.txt")
        network_qc <- try(fread(network_qc, header = TRUE, data.table = FALSE))

        if (class(network_qc) != "try-error") {
          # add network settings
          network_qc <- data.frame(meff_model = model, norm_method = norm,
                                   minsize = size, pamStage = pam,
                                   network_qc, stringsAsFactors = FALSE)
          # get accuracy results
          qc_results <- rbind(qc_results, network_qc)
        }

      }
    }
  }
}

# write summary
fwrite(qc_results, file = paste0(folder_base, "/summary_qc_networks.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# create output folder
output_folder <- paste0(folder_base, "/qc_networks")
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# number of markers per network
plot1 <- qc_results %>%
  filter(source == "TOM") %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>%
  summarize(n_markers = n()) %>%
  ungroup() %>%
  group_by(meff_model, norm_method) %>%
  summarize(n_markers = mean(n_markers)) %>%
  ungroup() %>%
  ggplot() +
  geom_col(aes(x = paste(meff_model, norm_method, sep = "_"), y = n_markers)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(x = "network", y = "number of markers")

ggsave(plot = plot1, filename = paste0(output_folder, "/n_markers.pdf"),
       device = "pdf", width = 8, height = 6)

# number of markers per module - with grey module
plot2 <- qc_results %>%
  filter(source == "TOM") %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n()) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_y", independent = "y") +
  geom_histogram(aes(x = n_markers), binwidth = 50) +
  labs(title = "Number of markers per module",
       subtitle = "(grey module included)")

ggsave(plot = plot2, filename = paste0(output_folder, "/dist_markers.pdf"),
       device = "pdf", width = 12, height = 6)

# number of markers per module - without grey module
plot3 <- qc_results %>%
  filter(source == "TOM", module != "grey") %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n()) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_y", independent = "y") +
  geom_histogram(aes(x = n_markers), binwidth = 50) +
  labs(title = "Number of markers per module",
       subtitle = "(grey module not included)")

ggsave(plot = plot3, filename = paste0(output_folder, "/dist_markers.no-grey-mod.pdf"),
       device = "pdf", width = 12, height = 6)


# plot stats distribution per network
plot4 <- ggplot(qc_results, aes(x = source, y = kDiff, fill = source)) +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
  geom_violin(draw_quantiles = 0.5, show.legend = FALSE)

ggsave(plot = plot4, filename = paste0(output_folder, "/dist_kDiff.pdf"),
       device = "pdf", width = 10, height = 6)

plot5 <- ggplot(qc_results, aes(x = source, y = kWithin, fill = source)) +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
  geom_violin(draw_quantiles = 0.5, show.legend = FALSE)

ggsave(plot = plot5, filename = paste0(output_folder, "/dist_kWithin.pdf"),
       device = "pdf", width = 10, height = 6)

plot6 <- ggplot(qc_results, aes(x = source, y = kOut, fill = source)) +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
  geom_violin(draw_quantiles = 0.5, show.legend = FALSE)

ggsave(plot = plot6, filename = paste0(output_folder, "/dist_kOut.pdf"),
       device = "pdf", width = 10, height = 6)

plot7 <- ggplot(qc_results, aes(x = source, y = clusterCoeff, fill = source)) +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
  geom_violin(draw_quantiles = 0.5, show.legend = FALSE)

ggsave(plot = plot7, filename = paste0(output_folder, "/dist_clusterCoeff.pdf"),
       device = "pdf", width = 10, height = 6)

# plot average median kDiff
plot8 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            median = median(kDiff)) %>%
  ungroup() %>%
  group_by(meff_model, norm_method, minsize, pamStage, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_median = mean(median, na.rm = TRUE),
            se_median = sd(median, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
  geom_col(aes(x = source, y = avg_median, fill = source), show.legend = FALSE) +
  geom_errorbar(aes(x = source, ymin = avg_median - se_median, ymax = avg_median + se_median),
                position = position_dodge(0.9), width = 0.2) +
  labs(y = "average median kDiff")

ggsave(plot = plot8, filename = paste0(output_folder, "/avg_median_kDiff.pdf"),
       device = "pdf", width = 10, height = 6)

plot9 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            median = median(clusterCoeff)) %>%
  ungroup() %>%
  group_by(meff_model, norm_method, minsize, pamStage, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_median = mean(median, na.rm = TRUE),
            se_median = sd(median, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage) +
  geom_col(aes(x = source, y = avg_median, fill = source), show.legend = FALSE) +
  geom_errorbar(aes(x = source, ymin = avg_median - se_median, ymax = avg_median + se_median),
                position = position_dodge(0.9), width = 0.2) +
  labs(y = "average median clusterCoeff")

ggsave(plot = plot9, filename = paste0(output_folder, "/avg_median_clusterCoeff.pdf"),
       device = "pdf", width = 10, height = 6)

plot10 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            prop_pos_kDiff = round(sum(kDiff > 0) / n_markers, digits = 2)) %>%
  ungroup() %>%
  group_by(meff_model, norm_method, minsize, pamStage, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_prop_pos_kDiff = mean(prop_pos_kDiff, na.rm = TRUE),
            se_prop_pos_kDiff = sd(prop_pos_kDiff, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_col(aes(x = source, y = avg_prop_pos_kDiff, fill = source), show.legend = FALSE) +
  geom_errorbar(aes(x = source, ymin = avg_prop_pos_kDiff - se_prop_pos_kDiff, ymax = avg_prop_pos_kDiff + se_prop_pos_kDiff),
                position = position_dodge(0.9), width = 0.2) +
  labs(y = "average proportion of positive kDiff")

ggsave(plot = plot10, filename = paste0(output_folder, "/avg_prop_pos_kDiff.pdf"),
       device = "pdf", width = 10, height = 6)

# average stats per type network setting
plot11 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            median = median(kDiff)) %>%
  ungroup() %>%
  group_by(pamStage, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_median = mean(median, na.rm = TRUE),
            se_median = sd(median, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(~ pamStage, scales = "free_x", independent = "x") +
  geom_col(aes(x = source, y = avg_median)) +
  geom_errorbar(aes(x = source, ymin = avg_median - se_median, ymax = avg_median + se_median),
                position = position_dodge(0.9), width = 0.2) +
  labs(title = "average median kDiff per pamStage")

ggsave(plot = plot11, filename = paste0(output_folder, "/avg_median_kDiff.per-pamStage.pdf"),
       device = "pdf", width = 8, height = 6)

plot12 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            median = median(kDiff)) %>%
  ungroup() %>%
  group_by(minsize, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_median = mean(median, na.rm = TRUE),
            se_median = sd(median, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(~ minsize, scales = "free_x", independent = "x") +
  geom_col(aes(x = source, y = avg_median)) +
  geom_errorbar(aes(x = source, ymin = avg_median - se_median, ymax = avg_median + se_median),
                position = position_dodge(0.9), width = 0.2) +
  labs(title = "average median kDiff per minsize")

ggsave(plot = plot12, filename = paste0(output_folder, "/avg_median_kDiff.per-minsize.pdf"),
       device = "pdf", width = 8, height = 6)

plot13 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            median = median(kDiff)) %>%
  ungroup() %>%
  group_by(norm_method, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_median = mean(median, na.rm = TRUE),
            se_median = sd(median, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(~ norm_method, scales = "free_x", independent = "x") +
  geom_col(aes(x = source, y = avg_median)) +
  geom_errorbar(aes(x = source, ymin = avg_median - se_median, ymax = avg_median + se_median),
                position = position_dodge(0.9), width = 0.2) +
  labs(title = "average median kDiff per norm_method")

ggsave(plot = plot13, filename = paste0(output_folder, "/avg_median_kDiff.per-norm_method.pdf"),
       device = "pdf", width = 8, height = 6)

plot14 <- qc_results %>%
  group_by(meff_model, norm_method, minsize, pamStage, source, module) %>%
  summarize(n_markers = n(),
            median = median(kDiff)) %>%
  ungroup() %>%
  group_by(meff_model, source) %>%
  summarize(n_mods = n(),
            n_markers = sum(n_markers),
            avg_median = mean(median, na.rm = TRUE),
            se_median = sd(median, na.rm = TRUE) / sqrt(n_mods)) %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(~ meff_model, scales = "free_x", independent = "x") +
  geom_col(aes(x = source, y = avg_median)) +
  geom_errorbar(aes(x = source, ymin = avg_median - se_median, ymax = avg_median + se_median),
                position = position_dodge(0.9), width = 0.2) +
  labs(title = "average median kDiff per meff_model")

ggsave(plot = plot14, filename = paste0(output_folder, "/avg_median_kDiff.per-meff_model.pdf"),
       device = "pdf", width = 8, height = 6)

# plot all kDiff and kRatio (kDiff/kTotal)
qc_results$kRatio <- round(qc_results$kDiff / qc_results$kTotal, digits = 2)
qc_results$module <- factor(qc_results$module)

plot15 <- qc_results %>%
  filter(source == "TOM") %>%
  mutate(module = factor(module)) %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_boxplot(aes(x = module, y = kDiff, fill = module), outlier.size = 0.25, show.legend = FALSE) +
  # geom_violin(aes(x = module, y = kDiff, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(qc_results$module)) +
  labs(title = "kDiff per network", subtitle = "(TOM matrix)") +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
ggsave(plot = plot15, filename = paste0(output_folder, "/kDiff_per_network.TOM.pdf"),
       device = "pdf", width = 28, height = 12)

plot16 <- qc_results %>%
  filter(source == "TOM") %>%
  mutate(module = factor(module)) %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_boxplot(aes(x = module, y = kRatio, fill = module), outlier.size = 0.25, show.legend = FALSE) +
  # geom_violin(aes(x = module, y = kRatio, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(qc_results$module)) +
  labs(title = "kRatio per network", subtitle = "(TOM matrix)") +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
ggsave(plot = plot16, filename = paste0(output_folder, "/kRatio_per_network.TOM.pdf"),
       device = "pdf", width = 28, height = 12)

plot17 <- qc_results %>%
  filter(source == "adjacency") %>%
  mutate(module = factor(module)) %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_boxplot(aes(x = module, y = kDiff, fill = module), outlier.size = 0.25, show.legend = FALSE) +
  # geom_violin(aes(x = module, y = kDiff, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(qc_results$module)) +
  labs(title = "kDiff per network", subtitle = "(adjacency matrix)") +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
ggsave(plot = plot17, filename = paste0(output_folder, "/kDiff_per_network.adjacency.pdf"),
       device = "pdf", width = 28, height = 12)

plot18 <- qc_results %>%
  filter(source == "adjacency") %>%
  mutate(module = factor(module)) %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_boxplot(aes(x = module, y = kRatio, fill = module), outlier.size = 0.25, show.legend = FALSE) +
  # geom_violin(aes(x = module, y = kRatio, fill = module), show.legend = FALSE) +
  scale_fill_manual(values = levels(qc_results$module)) +
  labs(title = "kRatio per network", subtitle = "(adjacency matrix)") +
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() +
  theme(panel.grid.minor =  element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5))
ggsave(plot = plot18, filename = paste0(output_folder, "/kRatio_per_network.adjacency.pdf"),
       device = "pdf", width = 28, height = 12)

# n modules per net
plot19 <- qc_results %>%
  filter(source == "adjacency") %>%
  group_by(meff_model, norm_method, minsize, pamStage, module) %>%
  summarize(n_markers = n()) %>%
  ungroup() %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>% 
  mutate(n_modules = n(),
         network = cur_group_id()) %>%
  ungroup() %>% 
  mutate(n_markers_no_grey = NA,
         n_markers_no_grey = case_when(module != "grey" ~ n_markers)) %>% 
  pivot_longer(c(n_markers, n_markers_no_grey, n_modules), names_to = "metric", values_to = "values") %>%
  mutate(network = factor(network),
         metric = factor(metric, levels = c("n_modules", "n_markers", "n_markers_no_grey"),
                         labels = c("Total modules", "Unassigned markers", "Markers per module"))) %>% 
  ggplot() +
  facet_nested(metric ~ meff_model + norm_method + minsize + pamStage, scales = "free",
               nest_line = element_line(linetype = 1, color = "gray80", lineend = "butt")) +
  geom_col(data = function(x) subset(x, metric == "Total modules"),
           aes(x = network, y = values), show.legend = FALSE) +
  geom_col(data = function(x) subset(x, metric == "Unassigned markers"  & module == "grey"),
           aes(x = network, y = values), fill = "grey", show.legend = FALSE) +
  geom_boxplot(data = function(x) subset(x, metric == "Markers per module"),
               aes(x = network, y = values, fill = metric), outlier.size = 0.25, show.legend = FALSE) +
  labs(x = "Networks", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0, "lines"))
ggsave(plot = plot19, filename = paste0(output_folder, "/n_modules_per_network.png"),
       device = "png", dpi = 300, width = 10, height = 6)
# need to add white spaces between strip groups

# connectivity + n markers
plot20 <- qc_results %>%
  # filter(source == "adjacency") %>%
  filter(source == "adjacency" & module != "grey") %>%
  group_by(meff_model, norm_method, minsize, pamStage, module) %>%
  mutate(n_markers = n()) %>%
  ungroup() %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>% 
  mutate(network = cur_group_id()) %>%
  ungroup() %>% 
  pivot_longer(c(kDiff, clusterCoeff, kRatio), names_to = "metric", values_to = "values") %>%
  mutate(module = factor(module),
         network = factor(network)) %>%
  ggplot() +
  facet_nested(metric ~ meff_model + norm_method + minsize + pamStage, scales = "free",
               nest_line = element_line(linetype = 1, color = "gray80", lineend = "butt")) +
  geom_boxplot(aes(x = network, y = values, fill = metric), outlier.size = 0.25, show.legend = FALSE) +
  labs(x = "Networks", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0, "lines"))
ggsave(plot = plot20, filename = paste0(output_folder, "/summary_metrics_per_network.png"),
       device = "png", dpi = 300, width = 10, height = 6)
# need to add white spaces between strip groups


# proportion of average positive kDiff
plot21 <- qc_results %>%
  # filter(source == "adjacency") %>%
  filter(source == "adjacency" & module != "grey") %>%
  group_by(meff_model, norm_method, minsize, pamStage, module) %>%
  summarize(n_markers = n(),
            median_kDiff = median(kDiff),
            mean_kDiff = mean(kDiff),
            se_kDiff = mean_kDiff / sqrt(n_markers)) %>%
  ungroup() %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>% 
  summarize(network = cur_group_id(),
            n_modules = n(),
            # positive_median_kDiff = round(sum(median_kDiff > 0) / n_modules, digits = 2),
            positive_mean_kDiff = round(sum(mean_kDiff > 0) / n_modules, digits = 2)) %>%
  ungroup() %>%
  # pivot_longer(c(positive_median_kDiff, positive_mean_kDiff), names_to = "metric", values_to = "values") %>%
  pivot_longer(positive_mean_kDiff, names_to = "metric", values_to = "values") %>%
  mutate(network = factor(network)) %>%
  ggplot() +
  facet_nested(metric ~ meff_model + norm_method + minsize + pamStage, scales = "free",
               nest_line = element_line(linetype = 1, color = "gray80", lineend = "butt")) +
  geom_col(aes(x = network, y = values, fill = metric), show.legend = FALSE) +
  labs(x = "Networks", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0, "lines"))
ggsave(plot = plot21, filename = paste0(output_folder, "/proportion_positive_kDiff_per_module.png"),
       device = "png", dpi = 300, width = 10, height = 6)
# need to add white spaces between strip groups


# qc_results %>%
#   filter(source == "adjacency") %>%
#   group_by(meff_model, norm_method, minsize, pamStage, module) %>%
#   summarize(n_markers = n()) %>%
#   ungroup() %>%
#   group_by(meff_model, norm_method, minsize, pamStage) %>%
#   summarize(n_modules = n(),
#             network = cur_group_id()) %>%
#   ungroup() %>% 
#   arrange(minsize, n_modules) %>% View()

# qc_results %>%
#   filter(source == "adjacency") %>%
#   group_by(meff_model, norm_method, minsize, pamStage, module) %>%
#   summarize(n_markers = n()) %>%
#   ungroup() %>%
#   group_by(meff_model, norm_method, minsize, pamStage) %>% 
#   summarize(n_modules = n(),
#             n_markers_median = median(n_markers, na.rm = TRUE), 
#             n_markers_mean = mean(n_markers, na.rm = TRUE),
#             n_markers_se = mean(n_markers, na.rm = TRUE) / sqrt(n_modules),
#             network = cur_group_id()) %>%
#   ungroup() %>%
#   arrange(minsize, n_modules) %>% View()

# qc_results %>%
#   filter(source == "adjacency") %>%
#   group_by(meff_model, norm_method, minsize, pamStage, module) %>%
#   summarize(n_markers = n()) %>%
#   ungroup() %>%
#   group_by(meff_model, norm_method, minsize, pamStage) %>%
#   summarize(n_modules = n(),
#             n_tests = n_modules * 9) %>%
#   ungroup() %>%
#   arrange(minsize, n_modules) %>% View()



#### debug ----

# folder_base <- "analysis/networks/YLD"
# meff_model <- c("rrblup", "gwas")
# norm_method <- c("minmax", "zscore")
# minsize <- c(25, 50, 100)
# pamStage <- c("on", "off")
