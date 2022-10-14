library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggh4x)

usage <- function() {
  cat("
description: summarize gwas enrichment results of marker effect networks.

usage: Rscript summarize_gwas_enrich.R [folder_base] [...]

positional arguments:
  folder_base                 path to folder with results of gwas enrichment

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



#### summarize gwas enrichment results ----

# create empty df to store results
enrichment_results <- data.frame(stringsAsFactors = FALSE)

for (model in meff_model) {
  for (norm in norm_method) {
    for (size in minsize) {
      for (pam in pamStage) {

        # get folder with enrichment results for a network setting
        folder_enrich <- paste0(folder_base, "/meff_", model, "/norm_", norm,
                                  "/min_mod_size_", size, "/pamStage_", pam)

        # get results
        network_enrich <- paste0(folder_enrich, "/gwas_enrichment_per_module.txt")
        network_enrich <- try(fread(network_enrich, header = TRUE, data.table = FALSE))

        if (class(network_enrich) != "try-error") {
          # add network settings
          network_enrich <- data.frame(meff_model = model, norm_method = norm,
                                       minsize = size, pamStage = pam,
                                       network_enrich, stringsAsFactors = FALSE)
          # get accuracy results
          enrichment_results <- rbind(enrichment_results, network_enrich)
        }

      }
    }
  }
}

# write summary
fwrite(enrichment_results, file = paste0(folder_base, "/gwas_enrichment_per_network.txt"))

# plot summary
# enrichment of gwas hits per signficant trait-module associations
plot_gwas_per_sig_mod <- enrichment_results %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>%
  filter(pval_trait_mod_cor < 0.05) %>%
  mutate(prop_gwas_hits_mod = round(n_gwas_hits_mod / n_markers_mod, digits = 2),
         prop_gwas_top_ns_mod = round(n_gwas_top_ns_mod / n_markers_mod, digits = 2),
         prop_not_gwas_hits_mod = 1 - prop_gwas_hits_mod - prop_gwas_top_ns_mod,
         sig_gwas_hyper = pval_gwas_hyper < 0.05) %>%
  select(mod, prop_gwas_hits_mod, prop_gwas_top_ns_mod, prop_not_gwas_hits_mod, sig_gwas_hyper) %>%
  arrange(prop_gwas_hits_mod, .by_group = TRUE) %>%
  pivot_longer(-c(meff_model, norm_method, minsize, pamStage, mod, sig_gwas_hyper),
               names_to = "type", values_to = "prop") %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_col(aes(x = mod, y = prop, fill = type, color = sig_gwas_hyper), position = "stack") +
  scale_fill_manual(values = c("grey20", "gray50", "grey80")) +
  scale_color_manual(values = c("white", "firebrick")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(title = "Significant trait-module associations",
       x = "Modules", y = "Proportion")

ggsave(filename = paste0(folder_base, "/gwas_enrichment_per_sig-mod.pdf"),
       plot = plot_gwas_per_sig_mod, device = "pdf", width = 15, height = 14)

# enrichment of gwas hits per not signficant trait-module associations
plot_gwas_per_not_sig_mod <- enrichment_results %>%
  group_by(meff_model, norm_method, minsize, pamStage) %>%
  filter(pval_trait_mod_cor > 0.05) %>%
  mutate(prop_gwas_hits_mod = round(n_gwas_hits_mod / n_markers_mod, digits = 2),
         prop_gwas_top_ns_mod = round(n_gwas_top_ns_mod / n_markers_mod, digits = 2),
         prop_not_gwas_hits_mod = 1 - prop_gwas_hits_mod - prop_gwas_top_ns_mod,
         sig_gwas_hyper = pval_gwas_hyper < 0.05) %>%
  select(mod, prop_gwas_hits_mod, prop_gwas_top_ns_mod, prop_not_gwas_hits_mod, sig_gwas_hyper) %>%
  arrange(prop_gwas_hits_mod, .by_group = TRUE) %>%
  pivot_longer(-c(meff_model, norm_method, minsize, pamStage, mod, sig_gwas_hyper),
               names_to = "type", values_to = "prop") %>%
  ungroup() %>%
  ggplot() +
  facet_grid2(meff_model + norm_method ~ minsize + pamStage, scales = "free_x", independent = "x") +
  geom_col(aes(x = mod, y = prop, fill = type, color = sig_gwas_hyper), position = "stack") +
  scale_fill_manual(values = c("grey20", "gray50", "grey80")) +
  scale_color_manual(values = c("white", "firebrick")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(title = "Not significant trait-module associations",
       x = "Modules", y = "Proportion")

ggsave(filename = paste0(folder_base, "/gwas_enrichment_per_not-sig-mod.pdf"),
       plot = plot_gwas_per_not_sig_mod, device = "pdf", width = 22, height = 14)



#### debug ----

# folder_base <- "analysis/networks/YLD"
# meff_model <- c("rrblup", "gwas")
# norm_method <- c("minmax", "zscore")
# minsize <- c(25, 50, 100)
# pamStage <- c("on", "off")
