library(data.table)
library(tidyr)
library(tibble)
library(dplyr)
library(ggplot2)

usage <- function() {
  cat("
description: plot estimated marker effects and GEBVs from different models.

usage: Rscript plot_marker_effects.R [results_folder] [...]

positional arguments:
  results_folder              folder with results of marker effects
  outfolder                   folder to save plots

optional argument:
  --help                      show this helpful message
  --traits=VALUE              comma-separated list of traits to plot (default: 'EHT,Moisture,PHT,TWT,YLD')
  --models=VALUE              comma-separated list of models to plot (default: 'rrblup,bayescpi,mrr,gwas')
  --no-missing-genotypes      exclude genotypes with missing data in any environment


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
results_folder <- args[1]
outfolder <- args[2]
if (!dir.exists(outfolder)) dir.create(outfolder, recursive = TRUE)

# set default of optional args
traits <- "EHT,Moisture,PHT,TWT,YLD"
models <- "rrblup,bayescpi,mrr,gwas"
no_missing_genotypes <- FALSE

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--traits", "--models", "--no-missing-genotypes")
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

# transform optional arguments to list
traits <- unlist(strsplit(traits, split = ","))
models <- unlist(strsplit(models, split = ","))



#### plot results ----

for (trait in traits) {
  
  cat(trait, "\n", sep = "")
  
  # create plots folder
  plots_folder <- paste0(outfolder, "/", trait)
  if (!dir.exists(plots_folder)) dir.create(plots_folder)
  
  # create empty data frames to store results across models
  marker_effects_all <- data.frame(stringsAsFactors = FALSE)
  gebvs_all <- data.frame(stringsAsFactors = FALSE)
  
  for (model in models) {
    
    # get filename
    if (no_missing_genotypes) {
      file_effects <- paste0(results_folder, "/", trait, "/marker_effects.", model, ".no-missing-genos.txt")
      file_gebvs <- paste0(results_folder, "/", trait, "/GEBVs.", model, ".no-missing-genos.txt")
    } else {
      file_effects <- paste0(results_folder, "/", trait, "/marker_effects.", model, ".txt")
      file_gebvs <- paste0(results_folder, "/", trait, "/GEBVs.", model, ".txt")
    }
    # load file
    marker_effects <- fread(file_effects, header = TRUE, data.table = FALSE)
    gebvs <- fread(file_gebvs, header = TRUE, data.table = FALSE)
    # transform to long format
    marker_effects <- pivot_longer(marker_effects, -marker, names_to = "env", values_to = "effect")
    gebvs <- pivot_longer(gebvs, -genotype, names_to = "env", values_to = "gebv")
    
    # add marker number to facilitate plotting of effects
    marker_effects <- merge(x = marker_effects,
                            y = data.frame(marker = unique(marker_effects$marker),
                                           marker_n = 1:length(unique(marker_effects$marker))),
                            by = "marker")
    marker_effects$marker_n <- as.numeric(marker_effects$marker_n)
    marker_effects <- marker_effects[order(marker_effects$marker_n, marker_effects$env), ]
    # plot effects of each marker
    plot_eff_per_marker <- ggplot(marker_effects) +
      geom_hex(aes(x = marker_n, y = effect)) +
      geom_hline(yintercept = 0, color = "firebrick") +
      facet_wrap(~ env) +
      scale_fill_viridis_c() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      labs(title = trait, subtitle = model)
    ggsave(filename = paste0(plots_folder, "/effects_per_marker.", model, ".pdf"),
           plot = plot_eff_per_marker, device = "pdf", width = 14, height = 10)

    # plot distibution of effects
    plot_hist_effects <- ggplot(marker_effects) +
      geom_histogram(aes(x = effect)) +
      facet_wrap(~ env) +
      labs(title = trait, subtitle = model)
    ggsave(filename = paste0(plots_folder, "/histogram_marker_effects.", model, ".pdf"),
           plot = plot_hist_effects, device = "pdf", width = 14, height = 10)

    # plot boxplot of effects
    plot_boxplot_effects <- ggplot(marker_effects) +
      geom_boxplot(aes(x = env, y = effect)) +
      labs(title = trait, subtitle = model)
    ggsave(filename = paste0(plots_folder, "/boxplot_marker_effects.", model, ".pdf"),
           plot = plot_boxplot_effects, device = "pdf", width = 14, height = 10)

    # append results to main df
    marker_effects_all <- rbind(marker_effects_all, data.frame(marker_effects, model = model))
    gebvs_all <- rbind(gebvs_all, data.frame(gebvs, model = model))
    
  }
  rm(marker_effects, gebvs, plot_eff_per_marker, plot_hist_effects, plot_boxplot_effects)
  
  # transform df with effects to wide format to make it easier to calculate correlations
  marker_effects_all <- pivot_wider(marker_effects_all, names_from = "model", values_from = "effect")
  # calculate correlations of effects among different models
  cor_effects <- round(cor(as.matrix(marker_effects_all[, models]), use = "complete.obs"), digits = 2)
  cor_effects <- reshape2::melt(cor_effects)
  # plot correlation among models
  plot_cor_effects <- ggplot(data = cor_effects, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() +
    geom_text(aes(Var1, Var2, label = value), color = "black", fontface = "bold") +
    labs(title = trait, x = "", y = "") +
    scale_fill_gradient2(low = "#66ADE5", mid = "#FFC465", high = "#BF1B0B", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Pearson\nCorrelation") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  ggsave(filename = paste0(plots_folder, "/cor_effects_models.pdf"),
         plot = plot_cor_effects, device = "pdf")
  
  # read file with genotype means for each environment
  blues_file <- paste0("data/1stStage_BLUEs.", trait, "-per-env.txt")
  means <- fread(blues_file, header = TRUE, data.table = FALSE)
  means <- data.frame(means, model = "real_pheno")
  colnames(means) <- colnames(gebvs_all)
  # correct env names for compatibility
  means$env <- gsub("-", ".", means$env)
  # append means to gebv df
  gebvs_all <- rbind(gebvs_all, means)
  rm(means)
  
  # compare distribution of gebvs with real data
  plot_boxplot_gebvs <- ggplot(gebvs_all) +
    geom_boxplot(aes(x = env, y = gebv, fill = model)) +
    labs(title = trait)
  ggsave(filename = paste0(plots_folder, "/distribution_gebvs_per_model.pdf"),
         plot = plot_boxplot_gebvs, device = "pdf", width = 14, height = 10)
  
  # create a temporary df with order number per group
  gebvs_tmp <- gebvs_all %>% 
    group_by(env, model) %>%
    arrange(gebv, .by_group = TRUE) %>%
    mutate(order = row_number()) %>%
    ungroup()
  
  for (environment in unique(gebvs_tmp$env)) {
    
    # subset by env
    gebvs_tmp_env <- as.data.frame(subset(gebvs_tmp, env == environment))
    # get order based on 'real_pheno' values only
    geno_order <- gebvs_tmp_env[gebvs_tmp_env$model == "real_pheno", "genotype"]
    # create empt df to stored reordered values
    gebvs_env_reordered <- data.frame(stringsAsFactors = FALSE)
    for (model_eff in unique(gebvs_tmp_env$model)) {
      # rename order number to genotype name in real phenotypes
      gebvs_tmp_env_model <- gebvs_tmp_env[gebvs_tmp_env$model == model_eff, ]
      gebvs_tmp_env_model <- gebvs_tmp_env_model[match(geno_order, gebvs_tmp_env_model[, "genotype"]), ]
      # add real phenotyps in df to appear in the same facet as each model
      gebvs_tmp_env_real <- gebvs_tmp_env[gebvs_tmp_env$model == "real_pheno", ]
      # add facet column
      gebvs_tmp_env_real$facet <- model_eff
      gebvs_tmp_env_model$facet <- model_eff
      # append to df with other models
      gebvs_env_reordered <- rbind(gebvs_env_reordered, rbind(gebvs_tmp_env_model, gebvs_tmp_env_real))
    }
    rm(gebvs_tmp_env, gebvs_tmp_env_model, gebvs_tmp_env_real)
    
    # reorder levels
    gebvs_env_reordered$genotype <- factor(gebvs_env_reordered$genotype, levels = geno_order)
    # remove facet with real pheno only
    gebvs_env_reordered <- subset(gebvs_env_reordered, facet != "real_pheno")
    # plot
    plot_gebv_reorder <- ggplot(gebvs_env_reordered) +
      geom_line(aes(x = genotype, y = gebv, color = model, group = model)) +
      facet_wrap(~ facet) +
      labs(title = trait, subtitle = environment) +
      scale_color_manual(values = c("rrblup" = "black", "bayescpi" = "black", "mrr" = "black",
                                    "gwas" = "black", "real_pheno" = "#BF1B0B"),
                         breaks = "real_pheno") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.title = element_blank())
    ggsave(filename = paste0(plots_folder, "/gebvs-vs-real_per-model.", environment, ".pdf"),
           plot = plot_gebv_reorder, device = "pdf")
    
  }
  rm(gebvs_tmp, environment, geno_order, gebvs_env_reordered, plot_gebv_reorder)
  
  # reformat df to plot correlations
  gebvs_all <- pivot_wider(gebvs_all, names_from = "model", values_from = "gebv")
  # calculate correlations among gebvs from different models
  cor_gebvs <- round(cor(as.matrix(gebvs_all[, c(models, "real_pheno")]), use = "complete.obs"), digits = 2)
  cor_gebvs <- reshape2::melt(cor_gebvs)
  # plot correlation among models
  plot_cor_gebvs <- ggplot(data = cor_gebvs, aes(x = Var1, y = Var2, fill = value)) + 
    geom_tile() +
    geom_text(aes(Var1, Var2, label = value), color = "black", fontface = "bold") +
    labs(title = trait, x = "", y = "") +
    scale_fill_gradient2(low = "#66ADE5", mid = "#FFC465", high = "#BF1B0B", 
                         midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Pearson\nCorrelation") +
    theme_minimal() +
    theme(panel.grid = element_blank())
  ggsave(filename = paste0(plots_folder, "/cor_gebvs_models.pdf"),
         plot = plot_cor_gebvs, device = "pdf")
}



#### debug ----

# results_folder <- "analysis/marker_effects"
# traits <- c("EHT", "Moisture", "PHT", "TWT", "YLD")
# models <- c("rrblup", "bayescpi", "mrr", "gwas")
# no_missing_genotypes <- FALSE
