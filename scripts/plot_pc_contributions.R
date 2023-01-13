library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)

usage <- function() {
  cat("
description: plot the top covariables that contributes the most to a principal component.

usage: Rscript plot_pc_contributions.R [pca_contrib_file] [output_folder] [...]

positional arguments:
  pca_contrib_file            PCA loadings of env covariables
  output_folder               folder to output plots

optional argument:
  --help                      show this helpful message
  --n-covariables=VALUE       the number of top covariables to select (default: 10)
  --prop-covariables=VALUE    the proportion of top covariables to select
                              (default is to use '--n-covariables')


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
pca_contrib_file <- args[1]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default
n_covariables <- NULL
prop_covariables <- NULL

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {
  
  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--n-covariables", "--prop-covariables")
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
if (!is.null(n_covariables)) {
  if (suppressWarnings(!is.na(as.integer(n_covariables)))) {
    n_covariables <- as.integer(n_covariables)
    prop_covariables <- NULL
  } else {
    stop("Optional argument '--n-covariables' should be an integer")
  }
}

if (!is.null(prop_covariables)) {
  if (suppressWarnings(!is.na(as.numeric(prop_covariables)))) {
    prop_covariables <- as.numeric(prop_covariables)
    n_covariables <- NULL
  } else {
    stop("Optional argument '--n-covariables' should be a number between 0 and 1")
  }
  if (prop_covariables < 0 | prop_covariables > 1) {
    stop("Optional argument '--n-covariables' should be a number between 0 and 1")
  }
}

if (is.null(n_covariables) & is.null(prop_covariables)) {
  n_covariables <- 10
  prop_covariables <- NULL
}



#### plot PC contributions ----

# load pc contributions
pca_contrib <- fread(pca_contrib_file, header = TRUE, data.table = FALSE)

for (pc in colnames(pca_contrib)[-1]) {

  # get contributions to that pc
  contrib <- pca_contrib[, c("env_idx", pc)]
  colnames(contrib)[2] <- "contrib"
  # keep only top covariables
  if (!is.null(n_covariables)) top_covariables <- 1:n_covariables
  if (!is.null(prop_covariables)) top_covariables <- 1:ceiling(nrow(contrib) * prop_covariables)
  contrib <- contrib[order(contrib$contrib, decreasing = TRUE)[top_covariables], ]
  contrib$env_idx <- factor(contrib$env_idx, levels = contrib$env_idx[order(contrib$contrib, decreasing = TRUE)])

  # plot contributions
  plot_pc_contrib <- ggplot(contrib) +
    geom_col(aes(x = env_idx, y = contrib)) +
    labs(title = bquote("Contributions to" ~ bold(.(pc))),
         x = "Environmental covariables",
         y = "Contributions (%)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  if (!is.null(n_covariables)) outfile <- paste0(output_folder, "/", pc, ".top-n-", n_covariables, ".pdf")
  if (!is.null(prop_covariables)) outfile <- paste0(output_folder, "/", pc, ".top-prop-", prop_covariables, ".pdf")
  ggsave(plot_pc_contrib, filename = outfile, device = "pdf", height = 8, width = 0.2 * nrow(contrib))

}

# summarize contributions to PC in a bar plot (frequency of certain env idx)
# across categories (beginning, mid, end season) - also add total % variance
# explained per category

# split idx and interval - careful with some idx that has more than one "_" to split
pca_contrib[, c("idx", "interval")] <- do.call(rbind, lapply(pca_contrib$env_idx, function(idx) {
  idx <- unlist(strsplit(idx, split = "_"))
  if (length(idx) > 2) idx <- c(paste0(idx[1:(length(idx) - 1)], collapse = "_"), idx[length(idx)])
  return(idx)
}))

# reformat table
pca_contrib <- pca_contrib %>% 
  select(-env_idx) %>% 
  pivot_longer(-c(idx, interval), names_to = "pc", values_to = "contrib") %>% 
  mutate(interval = as.numeric(interval),
         season = case_when(interval < max(interval) / 3 ~ "early",
                            interval >= max(interval) / 3 & interval < (max(interval) / 3) * 2 ~ "mid",
                            interval > (max(interval) / 3) * 2 ~ "late"),
         season = factor(season, levels = c("early", "mid", "late"))) %>% 
  group_by(idx, pc, season) %>% 
  summarize(contrib = sum(contrib))

for (PC in unique(pca_contrib$pc)) {
  
  # create plot
  plot_summary_pc <- pca_contrib %>% 
    filter(pc == PC) %>% 
    group_by(season) %>%
    arrange(desc(contrib), .by_group = TRUE) %>% 
    ggplot() +
    facet_wrap(~ idx, nrow = 3) +
    geom_line(aes(x = season, y = contrib, group = 1), show.legend = FALSE) +
    labs(x = "Season stage", y = "Total contributions (%)") +
    theme_bw()
  # save plot
  outfile <- paste0(output_folder, "/summary_contributions1.", PC, ".pdf")
  ggsave(plot_summary_pc, filename = outfile, device = "pdf", height = 5, width = 8)
  
  # pca_contrib %>% 
  #   filter(pc == PC) %>% 
  #   group_by(season) %>%
  #   arrange(desc(contrib), .by_group = TRUE) %>% 
  #   ggplot() +
  #   facet_wrap(~ season, nrow = 3) +
  #   geom_col(aes(x = idx, y = contrib))
  
}


#### alternative way to plot PC contributions ----

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

for (PC in unique(pca_contrib$pc)) {
  
  plot_summary_pc2 <- pca_contrib %>%
    filter(pc == PC) %>%
    arrange(desc(contrib)) %>%
    ggplot() +
    facet_wrap(~ idx, nrow = 3) +
    geom_line(aes(x = interval, y = contrib, group = 1), show.legend = FALSE) +
    labs(title = PC, x = "Time interval", y = "Total contributions (%)") +
    theme_bw()
  # save plot
  outfile <- paste0(output_folder, "/summary_contributions2.", PC, ".pdf")
  ggsave(plot_summary_pc2, filename = outfile, device = "pdf", height = 5, width = 8)
  
  
  plot_summary_pc3 <- pca_contrib %>%
    filter(pc == PC) %>%
    arrange(desc(contrib)) %>%
    ggplot() +
    facet_wrap(~ idx, nrow = 3) +
    geom_point(aes(x = interval, y = contrib, group = 1), show.legend = FALSE) +
    labs(title = PC, x = "Time interval", y = "Total contributions (%)") +
    theme_bw()
  # save plot
  outfile <- paste0(output_folder, "/summary_contributions3.", PC, ".pdf")
  ggsave(plot_summary_pc3, filename = outfile, device = "pdf", height = 5, width = 8)
  

  plot_summary_pc4 <- pca_contrib %>%
    filter(pc == PC) %>%
    arrange(desc(contrib)) %>%
    ggplot(aes(x = interval, y = contrib)) +
    facet_wrap(~ idx, nrow = 3) +
    geom_point(size = 1, show.legend = FALSE) +
    geom_smooth() +
    labs(title = PC, x = "Time interval", y = "Total contributions (%)") +
    theme_bw()
  # save plot
  outfile <- paste0(output_folder, "/summary_contributions4.", PC, ".pdf")
  ggsave(plot_summary_pc4, filename = outfile, device = "pdf", height = 5, width = 8)
  
  
  plot_summary_pc5 <- pca_contrib %>%
    filter(pc == PC) %>%
    arrange(desc(contrib)) %>%
    ggplot(aes(x = interval, y = contrib)) +
    facet_wrap(~ idx, nrow = 3) +
    # geom_point(size = 1, show.legend = FALSE) +
    geom_smooth(se = TRUE) +
    labs(title = PC, x = "Time interval", y = "Total contributions (%)") +
    theme_bw()
  # save plot
  outfile <- paste0(output_folder, "/summary_contributions5.", PC, ".pdf")
  ggsave(plot_summary_pc5, filename = outfile, device = "pdf", height = 5, width = 8)
  
}




#### debug ----

# pca_contrib_file <- "data/env_covariables/pca_contributions.txt"
# output_folder <- "analysis/networks/YLD/pc_contributions"
# n_covariables <- NULL
# prop_covariables <- 0.1