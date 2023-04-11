library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot LD decay from PLINK files

usage: Rscript plot_ld_decay.R [ld_filename] [outfile] [...]

positional arguments:
  ld_filename             input ld file from PLINK
  outfile                 name of ld decay plot

optional argument:
  --help                  show this helpful message
  --decay-plot=VALUE      type of plot to show LD decay. Options available are 'bins' (default), 'average' or 'raw'
                          - 'bins' will plot distribution of ld values across different window sizes (i.e. smaller windows for
                            closer markers and larger windows for distant markers)
                          - 'average' will plot average ld from 100bp windows (window size and window step can be changed with
                            the options '--window-size' and '--window-step')
                          - 'raw' will plot all ld values by physical distance
  --window-size=VALUE     control window size in bp (default: 100)
  --window-step=VALUE     control window step in bp (default: 100)
  --background-ld=VALUE   value for background LD (default: 0.2)

"
  )
}


getArgValue <- function(arg) {

  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])

}

avgR2byWindow <- function(data, col_r2, col_dist, window_size, window_step, window_start, all_points = FALSE) {

  window_stop <- window_start + (window_size - 1)

  # set df to store average of each window
  df_window_avgs <- data.frame(stringsAsFactors = FALSE)

  while (window_stop < max(data[, col_dist])) {

    window_avg <- data[data[, col_dist] > window_start & data[, col_dist] <= window_stop, c(col_r2, col_dist)]

    if (all_points) {
      # get all data points for each window
      window_avg <- data.frame(window_avg, n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
    } else {
      # get r2 average for each window
      window_avg <- data.frame(t(colMeans(window_avg)), n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
    }

    df_window_avgs <- rbind(df_window_avgs, window_avg)

    # set up the start of next window
    window_start <- window_start + window_step
    window_stop <- window_start + (window_size - 1)

  }

  # get all data points for last window
  window_avg <- data[data[, col_dist] > window_start & data[, col_dist] <= window_stop, c(col_r2, col_dist)]
  window_avg <- data.frame(window_avg, n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
  df_window_avgs <- rbind(df_window_avgs, window_avg)

  return(df_window_avgs)

}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# get positional arguments
ld_file <- args[1]
outfile <- args[2]

# set default of optional args
decay_plot <- "bins"
window_size <- "100"
window_step <- "100"
background_ld <- "0.2"

# assert to have the correct optional arguments
pos_args <- 2
if (length(args) < pos_args) stop(usage(), "missing positional argument(s)")

if (length(args) > pos_args) {

  opt_args <- args[-1:-pos_args]
  opt_args_allowed <- c("--decay-plot", "--window-size", "--window-step", "--background-ld")
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
if (!decay_plot %in% c("bins", "average", "raw")) {
  stop("Optional argument '--decay-plot' should be 'bins', 'average' or 'raw'")
}

if (suppressWarnings(!is.na(as.integer(window_size)))) {
  window_size <- as.integer(window_size)
} else {
  stop("Optional argument '--window-size' should be an integer")
}

if (suppressWarnings(!is.na(as.integer(window_step)))) {
  window_step <- as.integer(window_step)
} else {
  stop("Optional argument '--window-step' should be an integer")
}

if (suppressWarnings(!is.na(as.numeric(background_ld)))) {
  background_ld <- as.numeric(background_ld)
} else {
  stop("Optional argument '--background-ld' should be a number")
}



#### plot ld decay ----

# load data
ld_results <- fread(ld_file, header = TRUE, data.table = FALSE)

# make sure R2 is numeric
ld_results$R2 <- as.numeric(ld_results$R2)

# get distances between SNPs
ld_results$dist_markers <- abs(as.numeric(ld_results[, "BP_B"]) - as.numeric(ld_results[, "BP_A"]))

if (decay_plot == "bins") {

  # set default window parameters -- max distance possible is 100mb
  window_parameters <- data.frame(bp_start =    c(   1,  1001, 10001,  50001,  100001,  1000001,  10000000),
                                  bp_stop =     c(1000, 10000, 50000, 100000, 1000000, 10000000, 100000000),
                                  window_size = c( 100,  1000, 10000,  50000,  100000,  1000000,  10000000))

  # adjust parameters based on maximum distance between markers
  window_parameters <- window_parameters[which(window_parameters[, "bp_start"] < max(ld_results[, "dist_markers"])), ]

  df_plot <- data.frame(stringsAsFactors = FALSE)
  for (row in 1:NROW(window_parameters)) {

    window_start <- window_parameters[row , "bp_start"]
    window_stop <- window_parameters[row, "bp_stop"]
    window_size <- window_parameters[row , "window_size"]

    # get ld for window
    ld_window <- ld_results[which(ld_results[, "dist_markers"] >= window_start & ld_results[, "dist_markers"] < window_stop), c("R2", "dist_markers")]

    df_window <- avgR2byWindow(ld_window, "R2", "dist_markers",
                               window_size = window_size,
                               window_step = window_size,
                               window_start = window_start,
                               all_points = TRUE)

    df_plot <- rbind(df_plot, df_window)

  }

  # adjust scale for x axis
  old_levels <- as.integer(levels(as.factor(df_plot[, "bp_stop"])))
  new_levels <- c()
  for (i in 1:length(old_levels)) {

    if (i == 1) {
      new_levels <- c(new_levels, paste0("0\n|\n", old_levels[1] / 1000))
    } else {
      new_levels <- c(new_levels, paste0(old_levels[i - 1] / 1000, "\n|\n", old_levels[i] / 1000))
    }

  }

  n_per_bin <- function(x) {
    return(c(y = 1.02, label = length(x)))
  }

  # plot decay across unequal windows
  out_plot <- ggplot(df_plot, aes(x = as.factor(bp_stop), y = R2)) +
    geom_boxplot(width = 0.2) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(labels = new_levels) +
    stat_summary(fun.data = n_per_bin, geom = "text", size = 3) +
    geom_boxplot(fill = "gray70", width = 0.8) +
    geom_hline(yintercept = background_ld, color = "firebrick", linetype = "dashed") +
    labs(x = "Physical distance (kb)",
         y = bquote("LD"~(r^2))) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 18))

} else if (decay_plot == "average") {

  cat("Running rolling window")

  # get r2 average for each window
  df_window_avgs <- avgR2byWindow(ld_results, "R2", "dist_markers", window_size, window_step, window_start = 1)

  # plot average r2
  out_plot <- ggplot(df_window_avgs, aes(x = dist_markers, y = R2)) +
    geom_point(size = 0.2, alpha = 0.2)  +
    geom_smooth() +
    geom_hline(yintercept = background_ld, color = "firebrick", linetype = "dashed") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(labels = function(x) x/1000) +
    labs(title = "LD decay",
         x = "Physical distance (kb)",
         y = bquote("LD"~(r^2)))

} else if (decay_plot == "raw") {

  # plot ld decay without averaging r2 values
  out_plot <- ggplot(ld_results, aes(x = dist_markers, y = R2)) +
    geom_point(size = 0.2, alpha = 0.2)  +
    geom_smooth() +
    geom_hline(yintercept = background_ld, color = "firebrick", linetype = "dashed") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(labels = function(x) x/1000) +
    labs(title = "LD decay",
         x = "Physical distance (kb)",
         y = bquote("LD"~(r^2)))


}

# save plot
bitmap(file = outfile, type = "png16m", res = 600, width = 20)
print(out_plot)
dev.off()


#### debug ----

# ld_file <- "analysis/ld_decay/hybrids_geno.ld.gz"
# outfile <- "analysis/ld_decay/ld_decay.bins.png"
# decay_plot <- "bins"
# window_size <- 100
# window_step <- 100
# background_ld <- 0.16
