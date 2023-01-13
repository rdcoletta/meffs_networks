library(data.table)
library(ggplot2)
library(doParallel)

usage <- function() {
  cat("
description: plot how many markers are missing per marker and per RIL, and also
             plot the distribution of these markers along the chromosome

usage: Rscript qc_hybrid_hmp.R [infile] [sv_list] [outfile]

positional arguments:
  infile            hapmap file
  outfile           folder to save results

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 2) stop(usage(), "missing positional argument(s)")

# get arguments
infile <- args[1]
outfolder <- args[2]
if (!dir.exists(outfolder)) dir.create(outfolder)
# infile <- "tests/usda_hybrids.all_markers.adjusted-n-markers.hmp.txt"
# outfolder <- "analysis/qc/hybrid_hmp"




#### distribution missing data ----

# load hmp file
hmp <- fread(infile, header = TRUE, data.table = FALSE)

cat("Calculating missing data per marker\n")

missingMarker <- function(hmp_row) {

  # get marker name
  marker_info <- as.character(hmp_row[c(1, 3, 4)])
  # get genotypes for that marker
  marker_geno <- as.character(hmp_row[12:length(hmp_row)])
  # calculate how much missing data that marker has
  missing <- round(sum(marker_geno == "NN") / length(marker_geno), digits = 4)

  # return info about markers
  return(c(marker_info, missing))

}

num_cores <- ifelse(detectCores() > 10, yes = 10, no = detectCores())
registerDoParallel(num_cores)
summary_markers <- foreach(chr = unique(hmp[, "chrom"]), .combine = rbind) %dopar% {

  # subset by chromosome
  hmp_chr <- subset(hmp, chrom == chr)
  # calculate missing data
  summary_chr <- apply(hmp_chr, MARGIN = 1, FUN = missingMarker)
  summary_chr <- data.frame(t(summary_chr), stringsAsFactors = FALSE)
  colnames(summary_chr) <- c(colnames(hmp_chr)[c(1, 3, 4)], "missing")
  # return list of dups for a chromosome
  summary_chr

}
stopImplicitCluster()

# adjust column type
summary_markers[, 2:4] <- apply(summary_markers[, 2:4], MARGIN = 2, as.numeric)
# plot missing data
marker_plot <- ggplot(data = summary_markers, aes(x = missing)) +
  geom_histogram(binwidth = 0.01) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = paste0("Percent missing data"),
       y = paste0("Number of markers")) +
  theme(text = element_text(size = 20))
# save plot
outfile_ril <- paste0(outfolder, "/missing_data_per_marker.pdf")
ggsave(filename = outfile_ril, plot = marker_plot, device = "pdf")


cat("Calculating missing data per RIL at...\n")

for (threshold in c(0.2, 0.3, 0.4, 0.5)) {

  cat(" ", threshold, "missing data threshold\n")

  # get marker names for that missing data threshold
  markers_to_keep <- summary_markers[summary_markers$missing <= threshold, 1]

  registerDoParallel(num_cores)
  summary_rils <- foreach(chr = unique(hmp[, "chrom"]), .combine = rbind) %dopar% {

    # subset by chromosome
    hmp_chr <- subset(hmp, chrom == chr)
    markers_to_keep_chr <- markers_to_keep[markers_to_keep %in% hmp_chr[, 1]]
    # keep only markers that pass missing threshold
    hmp_chr <- hmp_chr[hmp_chr[, 1] %in% markers_to_keep_chr, ]
    # calculate missing data
    missing_per_ril <- apply(hmp_chr[, 12:NCOL(hmp_chr)], MARGIN = 2, function(ril) {
      return(round(sum(ril == "NN") / length(ril), digits = 4))
    })
    # return list of dups for a chromosome
    missing_per_ril

  }
  stopImplicitCluster()

  # take mean missing data per chromosome
  summary_rils <- data.frame(colMeans(summary_rils))
  colnames(summary_rils) <- "missing"

  # plot missing data
  ril_plot <- ggplot(data = summary_rils, aes(x = missing)) +
    geom_histogram(binwidth = 0.01) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(subtitle = paste0("Markers with <= ", threshold, " missing data"),
         x = paste0("Percent missing data"),
         y = paste0("Number of hybrids")) +
    theme(text = element_text(size = 20))
  # save plot
  outfile_ril <- paste0(outfolder, "/missing_data_per_hybrid.threshold-", threshold, ".pdf")
  ggsave(filename = outfile_ril, plot = ril_plot, device = "pdf")

}

# remove data to save memory
rm(hmp, summary_rils, marker_plot, ril_plot)




#### distribution along chromosome ----

for (threshold in c(0.2, 0.3, 0.4, 0.5)) {

  cat(" ", threshold, "missing data threshold\n")

  # add column identifying which marker pass threshold
  summary_markers$filter <- ifelse(summary_markers$missing <= threshold,
                                        yes = "after", no = "before")
  summary_markers$filter <- factor(summary_markers$filter, levels = c("before", "after"))

  # plot distribution of markers along each chomosome
  dist_chr_plot <- ggplot() +
    geom_density(data = summary_markers,
                 mapping = aes(x = pos, fill = "black"),
                 alpha = 0.5, position = "identity") +
    geom_density(data = subset(summary_markers, filter == "after"),
                 mapping = aes(x = pos, fill = "firebrick"),
                 alpha = 0.3, position = "identity") +
    facet_wrap(~ chrom, scales = "free_x") +
    labs(subtitle = paste0("markers with <= ", threshold, " missing data"),
         x = "Position (Mb)",
         y = "Density") +
    scale_x_continuous(labels = function(x) x/1000000) +
    scale_fill_identity(name = "Missing data filter",
                        breaks = c("black", "firebrick"),
                        labels = c("before", "after"),
                        guide = "legend")+
    theme(text = element_text(size = 20))
  # save plot
  outfile_dist <- paste0(outfolder, "/markers_dist.missing-threshold-", threshold, ".pdf")
  ggsave(filename = outfile_dist, plot = dist_chr_plot, device = "pdf", width = 15, height = 9)

  # write marker names to keep
  outfile_markers <- paste0(outfolder, "/markers_to_keep.missing-threshold-", threshold, ".txt")
  markers_after_filter <- summary_markers[summary_markers$filter == "after", 1, drop = FALSE]
  fwrite(markers_after_filter, file = outfile_markers, quote = FALSE,
         sep = "\t", na = NA, row.names = FALSE, col.names = FALSE)

}
