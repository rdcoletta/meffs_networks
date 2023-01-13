#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: This script applies a sliding window approach to correct possible miscalls on RIL data.

Usage: Rscript sliding_window_approach.R [...]")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 6 arguments
if (length(args) != 6) {
  stop("incorrect number of arguments provided.

       Usage: Rscript sliding_window_approach.R [...]
       ")
}
cross <- args[1]
hmp.rils.filename <- args[2]
hmp.parents.filename <- args[3]

if (grepl("--window_size=", args[4])) {
  window.size <- as.numeric(unlist(strsplit(args[4], split = "="))[2])
} else {
  stop("Invalid argument 4")
}

if (grepl("--window_step=", args[5])) {
  window.step <- as.numeric(unlist(strsplit(args[5], split = "="))[2])
} else {
  stop("Invalid argument 5")
}

if (grepl("--min_snps_per_window=", args[6])) {
  min.snps.per.window <- as.numeric(unlist(strsplit(args[6], split = "="))[2])
} else {
  stop("Invalid argument 6")
}

# cross <- "B73xLH82"
# hmp.rils.filename <- "data/merged_hapmaps_by_cross/usda_SNPs-SVs_B73xLH82_RILs.sorted.hmp.txt"
# hmp.parents.filename <- "data/merged_hapmaps_by_cross/usda_SNPs-SVs_B73xLH82_parents.sorted.hmp.txt"
# window.size <- 10
# window.step <- 1
# min.snps.per.window <- 3




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("foreach")) install.packages("foreach")
if(!require("doParallel")) install.packages("doParallel")
if(!require("stringr")) install.packages("stringr")

CountMarkers <- function(list_individuals, genotype_data) {

  df.results <- data.frame(line = as.character(), total = as.numeric(),
                           one_allele_missing = as.numeric(), both_alleles_missing = as.numeric(),
                           percent_missing = as.numeric(), available = as.numeric(),
                           homo = as.numeric(), het = as.numeric(),  percent_het = as.numeric(),
                           stringsAsFactors = FALSE)

  for (individual in list_individuals) {

    count.two.missing <- 0
    count.one.missing <- 0
    count.homo <- 0
    count.het <- 0

    for (marker in 1:NROW(genotype_data)) {

      alleles <- as.character(genotype_data[marker, individual])

      # if there are missing allele(s)...
      if (grepl("N", alleles)) {
        if (str_count(alleles, pattern = "N") == 2) {
          # ...check if both are missing...
          count.two.missing <- count.two.missing + 1
        } else {
          # ...or just one
          count.one.missing <- count.one.missing + 1
        }
      } else {
        # if there both alleles are present, split...
        alleles <- unlist(strsplit(alleles, split = ""))
        if (alleles[1] == alleles[2]) {
          # ...and check if they are homozygous...
          count.homo <- count.homo + 1
        } else {
          # ...or het
          count.het <- count.het + 1
        }
      }

    }

    total.missing <- count.one.missing + count.two.missing
    total.markers <- count.homo + count.het + total.missing

    df.results <- rbind(df.results,
                        list(line = individual,
                             total = total.markers,
                             one_allele_missing = count.one.missing,
                             both_alleles_missing = count.two.missing,
                             percent_missing = round(total.missing/total.markers, digits = 3),
                             available = total.markers - total.missing,
                             homo = count.homo,
                             het = count.het,
                             percent_het = round(count.het / (total.markers - total.missing),
                                                 digits = 3)),
                        stringsAsFactors = FALSE)

  }

  return(df.results)

}



#### sliding window approach ----

cat("\nLoading cross ", cross, "\n", sep = "")


# get available cores for paralellizing
if (detectCores() < 10) {
  num.cores <- detectCores()
} else {
  num.cores <- 10
}


# load data
hmp.rils <- fread(hmp.rils.filename, header = TRUE, data.table = FALSE)
hmp.parents <- fread(hmp.parents.filename, header = TRUE, data.table = FALSE)

# get parents names
parent1 <- unlist(strsplit(cross, "x"))[1]
parent2 <- unlist(strsplit(cross, "x"))[2]


# get parents column numbers in resequencing data
p1.col.reseq <- grep(parent1, colnames(hmp.parents))
p2.col.reseq <- grep(parent2, colnames(hmp.parents))


# make sure that gbs data has the snps as reseq
if (all(hmp.rils[,1] != hmp.parents[,1])) stop("Data have different length")


# remove SVs (because RILs don't have info yet!)
hmp.parents.noSV <- hmp.parents[grep("^del|^dup|^ins|^inv|^tra", hmp.parents[, 1], perl = TRUE, invert = TRUE), ]
hmp.rils.noSV <- hmp.rils[grep("^del|^dup|^ins|^inv|^tra", hmp.rils[, 1], perl = TRUE, invert = TRUE), ]


# select only polymorphic (keep missing and mono in another data frame)
marker.type <- apply(X = hmp.parents.noSV[, c(p1.col.reseq, p2.col.reseq)],
                         MARGIN = 1, FUN = function(snp) {
                           # get unique genotypes between parents
                           genotypes <- unique(snp)
                           # genotypes <- genotypes[genotypes != "NN"]

                           if (any(grepl("NN", genotypes))) {

                             # if there is no genotype, snp is missing
                             return("missing")

                           } else if (length(genotypes) == 1) {

                             # if there is one genotype, it's monomorphic
                             # but distinguish if SNP is het
                             alleles <- unlist(strsplit(genotypes, split = ""))
                             if (alleles[1] == alleles[2]) {
                               return("mono")
                             } else {
                               return("het")
                             }

                           } else {

                             # if there are two genotypes, it's polymorphic
                             # but distiguish if one of the genotypes is het
                             p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
                             p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
                             if (p1.alleles[1] == p1.alleles[2] & p2.alleles[1] == p2.alleles[2]) {
                               return("poly")
                             } else {
                               return("het")
                             }

                           }
                         })

hmp.parents.noSV.poly <- hmp.parents.noSV[which(marker.type == "poly"), ]
hmp.rils.noSV.poly <- hmp.rils.noSV[which(marker.type == "poly"), ]


# sliding window approach
cat("Sliding window correction\n")
hmp.rils.outfile <- data.frame(stringsAsFactors = FALSE)
for (chr in 1:10) {

  cat("  chromosome ", chr, "\n", sep = "")
  # subset by chr
  hmp.parents.chr <- subset(hmp.parents.noSV.poly, chrom == chr)
  hmp.rils.chr <- subset(hmp.rils.noSV.poly, chrom == chr)

  registerDoParallel(cores = num.cores)
  hmp.rils.window <- foreach(ril.col=12:NCOL(hmp.rils.chr), .combine = cbind) %dopar% {

    # set up first window
    window.start <- 1
    window.stop <- window.start + (window.size - 1)

    # create a vector to store consenus genotype for each window
    window.consensus <- c()

    # use slide window approach until end of the window reaches the last SNP
    while (window.stop <= NROW(hmp.rils.chr)) {

      # get genotypes from parents and ril for that window
      window <- cbind(hmp.parents.chr[window.start:window.stop, p1.col.reseq],
                      hmp.parents.chr[window.start:window.stop, p2.col.reseq],
                      hmp.rils.chr[window.start:window.stop, ril.col])

      # define from which parents the SNP in ril comes from
      window.calls <- apply(window, MARGIN = 1, FUN = function(genotypes) {
        if (genotypes[3] == "NN") {
          # if ril snp is NN
          return("missing")
        } else if (genotypes[1] == genotypes[3]) {
          # if ril snp is the same as parent1
          return("p1")
        } else if (genotypes[2] == genotypes[3]) {
          # if ril snp is the same as parent2
          return("p2")
        } else {
          # check if ril snp is a het or if it has a different allele from parents
          p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
          p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
          ril.alleles <- unlist(strsplit(genotypes[3], split = ""))
          if (unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles & ril.alleles[1] != ril.alleles[2]) {
            return("het")
          } else {
            return("missing")
          }
        }
      })

      # check if there is enough ril snps genotyped
      if (sum(window.calls != "missing") >= min.snps.per.window) {

        # get number of alleles for each parent (multiply by 2 because a homozygous has 2 alleles)
        n.p1.alleles <- (sum(window.calls == "p1") * 2) + (sum(window.calls == "het"))
        n.p2.alleles <- (sum(window.calls == "p2") * 2) + (sum(window.calls == "het"))
        total.alleles <- sum(window.calls != "missing") * 2
        # from those not missing, what proportion is p1, p2 and het?
        prop.p1 <- n.p1.alleles/total.alleles

        # define consensus based on threshold (p1: p1>0.7, p2: p1<0.3, het: 0.3<p1<0.7)
        if (prop.p1 >= 0.7) {
          # assign parent1 genotype of first snp on window to consensus
          window.consensus <- append(window.consensus, window[1:window.step, 1])
        } else if (prop.p1 <= 0.3) {
          # assign parent2 genotype of first snp on window to consensus
          window.consensus <- append(window.consensus, window[1:window.step, 2])
        } else {
          # assign het genotype to consensus
          p1.allele <- sapply(window[1:window.step, 1], function(x) return(unlist(strsplit(x, split = ""))[1]))
          p2.allele <- sapply(window[1:window.step, 2], function(x) return(unlist(strsplit(x, split = ""))[1]))
          window.consensus <- append(window.consensus, paste0(p1.allele, p2.allele))
        }
      } else {
        # if there is very few or none ril SNPs genotyped, consider consensus as missing data
        window.consensus <- append(window.consensus, rep("NN", times = window.step))
      }

      # set up the start of next window
      window.start <- window.start + window.step
      window.stop <- window.start + (window.size - 1)
    }

    # after that, add (window.size - 1) NNs in the consensus
    window.consensus <- append(window.consensus, rep("NN", times = NROW(hmp.rils.chr) - length(window.consensus)))

    window.consensus

  }
  stopImplicitCluster()

  # correct column names
  colnames(hmp.rils.window) <- colnames(hmp.rils.chr)[12:NCOL(hmp.rils.chr)]

  # create hapmap again
  hmp.rils.outfile.chr <- cbind(hmp.rils.chr[, 1:11], hmp.rils.window, stringsAsFactors = FALSE)

  # bind to final data frame
  hmp.rils.outfile <- rbind(hmp.rils.outfile, hmp.rils.outfile.chr)
}

# add SVs and rest of SNPs to outfile
hmp.rils.outfile <- rbind(hmp.rils.outfile,
                          hmp.rils.noSV[which(marker.type != "poly"), ],
                          hmp.rils[grep("^del|^dup|^ins|^inv|^tra", hmp.rils[, 1], perl = TRUE), ])
hmp.rils.outfile <- hmp.rils.outfile[order(hmp.rils.outfile$chrom, hmp.rils.outfile$pos), ]

if (!all(hmp.rils.outfile[, 1] == hmp.rils[, 1])) stop("Order of outfile is different from original ril file")

# write file
outfile <- gsub(".hmp.txt", ".sliding-window.hmp.txt", hmp.rils.filename, fixed = TRUE)
fwrite(hmp.rils.outfile, outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)



#### QC ---

cat("Performing QC\n")
qc.before <- CountMarkers(list_individuals = colnames(hmp.rils)[12:NCOL(hmp.rils)],
                          genotype_data = hmp.rils)
hets.before <- round(mean(qc.before$percent_het) * 100, digits = 2)
missing.before <- round(mean(qc.before$percent_missing) * 100, digits = 2)

qc.after <- CountMarkers(list_individuals = colnames(hmp.rils.outfile)[12:NCOL(hmp.rils.outfile)],
                         genotype_data = hmp.rils.outfile)
hets.after <- round(mean(qc.after$percent_het) * 100, digits = 2)
missing.after <- round(mean(qc.after$percent_missing) * 100, digits = 2)

cat("  Hets: ", hets.before, "% (before) vs ", hets.after, "% (after)\n", sep = "")
cat("  Missing: ", missing.before, "% (before) vs ", missing.after, "% (after)\n", sep = "")



#### karyoplots ----

# library(ggplot2)
#
# chr.info <- "data/B73_RefGen_V4_chrm_info.txt"
# cent.info <- "data/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xLH82"
# proj.folder <- "analysis/projection_reseq-snps"
# parents.folder <- "data/reseq_snps"
# out.folder <- "analysis/qc/karyotypes"
# random.rils <- FALSE
# set.seed(1626)
# RIL <- sample(colnames(hmp.rils.outfile)[12:NCOL(hmp.rils.outfile)], size = 1, replace = FALSE)
#
#
# # chromosomes
# chrms <- fread(chr.info, header = TRUE, data.table = FALSE)
# chrms <- data.frame(chr = chrms$chr, start_pos = 0, end_pos = chrms$length)
#
# # centromeres
# centros <- read.delim(cent.info, sep = "\t", header = TRUE)
# centros <- data.frame(chr = centros$chr, start_pos = centros$start_pos, end_pos = centros$end_pos)
#
#
# plot_karyotype <- function(hmp.rils, hmp.parents, RIL) {
#
#   cat("Plotting ", cross, "...\n", sep = "")
#
#   # remove SVs
#   geno.data.cross <- hmp.rils[grep("^del|^dup|^ins|^inv|^tra", hmp.rils[, 1], perl = TRUE, invert = TRUE), ]
#   parents.data.cross <- hmp.parents[grep("^del|^dup|^ins|^inv|^tra", hmp.parents[, 1], perl = TRUE, invert = TRUE), ]
#
#   # get parents names
#   parent1 <- unlist(strsplit(cross, "x"))[1]
#   parent2 <- unlist(strsplit(cross, "x"))[2]
#
#   # get parents column numbers in resequencing data
#   p1.col.reseq <- grep(parent1, colnames(parents.data.cross), fixed = TRUE)
#   p2.col.reseq <- grep(parent2, colnames(parents.data.cross), fixed = TRUE)
#
#   # get ril column number
#   ril.col <- grep(RIL, colnames(geno.data.cross))
#
#   cat("  RIL", RIL, "\n")
#
#   # merge information of RIL of interest with respective marker positions
#   geno.data <- cbind(geno.data.cross[, c(1,3,4)], geno.data.cross[, ril.col])
#   colnames(geno.data) <- c("marker", "chr", "pos", "geno")
#
#   # select only chip SNPs
#   parents.SNPs <- parents.data.cross[grep("^snp", parents.data.cross[, 1], perl = TRUE, invert = TRUE), ]
#   geno.data.SNPs <- geno.data[grep("^snp", geno.data$marker, perl = TRUE, invert = TRUE), ]
#   # add parental info
#   geno.data.SNPs <- cbind(geno.data.SNPs, parents.SNPs[, c(p1.col.reseq, p2.col.reseq)])
#   colnames(geno.data.SNPs) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")
#
#   # select missing data and non-missing data
#   geno.data.SNPs.not.missing <- subset(geno.data.SNPs, geno != "NN")
#   # get proportion of non-missing data to add in the plot
#   prop.SNPs.not.missing <- NROW(geno.data.SNPs.not.missing) / NROW(geno.data.SNPs)
#
#   # check from which parent an allele came
#   geno.data.SNPs.not.missing.p1 <- geno.data.SNPs.not.missing[which(geno.data.SNPs.not.missing[, "geno"] == geno.data.SNPs.not.missing[, "parent1"] &
#                                                                       geno.data.SNPs.not.missing[, "geno"] != geno.data.SNPs.not.missing[, "parent2"]), ]
#   geno.data.SNPs.not.missing.p2 <- geno.data.SNPs.not.missing[which(geno.data.SNPs.not.missing[, "geno"] != geno.data.SNPs.not.missing[, "parent1"] &
#                                                                       geno.data.SNPs.not.missing[, "geno"] == geno.data.SNPs.not.missing[, "parent2"]), ]
#   geno.data.SNPs.not.missing.rest <- geno.data.SNPs.not.missing[which(!geno.data.SNPs.not.missing[, "marker"] %in% geno.data.SNPs.not.missing.p1[, "marker"] &
#                                                                         !geno.data.SNPs.not.missing[, "marker"] %in% geno.data.SNPs.not.missing.p2[, "marker"]), ]
#
#   # find hets
#   geno.data.SNPs.not.missing.het <- data.frame(matrix(nrow = 0, ncol = NCOL(geno.data.SNPs.not.missing.rest)))
#   colnames(geno.data.SNPs.not.missing.het) <- colnames(geno.data.SNPs.not.missing.rest)
#   for (snp in 1:NROW(geno.data.SNPs.not.missing.rest)) {
#     # check if ril snp is a het or if it has a different allele from parents
#     p1.alleles <- unlist(strsplit(as.character(geno.data.SNPs.not.missing.rest[snp, "parent1"]), split = ""))
#     p2.alleles <- unlist(strsplit(as.character(geno.data.SNPs.not.missing.rest[snp, "parent2"]), split = ""))
#     ril.alleles <- unlist(strsplit(as.character(geno.data.SNPs.not.missing.rest[snp, "geno"]), split = ""))
#     # make sure to exclude hets in parental calls
#     if (length(unique(p1.alleles)) == 1 & length(unique(p2.alleles)) == 1) {
#       if (ril.alleles[1] != ril.alleles[2] & unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles) {
#         geno.data.SNPs.not.missing.het <- rbind(geno.data.SNPs.not.missing.het, geno.data.SNPs.not.missing.rest[snp, ])
#       }
#     }
#   }
#
#
#   # plot karyotypes before SV projection (i.e. only parental SNPs)
#   karyo.plot <- ggplot() +
#     geom_segment(data = chrms,
#                  aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
#                  lineend = "round", color = "Gainsboro", size = 5) +
#     geom_segment(data = centros,
#                  aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
#                  lineend = "round", color = "DimGray", size = 5) +
#     geom_segment(data = geno.data.SNPs.not.missing.p1,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
#                  lineend = "butt", color = "firebrick", size = 5, alpha = 0.3) +
#     geom_segment(data = geno.data.SNPs.not.missing.p2,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#386cb0", size = 5, alpha = 0.3) +
#     geom_segment(data = geno.data.SNPs.not.missing.het,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#fec44f", size = 5, alpha = 0.5) +
#     scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
#     scale_x_discrete(position = "top") +
#     theme_minimal() +
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           plot.caption = element_text(size = rel(1.1), color = "DimGray"),
#           axis.text = element_text(size=rel(2)),
#           axis.title = element_text(size=rel(2)),
#           strip.text.x = element_text(size=rel(2))) +
#     facet_grid(~chr, switch = "y") +
#     labs(caption = paste0(cross, " - ", gsub("RIL_", "RIL ", RIL), "\n\n",
#                           "Not missing: ", round(prop.SNPs.not.missing, digits = 3), "\n"),
#          x = "Chromosomes", y = "Genomic positions (Mb)")
#
#   return(karyo.plot)
#
# }
#
#
# karyo.plot.before <- plot_karyotype(hmp.rils, hmp.parents, RIL)
# karyo.plot.after <- plot_karyotype(hmp.rils.outfile, hmp.parents, RIL)
#
#
# # save plots
# karyo.name.before <- paste0(out.folder, "/", cross, "_", RIL,"_before-sliding-window.pdf")
# karyo.name.after <- paste0(out.folder, "/", cross, "_", RIL,"_after-sliding-window.pdf")
# ggsave(filename = karyo.name.before, plot = karyo.plot.before, device = "pdf")
# ggsave(filename = karyo.name.after, plot = karyo.plot.after, device = "pdf")




#### credits ----

# code for karyotype plot was adapted from Carles Hernandez-Ferrer's blog:
# https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
# assign arguments to variables
