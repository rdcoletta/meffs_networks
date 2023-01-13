#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script will filter the hapmap files of parents and RILs to have only data for each
                   biparental cross, and then run rqtl to estimate the recombination frequency and plot them
                   against the physical positions for each chromosome.

      Usage: Rscript recomb-freq_biparental-crosses.R [infile_parents] [infile_rils] [cross_info] [qc_folder] [out_dir]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 5) {
  stop("incorrect number of arguments provided.

       Usage: Rscript recomb-freq_biparental-crosses.R [infile_parents] [infile_rils] [cross_info] [qc_folder] [out_dir]
       ")
}

# assign arguments to variables
infile.parents <- args[1]
infile.rils <- args[2]
cross.info <- args[3]
qc.folder <- args[4]
out.dir <- args[5]


# infile.parents <- "data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.hmp.txt"
# infile.rils <- "data/usda_22kSNPs_rils.sorted.diploid.hmp.txt"
# cross.info <- "data/usda_biparental-crosses.txt"
# qc.folder <- "analysis/qc"
# out.dir <- "data/biparental-crosses"



#### libraries used ----

if(!require("data.table")) install.packages("data.table")
if(!require("tibble")) install.packages("tibble")
if(!require("qtl")) install.packages("qtl")
if(!require("dplyr")) install.packages("dplyr")
if(!require("tidyr")) install.packages("tidyr")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("grid")) install.packages("grid")
if(!require("gridExtra")) install.packages("gridExtra")



#### functions ----

Hmp2Rqtl <- function(parents, rils, out.dir) {

  # combine parent and ril data
  combined.hmp <- cbind(parents, rils[,12:NCOL(rils)])
  # keep only marker name, chromosome, position, parents and rils
  combined.hmp <- combined.hmp[,-c(2,5:11)]

  # remove monomorphic markers
  polymorphic.markers <- c()
  for (row in 1:NROW(combined.hmp)) {
    # return TRUE if the genotypes from parents are different
    if (combined.hmp[row, 4] != combined.hmp[row, 5]) {
      polymorphic.markers <- append(polymorphic.markers, row)
    }
  }
  combined.hmp.filter <- combined.hmp[polymorphic.markers,]

  # remove missing data on parents
  missing.geno.parents <- c()
  for (row in 1:NROW(combined.hmp.filter)) {
    # return TRUE if there is a "-" on at least one of the parents genotypes
    if (grepl("N", combined.hmp.filter[row,4]) || grepl("N", combined.hmp.filter[row,5])) {
      missing.geno.parents <- append(missing.geno.parents, row)
    }
  }
  combined.hmp.filter <- combined.hmp.filter[-missing.geno.parents,]

  # remove markers with at least one heterozygous parent
  het.parents <- c()
  for (row in 1:NROW(combined.hmp.filter)) {
    # get parents' genotype
    parentA.geno <- strsplit(combined.hmp.filter[row,4], split = "")[[1]]
    parentB.geno <- strsplit(combined.hmp.filter[row,5], split = "")[[1]]
    # if one of the parents is het, append to vector
    if (parentA.geno[1] != parentA.geno[2] || parentB.geno[1] != parentB.geno[2]) {
      het.parents <- append(het.parents, row)
    }
  }
  combined.hmp.filter <- combined.hmp.filter[-het.parents,]

  # code parent 1 geno as A, parent 2 geno as B, and het as H
  for (row in 1:NROW(combined.hmp.filter)) {
    for (ril in 6:NCOL(combined.hmp.filter)) {
      # get genotypes
      parentA.geno <- strsplit(combined.hmp.filter[row,4], split = "")[[1]]
      parentB.geno <- strsplit(combined.hmp.filter[row,5], split = "")[[1]]
      ril.geno <- strsplit(combined.hmp.filter[row,ril], split = "")[[1]]
      # since parents are polymorphic, parent A and B allele will always be different
      hetAB <- c(parentA.geno[1],parentB.geno[1])
      hetBA <- c(parentB.geno[1],parentA.geno[1])

      # transform genotype according to condition
      if (all(ril.geno == parentA.geno)) {
        combined.hmp.filter[row,ril] <- "A"
      }
      if (all(ril.geno == parentB.geno)) {
        combined.hmp.filter[row,ril] <- "B"
      }
      if (all(ril.geno == c("N","N"))) {
        combined.hmp.filter[row,ril] <- "N"
      }
      if (all(ril.geno == hetAB) || all(ril.geno == hetBA)) {
        combined.hmp.filter[row,ril] <- "H"
      }

    }
  }

  # remove parent data
  combined.hmp.filter <- combined.hmp.filter[,-c(4,5)]

  # format df to meet rqtl requirements for input
  input.rqtl <- t(combined.hmp.filter)
  input.rqtl <- rownames_to_column(as.data.frame(input.rqtl))
  input.rqtl[1,1] <- ""
  input.rqtl[2,1] <- ""
  input.rqtl[3,1] <- ""
  names(input.rqtl) <- NULL

  # write table
  write.csv(input.rqtl, paste0(out.dir, "/usda_22kSNPs_", cross, "_rqtl-format.csv"),
            quote = FALSE, row.names = FALSE)
}


EstimateRecombinationFreq <- function(cross, dir.csv, dir.qc, plot = FALSE) {

  # print message about the cross being analyzed
  cat("\n\n---------------------------\n",
      paste("Analyzing cross:", cross),
      "\n---------------------------\n\n")

  # load cross
  rqtl <- read.cross(format = "csv", dir = dir.csv,
                     file = paste0("usda_22kSNPs_", cross, "_rqtl-format.csv"),
                     crosstype = "riself", estimate.map = FALSE)

  # print initial summary of the cross
  summary(rqtl)
  # write initial summary of the cross
  summary.outfile.initial <- paste0(dir.qc, "/", cross, "/summary_", cross, "_rils.txt")
  capture.output(summary(rqtl), file = summary.outfile.initial, type = "output")

  # remove individual if half of the markers are missing
  rqtl <- subset(rqtl, ind = (ntyped(rqtl) > totmar(rqtl) / 2))
  cat(paste("\nIndividuals with lots of markers missing removed:",
              sum(ntyped(rqtl) <= totmar(rqtl) / 2)))

  # remove markers missing in half of individuals
  n.genotyped.markers <- ntyped(rqtl, "mar")
  to.drop <- names(n.genotyped.markers[n.genotyped.markers < length(ntyped(rqtl)) / 2])
  rqtl <- drop.markers(rqtl, to.drop)
  to.drop.filename <- paste0(dir.qc, "/", cross, "/markers-missing_", cross, ".txt")
  write.table(data.frame(to.drop), file = to.drop.filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(paste("\nMarkers missing in many individuals removed:", length(to.drop)))

  # identify duplicate individuals
  prop.matching.ind <- comparegeno(rqtl)
  # par(mfrow = c(1, 1), las = 1)
  # hist(prop.matching.ind[lower.tri(prop.matching.ind)], breaks = seq(0, 1, len = 101),
  #      xlab = "No. matching genotypes")
  # rug(prop.matching.ind[lower.tri(prop.matching.ind)])

  # remove duplicate individuals
  same.ind <- which(prop.matching.ind > 0.8, arr = TRUE)
  same.ind <- same.ind[same.ind[, 1] < same.ind[, 2], ]
  if (length(same.ind) > 0) {
    # make sure "same.ind" is a matrix, because if "length(same.ind) == 1" it would be a vector
    same.ind <- matrix(same.ind, ncol =  2)
    rqtl <- subset(rqtl, ind = -same.ind[, 2])
    cat(paste("\nIndividuals with same genotypes removed:", length(same.ind[, 2])))
  }

  # remove duplicate markers
  dup.markers <- findDupMarkers(rqtl, exact.only = TRUE, adjacent.only = TRUE)
  rqtl <- drop.markers(rqtl, unlist(dup.markers))
  cat(paste("\nDuplicated markers removed:", length(unlist(dup.markers))))

  # markers with distorted segregation
  geno.table.chi.sq <- geno.table(rqtl)
  # adjust p-values for multiple testing with FDR
  geno.table.chi.sq[, "p.val.adj"] <- p.adjust(geno.table.chi.sq$P.value, method = "fdr")
  to.drop <- rownames(geno.table.chi.sq[geno.table.chi.sq$p.val.adj < 0.05, ])
  rqtl <- drop.markers(rqtl, to.drop)
  to.drop.filename <- paste0(dir.qc, "/", cross, "/markers-seg-distortion_", cross, ".txt")
  write.table(data.frame(to.drop), file = to.drop.filename, quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(paste("\nMarkers with segregation distortion removed:", length(to.drop)))

  # plot genotype frequencies by individual
  geno.data <- pull.geno(rqtl)
  geno.freq <- apply(geno.data, 1, function(a) table(factor(a, levels = 1:2)))
  geno.freq <- as.data.frame(t(geno.freq) / colSums(geno.freq))
  # format before ploting
  colnames(geno.freq) <- c("AA", "BB")
  geno.freq.plot <- geno.freq %>%
    rownames_to_column(var = "RIL") %>%
    transform(RIL = as.numeric(RIL)) %>%
    gather(key = "genotype", value = "genotype_freq", -RIL) %>%
    ggplot(aes(x = RIL, y = genotype_freq, color = genotype)) +
    geom_point(size = 3, show.legend = FALSE) +
    facet_wrap(~ genotype) +
    ylim(c(0,1)) +
    labs(title = "Genotype frequency by individual RIL", x = "RIL", y = "Genotype Frequency") +
    scale_color_manual(values = c("#053061", "#b2182b"))
  # save figure
  geno.freq.figure <- paste0(dir.qc, "/", cross, "/geno-freq_", cross, "_by_individual_rils.pdf")
  ggsave(filename = geno.freq.figure, plot = geno.freq.plot, device = "pdf")

  # save the genotypic data after filtering
  geno.data.filtered <- rownames_to_column(as.data.frame(t(geno.data)))
  # change column names
  ril.names <- sprintf("RIL_%d", 1:length(ntyped(rqtl)))
  colnames(geno.data.filtered) <- c("marker", ril.names)
  # use AA and BB instead of 1 and 2
  geno.data.filtered[geno.data.filtered == 1] <- "AA"
  geno.data.filtered[geno.data.filtered == 2] <- "BB"
  # write table
  geno.filter.name <- paste0(dir.qc, "/", cross, "/geno-data_", cross, "_after_filtering.txt")
  write.table(geno.data.filtered, geno.filter.name, sep = "\t", quote = FALSE, row.names = FALSE)

  # estimate recombination frequency
  rqtl <- est.rf(rqtl)
  # check for switched alleles
  cat("\nChecking for switched alleles...")
  checkAlleles(rqtl, threshold = 5)

  # retrieve recombination frequency and physical positions and add them in the same df
  df.gen.phys.positions <- data.frame(chr = as.numeric(), pos = as.numeric(), rf = as.numeric(),
                                      pop_size = as.numeric())

  # chrm names are strings, and they have a white space for chrm <10
  for (chr in c(" 1", " 2", " 3", " 4", " 5", " 6", " 7", " 8", " 9", "10")) {
    # retrieve recombination frequency for a chromosome
    rf <- pull.rf(rqtl, what = "rf", chr = chr)
    # since markers are already ordered by physical positions, retrieve the recombination frequency
    # of each allele relatively to its neighbor (i.e., upper off-diagonal)
    rf <- rf[row(rf) == (col(rf) + 1)]
    # but code the above doesn't capture the first marker (i.e., recombination frequency of the
    # first marker with itself, which is 0), so i will do it manually
    rf <- c(0, rf)
    # finally, i need to calculate the cumulative sum of the markers to get the total genetic distance
    rf <- cumsum(rf)

    # retrieve information of the physical positions of markers in a chromosome
    pos.chrm <- pull.map(rqtl, chr = chr, as.table = TRUE)

    # add info about how many individual there in the population
    pop_size <- NROW(rqtl$pheno)

    # make sure the marker names match in both rf and physical positions
    if (all(names(rf) == rownames(pos.chrm))) {
      for (i in 1:length(rf)) {
        # add them into the same df
        df.gen.phys.per.chr <- cbind(pos.chrm, rf, pop_size)
      }
    } else {
      cat("\n'rf' and 'pos.chrm' have different marker names")
    }

    # combine information from multiple chromosomes
    df.gen.phys.positions <- rbind(df.gen.phys.positions, df.gen.phys.per.chr)
  }

  # transform row name into a column
  df.gen.phys.positions <- rownames_to_column(df.gen.phys.positions, var = "marker")
  # erase white space from chromosome name
  df.gen.phys.positions[, "chr"] <- gsub(" ", "", df.gen.phys.positions[, "chr"])

  # write final summary of cross after all filtering
  summary.outfile.final <- paste0(dir.qc, "/", cross, "/summary_", cross, "_rils_after-filters.txt")
  capture.output(summary(rqtl), file = summary.outfile.final, type = "output")

  # write table with both physical and genetic distances
  outfile.name <- paste0(dir.qc, "/", cross, "/recomb-freq_", cross, "_rils.txt")
  write.table(df.gen.phys.positions, outfile.name, sep = "\t", quote = FALSE, row.names = FALSE)


  # plot recombination frequency by physical position
  if (plot == TRUE) {
    plot.list <- list()
    for (i in 1:10) {
      plot.chrm <- df.gen.phys.positions[which(df.gen.phys.positions[,"chr"] == i),]
      plot.chrm[,"pos"] <- plot.chrm[,"pos"]/1000000  # adjust scale to Mb
      plot.chrm[,"rf"] <- plot.chrm[,"rf"]*100        # adjust scale to cM
      plot.list[[i]] <- ggplot(plot.chrm, aes(x = pos, y = rf)) +
        geom_point(size = 0.5) +
        # ggtitle(paste("chrm", i)) +
        labs(x = "Physical position (Mb)",
             y = "Genetic position (cM)",
             title = paste("chrm", i),
             subtitle = paste0("(", NROW(plot.chrm), " markers)")) +
        scale_x_continuous(labels = scales::comma)
    }

    fig.name <- paste0(dir.qc, "/", cross, "/plot_genetic_vs_physical_dist_", cross, ".pdf")
    ggsave(filename = fig.name,
           plot = grid.arrange(grobs = plot.list,
                               bottom = textGrob(paste("Population size:", NROW(rqtl$pheno), "RILs"),
                                                 gp = gpar(fontsize = 20, col = "#737373"),
                                                 hjust = 0.1, vjust = -4.5)),
           device = "pdf")

  }
}



#### filtering hmp according to biparental crosses ----

# read genotypic data for parents and RILs
geno.data.parents <- fread(infile.parents, header = TRUE, data.table = FALSE)
geno.data.rils <- fread(infile.rils, header = TRUE, data.table = FALSE)

# read table with biparental crosses
rils.per.cross <- fread(cross.info, header = TRUE, data.table = FALSE)
# change "*" to "x" in cross column (use "fixed = TRUE" to avoid special character)
rils.per.cross[, "cross"] <- gsub("*", "x", fixed = TRUE, rils.per.cross[,"cross"])

# get the names of crosses with genotypic data
crosses.with.geno.data <- which(rils.per.cross[, "cross"] %in% list.dirs(qc.folder,
                                                                         full.names = FALSE,
                                                                         recursive = FALSE))
crosses.with.geno.data <- rils.per.cross[crosses.with.geno.data, "cross"]

# create output folder
if (!dir.exists(out.dir)) dir.create(out.dir)

# separate hmp by cross
for (cross in crosses.with.geno.data) {
  # parents
  parent <- strsplit(cross, split = "x")[[1]]
  geno.data.parents.cross <- cbind(geno.data.parents[, c(1:11)],
                                   geno.data.parents[, c(parent[1], parent[2])])

  parents.filename <- paste0(out.dir, "/usda_22kSNPs_", cross, "_parents.hmp.txt")
  write.table(geno.data.parents.cross, parents.filename, sep = "\t", quote = FALSE,
              row.names = FALSE)

  # rils
  rils <- rils.per.cross[which(rils.per.cross[,"cross"] == cross),"RILs"]
  rils <- strsplit(rils, split = ",")[[1]]
  geno.data.rils.cross <- cbind(geno.data.rils[, c(1:11)],
                                geno.data.rils[, which(colnames(geno.data.rils) %in% rils)])

  rils.filename <- paste0(out.dir, "/usda_22kSNPs_", cross, "_rils.hmp.txt")
  write.table(geno.data.rils.cross, rils.filename, sep = "\t", quote = FALSE,
              row.names = FALSE)
}



#### prepare input file for rqtl ----

for (cross in crosses.with.geno.data) {

  cat(cross, "\n")

  # load data first
  geno.cross.parents <- fread(paste0(out.dir, "/usda_22kSNPs_", cross, "_parents.hmp.txt"),
                              header = TRUE, data.table = FALSE)
  geno.cross.rils <- fread(paste0(out.dir, "/usda_22kSNPs_", cross, "_rils.hmp.txt"),
                           header = TRUE, data.table = FALSE)

  # write csv files for input in rqtl
  Hmp2Rqtl(parents = geno.cross.parents, rils = geno.cross.rils, out.dir = out.dir)

}



#### run rqtl to get estimated recombination frequency ----

# write log of this program
sink(paste0(qc.folder, "/rqtl_log.txt"))

for (cross in crosses.with.geno.data) {
  EstimateRecombinationFreq(cross = cross, dir.csv = out.dir, dir.qc = qc.folder, plot = TRUE)
}

# close sink connection
sink()
