library(data.table)
library(gtools)

usage <- function() {
  cat("
description: creates hybrid genotypes based on RIL genotypes used for the cross.

usage: Rscript create_hybrid_genotypes.R [pheno_file] [hmp_rils_filename] [outfile]

positional arguments:
  pheno_file                file with phenotypic data
  hmp_rils_filename
  outfile

optional argument:
  --help                    show this helpful message
  --checks=VALUE            comma-separated list of checks to remove from data
  --het-parents=VALUE       determine what to do when parental marker is heterozygous.
                            If 'random' (default), one of the two parental alleles from
                            the het marker will be randomly sampled when creating a hybrid
                            marker. If 'missing', then the hybrid marker will be set
                            to missing if one of the parents have a het marker.

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

# set default
checks <- NULL
het_parents <- "random"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {

  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--checks", "--het-parents")
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
if (!is.null(checks)) checks <- unlist(strsplit(checks, split = ","))
if (!het_parents %in% c("random", "missing")) {
  stop("Optional argument '--het-parents' should be either 'random' or 'missing'")
}

# get positional arguments
pheno_file <- args[1]
hmp_rils_filename <- args[2]
outfile <- args[3]



#### get hybrid information ----

# load pheno data
hybrid_info <- fread(pheno_file, header = TRUE, data.table = FALSE)

# keep only pedigree information
hybrid_info <- hybrid_info[, c("Hybrid", "ParentA", "ParentB")]
hybrid_info <- hybrid_info[!duplicated(hybrid_info$Hybrid), ]
# remove checks
hybrid_info <- hybrid_info[!hybrid_info$Hybrid %in% checks, ]
# adjust format
hybrid_info <- as.data.frame(hybrid_info[mixedorder(hybrid_info$Hybrid), ])



#### create genotypes ----

# load ril data
hmp_rils <- fread(hmp_rils_filename, header = TRUE, data.table = FALSE)

# create empty hmp to store hybrid genotypes
hmp_hybrids <- hmp_rils[, 1:11]

cat("Creating hybrid genotypes ('--het-parents=", het_parents, "'):\n", sep = "")

for (hybrid in hybrid_info$Hybrid) {

  cat("  hybrid ", hybrid, "\n")

  # get RIL names used in for single cross name
  p1 <- hybrid_info[which(hybrid_info$Hybrid == hybrid), "ParentA"]
  p2 <- hybrid_info[which(hybrid_info$Hybrid == hybrid), "ParentB"]

  # make sure there is genotypic data for both parents
  if (p1 %in% colnames(hmp_rils) & p2 %in% colnames(hmp_rils)) {

    cat("  ", p1, " x ", p2, "\n", sep = "")

    # check parental marker type
    marker_type <- apply(X = hmp_rils[, c(p1, p2)], MARGIN = 1, FUN = function(marker) {

      # make sure alleles are in the same order (to avoid cases like "AT" and "TA")
      genotypes <- sapply(marker, function(geno) {
        geno <- sort(unlist(strsplit(geno, split = "")))
        geno <- paste0(geno, collapse = "")
        return(geno)
      })
      # get unique genotypes between parents
      genotypes <- unique(as.character(genotypes))

      if (any(grepl("NN", genotypes))) {

        # if NN, marker is missing
        return("missing")

      } else if (length(genotypes) == 1) {

        # if there is one genotype, it's monomorphic
        # but distinguish if marker is het
        alleles <- unlist(strsplit(genotypes, split = ""))
        if (alleles[1] == alleles[2]) {
          return("mono")
        } else {
          return("het")
        }

      } else {

        # if there are two genotypes, it's polymorphic
        # but distiguish if one of the genotypes is het
        p1_alleles <- sort(unlist(strsplit(genotypes[1], split = "")))
        p2_alleles <- sort(unlist(strsplit(genotypes[2], split = "")))
        if (p1_alleles[1] == p1_alleles[2] & p2_alleles[1] == p2_alleles[2]) {
          return("poly")
        } else if (p1_alleles[1] != p1_alleles[2] & p2_alleles[1] == p2_alleles[2]) {
          return("het_p1")
        } else if (p1_alleles[1] == p1_alleles[2] & p2_alleles[1] != p2_alleles[2]) {
          return("het_p2")
        }

      }
    })

    # add marker type to single cross info
    geno_parents <- cbind(hmp_rils[, c(p1, p2)], marker_type, stringsAsFactors = FALSE)

    # if sampling alleles from het parents...
    if (het_parents == "random") {

      # ...create hybrid genotype
      set.seed(67251)
      geno_hybrid <- apply(geno_parents, MARGIN = 1, function(marker) {

        if (marker[3] == "mono" | marker[3] == "poly") {
          # for markers that are poly or monomorphic, get one allele from each parent
          p1_allele <- unlist(strsplit(marker[1], split = ""))[1]
          p2_allele <- unlist(strsplit(marker[2], split = ""))[1]
          geno_hybrid <- paste0(p1_allele, p2_allele)
        } else if (marker[3] == "het") {
          # for markers that are het for both parents, randomly sample an allele from each
          p1_allele <- sample(unlist(strsplit(marker[1], split = "")), size = 1)
          p2_allele <- sample(unlist(strsplit(marker[2], split = "")), size = 1)
          geno_hybrid <- paste0(p1_allele, p2_allele)
        } else if (marker[3] == "het_p1") {
          # for markers that are het for p1, randomly sample an allele that parent
          p1_allele <- sample(unlist(strsplit(marker[1], split = "")), size = 1)
          p2_allele <- unlist(strsplit(marker[2], split = ""))[1]
          geno_hybrid <- paste0(p1_allele, p2_allele)
        } else if (marker[3] == "het_p2") {
          # for markers that are het for p2, randomly sample an allele that parent
          p1_allele <- unlist(strsplit(marker[1], split = ""))[1]
          p2_allele <- sample(unlist(strsplit(marker[2], split = "")), size = 1)
          geno_hybrid <- paste0(p1_allele, p2_allele)
        } else if (marker[3] == "missing") {
          # for markers that are missing in at least one parent, set marker to NN
          geno_hybrid <- "NN"
        }

        return(geno_hybrid)

      })
    }

    # if turning het parental markers into missing hybrid marker...
    if (het_parents == "missing") {

      # ...create hybrid genotype
      geno_hybrid <- apply(geno_parents, MARGIN = 1, function(marker) {

        if (marker[3] == "mono" | marker[3] == "poly") {
          # for markers that are poly or monomorphic, get one allele from each parent
          p1_allele <- unlist(strsplit(marker[1], split = ""))[1]
          p2_allele <- unlist(strsplit(marker[2], split = ""))[1]
          geno_hybrid <- paste0(p1_allele, p2_allele)
        } else {
          # for markers that are missing or het in at least one parent, set marker to NN
          geno_hybrid <- "NN"
        }

        return(geno_hybrid)

      })
    }

    # View(cbind(geno_parents, geno_hybrid))

    # append hybrid genotype to final hmp
    hmp_hybrids <- cbind(hmp_hybrids, geno_hybrid)
    # correct column name for that hybrid
    colnames(hmp_hybrids)[NCOL(hmp_hybrids)] <- hybrid

  }
  else {
    # debug
    cat("  missing genotypic data in at least one of the parents\n")
  }

}

# keep consistent allele order in het markers (i.e. major allele first, then minor)
corrected_markers <- apply(X = hmp_hybrids[, 12:NCOL(hmp_hybrids)], MARGIN = 1, FUN = function(marker) {

  # get number of alleles
  alleles <- unlist(strsplit(paste0(marker[marker != "NN"], collapse = ""), split = ""))
  alleles <- sort(table(alleles), decreasing = TRUE)

  # define major and minor alleles
  major <- names(alleles)[1]
  minor <- names(alleles)[2]

  # define hets
  het_correct <- paste0(major, minor)
  het_wrong <- paste0(minor, major)

  # correct wrong hets
  marker[marker == het_wrong] <- het_correct

  return(marker)

})
corrected_markers <- data.frame(t(corrected_markers))
hmp_hybrids[, 12:NCOL(hmp_hybrids)] <- corrected_markers

# write file
fwrite(hmp_hybrids, outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

# fix allele columns using TASSEL
commands.hapdip <- paste0("/home/hirschc1/della028/software/tassel-5-standalone/run_pipeline.pl",
                          " -Xmx100g -importGuess ", outfile,
                          " -export ", outfile,
                          " -exportType HapmapDiploid")
system(commands.hapdip)




#### debug ----

# pheno_file <- "data/NIFA_CompleteDataset.csv"
# hmp_rils_filename <- "data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.sliding-window.hmp.txt"
# outfile <- "data/usda_22kmarkers_hybrids.hmp.txt"
# checks <- c("UIUC-CHECK1", "UIUC-CHECK2", "UIUC-CHECK3", "UIUC-CHECK4", "6049V2P")
# het_parents <- "random"
