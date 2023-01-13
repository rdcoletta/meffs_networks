#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script correct marker alleles in the wrong strand from parental and RIL
                   datasets from the USDA project.

      Usage: Rscript correct_SNP_strands.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

       Usage: Rscript correct_SNP_strands.R [...]
       ")
}

# assign arguments to variables
v4.coord.file <- args[1]
hmp.parents.file <- args[2]
hmp.rils.file <- args[3]

# v4.coord.file <- "data/SNP_positions_v2-to-v4_probes-100bp.txt"
# hmp.parents.file <- "data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt"
# hmp.rils.file <- "data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt"



#### libraries used ----

if(!require("data.table")) install.packages("data.table")



#### functions ----

GetComp <- function(allele) {
  if (allele == "A") {
    complement <- "T"
  } else if (allele == "T") {
    complement <- "A"
  } else if (allele == "C") {
    complement <- "G"
  } else if (allele == "G") {
    complement <- "C"
  } else if (allele == "N") {
    complement <- "N"
  } else {
    stop("Not a nucleotide")
  }
  return(complement)
}


# load v4 coordinates
v4.coord <- fread(v4.coord.file, header = TRUE, stringsAsFactors = FALSE)
# make SNP in v4.coord diploid
v4.coord[, "SNP"] <- sapply(v4.coord[, "SNP"], function(snp) {
  snp.diploid <- paste0(snp, snp)
})

# load hapmap data
hmp.parents <- fread(hmp.parents.file, header = TRUE, data.table = FALSE)
hmp.rils <- fread(hmp.rils.file, header = TRUE, data.table = FALSE)

# sort all data data
v4.coord <- v4.coord[order(v4.coord$chr_v4, v4.coord$pos_v4), ]
hmp.parents <- hmp.parents[order(hmp.parents$chrom, hmp.parents$pos), ]
hmp.rils <- hmp.rils[order(hmp.rils$chrom, hmp.rils$pos), ]

# create empty dfs to store results
hmp.parents.corrected <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.parents)))
colnames(hmp.parents.corrected) <- colnames(hmp.parents)

hmp.rils.corrected <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.rils)))
colnames(hmp.rils.corrected) <- colnames(hmp.rils)

for (chr in unique(hmp.parents[, "chrom"])) {

  cat("Analyzing chr", chr, "\n")

  # filter v4 data to have only markers for that chromosome
  v4.coord.chr <- v4.coord[which(v4.coord[, "chr_v4"] == chr), ]
  hmp.parents.chr <- hmp.parents[which(hmp.parents[, "chrom"] == chr), ]
  hmp.rils.chr <- hmp.rils[which(hmp.rils[, "chrom"] == chr), ]

  # get mismatches
  hmp.parents.mismatches <- hmp.parents.chr[which(hmp.parents.chr[, "B73"] != v4.coord.chr[, "SNP"]), ]
  v4.coord.mismatches <- v4.coord.chr[which(v4.coord.chr[, "SNP"] != hmp.parents.chr[, "B73"]), ]

  count.missing <- 0
  count.hets <- 0
  count.mismatches <- 0
  count.complement <- 0

  SNPs.opposite.strand <- c()
  SNPs.disagree <- c()

  for (row in 1:NROW(hmp.parents.mismatches)) {

    b73.v2.alleles <- unlist(strsplit(as.character(hmp.parents.mismatches[row, "B73"]), split = ""))
    b73.v4.alleles <- unlist(strsplit(as.character(v4.coord.mismatches[row, "SNP"]), split = ""))

    if (paste0(b73.v2.alleles, collapse = "") == "NN") {

      # count missing
      count.missing <- count.missing + 1

    } else if (b73.v2.alleles[1] != b73.v2.alleles[2]) {

      # count hets
      count.hets <- count.hets + 1

    } else if (b73.v2.alleles[1] == b73.v2.alleles[2]) {

      # count mismatches between v2 and v4
      count.mismatches <- count.mismatches + 1

      # since the genotype is homozygous, check if they are all from opposite strands
      complement <- GetComp(b73.v2.alleles[1])
      if (complement == b73.v4.alleles[1]) {
        # count opposite strand
        count.complement <- count.complement + 1
        # get SNPs at opposite strand
        SNPs.opposite.strand <- append(SNPs.opposite.strand, hmp.parents.mismatches[row, "pos"])
      } else {
        cat("SNP at position", hmp.parents.mismatches[row, "pos"], "from hapmap doesn't match SNP",
            "at position", as.character(v4.coord.mismatches[row, "pos_v4"]), "from v4 coordinates",
            "on chr", chr, "\n")
        # get SNPs that disagree and are not from oppsotie strand
        SNPs.disagree <- append(SNPs.disagree, hmp.parents.mismatches[row, "pos"])
      }
    }

  }

  cat("  SNPs in opposite strand:", count.complement, "\n")

  cat("  Reversing strand for parents...\n")

  for (parent in 12:NCOL(hmp.parents.chr)) {
    for (pos in SNPs.opposite.strand) {
      alleles <- hmp.parents.chr[which(hmp.parents.chr[, "pos"] == pos), parent]
      alleles <- unlist(strsplit(alleles, split = ""))
      complement <- sapply(alleles, GetComp)
      complement <- paste(complement, collapse = "")
      hmp.parents.chr[which(hmp.parents.chr[, "pos"] == pos), parent] <- complement
    }
  }

  cat("  Done!\n  Reversing strands for RILs...\n")

  for (parent in 12:NCOL(hmp.rils.chr)) {
    for (pos in SNPs.opposite.strand) {
      alleles <- hmp.rils.chr[which(hmp.rils.chr[, "pos"] == pos), parent]
      alleles <- unlist(strsplit(alleles, split = ""))
      complement <- sapply(alleles, GetComp)
      complement <- paste(complement, collapse = "")
      hmp.rils.chr[which(hmp.rils.chr[, "pos"] == pos), parent] <- complement
    }
  }

  cat("  Done!\n")

  # append to final table
  hmp.parents.corrected <- rbind(hmp.parents.corrected, hmp.parents.chr)
  hmp.rils.corrected <- rbind(hmp.rils.corrected, hmp.rils.chr)
}

# overwrite files
fwrite(x = hmp.parents.corrected, file = hmp.parents.file, sep = "\t", na = NA, quote = FALSE)
fwrite(x = hmp.rils.corrected, file = hmp.rils.file, sep = "\t", na = NA, quote = FALSE)
