#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script reconstruct the PHG35 haplotype based on SNP information of
                   RILs that have this inbred line as parent.

      Usage: Rscript reconstruct_PHG35_from_RIL_data.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

       Usage: Rscript reconstruct_PHG35_from_RIL_data.R [...]
       ")
}

# assign arguments to variables
crosses.file <- args[1]
parents.filename <- args[2]
rils.filename <- args[3]

# crosses.file <- "data/usda_biparental-crosses.txt"
# parents.filename <- "data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt"
# rils.filename <- "data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### reconstruction ----

cat("Loading data...\n")

# load hapmap files
parents.hmp <- fread(parents.filename, header = TRUE, data.table = FALSE)
rils.hmp <- fread(rils.filename, header = TRUE, data.table = FALSE)

# open file with crosses information and select only crosses with PHG35
crosses <- fread(crosses.file, header = TRUE, data.table = FALSE)
crosses <- crosses[grep("PHG35", crosses[, "cross"]), ]
crosses[, "cross"] <- gsub("*", "x", crosses[, "cross"], fixed = TRUE)

# create empty list to store PHG35 reconstruction for each cross
reconstructed.PHG35.all.crosses <- list()


for (cross in crosses[, "cross"]) {

  # select the RILs that were actually genotyped for that cross
  cross.rils <- unlist(strsplit(crosses[which(crosses[, "cross"] == cross), "RILs"], split = ","))
  cross.rils <- cross.rils[cross.rils %in% colnames(rils.hmp)]

  # only parse PHG35 crosses that have genotyped rils
  if (length(cross.rils) > 0) {

    cat("  analyzing cross", cross, "\n")

    # filter parental data to have only parents of the cross
    parent2 <- unlist(strsplit(cross, split = "x"))
    parent2 <- parent2[grep("PHG35", parent2, invert = TRUE)]
    parents.hmp.cross <- cbind(parents.hmp[, 1:11], parents.hmp[, c("PHG35", parent2)])

    # filter ril data to have only genotyped rils for that cross
    rils.hmp.cross <- cbind(rils.hmp[, 1:11], rils.hmp[, cross.rils])

    # create empty vector to store reconstructed PHG35 alleles
    reconstructed.PHG35 <- rep("NN", times = NROW(parents.hmp.cross))

    for (row in 1:NROW(parents.hmp.cross)) {

      # get allele on resequencing data from non-PHG35 parent
      alleles.parent2 <- unlist(strsplit(parents.hmp.cross[row, parent2], split = ""))

      # only proceed if allele in non-PHG35 parent is not missing and is homozygous
      if (alleles.parent2[1] == alleles.parent2[2] & all(unique(alleles.parent2) != "N")) {

        # get alleles on ril data
        alleles.rils <- as.character(rils.hmp.cross[row, 12:NCOL(rils.hmp.cross)])
        alleles.rils <- unlist(strsplit(paste0(alleles.rils, collapse = ""), split = ""))
        # count allele frequency
        df.allele.count <- data.frame(table(alleles.rils), stringsAsFactors = FALSE)
        df.allele.count <- cbind(df.allele.count, percent = NA)
        # but first exclude missing allele ("N") if there are other two alleles
        if (NROW(df.allele.count) > 2) {
          df.allele.count <- df.allele.count[which(df.allele.count[, "alleles.rils"] != "N"), ]
        }
        # get percent allele frequency
        for (allele in 1:NROW(df.allele.count)) {
          allele.freq <- df.allele.count[allele, "Freq"] / sum(df.allele.count[, "Freq"])
          df.allele.count[allele, "percent"] <- round(allele.freq, digits = 2)
        }

        # keep only alleles with allele frequency > 0.05 (since it can be an artifact from SNP calling)
        alleles.rils <- df.allele.count[which(df.allele.count[, "percent"] > 0.05), "alleles.rils"]
        alleles.rils <- as.character(alleles.rils)

        # exclude allele from non-PHG35 parent to get PHG35 allele
        allele.PHG35 <- alleles.rils[!alleles.rils %in% alleles.parent2]

        # if there's only one allele in "alleles.rils" it means it is monomorphic between the two
        # parents, thus "allele.PHG35" should have the same allele as the other parent
        if (length(allele.PHG35) == 0) {
          allele.PHG35 <- unique(alleles.parent2)
        }
        # if there are two alleles remaining in "allele.PHG35" it means that there is either a
        # third allele for this locus or there is missing data on rils and I can't confirm whether
        # this missing allele is from PHG35 or parent2
        if (length(allele.PHG35) > 1) {
          allele.PHG35 <- "N"
        }

        # debug
        if (length(reconstructed.PHG35[row]) != length(paste0(allele.PHG35, allele.PHG35))) {
          print(row)
        }

        reconstructed.PHG35[row] <- paste0(allele.PHG35, allele.PHG35)
        # note: there won't be any heterozygote for reconstructed PHG35

      }
    }

    # append reconstructed PHG35 for this cross into list
    reconstructed.PHG35.all.crosses[[cross]] <- reconstructed.PHG35

  }
}

cat("Merging results\n")

# get final reconstructed PHG35 by merging results from all crosses and transforming disagreements
# in missing data
reconstructed.PHG35.df <- data.frame(do.call(cbind, reconstructed.PHG35.all.crosses), stringsAsFactors = FALSE)
reconstructed.PHG35.final <- apply(reconstructed.PHG35.df[, ], MARGIN = 1, function(genotypes) {
  # get unique genotypes among all crosses for that SNP
  genotypes <- unique(genotypes)
  # if there are more than one genotype
  if (length(genotypes) > 1) {
    # exclude missing data
    genotypes <- genotypes[genotypes != "NN"]
    # even if after excluding missing data, there is more than one unique genotype, make SNP as NN
    if (length(genotypes) > 1) {
      genotypes <- "NN"
    }
  }
  return(genotypes)
})


# sum(parents.hmp[, "PHG35"] != reconstructed.PHG35.final)


# fix parental data with the new reconstructed PHG35 data
parents.hmp[, "PHG35"] <- reconstructed.PHG35.final

# write fixed hmp file
fixed.outfile <- gsub(pattern = "hmp.txt",
                      replacement = "PHG35-reconstructed.hmp.txt",
                      parents.filename, fixed = TRUE)
fwrite(parents.hmp, file = fixed.outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)
