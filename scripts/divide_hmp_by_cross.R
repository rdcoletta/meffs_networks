#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script creates hapmap files for each USDA family for
                   both parental and RIL data

      Usage: Rscript divide_hmp_by_cross.R [...]
"
      )
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 4) {
  stop("incorrect number of arguments provided.

       Usage: Rscript divide_hmp_by_cross.R [...]
       ")
}

# assign arguments to variables
file.hmp <- args[1]
cross.info <- args[2]
out.dir <- args[3]

if (args[4] == "--parents") {
  data.type <- "parents"
} else if (args[4] == "--rils") {
  data.type <- "rils"
} else {
  stop("Last argument should be either '--parents' or '--rils'")
}

# file.hmp <- "data/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.hmp.txt"
# file.hmp <- "data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.hmp.txt"
# cross.info <- "data/usda_biparental-crosses.txt"
# out.dir <- "data/hapmap_by_cross"
# data.type <- "parents"


#### libraries used ----

if(!require("data.table")) install.packages("data.table")



#### subset by cross ----

# load data
hmp <- fread(file.hmp, header = TRUE, data.table = FALSE)

# load table with cross information
df.crosses <- fread(cross.info, header = TRUE, data.table = FALSE)

# make dir to store results
if (!dir.exists(out.dir)) dir.create(out.dir)

# write merged RILs hapmap by cross
for (row in 1:NROW(df.crosses)) {

  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)

  if (data.type == "parents") individuals <- unlist(strsplit(cross, split = "x"))
  if (data.type == "rils") individuals <- unlist(strsplit(df.crosses[row, "RILs"], split = ","))

  # subset hapmap
  subset.hmp <- cbind(hmp[, 1:11], hmp[, which(colnames(hmp) %in% individuals)])
  # since i'm subsetting, the alleles column might be wrong...
  subset.hmp$alleles <- NA

  # write results only if individuals were actually in the file
  if (NCOL(subset.hmp) > 12) {

    outfile.cross <- rev(unlist(strsplit(file.hmp, split = "/")))[1]
    outfile.cross <- gsub("hmp.txt", paste0(cross, ".hmp.txt"), outfile.cross)
    outfile.cross <- paste0(out.dir, "/", outfile.cross)
    fwrite(subset.hmp, file = outfile.cross, sep = "\t", na = "NN", quote = FALSE)

  }
}
