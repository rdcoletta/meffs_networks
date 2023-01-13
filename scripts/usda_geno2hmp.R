#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script tranforms genotypic data of 22k SNPs from both
             parental lines and RILs from the USDA project

Usage: Rscript usda_geno2hmp.R parents_file rils_file1 rils_file2 output_folder")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 4) {
  stop("incorrect number of arguments provided.

Usage: Rscript usda_geno2hmp.R parents_file rils_file1 rils_file2 output_folder
       ")
}

# assign arguments to variables
parents.file <- args[1]
rils1.file <- args[2]
rils2.file <- args[3]
output.folder <- args[4]


# parents.file <- "data/SNP_chip/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv"
# rils1.file <- "data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318.csv"
# rils2.file <- "data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318.csv"
# output.folder <- "data"



#### libraries ----

if(!require("data.table")) install.packages("data.table")




#### functions ----

geno2hmp <- function(geno.file, output.folder) {

  # read in table
  geno.data <- fread(geno.file, header = FALSE, data.table = FALSE)

  # create hapmap structrue
  hmp.data <- data.frame("rs" = geno.data[-(1:2),1],
                         "alleles" = NA,
                         "chrom" = geno.data[-(1:2),2],
                         "pos" = geno.data[-(1:2),3],
                         "strand" = NA,
                         "assembly" = NA,
                         "center" = NA,
                         "protLSID" = NA,
                         "assayLSID" = NA,
                         "panel" = NA,
                         "QCcode" = NA)

  # add genotypes
  hmp.data <- cbind(hmp.data, geno.data[-(1:2),4:NCOL(geno.data)])

  # include genotype names as column names
  colnames(hmp.data)[12:NCOL(hmp.data)] <- as.character(geno.data[1,4:NCOL(geno.data)])

  # convert missing data "--" to "NN"
  hmp.data[, 12:NCOL(hmp.data)] <- apply(X = hmp.data[, 12:NCOL(hmp.data)], MARGIN = 2,
                                         FUN = function(x) gsub(pattern = "--", replacement = "NN", x))

  # create output name
  outname <- rev(unlist(strsplit(geno.file, split = "/")))[1]
  outname <- gsub(".csv", ".hmp.txt", outname, fixed = TRUE)
  outname <- gsub("-Part[0-9]", "", outname , perl = TRUE)
  # hmp file
  outfile.path <- paste0(output.folder, "/", outname)
  # id table
  id.table.outfile.path <- gsub(".hmp.txt", ".txt", outname, fixed = TRUE)
  id.table.outfile.path <- paste0(output.folder, "/id_table_", id.table.outfile.path)

  # print hapmap converted genotypic data
  fwrite(hmp.data, file = outfile.path, sep = "\t", na = NA, quote = FALSE)

  # also print a table relating parent name and its ID when genotyping
  id.table <- data.frame("genotype_name" = as.character(geno.data[1,4:NCOL(geno.data)]),
                         "source_id" = as.character(geno.data[2,4:NCOL(geno.data)]))
  fwrite(id.table, file = id.table.outfile.path, sep = "\t")
}



#### parental data ----

cat("Converting parental data\n")
geno2hmp(geno.file = parents.file, output.folder = output.folder)



#### RIL data ----

# this dataset is divided in two separate files. First, i will use function above to convert part 1
# to hapmap, and then will append only genotypic information of part 2 into part 1

cat("Converting ril data\n")
cat("  file 1\n")
geno2hmp(geno.file = rils1.file, output.folder = output.folder)

cat("  file 2\n")
# read ril data again
ril1.hmp <- rev(unlist(strsplit(rils1.file, split = "/")))[1]
ril1.hmp <- gsub(".csv", ".hmp.txt", ril1.hmp, fixed = TRUE)
ril1.hmp <- gsub("-Part[0-9]", "", ril1.hmp , perl = TRUE)
ril1.hmp <- paste0(output.folder, "/", ril1.hmp)

ril.data1 <- fread(ril1.hmp, header = TRUE, data.table = FALSE)
ril.data2 <- fread(rils2.file, header = FALSE, data.table = FALSE)

# filter part 2 so it has only the genotypic information
ril.data2.filter <- ril.data2[-(1:2), 4:NCOL(ril.data2)]
colnames(ril.data2.filter) <- ril.data2[1, 4:NCOL(ril.data2)]

# convert missing data "--" to "NN" from part 2
ril.data2.filter <- data.frame(apply(X = ril.data2.filter, MARGIN = 2,
                                     FUN = function(x) gsub(pattern = "--", replacement = "NN", x)),
                               check.names = FALSE)

# append genotypic info of part 2 to part 1
ril.data.combined <- cbind(ril.data1, ril.data2.filter)
# View(ril.data.combined)

# overwrite RIL data so it has both parts now
fwrite(ril.data.combined, file = ril1.hmp, sep = "\t", na = NA, quote = FALSE)


# append the rest of ids of part 2 on table that already has the part 1

# id table
id.table.file <- gsub(".hmp.txt", ".txt", ril1.hmp, fixed = TRUE)
id.table.file <- gsub("/", "/id_table_", id.table.file)


id.table.rils <- fread(id.table.file, header = TRUE, data.table = F)

id.table.ril.part2 <- data.frame("genotype_name" = as.character(ril.data2[1,4:NCOL(ril.data2)]),
                                 "source_id" = as.character(ril.data2[2,4:NCOL(ril.data2)]))

id.table.rils.combined <- rbind(id.table.rils, id.table.ril.part2)

# overwrite previous table file
fwrite(id.table.rils.combined, file = id.table.file, sep = "\t")
cat("Done!\n")
