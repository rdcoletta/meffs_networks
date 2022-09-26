library(data.table)
library(FW)
library(tidyr)
library(tibble)
library(gtools)

usage <- function() {
  cat("
description: get environmental effects from a a Finlay-Wilkson regression.

usage: Rscript get_FW_env-idx.R [trait_filename] [outfolder] [...]

positional arguments:
  trait_filename         file with phenotypic data (3 columns: genotype, env, trait_value)
  outfolder              name of folder to save results

optional argument:
  --help                 show this helpful message


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
trait_filename <- args[1]
outfolder <- args[2]
if (!dir.exists(outfolder)) dir.create(outfolder, recursive = TRUE)

# assert to have the correct optional arguments
if (length(args) != 2) stop(usage(), "missing positional argument(s)")



#### get environmental indices ----

# load data
trait_data <- fread(trait_filename, header = TRUE, data.table = FALSE)
# make sure column names are correct
colnames(trait_data) <- c("genotype", "env", "trait_value")

# remove missing data
trait_data <- pivot_wider(trait_data, names_from = "env", values_from = "trait_value")
trait_data <- na.omit(trait_data)
trait_data <- pivot_longer(trait_data, -genotype, names_to = "env", values_to = "trait_value")

# reorder factor levels for plotting
trait_data$genotype <- factor(trait_data$genotype, levels = mixedsort(unique(trait_data$genotype)))
trait_data$env <- factor(trait_data$env, levels = mixedsort(unique(trait_data$env)))

# run FW regression
fw <- FW(y = trait_data$trait_value, VAR = trait_data$genotype, ENV = trait_data$env, method = "OLS")

# plot results
pdf(file = paste0(outfolder, "/FW_plot.pdf"))
plot.FW(fw, cex = 0.2, lwd = 0.2)
dev.off()

# save env idx, but make sure it's in the format for downstream analysis
fw_env_idx <- rownames_to_column(data.frame(value = fw$h), var = "env")
fw_env_idx <- data.frame(env = fw_env_idx$env, covariable = "FW", value = fw_env_idx$value)
fwrite(fw_env_idx, file = paste0(outfolder, "/FW_env_idx.txt"), quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# trait_filename <- "data/1stStage_BLUEs.YLD-per-env.txt"
# outfolder <- "analysis/FW/YLD"

