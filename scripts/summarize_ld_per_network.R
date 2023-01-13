library(data.table)


usage <- function() {
  cat("
description: summarize LD results for each module of each network.

usage: Rscript summarize_ld_per_network.R [net_ld_folder] [...]

positional arguments:
  net_ld_folder                path to folder with LD results for network

optional argument:
  --help                       show this helpful message
  --meff-model=[VALUE]         name of marker effect models
  --norm-method=[VALUE]        name of normalization methods
  --minsize=[VALUE]            name of minimum module size
  --pamStage=[VALUE]           name of pamStage values

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
net_ld_folder <- args[1]

# set default
meff_model <- "rrblup,gwas"
norm_method <- "minmax,zscore"
minsize <- "25,50,100"
pamStage <- "on,off"

# assert to have the correct optional arguments
if (length(args) < 1) stop(usage(), "missing positional argument(s)")

if (length(args) > 1) {
  
  opt_args <- args[-1]
  opt_args_allowed <- c("--meff-model", "--norm-method", "--minsize", "--pamStage")
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



#### summarize ld ----

# create empty dataframe to store results
net_ld_results <- data.frame(stringsAsFactors = FALSE)

# get list of LD files
net_modules <- list.files(net_ld_folder, pattern = ".ld.gz")

for (mod in net_modules) {
  
  # load ld file
  ld_module <- fread(paste0(net_ld_folder, "/", mod), header = TRUE, data.table = FALSE)
  # select columns of interest
  ld_module <- ld_module[, c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2")]
  # get module name
  mod_name <- gsub("ld_markers_", "", gsub(".ld.gz", "", mod, fixed = TRUE))
  # add network settings
  net_ld_results_mod <- data.frame(meff_model = meff_model, norm_method = norm_method,
                                   minsize = minsize, pamStage = pamStage, module = mod_name,
                                   ld_module, stringsAsFactors = FALSE)
  # get accuracy results
  net_ld_results <- rbind(net_ld_results, net_ld_results_mod)
  
}

# write summary
fwrite(net_ld_results, file = paste0(net_ld_folder, "/summary_ld.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


#### debug ----

# net_ld_folder <- "analysis/networks/YLD"
# meff_model <- c("rrblup", "gwas")
# norm_method <- c("minmax", "zscore")
# minsize <- c(25, 50, 100)
# pamStage <- c("on", "off")
