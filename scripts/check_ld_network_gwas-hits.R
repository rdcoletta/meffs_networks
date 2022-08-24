library(data.table)

usage <- function() {
  cat("
description: check which network markers are in LD with gwas hits or top non-significant hits.

usage: Rscript check_ld_network_gwas-hits.R [modules_filename] [ld_filename] [gwas_filename] [output_filename] [...]

positional arguments:
  modules_filename          file with markers assigned to modules and with networks stats (e.g. kDiff)
  ld_filename               file containing plink LD results (R2 > 0.9) between gwas hits or top 
                            non-significant hits and network markers.
  gwas_filename             file with gwas hits and top non-significant hits for each environment
  output_filename           name of file to save results

optional argument:
  --help                    show this helpful message


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
modules_filename <- args[1]
ld_filename <- args[2]
gwas_filename <- args[3]
output_filename <- args[4]

# assert to have the correct optional arguments
pos_args <- 4
if (length(args) != pos_args) stop(usage(), "missing positional argument(s)")



#### find markers in network in ld with gwas hits ----

# load network file
network <- fread(modules_filename, header = TRUE, data.table = FALSE)

# load gwas info
gwas <- fread(gwas_filename, header = TRUE, data.table = FALSE)
# get gwas hits
gwas_hits <- unique(gwas[gwas$signif == TRUE, "marker"])
# get top non significant gwas hits
gwas_top_ns <- unique(gwas[gwas$signif == FALSE, "marker"])
# remove any non-significant hit that may have been significant in a different env
gwas_top_ns <- gwas_top_ns[!gwas_top_ns %in% gwas_hits]

# load ld file
ld <- fread(ld_filename, header = TRUE, data.table = FALSE)
# compute distance between markers
ld$bp_dist <- apply(ld, MARGIN = 1, function(marker) abs(as.numeric(marker["BP_B"]) - as.numeric(marker["BP_A"])))

# get ld info between gwas and network markers
ld_gwas_net <- data.frame(stringsAsFactors = FALSE)
for (gwas_marker in unique(gwas[order(gwas$chr, gwas$pos), "marker"])) {
  
  # get info about ld between network and gwas markers
  gwas_marker_ld <- ld[ld$SNP_A == gwas_marker | ld$SNP_B == gwas_marker, ]
  
  # proceed only if there's LD info about this marker
  if (nrow(gwas_marker_ld) > 0) {
    
    # for each pair of markers
    for (row in 1:nrow(gwas_marker_ld)) {
      # identify which marker belongs to the network
      net_marker <- ifelse(gwas_marker_ld[row, "SNP_A"] == gwas_marker,
                                    yes = gwas_marker_ld[row, "SNP_B"],
                                    no = gwas_marker_ld[row, "SNP_A"])
      # if marker is indeed in the network and is in ld to gwas marker
      if (net_marker %in% network$marker) {
        # append info to df
        ld_gwas_net <- rbind(ld_gwas_net,
                             data.frame(gwas_marker = gwas_marker,
                                        network_marker = net_marker,
                                        R2 = gwas_marker_ld[row, "R2"],
                                        bp_dist = gwas_marker_ld[row, "bp_dist"],
                                        stringsAsFactors = FALSE))
      }
    }
  }
}
# remove duplicates -- including whem columns 1 and 2 are swapped
ld_gwas_net <- ld_gwas_net[!duplicated(data.frame(t(apply(ld_gwas_net[1:2], MARGIN = 1, sort)),
                                                  ld_gwas_net[, c("R2", "bp_dist")])), ]

# check which markers are either a gwas hit or linked to one
ld_gwas_net$linked_to_gwas_hit <- ld_gwas_net$network_marker %in% gwas_hits | ld_gwas_net$gwas_marker %in% gwas_hits
# # check which markers are either a top non-sig gwas hit or linked to one
ld_gwas_net$linked_to_gwas_top_ns <- ld_gwas_net$network_marker %in% gwas_top_ns | ld_gwas_net$gwas_marker %in% gwas_top_ns




# annotate whether a marker in network is a gwas hit, top non sig hit, or not linked
network$gwas_status <- NA
for (row in 1:nrow(network)) {
  
  # get marker name
  net_marker <- network[row, "marker"]
  
  if (any(ld_gwas_net[ld_gwas_net$network_marker == net_marker, "linked_to_gwas_hit"])) {
    # if network marker is linked to at least one gwas hit, annotate as a hit
    network[row, "gwas_status"] <- "gwas_hit"
  } else if (any(ld_gwas_net[ld_gwas_net$network_marker == net_marker, "linked_to_gwas_top_ns"])) {
    # if that's not the case, check if marker is linked to at least one top non significant hit
    network[row, "gwas_status"] <- "gwas_top_non_sig"
  } else {
    # and annotate as not a gwas hit if anything else
    network[row, "gwas_status"] <- "not_gwas_hit"
  }
  
}

# how many markers in network are gwas hit or linked to one
round(sum(network$gwas_status == "gwas_hit")/ nrow(network), 3)
# how many markers in network are top non significant gwas hits or linked to one
round(sum(network$gwas_status == "gwas_top_non_sig")/ nrow(network), 3)
# how many markers in network are neither a gwas hit (or linked to one) or a top non significant gwas hits (or linked to one)
round(sum(network$gwas_status == "not_gwas_hit")/ nrow(network), 3)

# write updated network file
fwrite(network, file = output_filename, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)




#### debug ----

# modules_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/kDiff_per_module.txt"
# ld_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/ld_markers_gwas_network.ld.gz"
# gwas_filename <- "../genomic_prediction/hybrids/analysis/gwas/gwas_top-peaks_YLD-per-env.txt"
# output_filename <- "analysis/networks/YLD/meff_rrblup/norm_zscore/min_mod_size_50/pamStage_off/kDiff_per_module.gwas-status.txt"