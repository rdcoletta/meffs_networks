library(data.table)
library(gtools)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)

usage <- function() {
  cat("
description: run ANOVA and calculate percent variance explained on a phenotypic trait.

usage: Rscript anova_sim_traits_hybrids.R [pheno_file] [envs] [checks] [trait] [outfolder]

positional arguments:
  pheno_file            file with phenotypes
  envs                  list of comma-separated names of environments to keep
  checks                list of comma-separated names of checks to remove
  trait                 trait name
  outfolder             name of folder to save results

optional argument:
  --help                  show this helpful message

"
  )
}




#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 5) stop(usage(), "missing positional argument(s)")

# get positional arguments
pheno_file <- args[1]
envs <- unlist(strsplit(args[2], split = ","))
checks <- unlist(strsplit(args[3], split = ","))
trait <- args[4]
outfolder <- args[5]
if (!dir.exists(outfolder)) dir.create(outfolder, recursive = TRUE)



#### QC traits ----

# read phenotypic file
data <- fread(pheno_file, header = TRUE, data.table = FALSE)

# format columns e remove unuseful info
data$Experiment <- gsub("NIFA-20", "", data$Experiment)
data$environment <- paste0(data$Location, data$Experiment)
data <- data[ c("Hybrid", "environment", "Replication", trait)]
# colnames(data)[c(1, 3)] <- c("genotype", "rep")
colnames(data)[c(1, 3, 4)] <- c("genotype", "rep", "trait_value")

# remove checks, if requested
if (!is.null(checks)) data <- data[!data$genotype %in% checks, ]

# make sure all genotypes have a rep number
data <- data %>%
  filter(environment %in% envs) %>% 
  group_by(genotype, environment) %>%
  mutate(rep = 1:n())

# # plot summary
# plot_summary <- ggplot(data, aes(x = environment, y = trait_value)) +
#   geom_violin(color = "gray70", fill = "gray70") +
#   stat_summary(fun.data = "mean_cl_normal", color = "black", show.legend = FALSE)
# 
# out_plot_summary <- paste0(outfolder, "/distribution-per-env.pdf")
# ggsave(filename = out_plot_summary, plot = plot_summary, device = "pdf")

# fit mixed linear model
# fixed effects -- genotypes
# random effects -- reps within envs, envs, GxE)
trait_lmer <- lmer(trait_value ~ genotype + (1 | environment/rep) + (1 | genotype:environment),
                   data = data, REML = TRUE)
trait_lmer_summary <- summary(trait_lmer)

# calculate percent variation explained
pve_fixed <- as.numeric(r.squaredGLMM(trait_lmer)[1, "R2m"])
pve_random <- as.numeric(r.squaredGLMM(trait_lmer)[1, "R2c"] - pve_fixed)
pve_residual <- 1 - pve_fixed - pve_random

# extract variation from random effects
var_random <- data.frame(trait_lmer_summary$varcor) %>%
  filter(grp != "Residual") %>%
  select(source = grp, variance = vcov) %>%
  mutate(pve = (variance / sum(variance)) * pve_random) %>%
  select(source, pve)

# run anova on fixed and random effects separately
trait_lmer_anova_fixed <- anova(trait_lmer)
trait_lmer_anova_random <- ranova(trait_lmer)

# format anova results and add percent variance explained
pval_fixed <- cbind(source = rownames(trait_lmer_anova_fixed),
                    data.frame(trait_lmer_anova_fixed, row.names = NULL, check.names = FALSE)) %>%
  select(source, pval = `Pr(>F)`) %>%
  mutate(pve = pve_fixed)

pval_random <- cbind(source = rownames(trait_lmer_anova_random),
                     data.frame(trait_lmer_anova_random, row.names = NULL, check.names = FALSE)) %>%
  select(source, pval = `Pr(>Chisq)`) %>%
  mutate(source = gsub("(1 | ", "", source, fixed = TRUE),
         source = gsub(")", "", source, fixed = TRUE)) %>%
  filter(source != "<none>") %>%
  full_join(var_random, by = "source")

mm_pve_pval <- rbind(pval_fixed, pval_random) %>%
  mutate(pve = round(pve, digits = 3),
         signif = case_when(pval > 0.1 ~ "NS",
                            pval > 0.05 & pval <= 0.1 ~ ".",
                            pval > 0.01 & pval <= 0.05 ~ "*",
                            pval <= 0.01 ~ "***")) %>%
  select(source, pve, pval, signif)

# add residual info into anova results
mm_pve_pval <- rbind(mm_pve_pval, data.frame(source = "residual",
                                             pve = round(pve_residual, digits = 3),
                                             pval = NA,
                                             signif = NA))

# write anova results
out_pve <- paste0(outfolder, "/", trait, "_pve.txt")
fwrite(x = mm_pve_pval, file = out_pve, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# also print full anova results
options(max.print=2000)
sink(paste0(outfolder, "/", trait, "_anova.txt"))
print(trait_lmer_summary)
cat("\n---- ANOVA fixed effects ----\n")
print(trait_lmer_anova_fixed)
cat("\n---- ANOVA random effects ----\n")
print(trait_lmer_anova_random)
sink()

# More on calculating PVE for mixed models:
#   https://stats.stackexchange.com/questions/7240/proportion-of-explained-variance-in-a-mixed-effects-model
#   https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/
#   https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
#   Marginal R2 is concerned with variance explained by fixed factors
#   Conditional R2 is concerned with variance explained by both fixed and random factors.


#### debug ----

# pheno_file <- "data/NIFA_CompleteDataset.csv"
# envs <- c("BEC-BL19", "BEC-BL20", "COR19", "COR20", "MIN19", "MIN20", "SYN19", "SYN20", "URB19")
# checks <- c("UIUC-CHECK1", "UIUC-CHECK2", "UIUC-CHECK3", "UIUC-CHECK4", "6049V2P")
# trait <- "YLD"
# outfolder <- "data"