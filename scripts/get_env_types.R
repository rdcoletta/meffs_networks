library(data.table)
library(EnvRtype)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)

usage <- function() {
  cat("
description: get environment similarity from a file with geographic coordinates.

usage: Rscript get_usda_env-types.R [coords_filename] [output_folder] [...]

positional arguments:
  coords_filename           file containing geographic coordinates (6 columns:
                            Location, Symbol, Latitude, Longitude, Start_date, End_date)
  output_folder             name of folder to save results

optional argument:
  --help                    show this helpful message
  --country=VALUE           3-letter ISO code for country (default: USA)
  --start=VALUE             planting date (default: 2020-04-01)
  --end=VALUE               harvesting date (default: 2020-10-31)
  --interval-window=VALUE   window size (in days) to capture environmental covariables


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
coords_filename <- args[1]
output_folder <- args[2]
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# set default
country <- "USA"
start <- "2020-04-01"
end <- "2020-10-31"
interval_window <- "30"

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {

  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--country", "--start", "--end", "--interval-window")
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

if (suppressWarnings(!is.na(as.numeric(interval_window)))) {
  interval_window <- as.numeric(interval_window)
} else {
  stop("Optional argument '--interval-window' should be a number")
}



#### get environmental data ----

# load file with location coordinates
sites_coord <- fread(coords_filename, header = TRUE, data.table = FALSE)

# get weather data
env_data <- get_weather(env.id = sites_coord$Symbol,
                        lat = sites_coord$Latitude,
                        lon = sites_coord$Longitude,
                        country = country,
                        start.day = sites_coord$Start_date,
                        end.day = sites_coord$End_date)
# remove ALT columns
env_data <- env_data[, colnames(env_data) != "ALT"]

# process data
env_data <- processWTH(env.data = env_data)

# T2M           Temperature at 2 Meters
# T2M_MAX       Maximum Temperature at 2 Meters
# T2M_MIN       Minimum Temperature at 2 Meters
# PRECTOT       Precipitation Corrected (mm/day)
# WS2M          Wind Speed at 2 Meters
# RH2M          Relative Humidity at 2 Meters
# T2MDEW        Dew/Frost Point at 2 Meters
# n             Actual duration of sunshine (hour)
# N             Daylight hours (hour)
# RTA           Extraterrestrial radiation (MJ/m^2/day)
# SRAD          Solar radiation (MJ/m^2/day)
# SPV           Slope of saturation vapour pressure curve (kPa.Celsius)
# VPD           Vapour pressure deficit (kPa)
# ETP           Potential Evapotranspiration (mm.day)
# PETP          Deficit by Precipitation (mm.day)
# GDD           Growing Degree Day (oC/day)
# FRUE          Effect of temperature on radiation use efficiency (from 0 to 1)
# T2M_RANGE     Daily Temperature Range (oC day)

# select environmental variables to use
env_vars <- c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOT", "WS2M", "RH2M", "T2MDEW",
              "n", "N", "RTA", "SPV", "VPD", "ETP", "PETP", "GDD","FRUE", "T2M_RANGE")

# summarize data
summary_envs <- summaryWTH(env.data = env_data, var.id = env_vars)
summary_envs <- pivot_longer(summary_envs, mean:prob_0.75, names_to = "stat", values_to = "value")
# plot summary
summary_envs_plot <- ggplot(subset(summary_envs, stat != "sum")) +
  facet_wrap(~ variable) +
  geom_line(aes(x = env, y = value, color = stat, group = stat)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(output_folder, "/env_weather_data_summary.pdf"),
       plot = summary_envs_plot, device = "pdf")

# # get environmental types using data-driven limits (instead of known cardinals)
# ET <- env_typing(env.data = env_data, env.id = 'env', var.id = env_vars, format = "wide")
# pdf(file = paste0(output_folder, "/envs_types.pdf"))
# plot_panel(ET, title = 'Panel of Environmental Types')
# dev.off()

# remove downloaded temp files
unlink("USA*_msk_alt.*")

# plot cumulative GDDs for each environment
gdd_plot <- env_data %>%
  select(env, daysFromStart, GDD) %>%
  group_by(env) %>%
  mutate(GDD = cumsum(GDD)) %>%
  ungroup() %>%
  pivot_longer(-c(env, daysFromStart), names_to = "idx", values_to = "values") %>%
  ggplot() +
  geom_line(aes(x = daysFromStart, y = values, color = env)) +
  labs(y = "cumulative_GDD")
ggsave(filename = paste0(output_folder, "/GDDs_per_env.pdf"),
       plot = gdd_plot, device = "pdf")



#### get environmental covariables means ----

# get EC
EC <- W_matrix(env.data = env_data, env.id = 'env', var.id = env_vars,
               statistic = "mean", by.interval = FALSE)
# pdf(file = paste0(output_folder, "/env_covariables_panel.pdf"))
# plot_panel(EC, title = 'Panel of Environmental Covariables')
# dev.off()

# format data
EC <- data.frame(EC, stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "env") %>%
  pivot_longer(contains("_mean"), names_to = "covariable", values_to = "value") %>%
  mutate(covariable = gsub("_mean", "", covariable))

# plot environmental covariables
env_cov_plot <- ggplot(EC, aes(x = env, y = value, group = covariable, color = covariable)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ covariable) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(output_folder, "/env_covariables_means.pdf"),
       plot = env_cov_plot, device = "pdf")

# write file with covariables
fwrite(EC, file = paste0(output_folder, "/env_covariables_means.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)



#### get environmental covariables per intervals ----

# define interval window based on latest day present in all environments
n_envs_per_day <- split(env_data[, c("env", "daysFromStart")], env_data[, "daysFromStart"])
n_envs_per_day <- sapply(n_envs_per_day, nrow)
max_date <- max(which(n_envs_per_day == max(n_envs_per_day)))

# get EC
EC_intervals <- W_matrix(env.data = env_data, env.id = 'env', var.id = env_vars,
                         statistic = "mean", by.interval = TRUE,
                         time.window = seq(0, max_date, by = interval_window))
# pdf(file = paste0(output_folder, "/env_covariables_panel_per_intervals.pdf"))
# plot_panel(EC_intervals, title = 'Panel of Environmental Covariables')
# dev.off()

# format data
EC_intervals <- data.frame(EC_intervals, stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "env") %>%
  pivot_longer(contains("_mean"), names_to = "covariable", values_to = "value") %>%
  mutate(covariable = gsub("_mean", "", covariable)) %>%
  separate(col = covariable, sep = "_Interval_", into = c("covariable", "intervals")) %>%
  mutate(intervals = as.numeric(intervals))

# plot environmental covariables
env_cov_plot_int <- ggplot(EC_intervals) +
  geom_line(aes(x = env, y = value, color = intervals, group = intervals), alpha = 0.6) +
  facet_wrap(~ covariable) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = paste0(output_folder, "/env_covariables_means_per_intervals.pdf"),
       plot = env_cov_plot_int, device = "pdf")

# write file with covariables
fwrite(EC_intervals, file = paste0(output_folder, "/env_covariables_means_per_intervals.txt"),
       sep = "\t", quote = FALSE, row.names = FALSE)



#### debug ----

# coords_filename <- "data/env_sites_coord.TEST.csv"
# output_folder <- "data/env_covariables"
# country <- "USA"
# interval_window <- 3
