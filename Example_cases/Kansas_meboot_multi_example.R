##############################################################################
# Kansas_meboot_multi_example.R
#
# Maximum entropy bootstrapping with meboot_multi(): Kansas example
#
# meboot::meboot() bootstraps one series at a time. meboot_multi() loops it
# over several columns, lets you set trim / xmin / xmax limits for each
# column, and returns a list of bootstrapped data frames in long or wide
# format. synth_lasso_boot() uses it internally; this script shows the
# bootstrap step on its own.
#
# Two examples with the `kansas` data from the augsynth package:
#   1. Bootstrapping GDP per capita (gdpcapita) only, as one stacked series.
#   2. Bootstrapping GDP per capita together with population (popestimate),
#      within each state, with a different limit for each column.
#
# Run from the project root (open Synthetic_Control_Claude.Rproj).
# For the full synthetic control analysis, see
# Example_cases/Synth_lasso_boot_example.R and the README.
##############################################################################


## ---- Setup ------------------------------------------------------------------
library(augsynth)
library(tidyverse)
library(meboot)

source("Synthetic_control_functions/synthetic_control_functions_claude.R")

n_reps <- 1000        # lower for a quick run
boot_seed <- 12112025


## ---- Data -------------------------------------------------------------------
# Kansas tax experiment panel from augsynth: 50 states, 1990 Q1 to 2016 Q1
data("kansas")

kansas |>
  select(state, year, qtr, gdpcapita, popestimate) |>
  head()


##############################################################################
## Example 1: GDP per capita only
##############################################################################

## ---- Bootstrapping ------------------------------------------------------------
# gdpcapita is bootstrapped as one stacked series (all states together), and
# each replicate is pivoted wider so every state becomes a column.
#
# The lower limit is set to 15,000 rather than 0. With xmin = 0, meboot
# stretches the lower tail all the way to 0, so some draws fall far below any
# observed GDP per capita (the observed minimum is about 15,000).

kansas_boot_gdp <- meboot_multi(
  data       = kansas,
  boot_cols  = "gdpcapita",
  id_cols    = c("state", "year", "qtr"),
  reps       = n_reps,
  trim       = list(trim = 0.10, xmin = 15000),
  output     = "wide",
  names_from = "state",
  seed       = boot_seed
)

length(kansas_boot_gdp)
dim(kansas_boot_gdp$rep_1)
kansas_boot_gdp$rep_1[1:5, 1:6]

# The same call with output = "long" keeps the original panel layout
kansas_boot_gdp_long <- meboot_multi(
  data      = kansas,
  boot_cols = "gdpcapita",
  id_cols   = c("state", "year", "qtr"),
  reps      = n_reps,
  trim      = list(trim = 0.10, xmin = 15000),
  output    = "long",
  seed      = boot_seed
)

head(kansas_boot_gdp_long$rep_1)

# All draws respect the lower limit
min(sapply(kansas_boot_gdp_long, function(d) min(d$gdpcapita)))


## ---- Check against calling meboot() directly --------------------------------
# With the same seed and settings, replicate 1 from meboot_multi() is
# identical to running meboot() by hand and pivoting the result.
set.seed(boot_seed)
manual_ens <- meboot(as.matrix(select(kansas, gdpcapita)), reps = n_reps,
                     trim = list(trim = 0.10, xmin = 15000, xmax = NULL))$ensemble

manual_rep1 <- cbind.data.frame(state = kansas$state, year = kansas$year,
                                qtr = kansas$qtr, gdp_capita = manual_ens[, 1]) |>
  pivot_wider(names_from = state, values_from = gdp_capita)

all.equal(as.data.frame(manual_rep1), kansas_boot_gdp$rep_1,
          check.attributes = FALSE)


## ---- Pre- and post-treatment matrices ---------------------------------------
# meboot_multi() returns the full period. For a synthetic control with
# treatment in 2012 Q2 (row 90), each replicate is split into a pre-period
# matrix (Kansas + donors) and a post-period matrix (donors only).
# synth_lasso_boot() does this split internally.

kansas_boot_pre <- kansas_boot_gdp |>
  map(~ .x |>
        arrange(year, qtr) |>
        slice(1:89) |>
        relocate(Kansas, .after = qtr) |>
        as.matrix())

kansas_boot_post <- kansas_boot_gdp |>
  map(~ .x |>
        arrange(year, qtr) |>
        slice(90:n()) |>
        select(-Kansas, -year, -qtr) |>
        as.matrix())

dim(kansas_boot_pre[[1]])    # 89 x 52
dim(kansas_boot_post[[1]])   # 16 x 49


## ---- Plot -------------------------------------------------------------------
# Observed Kansas GDP per capita against 25 bootstrap replicates
kansas_obs <- kansas |>
  filter(state == "Kansas") |>
  mutate(time = year + (qtr - 1) / 4)

ex1_plot_df <- kansas_boot_gdp[seq_len(min(25, n_reps))] |>
  imap(~ tibble(rep = .y, time = .x$year + (.x$qtr - 1) / 4,
                gdpcapita = .x$Kansas)) |>
  list_rbind()

ex1_plot <- ggplot() +
  geom_line(data = ex1_plot_df, aes(time, gdpcapita, group = rep),
            color = "lightblue", alpha = 0.6) +
  geom_line(data = kansas_obs, aes(time, gdpcapita),
            color = "black", linewidth = 1) +
  labs(title = "Kansas GDP per capita: observed (black) vs 25 meboot replicates",
       subtitle = "Example 1: stacked series, xmin = 15000",
       x = "Year", y = "GDP per capita") +
  theme_minimal()

ex1_plot


##############################################################################
## Example 2: GDP per capita and population
##############################################################################

## ---- Bootstrapping ------------------------------------------------------------
# Both gdpcapita and popestimate are bootstrapped, each with its own limits.
# group_col = "state" bootstraps each state's series separately, so every
# state keeps its own time dependence. Arrange by time within state first.
#
# Population gets a floor just under the smallest observed state population
# (about 450,000) so no draws fall to implausibly small values.

kansas_sorted <- kansas |> arrange(state, year, qtr)

ex2_trim <- list(
  gdpcapita   = list(trim = 0.10, xmin = 15000),
  popestimate = list(trim = 0.10, xmin = 450000)
)

kansas_boot_multi_long <- meboot_multi(
  data      = kansas_sorted,
  boot_cols = c("gdpcapita", "popestimate"),
  id_cols   = c("state", "year", "qtr"),
  reps      = n_reps,
  trim      = ex2_trim,
  group_col = "state",
  output    = "long",
  seed      = boot_seed
)

head(kansas_boot_multi_long$rep_1)

# In wide format each state gets one column per variable,
# named <variable>_<state>
kansas_boot_multi_wide <- meboot_multi(
  data       = kansas_sorted,
  boot_cols  = c("gdpcapita", "popestimate"),
  id_cols    = c("state", "year", "qtr"),
  reps       = n_reps,
  trim       = ex2_trim,
  group_col  = "state",
  output     = "wide",
  names_from = "state",
  seed       = boot_seed
)

dim(kansas_boot_multi_wide$rep_1)
names(kansas_boot_multi_wide$rep_1)[1:6]


## ---- Plot -------------------------------------------------------------------
# Observed Kansas values against 25 replicates for both variables
ex2_plot_df <- kansas_boot_multi_long[seq_len(min(25, n_reps))] |>
  imap(~ .x |>
         filter(state == "Kansas") |>
         mutate(rep = .y, time = year + (qtr - 1) / 4)) |>
  list_rbind() |>
  pivot_longer(c(gdpcapita, popestimate), names_to = "variable")

ex2_obs_df <- kansas_obs |>
  pivot_longer(c(gdpcapita, popestimate), names_to = "variable")

ex2_plot <- ggplot() +
  geom_line(data = ex2_plot_df, aes(time, value, group = rep),
            color = "lightblue", alpha = 0.6) +
  geom_line(data = ex2_obs_df, aes(time, value),
            color = "black", linewidth = 1) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  labs(title = "Kansas: observed (black) vs 25 meboot replicates",
       subtitle = "Example 2: bootstrapped within state",
       x = "Year", y = NULL) +
  theme_minimal()

ex2_plot


## ---- Summary of the bootstrap distribution ----------------------------------
# Mean and 95% percentile range of the replicates for Kansas in the last period
last_period <- kansas_boot_multi_long |>
  map(~ .x |> filter(state == "Kansas") |> slice_tail(n = 1)) |>
  list_rbind()

last_period |>
  summarise(across(c(gdpcapita, popestimate),
                   list(mean = mean,
                        lo = ~ quantile(.x, 0.025),
                        hi = ~ quantile(.x, 0.975))))

