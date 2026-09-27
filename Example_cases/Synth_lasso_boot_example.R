##############################################################################
# Synth_lasso_boot_example.R
#
# synth_lasso_boot(): LASSO synthetic control with maximum entropy bootstrap
# inference. Argument reference and a worked example on the Kansas tax
# experiment.
#
# Run from the project root (open Synthetic_Control_Claude.Rproj).
# The README has the full write-up of the methods.
#
# Disclaimer: the functions were written by Claude, an AI model made by
# Anthropic, at the direction of the repository owner. They have not been
# independently peer reviewed; check results for your own data.
##############################################################################

## ---- What synth_lasso_boot() does -------------------------------------------
# 1. Fits a LASSO (glmnet) of the treated unit's pre-treatment outcome on the
#    donor units' outcomes, plus any control variables.
# 2. Predicts the treated unit for every period. Observed minus predicted is
#    the ATT; its running total after treatment is the cumulative effect.
# 3. Adds sensitivity bounds and, optionally, jackknife conformal prediction
#    intervals.
# 4. Bootstraps the outcome (and controls) for every unit with maximum
#    entropy bootstrapping, refits the LASSO on each replicate's
#    pre-treatment period at the same penalty, and summarises the effects
#    into bootstrap effect and cumulative-effect intervals. This approach was
#    inspired by Wong et al. (2023), American Journal of Epidemiology
#    192(7): 1166-80, doi:10.1093/aje/kwad061.


## ---- Argument reference -----------------------------------------------------
# Data (long format: one row per unit and period; balanced panel)
#   data              long data frame
#   unit_col          unit name column, e.g. "state"
#   time_cols         column(s) that order time, e.g. c("year", "qtr");
#                     periods are numbered 1, 2, ... = `rowid` in the output
#   outcome_col       outcome column, e.g. "gdpcapita"
#   treated_unit      treated unit, e.g. "Kansas"
#   treatment_period  rowid of the FIRST post-treatment period (Kansas: 90)
#
# Control variables
#   controls          extra numeric columns, e.g. "popestimate"; enter the
#                     LASSO as <control>_<unit> columns and are bootstrapped
#                     with the outcome (default NULL)
#   control_units     "all" (default), "donors" or "treated". The treated
#                     unit's own post-treatment controls may be affected by
#                     the treatment; "donors" avoids that.
#
# LASSO fit
#   intercept         include an intercept (default TRUE)
#   lower_constr      lower limit(s): one number for all donors, or one
#                     unnamed default plus named columns, e.g.
#                     c(0, Texas = -0.2). Must be <= 0. (default 0)
#   upper_constr      upper limit(s), same format. Must be >= 0. (default 1)
#   control_lower     default lower limit for controls (default -Inf)
#   control_upper     default upper limit for controls (default Inf)
#   unpenalized       columns kept in the model without a penalty, e.g.
#                     "Colorado" or "popestimate_Kansas" (default NULL)
#   standardize       glmnet standardize; TRUE rescales EVERY predictor,
#                     donors included, and changes the donor weights
#                     (default FALSE)
#   scale_controls    rescale only the controls to the donors' pre-period
#                     spread; coefficients and limits stay in original
#                     units (default FALSE)
#   nlambda           lambda path length for cross-validation (default 1000)
#   lambda_choice     "min_path" (smallest lambda on the path; default),
#                     "lambda.min" or "lambda.1se"
#   seed              seed for the cross-validation folds (required)
#
# Inference on the observed fit
#   M                 sensitivity-bound multiplier (default 1)
#   jack_conform      add jackknife conformal prediction intervals
#                     (default FALSE)
#   conform_level     coverage of the conformal intervals (default 0.95)
#
# Bootstrap
#   boot              run the bootstrap (default TRUE)
#   reps              number of replicates (default 1000)
#   trim              meboot trim/xmin/xmax: one list for all variables or a
#                     list named by variable
#   boot_by_unit      FALSE (default) = one stacked series per variable;
#                     TRUE = each unit bootstrapped separately
#   boot_nlambda      lambda path length for the bootstrap refits
#                     (default 100)
#   conf_int          level of the bootstrap intervals (default 0.95)
#   boot_seed         seed for the bootstrap (default = seed)
#   keep_boot         also return the replicates (default FALSE)
#   ...               passed to meboot::meboot()
#
# Output (a list)
#   effects    one row per period: rowid, time columns, period,
#              treated_observed, treated_predicted, att, att_post, cum_att,
#              bounds_lo, bounds_hi, pred_jack_lo, pred_jack_hi (if
#              jack_conform), eff_low, eff_high, average_eff, cum_eff_low,
#              cum_eff_high, average_cum_effect, treated_unit
#   fit_stats  CV MSE, lambda used, pre-period MSE, M bound
#   weights    non-zero coefficients (original units), Type, Unpenalized
#   model      the cv.glmnet object
#   boot       bootstrap effect and cumulative-effect matrices
#   settings   settings used, incl. each predictor's limits, penalty factor
#              and scale factor
#
# Sensitivity bounds (M bounds)
#   B = M x the largest absolute gap between observed and predicted in the
#   pre-treatment period. Each post-treatment ATT is shown as ATT +/- B. If
#   the band excludes zero, the effect is larger than the worst
#   pre-treatment miss (scaled by M), so fit error of that size cannot
#   explain it. The M at which a quarter's band first touches zero is
#   |ATT| / max pre-treatment gap. The bounds are a benchmark, not a
#   confidence interval.


## ---- Setup ------------------------------------------------------------------
library(augsynth)
library(tidyverse)
library(glmnet)
library(meboot)

source("Synthetic_control_functions/synthetic_control_functions_claude.R")

n_reps <- 1000   # lower for a quick run


## ---- Data -------------------------------------------------------------------
# Kansas tax experiment: in 2012 Kansas passed large income tax cuts. The
# `kansas` dataset in the augsynth package is a quarterly panel of the 50
# US states from 1990 Q1 to 2016 Q1. Treatment starts in 2012 Q2 (rowid 90).
data("kansas")

kansas |>
  select(state, year, qtr, gdpcapita, popestimate) |>
  head()


##############################################################################
## Model 1: outcome only, with jackknife intervals
##############################################################################

kansas_m1 <- synth_lasso_boot(
  data             = kansas,
  unit_col         = "state",
  time_cols        = c("year", "qtr"),
  outcome_col      = "gdpcapita",
  treated_unit     = "Kansas",
  treatment_period = 90,
  seed             = 12112025,
  M                = 1,
  jack_conform     = TRUE,
  conform_level    = 0.95,
  reps             = n_reps,
  trim             = list(trim = 0.10, xmin = 15000),
  conf_int         = 0.95
)

kansas_m1$fit_stats
kansas_m1$weights


## ---- Jackknife conformal prediction intervals --------------------------------
# Observed Kansas against the synthetic prediction and its 95% jackknife
# prediction interval. outside_interval = TRUE when observed Kansas falls
# outside the interval.
kansas_m1$effects |>
  filter(rowid >= 90) |>
  mutate(outside_interval = treated_observed < pred_jack_lo |
                            treated_observed > pred_jack_hi) |>
  select(period, treated_observed, treated_predicted,
         pred_jack_lo, pred_jack_hi, outside_interval)


# NOTE: in this example the jackknife prediction intervals are far too
# small. They are much narrower than the bootstrap intervals and the
# sensitivity band (see the widths below), because:
#   1. consecutive quarters are highly correlated, so leaving one quarter
#      out barely changes the fit and the leave-one-out residuals are much
#      smaller than genuine forecast errors;
#   2. with 49 donor states and a small penalty, the LASSO fits the
#      pre-treatment quarters very closely;
#   3. the interval has the same width in every post-treatment quarter,
#      although forecast error should grow further past 2012.
# outside_interval = TRUE is therefore not reliable evidence of an effect
# here; rely on the bootstrap intervals and sensitivity bounds instead.

# Average width of each interval over the post-treatment quarters
kansas_m1$effects |>
  filter(rowid >= 90) |>
  summarise(
    jackknife_prediction_interval = mean(pred_jack_hi - pred_jack_lo),
    bootstrap_effect_interval     = mean(eff_high - eff_low),
    sensitivity_band_M1           = mean(bounds_hi - bounds_lo)
  )


## ---- Sensitivity bounds -------------------------------------------------------
# Largest absolute pre-treatment gap between observed and predicted
max_pre_gap <- kansas_m1$effects |>
  filter(rowid < 90) |>
  summarise(max(abs(att))) |>
  pull()

max_pre_gap

# ATT with its bounds at M = 1, and the M at which each band includes zero
kansas_m1$effects |>
  filter(rowid >= 90) |>
  mutate(M_breakdown = abs(att) / max_pre_gap) |>
  select(period, att, bounds_lo, bounds_hi, M_breakdown)


## ---- Bootstrap effect and cumulative-effect intervals ------------------------
kansas_m1$effects |>
  filter(rowid >= 90) |>
  select(period, att, eff_low, eff_high,
         cum_att, cum_eff_low, cum_eff_high)


##############################################################################
## Model 2: population as a control
##############################################################################
# Donor-state population is added as a control and rescaled to the donors'
# scale. Colorado is kept in the model without a penalty, cross-validation
# chooses the penalty, and each state is bootstrapped separately.

kansas_m2 <- synth_lasso_boot(
  data             = kansas,
  unit_col         = "state",
  time_cols        = c("year", "qtr"),
  outcome_col      = "gdpcapita",
  treated_unit     = "Kansas",
  treatment_period = 90,
  controls         = "popestimate",
  control_units    = "donors",
  scale_controls   = TRUE,
  lower_constr     = 0,
  upper_constr     = 1,
  unpenalized      = "Colorado",
  lambda_choice    = "lambda.min",
  seed             = 12112025,
  M                = 1,
  jack_conform     = TRUE,
  reps             = n_reps,
  trim             = list(gdpcapita   = list(trim = 0.10, xmin = 15000),
                          popestimate = list(trim = 0.10, xmin = 450000)),
  boot_by_unit     = TRUE,
  conf_int         = 0.95
)

kansas_m2$fit_stats
kansas_m2$weights

# Limits, penalty factors and scale factors for the unpenalized and
# rescaled predictors
kansas_m2$settings$predictors |>
  filter(penalty_factor == 0 | scale_factor != 1) |>
  head(10)

# Post-treatment rows: jackknife intervals (far too small here, as in
# Model 1), sensitivity bounds and bootstrap intervals
kansas_m2$effects |>
  filter(rowid >= 90) |>
  select(period, treated_observed, treated_predicted, pred_jack_lo, pred_jack_hi,
         att, bounds_lo, bounds_hi, eff_low, eff_high)


##############################################################################
## Plots (Model 1)
##############################################################################

plot_df <- kansas_m1$effects |> filter(year >= 2000)

# Label only the first quarter of each year (e.g. 2016_1) on the x axis
q1_breaks <- function(x) x[grepl("_1$", x)]


## ---- Observed vs synthetic Kansas, with jackknife prediction intervals -------
plot_fit <- ggplot(plot_df, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = pred_jack_lo, ymax = pred_jack_hi),
              fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = treated_predicted, color = "Synthetic")) +
  geom_line(aes(y = treated_observed, color = "Observed")) +
  geom_vline(xintercept = "2012_2", linetype = 2, color = "gray40") +
  scale_color_manual(values = c(Observed = "black", Synthetic = "blue")) +
  scale_x_discrete(breaks = q1_breaks) +
  labs(title = "Kansas GDP per capita: observed vs synthetic",
       subtitle = "Shaded: 95% jackknife conformal prediction interval (far too narrow in this example)",
       x = "Period", y = "GDP per capita", color = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_fit


## ---- ATT with sensitivity bounds ---------------------------------------------
plot_bounds <- ggplot(plot_df, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = bounds_lo, ymax = bounds_hi), fill = "tomato", alpha = 0.2) +
  geom_line(aes(y = att), color = "blue") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = "2012_2", linetype = 2, color = "gray40") +
  scale_x_discrete(breaks = q1_breaks) +
  labs(title = "ATT with sensitivity bounds (M = 1)",
       subtitle = "Band = ATT \u00b1 the largest absolute pre-treatment gap",
       x = "Period", y = "ATT (GDP per capita)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_bounds


## ---- Per-period effect with bootstrap intervals ------------------------------
plot_effects <- ggplot(plot_df, aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = eff_low, ymax = eff_high), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = att), color = "blue") +
  geom_point(aes(y = att), color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = "2012_2", linetype = 2, color = "gray40") +
  scale_x_discrete(breaks = q1_breaks) +
  labs(title = "ATT: observed minus synthetic Kansas (95% bootstrap intervals)",
       x = "Period", y = "ATT (GDP per capita)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_effects


## ---- Cumulative effect with bootstrap intervals ------------------------------
plot_cumulative <- kansas_m1$effects |>
  filter(rowid >= 90) |>
  ggplot(aes(x = period, group = 1)) +
  geom_ribbon(aes(ymin = cum_eff_low, ymax = cum_eff_high), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = cum_att), color = "blue") +
  geom_point(aes(y = cum_att), color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = "Cumulative effect on Kansas GDP per capita (95% bootstrap intervals)",
       x = "Period", y = "Cumulative effect") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

plot_cumulative


##############################################################################
## Comparing the two models: cumulative effect in the final quarter
##############################################################################

bind_rows(
  kansas_m1$effects |> slice_tail(n = 1) |> mutate(model = "Model 1: outcome only"),
  kansas_m2$effects |> slice_tail(n = 1) |> mutate(model = "Model 2: + population")
) |>
  select(model, period, cum_att, cum_eff_low, cum_eff_high)


## ---- Reference ----------------------------------------------------------------
# Wong, Anabelle, Sarah C Kramer, Marco Piccininni, Jessica L Rohmann,
# Tobias Kurth, Sylvie Escolano, Ulrike Grittner, and Matthieu Domenech De
# Celles. 2023. "Using LASSO Regression to Estimate the Population-Level
# Impact of Pneumococcal Conjugate Vaccines." American Journal of
# Epidemiology 192(7): 1166-80. doi:10.1093/aje/kwad061.

