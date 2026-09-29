##############################################################################
# synthetic_control_functions_claude.R
#
# LASSO synthetic control with maximum entropy bootstrap inference.
# Written by Claude (Anthropic) at the direction of the repository owner.
#
# Contents:
#   1. meboot_multi()      - helper: maximum entropy bootstrap of multiple
#                            columns (used inside synth_lasso_boot())
#   2. synth_lasso_boot()  - LASSO synthetic control with sensitivity bounds,
#                            jackknife conformal intervals, and bootstrap
#                            effect / cumulative-effect intervals
#   3. synth_lasso_placebo() - the same model with exact placebo (permutation)
#                            tests and CIs instead of the bootstrap
#
# Usage:
#   source("Synthetic_control_functions/synthetic_control_functions_claude.R")
##############################################################################

## ===========================================================================
## 1. meboot_multi()
## ===========================================================================

##############################################################################
# meboot_multi(): Maximum entropy bootstrap for multiple columns
#
# meboot::meboot() only bootstraps one series at a time. This function loops
# meboot over several columns, lets you set trim / limits (xmin, xmax) per
# column, and returns a list of `reps` bootstrapped data frames in either
# long or wide format.
#
# Arguments
#   data        data frame (long format, e.g. the `kansas` panel)
#   boot_cols   character vector of numeric columns to bootstrap
#   id_cols     character vector of columns carried along unchanged
#               (e.g. c("state", "year", "qtr")). Default: all columns not
#               in boot_cols.
#   reps        number of bootstrap replicates (passed to meboot)
#   trim        EITHER a single list applied to every column, e.g.
#                 list(trim = 0.10, xmin = 0, xmax = NULL)
#               OR a named list with one list per column, e.g.
#                 list(gdpcapita   = list(trim = 0.10, xmin = 0),
#                      lngdpcapita = list(trim = 0.10, xmin = NULL))
#               Settings within each list:
#                 trim      trimmed-mean proportion (default 0.10)
#                 xmin/xmax absolute lower / upper tail limits. One number
#                           for every series, or - with group_col - a
#                           vector with at most one unnamed default and
#                           values named by group, e.g.
#                           xmin = c(15000, Alaska = 30000)
#                 xmin_rel  lower limit as a fraction of EACH series' own
#                           minimum, e.g. 0.9 (positive data only)
#                 xmax_rel  upper limit as a multiple of each series' own
#                           maximum, e.g. 1.1 (positive data only)
#               Give xmin or xmin_rel (xmax or xmax_rel), not both. Unset
#               limits use the meboot defaults; columns not listed get the
#               defaults. A limit inside the observed range is an error.
#   group_col   optional column name. NULL (default) bootstraps each column
#               as one stacked series, exactly like the Kansas example.
#               If set (e.g. "state"), each group's series is bootstrapped
#               separately. Rows are used in their current order within each
#               group, so arrange the data by time first.
#   output      "long" (default) or "wide"
#   names_from  column to pivot into new columns when output = "wide"
#               (e.g. "state"). Required for wide output.
#   names_sep   separator for wide column names when more than one column
#               is bootstrapped (e.g. "gdpcapita_Kansas")
#   seed        optional seed, set once before any bootstrapping
#   shared_draws  with group_col only: TRUE gives every group the same
#               random numbers within each column and replicate, so groups
#               are nudged up or down together (by rank). FALSE (default)
#               draws independently for every group.
#   as_matrix   if TRUE, each list element is converted with as.matrix()
#               (useful for glmnet). All columns must be numeric for a
#               numeric matrix.
#   ...         further arguments passed to meboot::meboot()
#               (e.g. reachbnd, expand.sd, force.clt, sym)
#
# Value
#   A list of length `reps` (named rep_1, rep_2, ...). Each element is a
#   data frame (or matrix) holding id_cols plus the bootstrapped columns.
##############################################################################

meboot_multi <- function(data,
                         boot_cols,
                         id_cols = NULL,
                         reps = 999,
                         trim = list(trim = 0.10, xmin = NULL, xmax = NULL),
                         group_col = NULL,
                         output = c("long", "wide"),
                         names_from = NULL,
                         names_sep = "_",
                         seed = NULL,
                         as_matrix = FALSE,
                         shared_draws = FALSE,
                         ...) {
  
  output <- match.arg(output)
  data <- as.data.frame(data)
  
  ## ---- Checks ------------------------------------------------------------
  missing_cols <- setdiff(c(boot_cols, id_cols, group_col, names_from), names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }
  
  not_numeric <- boot_cols[!vapply(data[boot_cols], is.numeric, logical(1))]
  if (length(not_numeric) > 0) {
    stop("boot_cols must be numeric: ", paste(not_numeric, collapse = ", "))
  }
  
  has_na <- boot_cols[vapply(data[boot_cols], anyNA, logical(1))]
  if (length(has_na) > 0) {
    stop("meboot cannot handle missing values. NAs found in: ",
         paste(has_na, collapse = ", "))
  }
  
  if (output == "wide" && is.null(names_from)) {
    stop("output = 'wide' requires names_from (e.g. names_from = 'state').")
  }
  
  # Default id columns: everything that is not being bootstrapped
  if (is.null(id_cols)) {
    id_cols <- setdiff(names(data), boot_cols)
  }
  
  # Grouping and pivoting columns must be kept with the output
  id_cols <- unique(c(id_cols, group_col, names_from))
  
  ## ---- Build a trim list for each column ---------------------------------
  # Keys: trim, xmin, xmax (absolute limits) and xmin_rel, xmax_rel (limits
  # relative to each series' own minimum / maximum).
  trim_default <- list(trim = 0.10, xmin = NULL, xmax = NULL,
                       xmin_rel = NULL, xmax_rel = NULL)
  trim_keys <- names(trim_default)
  
  # One list for all columns if every name is a trim key
  is_global_trim <- !is.null(names(trim)) && all(names(trim) %in% trim_keys)
  
  if (is_global_trim) {
    trim_list <- rep(list(utils::modifyList(trim_default, trim, keep.null = TRUE)),
                     length(boot_cols))
    names(trim_list) <- boot_cols
  } else {
    bad_names <- setdiff(names(trim), boot_cols)
    if (is.null(names(trim)) || length(bad_names) > 0) {
      stop("trim must be a single list (", paste(trim_keys, collapse = "/"),
           ") or a list named by boot_cols. Unrecognised name(s): ",
           paste(bad_names, collapse = ", "))
    }
    trim_list <- lapply(boot_cols, function(col) {
      if (col %in% names(trim)) {
        bad_keys <- setdiff(names(trim[[col]]), trim_keys)
        if (length(bad_keys) > 0) {
          stop("trim for ", col, ": unknown setting(s) ", paste(bad_keys, collapse = ", "),
               ". Use ", paste(trim_keys, collapse = ", "), ".")
        }
        utils::modifyList(trim_default, trim[[col]], keep.null = TRUE)
      } else {
        trim_default
      }
    })
    names(trim_list) <- boot_cols
  }
  
  groups <- if (is.null(group_col)) NULL else unique(as.character(data[[group_col]]))
  
  # Check limits named by unit before any bootstrapping
  for (col in boot_cols) {
    for (key in c("xmin", "xmax")) {
      val <- trim_list[[col]][[key]]
      if (!is.null(trim_list[[col]][[paste0(key, "_rel")]]) && !is.null(val)) {
        stop("trim for ", col, ": give ", key, " or ", key, "_rel, not both.")
      }
      nm <- names(val)
      if (!is.null(nm) && any(nm != "")) {
        if (is.null(group_col)) {
          stop("trim for ", col, ": limits named by unit (", key, ") need group_col ",
               "(boot_by_unit = TRUE in synth_lasso_boot()).")
        }
        if (sum(nm == "") > 1) {
          stop("trim for ", col, ": ", key, " can have at most one unnamed default value.")
        }
        unknown <- setdiff(nm[nm != ""], groups)
        if (length(unknown) > 0) {
          stop("trim for ", col, ": ", key, " names unit(s) not in the data: ",
               paste(unknown, collapse = ", "))
        }
      }
    }
  }
  
  # The meboot trim list for one series (one column, and one group if grouped)
  series_trim <- function(col, x, group = NULL) {
    spec <- trim_list[[col]]
    where <- if (is.null(group)) col else paste0(col, " (", group, ")")
    
    pick <- function(val) {
      if (is.null(val)) return(NULL)
      nm <- names(val)
      if (is.null(nm) || all(nm == "")) return(unname(val[1]))
      if (!is.null(group) && group %in% nm) return(unname(val[group]))
      if (any(nm == "")) return(unname(val[nm == ""][1]))
      NULL
    }
    
    xmin <- pick(spec$xmin)
    xmax <- pick(spec$xmax)
    
    if (!is.null(spec$xmin_rel)) {
      if (any(x <= 0)) stop("xmin_rel needs strictly positive values; ", where, " has values <= 0. Use xmin instead.")
      if (spec$xmin_rel > 1) stop("xmin_rel must be <= 1 (a fraction of the series minimum).")
      xmin <- spec$xmin_rel * min(x)
    }
    if (!is.null(spec$xmax_rel)) {
      if (any(x <= 0)) stop("xmax_rel needs strictly positive values; ", where, " has values <= 0. Use xmax instead.")
      if (spec$xmax_rel < 1) stop("xmax_rel must be >= 1 (a multiple of the series maximum).")
      xmax <- spec$xmax_rel * max(x)
    }
    
    if (!is.null(xmin) && xmin > min(x)) {
      stop("xmin for ", where, " (", signif(xmin, 6), ") is above its lowest value (",
           signif(min(x), 6), ").")
    }
    if (!is.null(xmax) && xmax < max(x)) {
      stop("xmax for ", where, " (", signif(xmax, 6), ") is below its highest value (",
           signif(max(x), 6), ").")
    }
    list(trim = spec$trim, xmin = xmin, xmax = xmax)
  }
  
  ## ---- Run meboot for each column (and group) ----------------------------
  if (shared_draws && is.null(group_col)) {
    warning("shared_draws only applies when bootstrapping by group (group_col); ignored.")
    shared_draws <- FALSE
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # With shared draws, every group of a column restarts from the same seed, so
  # meboot draws the same random numbers for each group in each replicate.
  col_seeds <- if (shared_draws) sample.int(.Machine$integer.max, length(boot_cols)) else NULL
  
  run_meboot <- function(x, col_trim) {
    ens <- meboot::meboot(x, reps = reps, trim = col_trim, ...)$ensemble
    as.matrix(ens)  # rows = observations, columns = replicates
  }
  
  ensembles <- lapply(seq_along(boot_cols), function(j) {
    col <- boot_cols[j]
    x <- data[[col]]
    
    if (is.null(group_col)) {
      return(run_meboot(x, series_trim(col, x)))
    }
    
    # Bootstrap each group separately and put results back in original rows
    ens <- matrix(NA_real_, nrow = length(x), ncol = reps)
    row_index <- split(seq_along(x), as.character(data[[group_col]]))
    for (g in names(row_index)) {
      rows <- row_index[[g]]
      if (shared_draws) set.seed(col_seeds[j])
      ens[rows, ] <- run_meboot(x[rows], series_trim(col, x[rows], g))
    }
    ens
  })
  names(ensembles) <- boot_cols
  
  ## ---- Assemble one data frame per replicate -----------------------------
  base_df <- data[id_cols]
  
  if (output == "long") {
    boot_list <- lapply(seq_len(reps), function(r) {
      df <- base_df
      for (col in boot_cols) {
        df[[col]] <- ensembles[[col]][, r]
      }
      if (as_matrix) df <- as.matrix(df)
      df
    })
  }
  
  if (output == "wide") {
    wide_id_cols <- setdiff(id_cols, names_from)
    
    # Each id / names_from combination must identify a single row
    if (anyDuplicated(data[c(wide_id_cols, names_from)]) > 0) {
      stop("id_cols and names_from do not uniquely identify rows, so the ",
           "data cannot be pivoted wider. Check id_cols.")
    }
    
    # Pivot ONCE: a template whose cells hold the original row numbers.
    # Every replicate is then filled in by indexing, instead of calling
    # pivot_wider() reps times.
    template <- base_df[c(wide_id_cols, names_from)]
    template$.row_index <- seq_len(nrow(template))
    template <- tidyr::pivot_wider(template,
                                   id_cols = tidyselect::all_of(wide_id_cols),
                                   names_from = tidyselect::all_of(names_from),
                                   values_from = ".row_index")
    template <- as.data.frame(template)
    
    unit_names <- setdiff(names(template), wide_id_cols)
    row_index <- as.matrix(template[unit_names])  # NA where a unit has no row
    id_part <- template[wide_id_cols]
    
    # Column names follow pivot_wider(): unit names for one column,
    # <column><names_sep><unit> when several columns are bootstrapped
    wide_names <- if (length(boot_cols) == 1) {
      unit_names
    } else {
      unlist(lapply(boot_cols, function(col) paste(col, unit_names, sep = names_sep)))
    }
    
    boot_list <- lapply(seq_len(reps), function(r) {
      vals <- do.call(cbind, lapply(boot_cols, function(col) {
        v <- ensembles[[col]][, r]
        matrix(v[row_index], nrow = nrow(row_index))
      }))
      colnames(vals) <- wide_names
      
      if (as_matrix) {
        return(cbind(as.matrix(id_part), vals))
      }
      df <- cbind(id_part, as.data.frame(vals, stringsAsFactors = FALSE))
      rownames(df) <- NULL
      df
    })
  }
  
  names(boot_list) <- paste0("rep_", seq_len(reps))
  
  boot_list
}


## ===========================================================================
## 2. synth_lasso_boot()
## ===========================================================================

##############################################################################
# synth_lasso_boot(): LASSO synthetic control with maximum entropy bootstrap
#                     effect and cumulative-effect intervals
#
# Combines three steps in one call:
#   1. synth_lasso_function()        - LASSO synthetic control on the observed
#                                      data (weights, predictions, ATT,
#                                      sensitivity bounds, optional jackknife
#                                      conformal intervals)
#   2. meboot_multi()                - maximum entropy bootstrap of the outcome
#                                      (and any control variables) for every
#                                      unit, split into pre / post
#   3. bootstrap_cum_effects_fixed() - refit the LASSO on each bootstrap
#                                      replicate at the fixed lambda, and
#                                      summarise the effects and cumulative
#                                      effects
#
# Uses meboot_multi(), defined earlier in this file.
#
# Data: LONG format panel, one row per unit and time period
#   (e.g. the augsynth `kansas` data).
#
# ---- Data arguments --------------------------------------------------------
#   data              long data frame
#   unit_col          column holding unit names (e.g. "state")
#   time_cols         column(s) that order time (e.g. c("year", "qtr")).
#                     Periods are numbered 1, 2, ... in sorted order; this
#                     number is the `rowid` column in the output.
#   outcome_col       outcome column (e.g. "gdpcapita")
#   treated_unit      name of the treated unit (e.g. "Kansas")
#   treatment_period  the FIRST post-treatment period, given as the value(s)
#                     of time_cols in that period:
#                       - one time column:  c(year = 2013), or just 2013
#                       - several columns:  c(year = 2012, qtr = 2), or
#                                           unnamed in time_cols order,
#                                           c(2012, 2)
#                     Values are matched exactly against the data (a list
#                     can hold mixed types, e.g. a Date). Every period
#                     before it is pre-treatment; it and every later period
#                     are post-treatment. At least 2 pre-treatment periods
#                     are required.
#
# ---- Control variables ------------------------------------------------------
#   controls          optional character vector of extra numeric columns in
#                     `data` (e.g. c("popestimate", "avg_wkly_wage")). Each
#                     is pivoted wide and added to the LASSO as predictors
#                     named <control>_<unit> (e.g. "popestimate_Texas"), and
#                     is bootstrapped by meboot along with the outcome.
#   control_units     which units' control series enter the LASSO:
#                       "all" (default), "donors", or "treated".
#                     NOTE: the treated unit's own controls in the post
#                     period may themselves be affected by the treatment.
#                     "donors" avoids that.
#
# ---- LASSO arguments --------------------------------------------------------
#   intercept       include an intercept (default TRUE)
#   lower_constr    lower limit(s) for the coefficients. Either
#                     - one number: applies to every donor column, or
#                     - a vector with at most one unnamed value (the donor
#                       default) and named values for specific columns,
#                       e.g. c(0, Texas = -0.5, popestimate_Kansas = -Inf)
#                   glmnet requires every lower limit to be <= 0.
#                   Default 0.
#   upper_constr    upper limit(s), same format as lower_constr. glmnet
#                   requires every upper limit to be >= 0. Default 1.
#   control_lower   default lower limit for control columns (default -Inf)
#   control_upper   default upper limit for control columns (default Inf)
#                   Named entries in lower_constr / upper_constr override
#                   these for individual control columns.
#   unpenalized     optional character vector of predictor columns that are
#                   NOT penalized (penalty.factor = 0), i.e. always kept in
#                   the model: donor names (e.g. "Texas") and/or control
#                   columns (e.g. "popestimate_Kansas"). Their limits still
#                   apply.
#   standardize     passed to glmnet (default FALSE). TRUE standardizes
#                   EVERY predictor, donors included, which changes the
#                   donor weights even without controls. Coefficients and
#                   limits stay on the original scale.
#   scale_controls  TRUE rescales only the control columns so their
#                   pre-period standard deviation matches the average
#                   pre-period SD of the donor columns (default FALSE).
#                   Donors are untouched; control coefficients and limits
#                   are reported / applied in the original units. The same
#                   scale factors are used for every bootstrap refit.
#                   Ignored (with a warning) when standardize = TRUE.
#   nlambda         lambda path length for cv.glmnet (default 1000)
#   lambda_choice   which lambda is used for predictions, weights and the
#                   bootstrap refits:
#                     "min_path"   smallest lambda on the path,
#                                  min(cv.fit$lambda). This is what the
#                                  existing scripts use, and it does not
#                                  depend on the CV folds.
#                     "lambda.min" lambda with the lowest CV error
#                     "lambda.1se" largest lambda within 1 SE of the minimum
#   seed            seed for the cv.glmnet folds
#
# ---- Inference on the observed fit ------------------------------------------
#   M               sensitivity-bound multiplier on the largest absolute
#                   pre-period residual (default 1)
#   jack_conform    TRUE adds jackknife conformal prediction intervals for the
#                   post period (needs conformalInference)
#                   The intervals come from conformalInference::conformal.pred.jack().
#                   That function mishandles a single prediction column (it then
#                   uses only the first leave-one-out residual, whatever the level),
#                   so the prediction function here returns two identical columns
#                   and the first is used.
#   conform_level   coverage of the conformal intervals (default 0.95)
#
# ---- Bootstrap arguments ------------------------------------------------------
#   boot            TRUE (default) runs the bootstrap
#   reps            number of bootstrap replicates (default 1000)
#   trim            meboot limits, passed to meboot_multi(). One list for
#                   every bootstrapped column, or a list named by column.
#                   Settings: trim, xmin, xmax, xmin_rel, xmax_rel, e.g.
#                     list(gdpcapita   = list(xmin = 15000),
#                          popestimate = list(xmin = 450000))
#                   With boot_by_unit = TRUE, limits can differ by unit:
#                     xmin = c(15000, Alaska = 30000)   named by unit, or
#                     xmin_rel = 0.9   90% of EACH unit's own minimum
#                     xmax_rel = 1.1   110% of each unit's own maximum
#                   (relative limits need strictly positive data)
#   boot_by_unit    FALSE (default) bootstraps each variable as one stacked
#                   series in the row order of `data` (reproduces the Kansas
#                   script). TRUE bootstraps each unit's series separately.
#   shared_draws    with boot_by_unit = TRUE: TRUE uses the same random
#                   numbers for every unit (within each variable and
#                   replicate), so units are nudged up or down together by
#                   rank and keep their co-movement. FALSE (default) draws
#                   independently for each unit. Ignored, with a warning,
#                   when boot_by_unit = FALSE.
#   boot_nlambda    lambda path length for the bootstrap glmnet refits
#                   (default 100, as in bootstrap_cum_effects_fixed)
#   conf_int        interval level for bootstrap effects (default 0.95)
#   boot_seed       seed for meboot (default = seed)
#   keep_boot       TRUE also returns the bootstrap replicate matrices
#   ...             further arguments passed to meboot::meboot()
#
# ---- Value --------------------------------------------------------------------
# A list:
#   effects    one row per period: rowid, time columns, period label,
#              post (TRUE from the first treated period on),
#              treated_observed, treated_predicted, att, att_post, cum_att,
#              bounds_lo / bounds_hi, (pred_jack_lo / pred_jack_hi),
#              eff_low / eff_high / average_eff (bootstrap effect intervals),
#              cum_eff_low / cum_eff_high / average_cum_effect (bootstrap
#              cumulative-effect intervals). Pre-period inference columns
#              are NA.
#   fit_stats  CV MSE, lambda used, pre-period MSE of the fit, M bound
#   weights    non-zero coefficients at the lambda used, with their type
#              (intercept / donor / control) and whether they were
#              unpenalized
#   model      the cv.glmnet object
#   boot       list: effects and cum_effects (post periods x reps matrices),
#              and the replicates if keep_boot = TRUE (NULL if boot = FALSE)
#   settings   the main settings used, including the limits, penalty
#              factors and scale factors for every predictor
##############################################################################

synth_lasso_boot <- function(data,
                             unit_col,
                             time_cols,
                             outcome_col,
                             treated_unit,
                             treatment_period,
                             controls = NULL,
                             control_units = c("all", "donors", "treated"),
                             intercept = TRUE,
                             lower_constr = 0,
                             upper_constr = 1,
                             control_lower = -Inf,
                             control_upper = Inf,
                             unpenalized = NULL,
                             standardize = FALSE,
                             scale_controls = FALSE,
                             nlambda = 1000,
                             lambda_choice = c("min_path", "lambda.min", "lambda.1se"),
                             seed,
                             M = 1,
                             jack_conform = FALSE,
                             conform_level = 0.95,
                             boot = TRUE,
                             reps = 1000,
                             trim = list(trim = 0.10, xmin = NULL, xmax = NULL),
                             boot_by_unit = FALSE,
                             shared_draws = FALSE,
                             boot_nlambda = 100,
                             conf_int = 0.95,
                             boot_seed = seed,
                             keep_boot = FALSE,
                             ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  control_units <- match.arg(control_units)
  data <- as.data.frame(data)
  vars <- c(outcome_col, controls)
  
  ## ---- Checks ------------------------------------------------------------
  if (missing(seed)) stop("Please supply a seed for reproducibility.")
  
  missing_cols <- setdiff(c(unit_col, time_cols, vars), names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }
  if (outcome_col %in% controls) stop("outcome_col cannot also be a control.")
  if (anyDuplicated(controls)) stop("controls contains duplicates.")
  
  not_numeric <- vars[!vapply(data[vars], is.numeric, logical(1))]
  if (length(not_numeric) > 0) {
    stop("Outcome and controls must be numeric: ", paste(not_numeric, collapse = ", "))
  }
  has_na <- vars[vapply(data[vars], anyNA, logical(1))]
  if (length(has_na) > 0) {
    stop("Missing values in: ", paste(has_na, collapse = ", "))
  }
  if (!treated_unit %in% data[[unit_col]]) {
    stop("treated_unit '", treated_unit, "' not found in ", unit_col, ".")
  }
  if (boot && !exists("meboot_multi", mode = "function")) {
    stop("meboot_multi() is not loaded. Run: ",
         "source(\"Synthetic_control_functions/synthetic_control_functions_claude.R\")")
  }
  if (jack_conform && !requireNamespace("conformalInference", quietly = TRUE)) {
    stop("jack_conform = TRUE needs the conformalInference package.")
  }
  
  ## ---- Number the time periods (rowid) -----------------------------------
  time_key <- do.call(paste, c(data[time_cols], sep = "\r"))
  periods <- unique(data[time_cols])
  periods <- periods[do.call(order, unname(as.list(periods))), , drop = FALSE]
  rownames(periods) <- NULL
  period_key <- do.call(paste, c(periods, sep = "\r"))
  data$.rowid <- match(time_key, period_key)
  n_periods <- nrow(periods)
  
  ## ---- Locate the first treated period -----------------------------------
  period_labels <- do.call(paste, c(periods, sep = "_"))
  tp <- as.list(treatment_period)
  
  if (length(tp) != length(time_cols)) {
    stop("treatment_period needs one value for each time column (",
         paste(time_cols, collapse = ", "), "), e.g. ",
         if (length(time_cols) == 1) paste0("c(", time_cols, " = 2013)")
         else paste0("c(", paste(time_cols, "= ...", collapse = ", "), ")"), ".")
  }
  tp_names <- names(tp)
  if (is.null(tp_names) || all(tp_names == "")) {
    names(tp) <- time_cols
  } else if (any(tp_names == "") || !setequal(tp_names, time_cols)) {
    stop("treatment_period names must match time_cols (",
         paste(time_cols, collapse = ", "), "), or leave them all unnamed.")
  } else {
    tp <- tp[time_cols]
  }
  
  is_tp <- Reduce(`&`, lapply(time_cols, function(tc) {
    as.character(periods[[tc]]) == as.character(tp[[tc]])
  }))
  if (!any(is_tp)) {
    stop("treatment_period (",
         paste(names(tp), vapply(tp, as.character, character(1)), sep = " = ", collapse = ", "),
         ") does not match any period in the data. Periods run from ",
         period_labels[1], " to ", period_labels[n_periods], ".")
  }
  treatment_rowid <- which(is_tp)
  if (treatment_rowid < 3) {
    stop("At least 2 pre-treatment periods are needed; the first treated period is ",
         period_labels[treatment_rowid], ".")
  }
  pre  <- seq_len(n_periods) < treatment_rowid
  post <- !pre
  
  ## ---- Units and predictor layout ----------------------------------------
  units <- unique(as.character(data[[unit_col]]))
  if (anyDuplicated(data[c(".rowid", unit_col)]) > 0) {
    stop("More than one row per unit and period. Check unit_col / time_cols.")
  }
  donors <- setdiff(units, treated_unit)
  
  ctrl_units <- switch(control_units,
                       all     = units,
                       donors  = donors,
                       treated = treated_unit)
  control_names <- unlist(lapply(controls, function(v) paste(v, ctrl_units, sep = "_")))
  x_names <- c(donors, control_names)
  
  if (anyDuplicated(x_names)) {
    stop("Predictor names clash (a unit name matches a <control>_<unit> name).")
  }
  
  # Build the X matrix from a list of period x unit matrices, one per variable
  build_X <- function(var_mats, scaled = TRUE) {
    X <- var_mats[[outcome_col]][, donors, drop = FALSE]
    for (v in controls) {
      Xc <- var_mats[[v]][, ctrl_units, drop = FALSE]
      colnames(Xc) <- paste(v, ctrl_units, sep = "_")
      X <- cbind(X, Xc)
    }
    if (scaled) X <- sweep(X, 2, x_scale, `*`)
    X
  }
  
  ## ---- Coefficient limits and penalty factors ----------------------------
  limit_vector <- function(spec, donor_fallback, control_default, arg) {
    nm <- names(spec)
    if (is.null(nm)) nm <- rep("", length(spec))
    unnamed <- spec[nm == ""]
    if (length(unnamed) > 1) {
      stop(arg, ": give at most one unnamed value (the default for donor columns).")
    }
    donor_default <- if (length(unnamed) == 1) unname(unnamed) else donor_fallback
    
    out <- c(stats::setNames(rep(donor_default, length(donors)), donors),
             stats::setNames(rep(control_default, length(control_names)), control_names))
    named <- spec[nm != ""]
    bad <- setdiff(names(named), x_names)
    if (length(bad) > 0) {
      stop(arg, ": unknown predictor name(s): ", paste(bad, collapse = ", "),
           ". Use donor names or <control>_<unit> names.")
    }
    out[names(named)] <- named
    out[x_names]
  }
  
  lower_vec <- limit_vector(lower_constr, 0, control_lower, "lower_constr")
  upper_vec <- limit_vector(upper_constr, 1, control_upper, "upper_constr")
  
  if (any(lower_vec > 0)) {
    stop("glmnet requires every lower limit to be <= 0. Check: ",
         paste(names(lower_vec)[lower_vec > 0], collapse = ", "))
  }
  if (any(upper_vec < 0)) {
    stop("glmnet requires every upper limit to be >= 0. Check: ",
         paste(names(upper_vec)[upper_vec < 0], collapse = ", "))
  }
  
  bad_unpen <- setdiff(unpenalized, x_names)
  if (length(bad_unpen) > 0) {
    stop("unpenalized: unknown predictor name(s): ", paste(bad_unpen, collapse = ", "),
         ". Use donor names or <control>_<unit> names.")
  }
  penalty_vec <- stats::setNames(ifelse(x_names %in% unpenalized, 0, 1), x_names)
  if (all(penalty_vec == 0)) stop("At least one predictor must be penalized.")
  
  # glmnet ignores the vectors' names; order matches x_names.
  # lower_fit / upper_fit are the limits on the (possibly rescaled) columns.
  fit_glmnet <- function(x, y, nl) {
    glmnet::glmnet(x, y, family = "gaussian", alpha = 1, standardize = standardize,
                   intercept = intercept, nlambda = nl,
                   lower.limits = unname(lower_fit), upper.limits = unname(upper_fit),
                   penalty.factor = unname(penalty_vec))
  }
  
  ## ---- Observed wide matrices: periods x units, one per variable ---------
  unit_index <- match(as.character(data[[unit_col]]), units)
  obs_mats <- lapply(vars, function(v) {
    m <- matrix(NA_real_, n_periods, length(units), dimnames = list(NULL, units))
    m[cbind(data$.rowid, unit_index)] <- data[[v]]
    m
  })
  names(obs_mats) <- vars
  if (anyNA(obs_mats[[outcome_col]])) {
    stop("The panel is unbalanced: some units are missing periods.")
  }
  
  y_obs <- obs_mats[[outcome_col]][, treated_unit]
  
  ## ---- Optional rescaling of the control columns -------------------------
  if (standardize && scale_controls) {
    warning("scale_controls is ignored when standardize = TRUE ",
            "(glmnet already standardizes every column).")
    scale_controls <- FALSE
  }
  
  x_scale <- stats::setNames(rep(1, length(x_names)), x_names)
  if (scale_controls && length(controls) > 0) {
    X_raw <- build_X(obs_mats, scaled = FALSE)
    pre_sd <- apply(X_raw[pre, , drop = FALSE], 2, stats::sd)
    if (any(pre_sd[control_names] == 0)) {
      stop("scale_controls: these control columns are constant in the pre-period: ",
           paste(control_names[pre_sd[control_names] == 0], collapse = ", "))
    }
    donor_sd <- mean(pre_sd[donors])
    x_scale[control_names] <- donor_sd / pre_sd[control_names]
  }
  
  # Limits on the scale glmnet sees: coefficient on (x * k) is b / k
  lower_fit <- lower_vec / x_scale
  upper_fit <- upper_vec / x_scale
  
  X_obs <- build_X(obs_mats)
  
  ## ---- 1. LASSO synthetic control on the observed data -------------------
  set.seed(seed)
  cv.fit <- glmnet::cv.glmnet(x = X_obs[pre, , drop = FALSE], y = y_obs[pre],
                              family = "gaussian", alpha = 1, standardize = standardize,
                              intercept = intercept, nlambda = nlambda,
                              lower.limits = unname(lower_fit),
                              upper.limits = unname(upper_fit),
                              penalty.factor = unname(penalty_vec))
  
  lambda_used <- switch(lambda_choice,
                        min_path   = min(cv.fit$lambda),
                        lambda.min = cv.fit$lambda.min,
                        lambda.1se = cv.fit$lambda.1se)
  
  # Coefficients at the lambda actually used for the predictions
  w <- as.matrix(stats::coef(cv.fit, s = lambda_used))
  # Convert coefficients on rescaled control columns back to original units
  w[x_names, 1] <- w[x_names, 1] * x_scale
  weights <- data.frame(Unit = rownames(w), Weights = w[, 1], row.names = NULL)
  weights$Type <- ifelse(weights$Unit == "(Intercept)", "intercept",
                         ifelse(weights$Unit %in% donors, "donor", "control"))
  weights$Unpenalized <- weights$Unit %in% unpenalized
  weights <- weights[weights$Weights != 0 | weights$Unpenalized, ]
  weights <- weights[order(-weights$Weights), ]
  rownames(weights) <- NULL
  
  treated_predicted <- as.numeric(stats::predict(cv.fit, newx = X_obs, s = lambda_used))
  att <- y_obs - treated_predicted
  
  effects <- data.frame(rowid = seq_len(n_periods), periods,
                        period = period_labels,
                        post = post,
                        treated_observed = y_obs,
                        treated_predicted = treated_predicted,
                        att = att,
                        att_post = ifelse(post, att, NA_real_),
                        cum_att = NA_real_)
  effects$cum_att[post] <- cumsum(att[post])
  
  # Sensitivity bounds: M x largest absolute pre-period residual
  bound <- M * max(abs(att[pre]))
  effects$bounds_lo <- ifelse(post, y_obs - (treated_predicted + bound), NA_real_)
  effects$bounds_hi <- ifelse(post, y_obs - (treated_predicted - bound), NA_real_)
  
  fit_stats <- data.frame(
    Name  = c("CV MSE (min)", "Lambda used", "Pre-period MSE", "M bound"),
    Stats = c(min(cv.fit$cvm), lambda_used, mean(att[pre]^2), bound)
  )
  
  ## ---- Optional jackknife conformal intervals ----------------------------
  if (jack_conform) {
    train_fun <- function(x, y, out = NULL) fit_glmnet(x, y, nlambda)
    # conformal.pred.jack() mishandles a single prediction column (t() of a vector makes res 1 x n, so only
    # the first leave-one-out residual is used). Returning two identical columns avoids this; column 1 is used.
    pred_fun <- function(out, newx) { p <- stats::predict(out, newx, s = lambda_used); cbind(p, p) }
    
    limits_jack <- conformalInference::conformal.pred.jack(
      x = X_obs[pre, , drop = FALSE], y = y_obs[pre],
      x0 = X_obs[post, , drop = FALSE],
      train.fun = train_fun, predict.fun = pred_fun,
      alpha = 1 - conform_level, verbose = FALSE, plus = FALSE)
    
    effects$pred_jack_lo <- NA_real_
    effects$pred_jack_hi <- NA_real_
    effects$pred_jack_lo[post] <- as.numeric(as.matrix(limits_jack$lo)[, 1])
    effects$pred_jack_hi[post] <- as.numeric(as.matrix(limits_jack$up)[, 1])
  }
  
  ## ---- 2 + 3. Bootstrap and refit ----------------------------------------
  boot_out <- NULL
  
  if (boot) {
    boot_df <- data[c(unit_col, ".rowid", vars)]
    boot_df[[unit_col]] <- as.character(boot_df[[unit_col]])
    if (boot_by_unit) {
      boot_df <- boot_df[order(boot_df[[unit_col]], boot_df$.rowid), ]
    }
    
    replicates <- meboot_multi(boot_df,
                               boot_cols  = vars,
                               id_cols    = c(unit_col, ".rowid"),
                               reps       = reps,
                               trim       = trim,
                               group_col  = if (boot_by_unit) unit_col else NULL,
                               shared_draws = shared_draws,
                               output     = "wide",
                               names_from = unit_col,
                               seed       = boot_seed,
                               as_matrix  = TRUE,
                               ...)
    
    # Split each replicate back into period-ordered matrices, one per variable
    replicates <- lapply(replicates, function(m) {
      m <- m[order(m[, ".rowid"]), , drop = FALSE]
      mats <- lapply(vars, function(v) {
        cols <- if (length(vars) == 1) units else paste(v, units, sep = "_")
        out <- m[, cols, drop = FALSE]
        colnames(out) <- units
        out
      })
      names(mats) <- vars
      mats
    })
    
    # Refit the LASSO on each bootstrapped pre period, predict the
    # bootstrapped post period at the fixed lambda, and compare with the
    # OBSERVED treated outcome (as in bootstrap_cum_effects_fixed)
    eff_sim <- vapply(replicates, function(mats) {
      X_b <- build_X(mats)
      y_b <- mats[[outcome_col]][, treated_unit]
      fit <- fit_glmnet(X_b[pre, , drop = FALSE], y_b[pre], boot_nlambda)
      pred <- stats::predict(fit, newx = X_b[post, , drop = FALSE], s = lambda_used)
      y_obs[post] - as.numeric(pred)
    }, numeric(sum(post)))
    eff_sim <- matrix(eff_sim, nrow = sum(post))  # post periods x reps
    
    cum_eff_sim <- apply(eff_sim, 2, cumsum)
    cum_eff_sim <- matrix(cum_eff_sim, nrow = sum(post))
    
    a <- (1 - conf_int) / 2
    row_q <- function(mat, p) apply(mat, 1, stats::quantile, probs = p, names = FALSE)
    
    for (nm in c("eff_low", "eff_high", "average_eff",
                 "cum_eff_low", "cum_eff_high", "average_cum_effect")) {
      effects[[nm]] <- NA_real_
    }
    
    effects$eff_low[post]            <- row_q(eff_sim, a)
    effects$eff_high[post]           <- row_q(eff_sim, 1 - a)
    effects$average_eff[post]        <- rowMeans(eff_sim)
    effects$cum_eff_low[post]       <- row_q(cum_eff_sim, a)
    effects$cum_eff_high[post]      <- row_q(cum_eff_sim, 1 - a)
    effects$average_cum_effect[post] <- rowMeans(cum_eff_sim)
    
    boot_out <- list(effects = eff_sim, cum_effects = cum_eff_sim)
    if (keep_boot) boot_out$replicates <- replicates
  }
  
  effects$treated_unit <- treated_unit
  
  settings <- list(treated_unit = treated_unit, treatment_period = tp,
                   treatment_rowid = treatment_rowid,
                   treatment_label = period_labels[treatment_rowid],
                   controls = controls, control_units = control_units,
                   intercept = intercept,
                   standardize = standardize, scale_controls = scale_controls,
                   predictors = data.frame(predictor = x_names,
                                           lower = unname(lower_vec),
                                           upper = unname(upper_vec),
                                           penalty_factor = unname(penalty_vec),
                                           scale_factor = unname(x_scale)),
                   lambda_choice = lambda_choice, lambda_used = lambda_used,
                   n_pre = sum(pre), n_post = sum(post),
                   boot = boot, reps = if (boot) reps else NA,
                   boot_by_unit = boot_by_unit, shared_draws = shared_draws,
                   trim = trim, conf_int = conf_int,
                   seed = seed, boot_seed = boot_seed)
  
  list(effects = effects, fit_stats = fit_stats, weights = weights,
       model = cv.fit, boot = boot_out, settings = settings)
}


## ===========================================================================
## 3. synth_lasso_placebo()
## ===========================================================================

##############################################################################
# synth_lasso_placebo(): LASSO synthetic control with exact placebo
#                        (permutation) tests and confidence intervals
#
# Same model and arguments as synth_lasso_boot(), but the bootstrap intervals
# are replaced by exact placebo tests:
#   1. Fit the treated unit (identical to synth_lasso_boot(), boot = FALSE).
#   2. Drop the treated unit. Fit every donor in turn as a placebo "treated"
#      unit, from the remaining donors, with the same settings.
#   3. Optionally keep only placebos with a good pre-treatment fit (Abadie,
#      Diamond and Hainmueller 2010): pre-period MSPE <= abadie_cutoff x the
#      treated unit's pre-period MSPE (typically 5 or 2).
#   4. For every post-treatment period, test the treated unit's gap against the
#      placebo gaps, for the per-period effect and for the cumulative effect
#      through that period. Gaps are relative: (observed - synthetic) /
#      synthetic; cumulative gaps are sum(gap) / sum(synthetic).
#
# Exact tests with N placebos (N + 1 units), alpha = 1 - conf_level:
#   absolute  p = (1 + #{|placebo gap| >= |treated gap|}) / (N + 1)
#             CI = treated gap +/- the k-th largest |placebo gap|,
#             k = floor(alpha x (N + 1))
#   signed    p = min(1, 2 x min(lower-tail p, upper-tail p))
#             CI = [treated gap - k2-th largest placebo gap,
#                   treated gap - k2-th smallest placebo gap],
#             k2 = floor(alpha / 2 x (N + 1))
#   The CIs are the effect sizes the test does not reject (test inversion);
#   a CI excludes zero exactly when p <= alpha. If k (or k2) is 0 there are
#   too few placebos for conf_level and the CI is infinite (with a warning):
#   the absolute test needs N >= 1/alpha - 1 placebos, the signed test
#   N >= 2/alpha - 1 (19 and 39 at 95%).
#   CIs are computed on the relative scale and converted to outcome units by
#   multiplying by the treated unit's synthetic value (or its running total).
#
# Arguments: as synth_lasso_boot() for the data, controls, LASSO fit, M bounds
# and jackknife (see its documentation). Bootstrap arguments are replaced by:
#   conf_level     confidence level of the placebo CIs (default 0.95)
#   abadie_cutoff  NULL (default) keeps every placebo; a number (e.g. 5 or 2)
#                  keeps placebos whose pre-period MSPE is at most that many
#                  times the treated unit's
#   n_cores        number of cores for the placebo fits (default 1; > 1 runs
#                  them in parallel)
#
# Named settings (lower_constr / upper_constr / unpenalized entries naming a
# donor or <control>_<unit> column) are applied in a placebo run only when that
# column exists in it; e.g. unpenalized = "Colorado" is dropped in the run
# where Colorado is the placebo.
#
# Value: a list
#   effects    one row per period, as in synth_lasso_boot() (without bootstrap
#              columns), plus for post-treatment periods:
#              att_rel, cum_att_rel         relative gap and cumulative gap
#              p_abs, eff_low_abs, eff_high_abs              per-period, absolute
#              p_signed, eff_low_signed, eff_high_signed     per-period, signed
#              cum_p_abs, cum_eff_low_abs, cum_eff_high_abs  cumulative, absolute
#              cum_p_signed, cum_eff_low_signed, cum_eff_high_signed
#                                                            cumulative, signed
#              (intervals in outcome units)
#   fit_stats, weights, model   as in synth_lasso_boot()
#   placebo    list: gaps and cum_gaps (post periods x placebos, relative),
#              fit (each placebo's pre-period MSPE, its ratio to the treated
#              unit's, and whether it was kept), n_kept, n_failed
#   settings   the settings used
##############################################################################

synth_lasso_placebo <- function(data,
                                unit_col,
                                time_cols,
                                outcome_col,
                                treated_unit,
                                treatment_period,
                                controls = NULL,
                                control_units = c("all", "donors", "treated"),
                                intercept = TRUE,
                                lower_constr = 0,
                                upper_constr = 1,
                                control_lower = -Inf,
                                control_upper = Inf,
                                unpenalized = NULL,
                                standardize = FALSE,
                                scale_controls = FALSE,
                                nlambda = 1000,
                                lambda_choice = c("min_path", "lambda.min", "lambda.1se"),
                                seed,
                                M = 1,
                                jack_conform = FALSE,
                                conform_level = 0.95,
                                conf_level = 0.95,
                                abadie_cutoff = NULL,
                                n_cores = 1) {
  
  control_units <- match.arg(control_units)
  lambda_choice <- match.arg(lambda_choice)
  data <- as.data.frame(data)
  
  ## ---- Checks ------------------------------------------------------------
  if (missing(seed)) stop("Please supply a seed for reproducibility.")
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number between 0 and 1, e.g. 0.95.")
  }
  if (!is.null(abadie_cutoff) &&
      (!is.numeric(abadie_cutoff) || length(abadie_cutoff) != 1 || abadie_cutoff <= 0)) {
    stop("abadie_cutoff must be NULL or a single positive number, e.g. 5 or 2.")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
    stop("n_cores must be a single number >= 1.")
  }
  alpha <- round(1 - conf_level, 10)   # avoid floating-point error in the rank (e.g. 1 - 0.9)
  
  # Arguments shared by every fit (treated and placebo)
  common <- list(unit_col = unit_col, time_cols = time_cols, outcome_col = outcome_col,
                 treatment_period = treatment_period, controls = controls,
                 control_units = control_units, intercept = intercept,
                 control_lower = control_lower, control_upper = control_upper,
                 standardize = standardize, scale_controls = scale_controls,
                 nlambda = nlambda, lambda_choice = lambda_choice, seed = seed, boot = FALSE)
  
  ## ---- 1. Treated unit (same as synth_lasso_boot, no bootstrap) ----------
  main <- do.call(synth_lasso_boot, c(common, list(
    data = data, treated_unit = treated_unit, lower_constr = lower_constr,
    upper_constr = upper_constr, unpenalized = unpenalized, M = M,
    jack_conform = jack_conform, conform_level = conform_level)))
  eff <- main$effects
  post <- eff$post
  treated_pre_mspe <- mean(eff$att[!post]^2)
  
  ## ---- 2. Placebo fits: every donor, treated unit removed ----------------
  pl_data <- data[as.character(data[[unit_col]]) != treated_unit, ]
  pl_units <- unique(as.character(pl_data[[unit_col]]))
  if (length(pl_units) < 3) stop("At least 3 donor units are needed for placebo runs.")
  
  # Keep only named settings that exist as predictors in a given placebo run
  predictor_names <- function(p) {
    donors_p <- setdiff(pl_units, p)
    ctrl_p <- switch(control_units, all = pl_units, donors = donors_p, treated = p)
    c(donors_p, unlist(lapply(controls, function(v) paste(v, ctrl_p, sep = "_"))))
  }
  keep_named <- function(spec, ok) {
    nm <- names(spec)
    if (is.null(nm)) return(spec)
    spec[nm == "" | nm %in% ok]
  }
  placebo_fit <- function(p) {
    ok <- predictor_names(p)
    unpen_p <- intersect(unpenalized, ok)
    m <- tryCatch(
      do.call(synth_lasso_boot, c(common, list(
        data = pl_data, treated_unit = p,
        lower_constr = keep_named(lower_constr, ok), upper_constr = keep_named(upper_constr, ok),
        unpenalized = if (length(unpen_p) > 0) unpen_p else NULL))),
      error = function(e) conditionMessage(e))
    if (is.character(m)) return(list(error = m))
    e <- m$effects
    list(gap_rel = e$att[e$post] / e$treated_predicted[e$post],
         cum_gap_rel = cumsum(e$att[e$post]) / cumsum(e$treated_predicted[e$post]),
         pre_mspe = mean(e$att[!e$post]^2))
  }
  
  if (n_cores > 1) {
    cl <- parallel::makeCluster(min(n_cores, length(pl_units)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    fn_env <- environment(sys.function())
    parallel::clusterExport(cl, intersect(c("synth_lasso_boot", "meboot_multi"), ls(fn_env)), envir = fn_env)
    fits <- parallel::parLapplyLB(cl, pl_units, placebo_fit)
  } else {
    fits <- lapply(pl_units, placebo_fit)
  }
  names(fits) <- pl_units
  
  failed <- vapply(fits, function(f) !is.null(f$error), logical(1))
  if (any(failed)) {
    warning(sum(failed), " placebo fit(s) failed and were left out: ",
            paste(names(fits)[failed], collapse = ", "), ". First error: ", fits[[which(failed)[1]]]$error)
  }
  fits <- fits[!failed]
  if (length(fits) == 0) stop("All placebo fits failed.")
  
  G  <- sapply(fits, `[[`, "gap_rel");     G  <- matrix(G,  nrow = sum(post), dimnames = list(NULL, names(fits)))
  Gc <- sapply(fits, `[[`, "cum_gap_rel"); Gc <- matrix(Gc, nrow = sum(post), dimnames = list(NULL, names(fits)))
  pl_mspe <- vapply(fits, `[[`, numeric(1), "pre_mspe")
  
  ## ---- 3. Abadie filter on pre-treatment fit -------------------------------
  kept <- if (is.null(abadie_cutoff)) rep(TRUE, length(fits)) else pl_mspe <= abadie_cutoff * treated_pre_mspe
  n_kept <- sum(kept)
  fit_table <- data.frame(unit = names(fits), pre_mspe = pl_mspe,
                          mspe_ratio_to_treated = pl_mspe / treated_pre_mspe, kept = kept, row.names = NULL)
  
  if (floor(alpha * (n_kept + 1)) == 0) {
    warning(n_kept, " placebos kept: too few for a ", 100 * conf_level,
            "% absolute CI (needs ", ceiling(1 / alpha - 1), "); absolute CIs are infinite.")
  }
  if (floor(alpha / 2 * (n_kept + 1)) == 0) {
    warning(n_kept, " placebos kept: too few for a ", 100 * conf_level,
            "% signed CI (needs ", ceiling(2 / alpha - 1), "); signed CIs are infinite.")
  }
  
  ## ---- 4. Exact tests and inverted CIs, every post period ----------------
  exact_tests <- function(t, g) {
    n1 <- length(g) + 1
    p_abs <- (1 + sum(abs(g) >= abs(t))) / n1
    k <- floor(alpha * n1)
    a_k <- if (k == 0) Inf else sort(abs(g), decreasing = TRUE)[k]
    p_sgn <- min(1, 2 * min((1 + sum(g <= t)) / n1, (1 + sum(g >= t)) / n1))
    k2 <- floor(alpha / 2 * n1)
    gs <- sort(g)
    c(p_abs = p_abs, abs_lo = t - a_k, abs_hi = t + a_k, p_signed = p_sgn,
      signed_lo = if (k2 == 0) -Inf else t - gs[length(g) + 1 - k2],
      signed_hi = if (k2 == 0)  Inf else t - gs[k2])
  }
  
  pred_post <- eff$treated_predicted[post]
  att_rel <- eff$att[post] / pred_post
  cum_att_rel <- cumsum(eff$att[post]) / cumsum(pred_post)
  test_names <- c("p_abs", "abs_lo", "abs_hi", "p_signed", "signed_lo", "signed_hi")
  per <- t(vapply(seq_along(att_rel), function(h) exact_tests(att_rel[h], G[h, kept]),
                  setNames(numeric(6), test_names)))
  cum <- t(vapply(seq_along(cum_att_rel), function(h) exact_tests(cum_att_rel[h], Gc[h, kept]),
                  setNames(numeric(6), test_names)))
  
  add_col <- function(name, values) {
    eff[[name]] <<- NA_real_
    eff[[name]][post] <<- values
  }
  add_col("att_rel", att_rel)
  add_col("cum_att_rel", cum_att_rel)
  add_col("p_abs", per[, "p_abs"])
  add_col("eff_low_abs", per[, "abs_lo"] * pred_post)
  add_col("eff_high_abs", per[, "abs_hi"] * pred_post)
  add_col("p_signed", per[, "p_signed"])
  add_col("eff_low_signed", per[, "signed_lo"] * pred_post)
  add_col("eff_high_signed", per[, "signed_hi"] * pred_post)
  add_col("cum_p_abs", cum[, "p_abs"])
  add_col("cum_eff_low_abs", cum[, "abs_lo"] * cumsum(pred_post))
  add_col("cum_eff_high_abs", cum[, "abs_hi"] * cumsum(pred_post))
  add_col("cum_p_signed", cum[, "p_signed"])
  add_col("cum_eff_low_signed", cum[, "signed_lo"] * cumsum(pred_post))
  add_col("cum_eff_high_signed", cum[, "signed_hi"] * cumsum(pred_post))
  eff$treated_unit <- NULL
  eff$treated_unit <- treated_unit
  
  settings <- main$settings
  settings$boot <- NULL; settings$reps <- NULL; settings$boot_by_unit <- NULL
  settings$shared_draws <- NULL; settings$trim <- NULL; settings$conf_int <- NULL
  settings$boot_seed <- NULL
  settings <- c(settings, list(inference = "placebo", conf_level = conf_level,
                               abadie_cutoff = abadie_cutoff, treated_pre_mspe = treated_pre_mspe,
                               n_placebos = length(fits), n_placebos_kept = n_kept))
  
  list(effects = eff, fit_stats = main$fit_stats, weights = main$weights, model = main$model,
       placebo = list(gaps = G, cum_gaps = Gc, fit = fit_table, n_kept = n_kept,
                      n_failed = sum(failed)),
       settings = settings)
}
