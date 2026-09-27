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
#               Missing elements (trim / xmin / xmax) are filled with the
#               meboot defaults (trim = 0.10, xmin = NULL, xmax = NULL);
#               columns not listed get the defaults.
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
  trim_default <- list(trim = 0.10, xmin = NULL, xmax = NULL)
  trim_keys <- c("trim", "xmin", "xmax")
  
  # One list for all columns if every name is trim / xmin / xmax
  is_global_trim <- !is.null(names(trim)) && all(names(trim) %in% trim_keys)
  
  if (is_global_trim) {
    trim_list <- rep(list(utils::modifyList(trim_default, trim, keep.null = TRUE)),
                     length(boot_cols))
    names(trim_list) <- boot_cols
  } else {
    bad_names <- setdiff(names(trim), boot_cols)
    if (is.null(names(trim)) || length(bad_names) > 0) {
      stop("trim must be a single list (trim/xmin/xmax) or a list named by ",
           "boot_cols. Unrecognised name(s): ", paste(bad_names, collapse = ", "))
    }
    trim_list <- lapply(boot_cols, function(col) {
      if (col %in% names(trim)) {
        utils::modifyList(trim_default, trim[[col]], keep.null = TRUE)
      } else {
        trim_default
      }
    })
    names(trim_list) <- boot_cols
  }
  
  ## ---- Run meboot for each column (and group) ----------------------------
  if (!is.null(seed)) set.seed(seed)
  
  run_meboot <- function(x, col_trim) {
    ens <- meboot::meboot(x, reps = reps, trim = col_trim, ...)$ensemble
    as.matrix(ens)  # rows = observations, columns = replicates
  }
  
  ensembles <- lapply(boot_cols, function(col) {
    x <- data[[col]]
    
    if (is.null(group_col)) {
      return(run_meboot(x, trim_list[[col]]))
    }
    
    # Bootstrap each group separately and put results back in original rows
    ens <- matrix(NA_real_, nrow = length(x), ncol = reps)
    row_index <- split(seq_along(x), data[[group_col]])
    for (g in names(row_index)) {
      rows <- row_index[[g]]
      ens[rows, ] <- run_meboot(x[rows], trim_list[[col]])
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
#   treatment_period  rowid of the FIRST post-treatment period
#                     (Kansas: 90 = 2012 Q2). Pre = rowid < treatment_period.
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
#   conform_level   coverage of the conformal intervals (default 0.95)
#
# ---- Bootstrap arguments ------------------------------------------------------
#   boot            TRUE (default) runs the bootstrap
#   reps            number of bootstrap replicates (default 1000)
#   trim            meboot trim / xmin / xmax, passed to meboot_multi().
#                   One list for every bootstrapped column, or a list named
#                   by column, e.g.
#                     list(gdpcapita   = list(xmin = 15000),
#                          popestimate = list(xmin = 450000))
#   boot_by_unit    FALSE (default) bootstraps each variable as one stacked
#                   series in the row order of `data` (reproduces the Kansas
#                   script). TRUE bootstraps each unit's series separately.
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
  
  if (treatment_period < 3 || treatment_period > n_periods) {
    stop("treatment_period must be between 3 and ", n_periods, ".")
  }
  pre  <- seq_len(n_periods) < treatment_period
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
                        period = do.call(paste, c(periods, sep = "_")),
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
    pred_fun <- function(out, newx) stats::predict(out, newx, s = lambda_used)
    
    limits_jack <- conformalInference::conformal.pred.jack(
      x = X_obs[pre, , drop = FALSE], y = y_obs[pre],
      x0 = X_obs[post, , drop = FALSE],
      train.fun = train_fun, predict.fun = pred_fun,
      alpha = 1 - conform_level, verbose = FALSE, plus = FALSE)
    
    effects$pred_jack_lo <- NA_real_
    effects$pred_jack_hi <- NA_real_
    effects$pred_jack_lo[post] <- as.numeric(limits_jack$lo)
    effects$pred_jack_hi[post] <- as.numeric(limits_jack$up)
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
  
  settings <- list(treated_unit = treated_unit, treatment_period = treatment_period,
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
                   boot_by_unit = boot_by_unit, trim = trim, conf_int = conf_int,
                   seed = seed, boot_seed = boot_seed)
  
  list(effects = effects, fit_stats = fit_stats, weights = weights,
       model = cv.fit, boot = boot_out, settings = settings)
}
