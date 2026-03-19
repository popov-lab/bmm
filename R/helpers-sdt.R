############################################################################# !
# SDT SHARED PIPELINE HELPERS                                            ####
# Used by model_sdt_*.R files for check_data, bmf2bf, configure_model    #
############################################################################# !


############################################################################# !
# SHARED S3 METHODS                                                      ####
############################################################################# !

#' @export
check_data.sdt <- function(model, data, formula) {
  NextMethod("check_data")
}


############################################################################# !
# DISTRIBUTION ID & NAME MAPPING                                          ####
############################################################################# !

# Map distribution name to Stan integer ID
.sdt_dist_id <- function(dist) {
  .SDT_DISTS[[dist]]$id
}

# Map dist_type integer (from Stan vint) back to distribution name
.sdt_dist_names <- names(.SDT_DISTS)[order(vapply(.SDT_DISTS, `[[`, 0L, "id"))]

# Stan CDF expression builder — returns function(x) producing Stan expression
.sdt_cdf_expr <- function(dist) {
  .SDT_DISTS[[dist]]$stan_expr
}

# Stan log-CDF expression builder
.sdt_log_cdf_expr <- function(dist) {
  .SDT_DISTS[[dist]]$log_stan_expr
}

# Stan log complementary CDF expression builder: log(1 - CDF(x))
.sdt_log1m_cdf_expr <- function(dist) {
  .SDT_DISTS[[dist]]$log1m_stan_expr
}


############################################################################# !
# DATA VALIDATION HELPERS                                                 ####
############################################################################# !

# Shared validation for count + trial columns
# Checks existence, non-negativity, integerishness, and count <= trials
.validate_sdt_counts <- function(data, resp_var, n_trials_var) {
  required <- c(resp_var, n_trials_var)
  missing <- setdiff(required, colnames(data))
  stopif(length(missing) > 0,
    "Variables {collapse_comma(missing)} missing in the data")

  resp_vals <- data[[resp_var]]
  stopif(any(resp_vals < 0, na.rm = TRUE),
    "Response variable '{resp_var}' must contain non-negative counts")
  warnif(any(resp_vals != round(resp_vals), na.rm = TRUE),
    "Response variable '{resp_var}' should contain integer counts")

  trial_vals <- data[[n_trials_var]]
  stopif(any(trial_vals <= 0, na.rm = TRUE),
    "Variable '{n_trials_var}' must contain positive values")

  stopif(any(resp_vals > trial_vals, na.rm = TRUE),
    "Response counts in '{resp_var}' must not exceed '{n_trials_var}'")
}

# Rating-specific data checks (shared by sdt_rating, sdt_dp, sdt_metad)
.check_data_sdt_rating <- function(model, data) {
  resp_cols <- model$resp_vars$response

  missing <- setdiff(resp_cols, colnames(data))
  stopif(length(missing) > 0,
    "Response columns {collapse_comma(missing)} missing in the data")

  for (col in resp_cols) {
    vals <- data[[col]]
    stopif(any(vals < 0, na.rm = TRUE),
      "Response column '{col}' must contain non-negative counts")
    warnif(any(vals != round(vals), na.rm = TRUE),
      "Response column '{col}' should contain integer counts")
  }

  Y <- as.matrix(data[resp_cols])
  data$nTrials <- rowSums(Y)
  data$Y <- Y

  stopif(any(data$nTrials <= 0, na.rm = TRUE),
    "Row sums of response columns must be positive (no empty rows)")

  data
}

# Check if model is a rating model (used by old sdt() code path)
.is_sdt_rating <- function(model) {
  n_ratings <- model$other_vars$n_ratings
  !is.null(n_ratings) && n_ratings > 2
}


############################################################################# !
# THRESHOLD EXPRESSIONS                                                   ####
############################################################################# !

# Equidistant: theta[k] = criterion + (k - mid) * exp(spacing)
.sdt_thresholds_equidistant <- function(n_ratings, mid) {
  function(k) {
    offset <- k - mid
    if (offset == 0L) return("criterion")
    if (offset > 0L) {
      return(paste0("criterion + ", offset, " * exp(spacing)"))
    }
    paste0("criterion - ", abs(offset), " * exp(spacing)")
  }
}

# Log-distance: theta[mid] = criterion, theta[k] = theta[k-1] + exp(delta_k)
.sdt_thresholds_log_distance <- function(n_ratings, mid) {
  K1 <- n_ratings - 1L
  theta_exprs <- character(K1)
  theta_exprs[mid] <- "criterion"

  if (mid < K1) {
    for (k in (mid + 1L):K1) {
      delta_name <- paste0("delta", k)
      theta_exprs[k] <- paste0(theta_exprs[k - 1L], " + exp(", delta_name, ")")
    }
  }

  if (mid > 1L) {
    for (k in (mid - 1L):1L) {
      delta_name <- paste0("delta", k)
      theta_exprs[k] <- paste0(theta_exprs[k + 1L], " - exp(", delta_name, ")")
    }
  }

  function(k) theta_exprs[k]
}

# Parsimonious (Selker et al., 2019):
# theta[k] = criterion + exp(spacing) * logit(k/K)
.sdt_thresholds_parsimonious <- function(n_ratings) {
  K <- n_ratings
  K1 <- K - 1L
  canonical <- log(seq_len(K1) / (K - seq_len(K1)))
  function(k) {
    gk <- canonical[k]
    if (gk == 0) return("criterion")
    if (gk > 0) {
      sprintf("criterion + exp(spacing) * %.10f", gk)
    } else {
      sprintf("criterion - exp(spacing) * %.10f", abs(gk))
    }
  }
}

# Log-ratio (Paulewicz & Blaut, 2020):
# theta[mid] = criterion, distances as ratios relative to reference spread
.sdt_thresholds_log_ratio <- function(n_ratings, mid) {
  K1 <- n_ratings - 1L
  theta_exprs <- character(K1)
  theta_exprs[mid] <- "criterion"

  spread_name <- paste0("delta", mid + 1L)
  spread_above <- paste0("exp(", spread_name, ")")
  theta_exprs[mid + 1L] <- paste0("criterion + ", spread_above)

  if (mid > 1L) {
    ratio_name <- paste0("delta", mid - 1L)
    spread_below <- paste0("exp(", ratio_name, ") * ", spread_above)
    theta_exprs[mid - 1L] <- paste0("criterion - ", spread_below)
  }

  if (mid + 2L <= K1) {
    for (k in (mid + 2L):K1) {
      delta_name <- paste0("delta", k)
      theta_exprs[k] <- paste0(
        theta_exprs[k - 1L], " + exp(", delta_name, ") * ", spread_above
      )
    }
  }

  if (mid > 1L && mid - 2L >= 1L) {
    for (k in (mid - 2L):1L) {
      delta_name <- paste0("delta", k)
      theta_exprs[k] <- paste0(
        theta_exprs[k + 1L], " - exp(", delta_name, ") * ", spread_below
      )
    }
  }

  function(k) theta_exprs[k]
}


############################################################################# !
# RATING FORMULA CONSTRUCTION                                             ####
############################################################################# !

# Shared infrastructure for rating SDT formula construction.
# Extracts model properties and builds the base formula, CDF expression builder,
# threshold expressions, and cumulative probability function.
# shift_expr: the Stan expression for the decision variable shift
# Returns a list with: base_formula, cum_prob, cdf_fn, theta_expr,
#   n_ratings, mid, resp_cols, raw_arg, dist
.sdt_rating_formula_parts <- function(model,
                                      shift_expr = "dprime / 2 * (2 * stimulus - 1)") {
  resp_cols <- model$resp_vars$response
  n_ratings <- model$other_vars$n_ratings
  threshold_type <- model$other_vars$threshold_type
  dist <- model$other_vars$dist
  mid <- n_ratings %/% 2L

  resp_str <- paste(resp_cols, collapse = ", ")
  base_formula <- brms::bf(
    paste0("cbind(", resp_str, ") | trials(nTrials) ~ 1"),
    nl = TRUE
  )

  cdf_fn <- .sdt_cdf_expr(dist)

  theta_expr <- if (threshold_type == "equidistant") {
    .sdt_thresholds_equidistant(n_ratings, mid)
  } else if (threshold_type == "parsimonious") {
    .sdt_thresholds_parsimonious(n_ratings)
  } else if (threshold_type == "log_ratio") {
    .sdt_thresholds_log_ratio(n_ratings, mid)
  } else {
    .sdt_thresholds_log_distance(n_ratings, mid)
  }

  # Always include UV-SDT scale factor: when sdratio=0, scale=1 (EV-SDT)
  scale <- "(stimulus * exp(sdratio) + (1 - stimulus))"
  raw_arg <- function(k) paste0("(", theta_expr(k), " - ", shift_expr, ") / ", scale)

  cum_prob <- function(k) cdf_fn(raw_arg(k))

  list(base_formula = base_formula, cum_prob = cum_prob, cdf_fn = cdf_fn,
       theta_expr = theta_expr, n_ratings = n_ratings, mid = mid,
       resp_cols = resp_cols, raw_arg = raw_arg, dist = dist)
}


############################################################################# !
# NLF FORMULA BUILDERS                                                    ####
############################################################################# !

# Standard SDT nlf formulas (log of raw category probabilities)
.sdt_build_nlf_standard <- function(parts) {
  nlf_formulas <- list()
  for (k in seq_len(parts$n_ratings)) {
    dpar <- paste0("mu", parts$resp_cols[k])
    if (k == 1L) {
      expr <- paste0("log(", parts$cum_prob(1), ")")
    } else if (k == parts$n_ratings) {
      expr <- paste0("log(1 - ", parts$cum_prob(parts$n_ratings - 1L), ")")
    } else {
      expr <- paste0("log(", parts$cum_prob(k), " - ", parts$cum_prob(k - 1L), ")")
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}

# Log-scale SDT nlf formulas (numerically stable)
# Uses log_CDF, log(1-CDF), and log_diff_exp for interior categories
.sdt_build_nlf_logscale <- function(parts) {
  log_cdf_fn <- .sdt_log_cdf_expr(parts$dist)
  log1m_cdf_fn <- .sdt_log1m_cdf_expr(parts$dist)

  nlf_formulas <- list()
  for (k in seq_len(parts$n_ratings)) {
    dpar <- paste0("mu", parts$resp_cols[k])
    if (k == 1L) {
      expr <- log_cdf_fn(parts$raw_arg(1))
    } else if (k == parts$n_ratings) {
      expr <- log1m_cdf_fn(parts$raw_arg(parts$n_ratings - 1L))
    } else {
      log_upper <- log_cdf_fn(parts$raw_arg(k))
      log_lower <- log_cdf_fn(parts$raw_arg(k - 1L))
      expr <- paste0("log_diff_exp(", log_upper, ", ", log_lower, ")")
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}

# DPSDT nlf formulas (familiarity-scaled probs + recollection)
# For old items (stimulus=1): p[k] = (1-Ro)*sdt_prob[k], p[K] += Ro
# For new items (stimulus=0): p[k] = (1-Rn)*sdt_prob[k], p[1] += Rn
.sdt_build_nlf_dpsdt <- function(parts) {
  fam_scale <- "(1 - inv_logit(Ro) * stimulus - inv_logit(Rn) * (1 - stimulus))"
  rec_add_K <- "inv_logit(Ro) * stimulus"
  rec_add_1 <- "inv_logit(Rn) * (1 - stimulus)"

  nlf_formulas <- list()
  for (k in seq_len(parts$n_ratings)) {
    dpar <- paste0("mu", parts$resp_cols[k])
    if (k == 1L) {
      expr <- paste0(
        "log(", fam_scale, " * ", parts$cum_prob(1),
        " + ", rec_add_1, ")"
      )
    } else if (k == parts$n_ratings) {
      expr <- paste0(
        "log(", fam_scale, " * (1 - ", parts$cum_prob(parts$n_ratings - 1L), ")",
        " + ", rec_add_K, ")"
      )
    } else {
      expr <- paste0(
        "log(", fam_scale, " * (",
        parts$cum_prob(k), " - ", parts$cum_prob(k - 1L), "))"
      )
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}

# Meta-d' nlf formulas (normalized meta-d' category probabilities)
# Uses meta-d' for confidence threshold placement, but normalizes so that the
# total "old"/"new" response rates match what type-1 d' predicts.
.sdt_build_nlf_metad <- function(parts_metad, parts_dprime) {
  crit_expr <- parts_metad$theta_expr(parts_metad$mid)
  d_shift <- "dprime / 2 * (2 * stimulus - 1)"
  metad_shift <- "metad / 2 * (2 * stimulus - 1)"
  cdf_fn <- parts_metad$cdf_fn

  norm_below <- paste0(
    cdf_fn(paste0(crit_expr, " - ", d_shift)),
    " / ",
    cdf_fn(paste0(crit_expr, " - ", metad_shift))
  )
  norm_above <- paste0(
    "(1 - ", cdf_fn(paste0(crit_expr, " - ", d_shift)), ")",
    " / ",
    "(1 - ", cdf_fn(paste0(crit_expr, " - ", metad_shift)), ")"
  )

  nlf_formulas <- list()
  for (k in seq_len(parts_metad$n_ratings)) {
    dpar <- paste0("mu", parts_metad$resp_cols[k])
    norm <- if (k <= parts_metad$mid) norm_below else norm_above

    if (k == 1L) {
      expr <- paste0("log(", parts_metad$cum_prob(1), " * ", norm, ")")
    } else if (k == parts_metad$n_ratings) {
      expr <- paste0("log((1 - ", parts_metad$cum_prob(parts_metad$n_ratings - 1L),
                      ") * ", norm, ")")
    } else {
      expr <- paste0("log((", parts_metad$cum_prob(k), " - ",
                      parts_metad$cum_prob(k - 1L), ") * ", norm, ")")
    }
    nlf_formulas[[k]] <- glue_nlf("{dpar} ~ {expr}")
  }
  nlf_formulas
}


############################################################################# !
# CONFIGURE_MODEL HELPER                                                  ####
############################################################################# !

# Shared rating model configuration (multinomial + nlf)
# Used by sdt_rating, sdt_dp, sdt_metad
.configure_sdt_rating <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  resp_cols <- model$resp_vars$response
  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- resp_cols
  formula$family$dpars <- paste0("mu", resp_cols)

  nlist(formula, data)
}


############################################################################# !
# PUBLIC HELPER FUNCTIONS                                                 ####
############################################################################# !

#' Combine stimulus type and confidence into SDT response categories
#'
#' Creates a combined response variable from separate stimulus and confidence
#' columns, suitable for use with SDT rating models. The combined response
#' maps "noise" trials to categories 1..K/2 (from highest to lowest confidence)
#' and "signal" trials to categories K/2+1..K (from lowest to highest).
#'
#' @param stimulus Integer vector (0/1). Stimulus type: 0 = noise, 1 = signal.
#' @param confidence Integer vector. Confidence rating (1 = lowest confidence,
#'   n_levels = highest confidence).
#' @param n_levels Integer. Number of confidence levels per stimulus type.
#'
#' @return Integer vector of combined response categories (1 to 2*n_levels).
#' @export
#' @examples
#' # 3 confidence levels per stimulus type -> K=6 combined categories
#' stim <- c(0, 0, 0, 1, 1, 1)
#' conf <- c(3, 2, 1, 1, 2, 3)
#' combine_sdt_response(stim, conf, n_levels = 3)
#' # Returns: 1, 2, 3, 4, 5, 6
combine_sdt_response <- function(stimulus, confidence, n_levels) {
  stopif(any(!stimulus %in% c(0, 1)),
         "stimulus must be 0 (noise) or 1 (signal)")
  stopif(any(confidence < 1 | confidence > n_levels),
         "confidence must be between 1 and {n_levels}")

  ifelse(stimulus == 0,
         n_levels - confidence + 1L,
         n_levels + confidence)
}
