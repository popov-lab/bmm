############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdt_cdp <- function(new_response = NULL, old_know = NULL,
                           old_remember = NULL, old_guess = NULL,
                           stimulus = NULL, dist = "normal",
                           n_ratings = NULL,
                           threshold_type = "parsimonious",
                           links = NULL, call = NULL, ...) {
  dist_int <- .sdt_dist_id(dist)
  n_new <- length(new_response)
  n_old <- length(old_know)
  has_guess <- !is.null(old_guess)

  if (is.null(n_ratings)) {
    n_ratings <- n_new + n_old
  }

  thresh_type_int <- if (threshold_type == "parsimonious") 1L else 2L

  parameters <- list(
    dprimef = glue(
      "Familiarity sensitivity: target mean on familiarity axis"
    ),
    dprimer = glue(
      "Recollection sensitivity: target mean on recollection axis"
    ),
    criterion = glue(
      "Response bias: central confidence criterion on F+R axis"
    ),
    spacing = glue(
      "Threshold spacing: controls distance between adjacent ",
      "confidence thresholds (exp(spacing) ensures positive spacing)"
    ),
    rcrit = glue(
      "Remember criterion: threshold on recollection axis"
    )
  )
  default_priors <- list(
    dprimef = list(main = "normal(1, 1)",
                   effects = "normal(0, 0.5)"),
    dprimer = list(main = "normal(1, 1)",
                   effects = "normal(0, 0.5)"),
    criterion = list(main = "normal(0, 1.5)",
                     effects = "normal(0, 0.5)"),
    spacing = list(main = "normal(0, 0.5)",
                   effects = "normal(0, 0.3)"),
    rcrit = list(main = "normal(0, 1)",
                 effects = "normal(0, 0.5)")
  )
  param_links <- list(
    dprimef = "identity", dprimer = "identity",
    criterion = "identity", spacing = "identity",
    rcrit = "identity"
  )

  parameters$sigmar <- glue(
    "SD ratio: log ratio of recollection target to noise SD ",
    "(exp(sigmar) ensures positivity, 0 = equal variance)"
  )
  default_priors$sigmar <- list(
    main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
  )
  param_links$sigmar <- "identity"

  fixed_parameters <- list(mu = 0, sigmar = 0)

  if (has_guess) {
    parameters$kcrit <- glue(
      "Know criterion: threshold on familiarity axis ",
      "(separates know from guess judgments)"
    )
    default_priors$kcrit <- list(
      main = "normal(0, 1)", effects = "normal(0, 0.5)"
    )
    param_links$kcrit <- "identity"
    fixed_parameters$kcrit <- -100
  }

  requirements <- glue(
    "Provide pre-aggregated data with the following columns:",
    "\n\n",
    "  - 'New' response counts: {n_new} columns", "\n",
    "  - 'Old-Know' counts: {n_old} columns", "\n",
    "  - 'Old-Remember' counts: {n_old} columns", "\n",
    if (has_guess) {
      glue("  - 'Old-Guess' counts: {n_old} columns", "\n")
    } else {
      ""
    },
    "  - Stimulus type (stimulus): 0 = noise, 1 = signal"
  )

  out <- structure(
    list(
      resp_vars = nlist(new_response, old_know,
                        old_remember, old_guess),
      other_vars = nlist(stimulus, dist, dist_int,
                         n_ratings, n_new, n_old,
                         threshold_type, thresh_type_int,
                         has_guess),
      domain = "Recognition Memory",
      task = "Old/New Recognition with Remember/Know Judgments",
      name = "Continuous Dual-Process Signal Detection Theory (CDP)",
      citation = glue(
        "Wixted, J. T., & Mickes, L. (2010). ",
        "A continuous dual-process model of remember/know ",
        "judgments. Psychological Review, 117(4), 1025-1054."
      ),
      version = "cdp",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = fixed_parameters,
      default_priors = default_priors,
      init_ranges = c(
        list(
          mu = c(0, 0),
          dprimef = c(0.3, 1.0),
          dprimer = c(0.3, 1.0),
          criterion = c(-0.3, 0.3),
          spacing = c(0.2, 0.8),
          rcrit = c(-0.5, 0.5),
          sigmar = c(-0.1, 0.1)
        ),
        if (has_guess) list(kcrit = c(-1.0, 0.0))
      ),
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdt", "sdt_cdp"),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title Continuous Dual-Process Signal Detection Theory Model
#' @name sdt_cdp
#' @details `r model_info(.model_sdt_cdp())`
#'
#' The CDP model (Wixted & Mickes, 2010) assumes two independent continuous
#' signal-detection dimensions: Familiarity (F) and Recollection (R). Old/new
#' confidence ratings are based on the aggregated signal F + R. Remember/Know
#' judgments are based on R exceeding a criterion on the recollection axis.
#'
#' The number of "new" and "old" confidence levels need not be equal. For
#' example, a 6-point scale with 1 "new" level and 5 "old" levels (each with
#' R/K split) is supported — as used in Rotello et al. (2005).
#'
#' Optionally, a Know/Guess split uses F exceeding a second criterion on the
#' familiarity axis. To enable the 3-way R/K/G model, provide `old_guess`
#' columns and add `kcrit ~ 1` to the formula.
#'
#' By default, the model assumes equal variance for recollection (sigmar
#' fixed to 0). To estimate unequal variance, add `sigmar ~ 1` to the formula.
#'
#' When only confidence rating data is available (no R/K split), the CDP model
#' reduces to unequal-variance SDT on the aggregated F + R axis. Use
#' [sdt_rating()] for that case.
#'
#' @param new_response A character vector of column names containing response
#'   counts for "new" confidence levels, ordered from "definitely new" to
#'   "guess new". Can have a different length than `old_know`/`old_remember`.
#' @param old_know A character vector of column names containing "old-know"
#'   response counts, ordered from "guess old" to "definitely old".
#' @param old_remember A character vector of column names containing
#'   "old-remember" response counts, same ordering as `old_know`.
#' @param old_guess Optional character vector of column names for "old-guess"
#'   response counts. If provided, enables the 3-way R/K/G model with the
#'   `kcrit` parameter. Must have the same length as `old_know`.
#' @param stimulus The name of the variable coding the stimulus type.
#'   Must be coded as 0 (noise/new) and 1 (signal/old).
#' @param n_ratings Integer. Total number of confidence levels. Defaults to
#'   `length(new_response) + length(old_know)`. Must be >= 3.
#' @param dist Character. The noise distribution to use:
#'   \itemize{
#'     \item "normal" (default): Gaussian SDT
#'     \item "logistic": Logistic SDT
#'     \item "gumbel_min": Extreme-value (Gumbel minimum) SDT
#'     \item "gumbel_max": Gumbel maximum SDT
#'   }
#' @param threshold_type Character. Threshold parameterization:
#'   "parsimonious" (default) or "equidistant".
#' @param links A named list of link functions for the parameters.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Wixted, J. T., & Mickes, L. (2010). A continuous dual-process model of
#'   remember/know judgments. \emph{Psychological Review}, \emph{117}(4),
#'   1025--1054. \doi{10.1037/a0020874}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # Symmetric 6-point scale (3 new + 3 old)
#' model <- sdt_cdp(
#'   new_response = c("n1", "n2", "n3"),
#'   old_know = c("k4", "k5", "k6"),
#'   old_remember = c("r4", "r5", "r6"),
#'   stimulus = "stimulus"
#' )
#'
#' # Asymmetric scale (1 new + 5 old, as in Rotello et al. 2005)
#' model <- sdt_cdp(
#'   new_response = c("new"),
#'   old_know = c("k2", "k3", "k4", "k5", "k6"),
#'   old_remember = c("r2", "r3", "r4", "r5", "r6"),
#'   stimulus = "stimulus"
#' )
#' }
sdt_cdp <- function(new_response, old_know, old_remember,
                    old_guess = NULL, stimulus,
                    n_ratings = NULL,
                    dist = c("normal", "logistic",
                             "gumbel_min", "gumbel_max"),
                    threshold_type = c("parsimonious",
                                       "equidistant"),
                    links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  dist <- match.arg(dist)
  threshold_type <- match.arg(threshold_type)

  n_new <- length(new_response)
  n_old <- length(old_know)
  stopif(n_new < 1,
         "CDP requires at least 1 new confidence level")
  stopif(n_old < 1,
         "CDP requires at least 1 old confidence level")
  stopif(n_new + n_old < 3,
         "CDP requires at least 3 total confidence levels")
  stopif(length(old_remember) != n_old,
         "old_remember must have the same length as old_know ({n_old})")
  if (!is.null(old_guess)) {
    stopif(length(old_guess) != n_old,
           "old_guess must have the same length as old_know ({n_old})")
  }

  if (is.null(n_ratings)) {
    n_ratings <- n_new + n_old
  }
  stopif(
    n_ratings != n_new + n_old,
    "n_ratings ({n_ratings}) must equal length(new_response) + length(old_know) ({n_new + n_old})"
  )

  .model_sdt_cdp(
    new_response = new_response, old_know = old_know,
    old_remember = old_remember, old_guess = old_guess,
    stimulus = stimulus, dist = dist,
    n_ratings = n_ratings,
    threshold_type = threshold_type,
    links = links, call = call, ...
  )
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdt_cdp <- function(model, data, formula) {
  stim_var <- model$other_vars$stimulus
  stopif(!stim_var %in% colnames(data),
         "Stimulus variable '{stim_var}' missing in the data")
  stim_vals <- unique(data[[stim_var]])
  stopif(
    !all(stim_vals %in% c(0, 1)),
    "Stimulus variable '{stim_var}' must be coded as 0 and 1"
  )

  new_cols <- model$resp_vars$new_response
  know_cols <- model$resp_vars$old_know
  rem_cols <- model$resp_vars$old_remember
  guess_cols <- model$resp_vars$old_guess
  has_guess <- model$other_vars$has_guess
  n_new <- model$other_vars$n_new
  n_old <- model$other_vars$n_old

  all_resp_cols <- if (has_guess) {
    c(new_cols, guess_cols, know_cols, rem_cols)
  } else {
    c(new_cols, know_cols, rem_cols)
  }

  for (col in all_resp_cols) {
    stopif(!col %in% colnames(data),
           "Response column '{col}' missing in the data")
    vals <- data[[col]]
    stopif(any(vals < 0, na.rm = TRUE),
           "Response column '{col}' must be non-negative")
    warnif(any(vals != round(vals), na.rm = TRUE),
           "Response column '{col}' should be integer counts")
  }

  n_total <- rowSums(data[all_resp_cols])
  other_cols <- setdiff(names(data), all_resp_cols)

  n_cats <- length(all_resp_cols)
  orig_nrow <- nrow(data)

  long_data <- data[rep(seq_len(orig_nrow), each = n_cats),
                    other_cols, drop = FALSE]
  rownames(long_data) <- NULL

  # Build category metadata with asymmetric block sizes
  # Ordering: new (type=1, n_new items),
  #           [guess (type=2, n_old items)],
  #           know (type=3, n_old items),
  #           remember (type=4, n_old items)
  if (has_guess) {
    cat_types <- c(rep(1L, n_new), rep(2L, n_old),
                   rep(3L, n_old), rep(4L, n_old))
    cat_confs <- c(seq_len(n_new), rep(seq_len(n_old), 3L))
  } else {
    cat_types <- c(rep(1L, n_new), rep(3L, n_old),
                   rep(4L, n_old))
    cat_confs <- c(seq_len(n_new), rep(seq_len(n_old), 2L))
  }

  long_data$cat_type <- rep(cat_types, orig_nrow)
  long_data$cat_conf <- rep(cat_confs, orig_nrow)
  long_data$stim_val <- as.integer(long_data[[stim_var]])
  long_data$n_new <- n_new
  long_data$n_old <- n_old
  long_data$has_guess <- as.integer(has_guess)
  long_data$dist_type <- model$other_vars$dist_int
  long_data$thresh_type <- model$other_vars$thresh_type_int
  long_data$n_total <- rep(n_total, each = n_cats)

  counts <- as.matrix(data[all_resp_cols])
  long_data$count <- as.vector(t(counts))

  data <- long_data
  NextMethod("check_data")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.sdt_cdp <- function(model, formula) {
  brms::bf(paste0(
    "count | vint(cat_type, cat_conf, stim_val, n_new, ",
    "n_old, has_guess, dist_type, thresh_type) ~ 1"
  ))
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.sdt_cdp <- function(model, data, formula) {
  sigmar_estimated <- !is.null(formula[["sigmar"]]) &&
    !is_constant(formula[["sigmar"]])
  kcrit_estimated <- model$other_vars$has_guess &&
    !is.null(formula[["kcrit"]]) &&
    !is_constant(formula[["kcrit"]])

  formula <- bmf2bf(model, formula)

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_binary <- read_lines2(
    paste0(sc_path, "/sdt_binary_funs.stan")
  )
  stan_cdp <- read_lines2(
    paste0(sc_path, "/sdt_cdp_funs.stan")
  )
  stanvars <- brms::stanvar(
    scode = paste(stan_binary, stan_cdp, sep = "\n"),
    block = "functions"
  )

  dpars <- c("mu", "dprimef", "dprimer", "criterion",
             "spacing", "rcrit", "sigmar")
  links <- c("identity", model$links$dprimef,
             model$links$dprimer, model$links$criterion,
             model$links$spacing, model$links$rcrit,
             model$links$sigmar)

  if (kcrit_estimated) {
    dpars <- c(dpars, "kcrit")
    links <- c(links, model$links$kcrit)
  } else if (model$other_vars$has_guess) {
    dpars <- c(dpars, "kcrit")
    links <- c(links, model$links$kcrit)
  }

  has_guess <- model$other_vars$has_guess
  family_name <- if (has_guess) "sdt_cdp_uv3" else "sdt_cdp_uv"

  vint_vars <- paste0("vint", 1:8, "[n]")

  formula$family <- brms::custom_family(
    family_name,
    dpars = dpars,
    links = links,
    type = "int",
    loop = TRUE,
    log_lik = log_lik_sdt_cdp_uv,
    posterior_predict = posterior_predict_sdt_cdp,
    vars = vint_vars
  )

  nlist(formula, data, stanvars)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

# R-side: compute log-probability for a single CDP observation
.cdp_log_prob_r <- function(dprimef, dprimer, sigmar,
                            criterion, spacing, rcrit, kcrit,
                            cat_type, cat_conf, stim,
                            n_new, n_old, dist, thresh_type) {
  n_ratings <- n_new + n_old
  threshold_type_str <- if (thresh_type == 1L) {
    "parsimonious"
  } else {
    "equidistant"
  }
  thresholds <- .sdt_make_thresholds(
    criterion, n_ratings, threshold_type_str, spacing
  )

  has_guess_val <- !is.null(kcrit) && is.finite(kcrit)
  probs <- .sdt_cdp_category_probs(
    thresholds, dprimef, dprimer, sigmar,
    rcrit, kcrit, stim, n_new, n_old,
    .sdt_dist_names[dist]
  )

  # Map (cat_type, cat_conf) to index in probs
  # Probs layout: new(n_new), [guess(n_old)], know(n_old),
  #               remember(n_old)
  if (cat_type == 1L) {
    idx <- cat_conf
  } else if (has_guess_val) {
    block <- cat_type - 1L  # 2→1, 3→2, 4→3
    idx <- n_new + (block - 1L) * n_old + cat_conf
  } else {
    block <- if (cat_type == 3L) 1L else 2L
    idx <- n_new + (block - 1L) * n_old + cat_conf
  }
  log(max(probs[idx], .Machine$double.eps))
}

# Reverse map from dist_int to dist name
.sdt_dist_names <- c(
  "normal", "gumbel_min", "gumbel_max", "logistic"
)

log_lik_sdt_cdp_uv <- function(i, prep) {
  dprimef <- brms::get_dpar(prep, "dprimef", i = i)
  dprimer <- brms::get_dpar(prep, "dprimer", i = i)
  sigmar <- brms::get_dpar(prep, "sigmar", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  spacing <- brms::get_dpar(prep, "spacing", i = i)
  rcrit <- brms::get_dpar(prep, "rcrit", i = i)

  cat_type <- prep$data$vint1[i]
  cat_conf <- prep$data$vint2[i]
  stim <- prep$data$vint3[i]
  n_new <- prep$data$vint4[i]
  n_old <- prep$data$vint5[i]
  has_guess <- prep$data$vint6[i]
  dist_type <- prep$data$vint7[i]
  thresh_type <- prep$data$vint8[i]
  y <- prep$data$Y[i]

  has_kcrit <- "kcrit" %in% names(prep$dpars)
  kcrit_vec <- if (has_kcrit) {
    brms::get_dpar(prep, "kcrit", i = i)
  } else if (has_guess == 1L) {
    rep(-100, length(dprimef))
  } else {
    rep(NULL, length(dprimef))
  }

  log_p <- mapply(
    function(df, dr, sr, cr, sp, rc, kc) {
      .cdp_log_prob_r(df, dr, sr, cr, sp, rc, kc,
                      cat_type, cat_conf, stim,
                      n_new, n_old, dist_type, thresh_type)
    },
    dprimef, dprimer, sigmar,
    criterion, spacing, rcrit, kcrit_vec
  )

  y * log_p
}

posterior_predict_sdt_cdp <- function(i, prep, ...) {
  dprimef <- brms::get_dpar(prep, "dprimef", i = i)
  dprimer <- brms::get_dpar(prep, "dprimer", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)
  spacing <- brms::get_dpar(prep, "spacing", i = i)
  rcrit <- brms::get_dpar(prep, "rcrit", i = i)

  cat_type <- prep$data$vint1[i]
  cat_conf <- prep$data$vint2[i]
  stim <- prep$data$vint3[i]
  n_new <- prep$data$vint4[i]
  n_old <- prep$data$vint5[i]
  has_guess <- prep$data$vint6[i]
  dist_type <- prep$data$vint7[i]
  thresh_type <- prep$data$vint8[i]

  has_sigmar <- "sigmar" %in% names(prep$dpars)
  sigmar_vec <- if (has_sigmar) {
    brms::get_dpar(prep, "sigmar", i = i)
  } else {
    rep(0, length(dprimef))
  }

  has_kcrit <- "kcrit" %in% names(prep$dpars)
  kcrit_vec <- if (has_kcrit) {
    brms::get_dpar(prep, "kcrit", i = i)
  } else if (has_guess == 1L) {
    rep(-100, length(dprimef))
  } else {
    rep(NULL, length(dprimef))
  }

  mapply(
    function(df, dr, sr, cr, sp, rc, kc) {
      n_ratings <- n_new + n_old
      ttype <- if (thresh_type == 1L) "parsimonious" else "equidistant"
      thresholds <- .sdt_make_thresholds(cr, n_ratings, ttype, sp)
      dist_name <- .sdt_dist_names[dist_type]
      probs <- .sdt_cdp_category_probs(
        thresholds, df, dr, sr, rc, kc,
        stim, n_new, n_old, dist_name
      )
      has_guess_val <- !is.null(kc) && is.finite(kc)
      if (cat_type == 1L) {
        idx <- cat_conf
      } else if (has_guess_val) {
        block <- cat_type - 1L
        idx <- n_new + (block - 1L) * n_old + cat_conf
      } else {
        block <- if (cat_type == 3L) 1L else 2L
        idx <- n_new + (block - 1L) * n_old + cat_conf
      }
      probs[idx]
    },
    dprimef, dprimer, sigmar_vec,
    criterion, spacing, rcrit, kcrit_vec
  )
}
