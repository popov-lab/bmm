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

.sdt_dist_id <- function(dist) {
  .SDT_DISTS[[dist]]$id
}

# ordered by id so vint2[i] in log_lik/posterior_predict indexes the right name
.sdt_dist_names <- names(.SDT_DISTS)[order(vapply(.SDT_DISTS, `[[`, 0L, "id"))]


############################################################################# !
# DATA VALIDATION HELPERS                                                 ####
############################################################################# !

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
  n_trials <- rowSums(Y)

  stopif(any(n_trials <= 0, na.rm = TRUE),
    "Row sums of response columns must be positive (no empty rows)")

  data
}

.sdt_threshold_type_ids <- c(
  parsimonious = 1L,
  equidistant = 2L,
  log_distance = 3L,
  log_ratio = 4L,
  softmax = 5L
)

.sdt_threshold_type_id <- function(threshold_type) {
  .sdt_threshold_type_ids[[threshold_type]]
}

.sdt_threshold_type_name <- function(thresh_type) {
  names(.sdt_threshold_type_ids)[match(thresh_type, .sdt_threshold_type_ids)]
}

# Shared parameter specification for confidence-threshold parameterizations.
.sdt_threshold_parameter_parts <- function(n_ratings, threshold_type) {
  parameters <- list()
  default_priors <- list()
  param_links <- list()

  if (threshold_type %in% c("equidistant", "parsimonious")) {
    parameters$spacing <- glue(
      "Threshold spacing: controls distance between adjacent thresholds ",
      "(exp(spacing) ensures positive spacing)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$spacing <- "identity"
  } else if (threshold_type %in% c("log_distance", "log_ratio")) {
    n_deltas <- n_ratings - 2L
    mid <- n_ratings %/% 2L
    for (i in seq_len(n_deltas)) {
      idx <- if (i < mid) i else i + 1L
      pname <- paste0("delta", idx)
      parameters[[pname]] <- glue("Threshold parameter for threshold {idx}")
      default_priors[[pname]] <- list(
        main = "normal(0, 1)", effects = "normal(0, 0.5)"
      )
      param_links[[pname]] <- "identity"
    }
  } else if (threshold_type == "softmax") {
    parameters$spacing <- glue(
      "Average log spacing across threshold intervals ",
      "(exp(spacing) is the mean interval size)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)"
    )
    param_links$spacing <- "identity"

    n_deltas <- max(0L, n_ratings - 3L)
    if (n_deltas > 0L) {
      for (i in seq_len(n_deltas)) {
        pname <- paste0("delta", i)
        parameters[[pname]] <- glue(
          "Softmax threshold allocation parameter for interval {i}"
        )
        default_priors[[pname]] <- list(
          main = "normal(0, 1)", effects = "normal(0, 0.5)"
        )
        param_links[[pname]] <- "identity"
      }
    }
  }

  nlist(parameters, default_priors, param_links)
}

.sdt_threshold_delta_names <- function(object) {
  grep("^delta", names(object$parameters), value = TRUE)
}

.sdt_brms_dpar_name <- function(name) {
  if (grepl("^delta[0-9]+$", name)) paste0(name, "par") else name
}

.sdt_brms_dpar_map <- function(model) {
  stats::setNames(
    vapply(names(model$parameters), .sdt_brms_dpar_name, character(1)),
    names(model$parameters)
  )
}

.sdt_remap_formula_dpars <- function(formula, model) {
  dpar_map <- .sdt_brms_dpar_map(model)
  out <- formula

  for (idx in seq_along(out)) {
    old_name <- names(out)[idx]
    new_name <- if (old_name %in% names(dpar_map)) dpar_map[[old_name]] else old_name
    names(out)[idx] <- new_name
    if (inherits(out[[idx]], "formula") && length(out[[idx]]) >= 2L) {
      out[[idx]][[2]] <- as.name(new_name)
    }
  }

  out
}

.sdt_rating_long_data <- function(model, data) {
  resp_cols <- model$resp_vars$response
  stim_var <- model$other_vars$stimulus
  n_ratings <- model$other_vars$n_ratings

  Y <- as.matrix(data[resp_cols])
  n_trials <- rowSums(Y)
  stopif(any(n_trials <= 0, na.rm = TRUE),
    "Row sums of response columns must be positive (no empty rows)")

  other_cols <- setdiff(names(data), resp_cols)
  orig_nrow <- nrow(data)
  long_data <- data[rep(seq_len(orig_nrow), each = n_ratings),
                    other_cols, drop = FALSE]
  rownames(long_data) <- NULL

  long_data$obs_id <- rep(seq_len(orig_nrow), each = n_ratings)
  long_data$category <- rep(seq_len(n_ratings), times = orig_nrow)
  long_data$count <- as.vector(t(Y))
  long_data$stim_val <- as.integer(long_data[[stim_var]])
  long_data$n_ratings <- n_ratings
  long_data$n_trials_total <- rep(n_trials, each = n_ratings)
  long_data$dist_type <- model$other_vars$dist_int
  long_data$thresh_type <- model$other_vars$thresh_type_int
  long_data
}

.sdt_rating_base_formula <- function() {
  brms::bf(
    "count | vint(category, stim_val, n_ratings, dist_type, thresh_type) ~ 0"
  )
}

.sdt_rating_custom_family <- function(model, family_name,
                                      log_lik, posterior_predict) {
  dpar_map <- .sdt_brms_dpar_map(model)
  model_dpars <- names(model$parameters)
  dpars <- c("mu", unname(dpar_map[model_dpars]))
  links <- c("identity", unname(unlist(model$links[model_dpars])))
  brms::custom_family(
    family_name,
    dpars = dpars,
    links = links,
    type = "int",
    loop = TRUE,
    log_lik = log_lik,
    posterior_predict = posterior_predict,
    vars = paste0("vint", 1:5, "[n]")
  )
}

.sdt_rating_stan_wrapper <- function(model, family_name, helper_fun) {
  dpar_map <- .sdt_brms_dpar_map(model)
  dpars <- unname(dpar_map[names(model$parameters)])
  sig <- paste(
    c("int y", "real mu", paste("real", dpars),
      "int category", "int stimulus", "int n_ratings",
      "int dist_type", "int thresh_type"),
    collapse = ", "
  )

  has_spacing <- "spacing" %in% names(model$parameters)
  delta_names <- .sdt_threshold_delta_names(model)
  delta_dpars <- unname(dpar_map[delta_names])
  spacing_value <- if (has_spacing) dpar_map[["spacing"]] else "0"
  n_deltas <- length(delta_dpars)

  delta_code <- if (n_deltas == 0L) {
    "  array[0] real deltas;\n"
  } else {
    paste0(
      "  array[", n_deltas, "] real deltas;\n",
      paste(
        sprintf("  deltas[%d] = %s;", seq_len(n_deltas), delta_dpars),
        collapse = "\n"
      ),
      "\n"
    )
  }

  # Stan 2.33+ requires `|` between first arg and the rest for _lpmf functions
  helper_args <- switch(
    family_name,
    sdt_rating = paste(
      "y | category, stimulus, n_ratings, dist_type, thresh_type,",
      "dprime, criterion,", spacing_value, ", sdratio, deltas"
    ),
    sdt_dp = paste(
      "y | category, stimulus, n_ratings, dist_type, thresh_type,",
      "dprime, criterion,", spacing_value, ", Ro, Rn, sdratio, deltas"
    ),
    sdt_metad = paste(
      "y | category, stimulus, n_ratings, dist_type, thresh_type,",
      "dprime, criterion,", spacing_value, ", metad, sdratio, deltas"
    ),
    stop2("Unsupported SDT family wrapper: '{family_name}'")
  )

  paste0(
    "real ", family_name, "_lpmf(", sig, ") {\n",
    delta_code,
    "  return ", helper_fun, "(", helper_args, ");\n",
    "}\n"
  )
}

.sdt_rating_stanvars <- function(model, family_name, helper_fun) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_binary <- read_lines2(paste0(sc_path, "/sdt_binary_funs.stan"))
  stan_rating <- read_lines2(paste0(sc_path, "/sdt_rating_funs.stan"))
  stan_wrapper <- .sdt_rating_stan_wrapper(model, family_name, helper_fun)
  brms::stanvar(
    scode = paste(stan_binary, stan_rating, stan_wrapper, sep = "\n"),
    block = "functions"
  )
}

.sdt_rating_thresholds_from_draw <- function(prep, i, draw_idx) {
  threshold_type <- .sdt_threshold_type_name(prep$data$vint5[i])
  n_ratings <- prep$data$vint3[i]
  criterion <- brms::get_dpar(prep, .sdt_brms_dpar_name("criterion"), i = i)[draw_idx]

  spacing <- if ("spacing" %in% names(prep$dpars)) {
    brms::get_dpar(prep, .sdt_brms_dpar_name("spacing"), i = i)[draw_idx]
  } else {
    NULL
  }

  delta_names <- grep("^delta", names(prep$dpars), value = TRUE)
  deltas <- if (length(delta_names) == 0L) {
    NULL
  } else {
    vapply(delta_names, function(x) brms::get_dpar(prep, x, i = i)[draw_idx],
           numeric(1))
  }

  .sdt_make_thresholds(
    criterion = criterion,
    n_ratings = n_ratings,
    threshold_type = threshold_type,
    spacing = spacing,
    deltas = deltas
  )
}

.sdt_loglik_rating_common <- function(i, prep, variant) {
  draws <- length(brms::get_dpar(prep, .sdt_brms_dpar_name("criterion"), i = i))
  category <- prep$data$vint1[i]
  stimulus <- prep$data$vint2[i]
  dist <- .sdt_dist_names[prep$data$vint4[i]]
  y <- prep$data$Y[i]

  dprime <- brms::get_dpar(prep, .sdt_brms_dpar_name("dprime"), i = i)
  sdratio <- if ("sdratio" %in% names(prep$dpars)) {
    exp(brms::get_dpar(prep, .sdt_brms_dpar_name("sdratio"), i = i))
  } else {
    rep(1, draws)
  }

  if (variant == "rating") {
    vapply(seq_len(draws), function(draw_idx) {
      thresholds <- .sdt_rating_thresholds_from_draw(prep, i, draw_idx)
      probs <- .sdt_category_probs(
        thresholds = thresholds,
        shift = dprime[draw_idx] / 2 * (2 * stimulus - 1),
        sdratio = sdratio[draw_idx],
        stimulus = stimulus,
        dist = dist
      )
      y * log(probs[category])
    }, numeric(1))
  } else if (variant == "dp") {
    Ro <- brms::get_dpar(prep, .sdt_brms_dpar_name("Ro"), i = i)
    Rn <- brms::get_dpar(prep, .sdt_brms_dpar_name("Rn"), i = i)
    vapply(seq_len(draws), function(draw_idx) {
      thresholds <- .sdt_rating_thresholds_from_draw(prep, i, draw_idx)
      probs <- .sdt_dpsdt_category_probs(
        thresholds = thresholds,
        shift = dprime[draw_idx] / 2 * (2 * stimulus - 1),
        sdratio = sdratio[draw_idx],
        stimulus = stimulus,
        dist = dist,
        Ro = plogis(Ro[draw_idx]),
        Rn = plogis(Rn[draw_idx])
      )
      y * log(probs[category])
    }, numeric(1))
  } else {
    metad <- brms::get_dpar(prep, .sdt_brms_dpar_name("metad"), i = i)
    vapply(seq_len(draws), function(draw_idx) {
      thresholds <- .sdt_rating_thresholds_from_draw(prep, i, draw_idx)
      probs <- .sdt_metad_category_probs(
        thresholds = thresholds,
        dprime = dprime[draw_idx],
        metad = metad[draw_idx],
        stimulus = stimulus,
        sdratio = sdratio[draw_idx],
        dist = dist
      )
      y * log(probs[category])
    }, numeric(1))
  }
}

.sdt_posterior_predict_rating_common <- function(i, prep, variant) {
  draws <- length(brms::get_dpar(prep, .sdt_brms_dpar_name("criterion"), i = i))
  category <- prep$data$vint1[i]
  stimulus <- prep$data$vint2[i]
  dist <- .sdt_dist_names[prep$data$vint4[i]]

  dprime <- brms::get_dpar(prep, .sdt_brms_dpar_name("dprime"), i = i)
  sdratio <- if ("sdratio" %in% names(prep$dpars)) {
    exp(brms::get_dpar(prep, .sdt_brms_dpar_name("sdratio"), i = i))
  } else {
    rep(1, draws)
  }

  if (variant == "rating") {
    vapply(seq_len(draws), function(draw_idx) {
      thresholds <- .sdt_rating_thresholds_from_draw(prep, i, draw_idx)
      probs <- .sdt_category_probs(
        thresholds = thresholds,
        shift = dprime[draw_idx] / 2 * (2 * stimulus - 1),
        sdratio = sdratio[draw_idx],
        stimulus = stimulus,
        dist = dist
      )
      probs[category]
    }, numeric(1))
  } else if (variant == "dp") {
    Ro <- brms::get_dpar(prep, .sdt_brms_dpar_name("Ro"), i = i)
    Rn <- brms::get_dpar(prep, .sdt_brms_dpar_name("Rn"), i = i)
    vapply(seq_len(draws), function(draw_idx) {
      thresholds <- .sdt_rating_thresholds_from_draw(prep, i, draw_idx)
      probs <- .sdt_dpsdt_category_probs(
        thresholds = thresholds,
        shift = dprime[draw_idx] / 2 * (2 * stimulus - 1),
        sdratio = sdratio[draw_idx],
        stimulus = stimulus,
        dist = dist,
        Ro = plogis(Ro[draw_idx]),
        Rn = plogis(Rn[draw_idx])
      )
      probs[category]
    }, numeric(1))
  } else {
    metad <- brms::get_dpar(prep, .sdt_brms_dpar_name("metad"), i = i)
    vapply(seq_len(draws), function(draw_idx) {
      thresholds <- .sdt_rating_thresholds_from_draw(prep, i, draw_idx)
      probs <- .sdt_metad_category_probs(
        thresholds = thresholds,
        dprime = dprime[draw_idx],
        metad = metad[draw_idx],
        stimulus = stimulus,
        sdratio = sdratio[draw_idx],
        dist = dist
      )
      probs[category]
    }, numeric(1))
  }
}


