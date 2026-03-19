############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdai <- function(resp = NULL, n_sources = 2, n_confidence = 4,
                        dist = "gumbel_min", links = NULL, call = NULL, ...) {
  dist <- match.arg(dist, c("gumbel_min", "normal"))
  stopif(n_sources < 2, "n_sources must be at least 2")
  stopif(n_confidence < 2, "n_confidence must be at least 2")

  n_thresholds <- n_confidence - 1

  parameters <- list(
    mu = "Internal parameter (fixed to 0)",
    dprime = glue(
      "Source discrimination: distance between source distributions. ",
      "Higher values indicate better source discrimination."
    ),
    criterion = "Central response criterion (anchor threshold)"
  )

  default_priors <- list(
    dprime = list(main = "student_t(3, 1, 2)", effects = "normal(0, 1)"),
    criterion = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
  )

  param_links <- list(
    mu = "identity",
    dprime = "identity",
    criterion = "identity"
  )

  # threshold margin parameters use letter suffixes (brms rejects trailing digits)
  n_margins <- n_thresholds - 1
  margin_suffixes <- letters[seq_len(n_margins)]
  if (n_margins > 0) {
    for (i in seq_len(n_margins)) {
      lname <- paste0("crl", margin_suffixes[i])
      hname <- paste0("crh", margin_suffixes[i])
      parameters[[lname]] <- glue(
        "Left threshold margin {i}: log-distance below criterion"
      )
      parameters[[hname]] <- glue(
        "Right threshold margin {i}: log-distance above criterion"
      )
      default_priors[[lname]] <- list(
        main = "normal(-0.5, 0.5)", effects = "normal(0, 0.3)"
      )
      default_priors[[hname]] <- list(
        main = "normal(-0.5, 0.5)", effects = "normal(0, 0.3)"
      )
      param_links[[lname]] <- "log"
      param_links[[hname]] <- "log"
    }
  }

  # Gaussian UV-SDT adds sdratio
  if (dist == "normal") {
    parameters$sdratio <- glue(
      "SD ratio: log ratio of signal to noise standard deviations. ",
      "exp(sdratio) is the actual ratio (0 = equal variance)."
    )
    default_priors$sdratio <- list(
      main = "student_t(3, 0.5, 1)", effects = "normal(0, 0.5)"
    )
    param_links$sdratio <- "log"
  }

  version <- ifelse(dist == "gumbel_min", "gumbel", "gaussian")
  resp_list <- if (is.character(resp)) list(resp = resp) else nlist(resp)

  out <- structure(
    list(
      resp_vars = resp_list,
      other_vars = nlist(n_sources, n_confidence, n_thresholds, n_margins,
                         margin_suffixes, dist),
      domain = "Signal detection",
      task = "Source memory",
      name = glue(
        "Source-Dependent Asymmetric Information model (SDAI; ",
        "Meyer-Grant et al., 2025)"
      ),
      citation = glue(
        "Meyer-Grant, C. G., Kellen, D., & Klauer, K. C. (2025). ",
        "A Gumbel-distribution model of recognition memory."
      ),
      version = version,
      requirements = glue(
        "- Aggregated count data: columns for each identification x confidence ",
        "combination\n",
        "- Response column names should follow the pattern ",
        "'i{{source}}_c{{confidence}}'\n",
        "- n_sources x n_confidence response columns required"
      ),
      parameters = parameters,
      links = param_links,
      fixed_parameters = list(mu = 0),
      default_priors = default_priors,
      void_mu = FALSE
    ),
    class = c("bmmodel", "sdai", paste0("sdai_", version)),
    call = call
  )
  out$links[names(links)] <- links
  out
}


#' @title `r .model_sdai()$name`
#' @name sdai
#' @details
#' The SDAI model is a signal detection theory model for source memory data
#' where participants make both identification decisions (which source?) and
#' confidence judgments.
#'
#' The response is a grid of identification categories x confidence levels,
#' producing two independent multinomial processes:
#' \itemize{
#'   \item Table process (p_tabs): discrimination between sources
#'   \item Test-present process (p_tpres): source confusion probabilities
#' }
#'
#' `r model_info(.model_sdai())`
#'
#' @param resp Character vector of response column names. Should be ordered
#'   as identification x confidence (e.g., c("i1_c1", "i1_c2", ..., "i3_c4")).
#'   The first `n_confidence` columns are for the first identification category,
#'   the next `n_confidence` for the second, etc.
#' @param n_sources Integer. Number of source identification categories
#'   (default 2). Standard source memory uses 3 (lure/neutral/target).
#' @param n_confidence Integer. Number of confidence levels (default 4).
#' @param dist Character. Noise distribution: "gumbel_min" (default) for
#'   Gumbel-min SDT or "normal" for Gaussian UV-SDT.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @export
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # Create an SDAI model for 3 identification x 4 confidence
#' resp_cols <- c("i1_c1", "i1_c2", "i1_c3", "i1_c4",
#'                "i2_c1", "i2_c2", "i2_c3", "i2_c4",
#'                "i3_c1", "i3_c2", "i3_c3", "i3_c4")
#' model <- sdai(resp = resp_cols, n_sources = 3, n_confidence = 4)
sdai <- function(resp, n_sources = 2, n_confidence = 4,
                 dist = "gumbel_min", ...) {
  call <- match.call()
  stop_missing_args()
  .model_sdai(resp = resp, n_sources = n_sources, n_confidence = n_confidence,
              dist = dist, call = call, ...)
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdai <- function(model, data, formula) {
  resp_cols <- model$resp_vars$resp
  n_sources <- model$other_vars$n_sources
  n_confidence <- model$other_vars$n_confidence
  n_expected <- n_sources * n_confidence

  stopif(length(resp_cols) != n_expected,
         "Expected {n_expected} response columns ({n_sources} sources x {n_confidence} confidence levels), but got {length(resp_cols)}")

  missing_cols <- setdiff(resp_cols, names(data))
  stopif(length(missing_cols) > 0,
         "Response columns not found in data: {collapse_comma(missing_cols)}")

  for (col in resp_cols) {
    stopif(any(data[[col]] < 0, na.rm = TRUE),
           "Response column '{col}' contains negative values")
  }

  # Store response matrix and row totals as attributes
  resp_mat <- as.matrix(data[resp_cols])
  attr(data, "resp_cols") <- resp_cols
  attr(data, "resp_mat") <- resp_mat

  # Compute n_trials for tabs and tpres
  tabs_cols <- resp_cols[seq_len(n_confidence)]
  tpres_cols <- resp_cols[(n_confidence + 1):n_expected]
  attr(data, "tabs_cols") <- tabs_cols
  attr(data, "tpres_cols") <- tpres_cols

  data$nTrials_tabs <- rowSums(data[tabs_cols])
  data$nTrials_tpres <- rowSums(data[tpres_cols])

  NextMethod("check_data")
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

# Internal: build dpars, links, and lower bounds for SDAI custom family
.sdai_build_dpars <- function(model) {
  n_margins <- model$other_vars$n_margins
  dist <- model$other_vars$dist

  dpars <- c("mu", "dprime", "criterion")
  dpar_links <- c("identity", model$links$dprime, model$links$criterion)
  dpar_lb <- c(NA, NA, NA)

  if (dist == "normal") {
    dpars <- c(dpars, "sdratio")
    dpar_links <- c(dpar_links, model$links$sdratio)
    dpar_lb <- c(dpar_lb, 0)
  }

  if (n_margins > 0) {
    suffixes <- model$other_vars$margin_suffixes
    for (i in seq_len(n_margins)) {
      dpars <- c(dpars, paste0("crl", suffixes[i]), paste0("crh", suffixes[i]))
      dpar_links <- c(dpar_links, "log", "log")
      dpar_lb <- c(dpar_lb, 0, 0)
    }
  }

  nlist(dpars, dpar_links, dpar_lb)
}

# Internal: generate Stan lpmf code for SDAI Gumbel model
.sdai_generate_gumbel_lpmf <- function(n_confidence, n_sources, n_margins,
                                       margin_suffixes = letters[seq_len(n_margins)]) {
  n_thres <- n_confidence - 1
  n_tabs <- n_confidence
  n_tpres <- (n_sources - 1) * n_confidence
  n_vint <- n_sources * n_confidence - 1

  dpar_args <- "real mu, real dprime, real criterion"
  if (n_margins > 0) {
    for (i in seq_len(n_margins)) {
      s <- margin_suffixes[i]
      dpar_args <- paste0(dpar_args, ", real crl", s, ", real crh", s)
    }
  }

  vint_args <- paste(paste0("int y", seq_len(n_vint)), collapse = ", ")

  thres_code <- "    vector[N_THRES] thres;\n"
  if (n_thres == 1) {
    thres_code <- paste0(thres_code, "    thres[1] = criterion;\n")
  } else {
    mid <- ceiling(n_thres / 2)
    thres_code <- paste0(thres_code, "    thres[", mid, "] = criterion;\n")
    for (i in seq_len(n_margins)) {
      if (mid - i >= 1) {
        cum_terms <- paste0("crl", margin_suffixes[seq_len(i)], collapse = " + ")
        thres_code <- paste0(thres_code, "    thres[", mid - i,
                             "] = criterion - (", cum_terms, ");\n")
      }
      if (mid + i <= n_thres) {
        cum_terms <- paste0("crh", margin_suffixes[seq_len(i)], collapse = " + ")
        thres_code <- paste0(thres_code, "    thres[", mid + i,
                             "] = criterion + (", cum_terms, ");\n")
      }
    }
  }
  thres_code <- gsub("N_THRES", n_thres, thres_code)

  # Build response count arrays
  count_code <- paste0("    array[", n_tabs, "] int res_tabs;\n")
  count_code <- paste0(count_code, "    res_tabs[1] = y;\n")
  for (j in 2:n_tabs) {
    count_code <- paste0(count_code, "    res_tabs[", j, "] = y", j - 1, ";\n")
  }
  count_code <- paste0(count_code, "    array[", n_tpres, "] int res_tpres;\n")
  for (j in seq_len(n_tpres)) {
    count_code <- paste0(count_code, "    res_tpres[", j, "] = y",
                         n_tabs - 1 + j, ";\n")
  }

  # Build tabs probability computation
  tabs_code <- paste0("    vector[", n_tabs, "] p_tabs;\n")
  tabs_code <- paste0(tabs_code,
                       "    p_tabs[1] = sdai_p_hit_gumbel(-100, thres[1]);\n")
  if (n_tabs > 2) {
    for (k in 2:(n_tabs - 1)) {
      tabs_code <- paste0(tabs_code, "    p_tabs[", k,
                           "] = sdai_p_hit_gumbel(thres[", k - 1,
                           "], thres[", k, "]);\n")
    }
  }
  tabs_code <- paste0(tabs_code, "    p_tabs[", n_tabs,
                       "] = sdai_p_hit_gumbel(thres[", n_thres, "], 100);\n")

  # Build tpres probability computation (2 sources)
  tpres_code <- paste0("    vector[", n_tpres, "] p_tpres;\n")
  # Source 1: g1=0, g2=dprime
  tpres_code <- paste0(tpres_code,
                         "    p_tpres[1] = sdai_pfun(-100, thres[1], 0, dprime);\n")
  if (n_confidence > 2) {
    for (k in 2:(n_confidence - 1)) {
      tpres_code <- paste0(tpres_code, "    p_tpres[", k,
                             "] = sdai_pfun(thres[", k - 1,
                             "], thres[", k, "], 0, dprime);\n")
    }
  }
  tpres_code <- paste0(tpres_code, "    p_tpres[", n_confidence,
                         "] = sdai_pfun(thres[", n_thres, "], 100, 0, dprime);\n")

  # Source 2: g1=dprime, g2=0
  offset <- n_confidence
  tpres_code <- paste0(tpres_code, "    p_tpres[", offset + 1,
                         "] = sdai_pfun(-100, thres[1], dprime, 0);\n")
  if (n_confidence > 2) {
    for (k in 2:(n_confidence - 1)) {
      tpres_code <- paste0(tpres_code, "    p_tpres[", offset + k,
                             "] = sdai_pfun(thres[", k - 1,
                             "], thres[", k, "], dprime, 0);\n")
    }
  }
  tpres_code <- paste0(tpres_code, "    p_tpres[", offset + n_confidence,
                         "] = sdai_pfun(thres[", n_thres,
                         "], 100, dprime, 0);\n")

  # Assemble full lpmf
  paste0(
    "  real sdai_gumbel_lpmf(int y, ", dpar_args, ", ", vint_args, ") {\n",
    thres_code,
    count_code,
    tabs_code,
    tpres_code,
    "    return multinomial_lpmf(res_tabs | p_tabs) + ",
    "multinomial_lpmf(res_tpres | p_tpres);\n",
    "  }\n"
  )
}

# Internal: generate Stan lpmf code for SDAI Gaussian model
.sdai_generate_gaussian_lpmf <- function(n_confidence, n_sources, n_margins,
                                         margin_suffixes = letters[seq_len(n_margins)]) {
  n_thres <- n_confidence - 1
  n_tabs <- n_confidence
  n_tpres <- (n_sources - 1) * n_confidence
  n_vint <- n_sources * n_confidence - 1

  dpar_args <- "real mu, real dprime, real criterion, real sdratio"
  if (n_margins > 0) {
    for (i in seq_len(n_margins)) {
      s <- margin_suffixes[i]
      dpar_args <- paste0(dpar_args, ", real crl", s, ", real crh", s)
    }
  }

  vint_args <- paste(paste0("int y", seq_len(n_vint)), collapse = ", ")

  thres_code <- "    vector[N_THRES] thres;\n    real disc = 1.0 / sdratio;\n"
  if (n_thres == 1) {
    thres_code <- paste0(thres_code, "    thres[1] = criterion;\n")
  } else {
    mid <- ceiling(n_thres / 2)
    thres_code <- paste0(thres_code, "    thres[", mid, "] = criterion;\n")
    for (i in seq_len(n_margins)) {
      if (mid - i >= 1) {
        cum_terms <- paste0("crl", margin_suffixes[seq_len(i)], collapse = " + ")
        thres_code <- paste0(thres_code, "    thres[", mid - i,
                             "] = criterion - (", cum_terms, ");\n")
      }
      if (mid + i <= n_thres) {
        cum_terms <- paste0("crh", margin_suffixes[seq_len(i)], collapse = " + ")
        thres_code <- paste0(thres_code, "    thres[", mid + i,
                             "] = criterion + (", cum_terms, ");\n")
      }
    }
  }
  thres_code <- gsub("N_THRES", n_thres, thres_code)

  # Count arrays (same as Gumbel)
  count_code <- paste0("    array[", n_tabs, "] int res_tabs;\n")
  count_code <- paste0(count_code, "    res_tabs[1] = y;\n")
  for (j in 2:n_tabs) {
    count_code <- paste0(count_code, "    res_tabs[", j, "] = y", j - 1, ";\n")
  }
  count_code <- paste0(count_code, "    array[", n_tpres, "] int res_tpres;\n")
  for (j in seq_len(n_tpres)) {
    count_code <- paste0(count_code, "    res_tpres[", j, "] = y",
                         n_tabs - 1 + j, ";\n")
  }

  # Tabs: Gaussian discrimination
  tabs_code <- paste0("    vector[", n_tabs, "] p_tabs;\n")
  tabs_code <- paste0(tabs_code,
                       "    p_tabs[1] = sdai_p_hit_gaussian(negative_infinity(), thres[1]);\n")
  if (n_tabs > 2) {
    for (k in 2:(n_tabs - 1)) {
      tabs_code <- paste0(tabs_code, "    p_tabs[", k,
                           "] = sdai_p_hit_gaussian(thres[", k - 1,
                           "], thres[", k, "]);\n")
    }
  }
  tabs_code <- paste0(tabs_code, "    p_tabs[", n_tabs,
                       "] = sdai_p_hit_gaussian(thres[", n_thres,
                       "], positive_infinity());\n")

  # Tpres: numerical integration for source confusion
  tpres_code <- paste0("    vector[", n_tpres, "] p_tpres;\n")

  # Source 1: d1=0, s1=1, d2=dprime, s2=disc
  bounds <- list()
  for (k in seq_len(n_confidence)) {
    lower <- if (k == 1) "negative_infinity()" else paste0("thres[", k - 1, "]")
    upper <- if (k == n_confidence) "positive_infinity()" else paste0("thres[", k, "]")
    tpres_code <- paste0(tpres_code,
      "    p_tpres[", k, "] = integrate_1d(sdai_int_inst_uvg, ",
      lower, ", ", upper,
      ", { 0, 1, dprime, disc }, x_r, x_i);\n")
  }

  # Source 2: d1=dprime, s1=disc, d2=0, s2=1
  offset <- n_confidence
  for (k in seq_len(n_confidence)) {
    lower <- if (k == 1) "negative_infinity()" else paste0("thres[", k - 1, "]")
    upper <- if (k == n_confidence) "positive_infinity()" else paste0("thres[", k, "]")
    tpres_code <- paste0(tpres_code,
      "    p_tpres[", offset + k, "] = integrate_1d(sdai_int_inst_uvg, ",
      lower, ", ", upper,
      ", { dprime, disc, 0, 1 }, x_r, x_i);\n")
  }

  paste0(
    "  real sdai_gaussian_lpmf(int y, ", dpar_args, ", ", vint_args,
    ", data array[] real x_r, data array[] int x_i) {\n",
    thres_code,
    count_code,
    tabs_code,
    tpres_code,
    "    return multinomial_lpmf(res_tabs | p_tabs) + ",
    "multinomial_lpmf(res_tpres | p_tpres);\n",
    "  }\n"
  )
}

#' @export
configure_model.sdai_gumbel <- function(model, data, formula) {
  resp_cols <- attr(data, "resp_cols")
  n_confidence <- model$other_vars$n_confidence
  n_sources <- model$other_vars$n_sources
  n_margins <- model$other_vars$n_margins
  margin_suffixes <- model$other_vars$margin_suffixes
  n_vint <- length(resp_cols) - 1

  dpar_info <- .sdai_build_dpars(model)

  ll_fn <- .sdai_make_log_lik(n_confidence, n_sources, n_margins,
                               margin_suffixes, "gumbel")
  pp_fn <- .sdai_make_posterior_predict(n_confidence, n_sources, n_margins,
                                        margin_suffixes, "gumbel")

  family <- brms::custom_family(
    "sdai_gumbel",
    dpars = dpar_info$dpars,
    links = dpar_info$dpar_links,
    lb = dpar_info$dpar_lb,
    type = "int",
    loop = TRUE,
    log_lik = ll_fn,
    posterior_predict = pp_fn,
    vars = paste0("vint", seq_len(n_vint), "[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_helpers <- read_lines2(paste0(sc_path, "/sdai_gumbel_funs.stan"))
  stan_lpmf <- .sdai_generate_gumbel_lpmf(n_confidence, n_sources, n_margins,
                                            margin_suffixes)

  stanvars <- brms::stanvar(scode = stan_helpers, block = "functions") +
    brms::stanvar(scode = stan_lpmf, block = "functions")

  formula <- bmf2bf(model, formula)
  formula$family <- family

  nlist(formula, data, stanvars)
}

#' @export
configure_model.sdai_gaussian <- function(model, data, formula) {
  resp_cols <- attr(data, "resp_cols")
  n_confidence <- model$other_vars$n_confidence
  n_sources <- model$other_vars$n_sources
  n_margins <- model$other_vars$n_margins
  margin_suffixes <- model$other_vars$margin_suffixes
  n_vint <- length(resp_cols) - 1

  dpar_info <- .sdai_build_dpars(model)

  ll_fn <- .sdai_make_log_lik(n_confidence, n_sources, n_margins,
                               margin_suffixes, "gaussian")
  pp_fn <- .sdai_make_posterior_predict(n_confidence, n_sources, n_margins,
                                        margin_suffixes, "gaussian")

  family <- brms::custom_family(
    "sdai_gaussian",
    dpars = dpar_info$dpars,
    links = dpar_info$dpar_links,
    lb = dpar_info$dpar_lb,
    type = "int",
    loop = TRUE,
    log_lik = ll_fn,
    posterior_predict = pp_fn,
    vars = c(paste0("vint", seq_len(n_vint), "[n]"), "x_r", "x_i")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_helpers <- read_lines2(paste0(sc_path, "/sdai_gaussian_funs.stan"))
  stan_lpmf <- .sdai_generate_gaussian_lpmf(n_confidence, n_sources, n_margins,
                                              margin_suffixes)

  stanvars <- brms::stanvar(scode = stan_helpers, block = "functions") +
    brms::stanvar(scode = stan_lpmf, block = "functions")

  tdata_code <- "  array[0] real x_r;\n  array[0] int x_i;"
  stanvars <- stanvars +
    brms::stanvar(scode = tdata_code, block = "tdata")

  formula <- bmf2bf(model, formula)
  formula$family <- family

  nlist(formula, data, stanvars)
}


############################################################################# !
# BMF2BF (Formula Conversion)                                            ####
############################################################################# !

#' @export
bmf2bf.sdai <- function(model, formula) {
  resp_cols <- model$resp_vars$resp

  # Build base formula: first column | vint(rest) ~ 1
  vint_str <- paste(resp_cols[-1], collapse = ", ")
  base_str <- paste0(resp_cols[1], " | vint(", vint_str, ") ~ 1")
  brms::bf(as.formula(base_str))
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

# Internal: compute SDAI thresholds from criterion and margin parameters
.sdai_thresholds <- function(criterion, margins_l, margins_h) {
  n_margins <- length(margins_l)
  n_thresholds <- n_margins + 1

  thresholds <- matrix(NA_real_, nrow = length(criterion), ncol = n_thresholds)
  # Central threshold is the criterion
  mid <- ceiling(n_thresholds / 2)

  if (n_thresholds == 1) {
    thresholds[, 1] <- criterion
    return(thresholds)
  }

  if (n_thresholds == 3) {
    thresholds[, 1] <- criterion - margins_l[, 1]
    thresholds[, 2] <- criterion
    thresholds[, 3] <- criterion + margins_h[, 1]
    return(thresholds)
  }

  # General case: build outward from criterion
  thresholds[, mid] <- criterion
  cum_l <- 0
  cum_h <- 0
  for (i in seq_len(n_margins)) {
    if (mid - i >= 1) {
      cum_l <- cum_l + margins_l[, i]
      thresholds[, mid - i] <- criterion - cum_l
    }
    if (mid + i <= n_thresholds) {
      cum_h <- cum_h + margins_h[, i]
      thresholds[, mid + i] <- criterion + cum_h
    }
  }
  thresholds
}

# Internal: Gumbel-min CDF
.gumbel_min_cdf <- function(x) {
  1 - exp(-exp(x))
}

# Internal: compute SDAI table probabilities (discrimination) for Gumbel-min
.sdai_tabs_probs_gumbel <- function(thresholds) {
  n_thresholds <- ncol(thresholds)
  n_cats <- n_thresholds + 1
  n_draws <- nrow(thresholds)

  probs <- matrix(NA_real_, nrow = n_draws, ncol = n_cats)

  for (k in seq_len(n_cats)) {
    if (k == 1L) {
      probs[, k] <- .gumbel_min_cdf(thresholds[, 1])^2
    } else if (k == n_cats) {
      probs[, k] <- 1 - .gumbel_min_cdf(thresholds[, n_thresholds])^2
    } else {
      probs[, k] <- .gumbel_min_cdf(thresholds[, k])^2 -
        .gumbel_min_cdf(thresholds[, k - 1])^2
    }
  }
  pmax(probs, .Machine$double.eps)
}

# Internal: compute SDAI test-present probabilities for Gumbel-min
# Source 1 has location 0, Source 2 has location mu (dprime)
.sdai_tpres_probs_gumbel <- function(thresholds, dprime) {
  n_thresholds <- ncol(thresholds)
  n_cats <- n_thresholds + 1
  n_draws <- nrow(thresholds)

  # 2 sources x n_cats categories
  probs <- matrix(NA_real_, nrow = n_draws, ncol = 2 * n_cats)

  for (draw in seq_len(n_draws)) {
    g1_src1 <- 0
    g2_src1 <- dprime[draw]
    g1_src2 <- dprime[draw]
    g2_src2 <- 0

    for (k in seq_len(n_cats)) {
      lower <- if (k == 1L) -100 else thresholds[draw, k - 1]
      upper <- if (k == n_cats) 100 else thresholds[draw, k]

      probs[draw, k] <- .sdai_pfun_gumbel(lower, upper, g1_src1, g2_src1)
      probs[draw, n_cats + k] <- .sdai_pfun_gumbel(lower, upper, g1_src2, g2_src2)
    }
  }
  pmax(probs, .Machine$double.eps)
}

# Internal: Gumbel SDAI pfun (source confusion probability in interval)
.sdai_pfun_gumbel <- function(l, u, g1, g2) {
  denom <- exp(g1) + exp(g2)
  p1_l <- exp(-exp(-g1 + l) + g1)
  p1_u <- exp(-exp(-g1 + u) + g1)
  p2_l <- exp(-exp(-g1 + l) + g2)
  p2_u <- exp(-exp(-g1 + u) + g2)
  p3_u <- exp(-exp(-g1 + u) - exp(-g2 + u) + g2)
  p3_l <- exp(-exp(-g1 + l) - exp(-g2 + l) + g2)

  (p1_l - p1_u + p2_l - p2_u + p3_u - p3_l) / denom
}

# Internal: compute SDAI table probabilities for Gaussian
.sdai_tabs_probs_gaussian <- function(thresholds) {
  n_thresholds <- ncol(thresholds)
  n_cats <- n_thresholds + 1
  n_draws <- nrow(thresholds)

  probs <- matrix(NA_real_, nrow = n_draws, ncol = n_cats)
  for (k in seq_len(n_cats)) {
    if (k == 1L) {
      probs[, k] <- stats::pnorm(thresholds[, 1])^2
    } else if (k == n_cats) {
      probs[, k] <- 1 - stats::pnorm(thresholds[, n_thresholds])^2
    } else {
      probs[, k] <- stats::pnorm(thresholds[, k])^2 -
        stats::pnorm(thresholds[, k - 1])^2
    }
  }
  pmax(probs, .Machine$double.eps)
}

# Internal: compute SDAI test-present probabilities for Gaussian
.sdai_tpres_probs_gaussian <- function(thresholds, dprime, sdratio) {
  n_thresholds <- ncol(thresholds)
  n_cats <- n_thresholds + 1
  n_draws <- nrow(thresholds)

  probs <- matrix(NA_real_, nrow = n_draws, ncol = 2 * n_cats)

  for (draw in seq_len(n_draws)) {
    disc <- sdratio[draw]
    for (k in seq_len(n_cats)) {
      lower <- if (k == 1L) -Inf else thresholds[draw, k - 1]
      upper <- if (k == n_cats) Inf else thresholds[draw, k]

      # Source 1: d1=0, s1=1; Source 2: d2=dprime, s2=disc
      probs[draw, k] <- stats::integrate(
        function(x) stats::pnorm(x, dprime[draw], disc)^1 *
          stats::dnorm(x, 0, 1),
        lower, upper, stop.on.error = FALSE
      )$value

      # Source 2: d1=dprime, s1=disc; Source 1: d2=0, s2=1
      probs[draw, n_cats + k] <- stats::integrate(
        function(x) stats::pnorm(x, 0, 1)^1 *
          stats::dnorm(x, dprime[draw], disc),
        lower, upper, stop.on.error = FALSE
      )$value
    }
  }
  pmax(probs, .Machine$double.eps)
}

# Internal: extract SDAI parameters from brms prep object
.sdai_extract_params <- function(i, prep, n_conf, n_margins,
                                 margin_suffixes, model_type) {
  dprime <- brms::get_dpar(prep, "dprime", i = i)
  criterion <- brms::get_dpar(prep, "criterion", i = i)

  margins_l <- margins_h <- NULL
  if (n_margins > 0) {
    margins_l <- matrix(NA_real_, nrow = length(dprime), ncol = n_margins)
    margins_h <- matrix(NA_real_, nrow = length(dprime), ncol = n_margins)
    for (j in seq_len(n_margins)) {
      margins_l[, j] <- brms::get_dpar(prep, paste0("crl", margin_suffixes[j]),
                                         i = i)
      margins_h[, j] <- brms::get_dpar(prep, paste0("crh", margin_suffixes[j]),
                                         i = i)
    }
  }

  thresholds <- .sdai_thresholds(criterion, margins_l, margins_h)

  sdratio <- NULL
  if (model_type == "gaussian") {
    sdratio <- brms::get_dpar(prep, "sdratio", i = i)
  }

  nlist(dprime, criterion, thresholds, sdratio)
}

# Internal: get multinomial counts from prep data
.sdai_get_counts <- function(i, prep, n_conf, n_src) {
  n_total <- n_src * n_conf
  counts <- integer(n_total)
  counts[1] <- prep$data$Y[i]
  for (j in seq_len(n_total - 1)) {
    counts[j + 1] <- prep$data[[paste0("vint", j)]][i]
  }
  tabs_counts <- counts[seq_len(n_conf)]
  tpres_counts <- counts[(n_conf + 1):n_total]
  nlist(tabs_counts, tpres_counts)
}


# Closure factory: create log_lik function with captured model constants
.sdai_make_log_lik <- function(n_conf, n_src, n_margins, margin_suffixes,
                               model_type) {
  force(n_conf)
  force(n_src)
  force(n_margins)
  force(margin_suffixes)
  force(model_type)

  function(i, prep) {
    params <- .sdai_extract_params(i, prep, n_conf, n_margins,
                                    margin_suffixes, model_type)
    counts <- .sdai_get_counts(i, prep, n_conf, n_src)

    if (model_type == "gumbel") {
      p_tabs <- .sdai_tabs_probs_gumbel(params$thresholds)
      p_tpres <- .sdai_tpres_probs_gumbel(params$thresholds, params$dprime)
    } else {
      p_tabs <- .sdai_tabs_probs_gaussian(params$thresholds)
      p_tpres <- .sdai_tpres_probs_gaussian(params$thresholds, params$dprime,
                                             params$sdratio)
    }

    ll_tabs <- extraDistr::dmnom(
      x = matrix(counts$tabs_counts, nrow = nrow(p_tabs),
                 ncol = length(counts$tabs_counts), byrow = TRUE),
      size = sum(counts$tabs_counts),
      prob = p_tabs, log = TRUE
    )
    ll_tpres <- extraDistr::dmnom(
      x = matrix(counts$tpres_counts, nrow = nrow(p_tpres),
                 ncol = length(counts$tpres_counts), byrow = TRUE),
      size = sum(counts$tpres_counts),
      prob = p_tpres, log = TRUE
    )

    ll_tabs + ll_tpres
  }
}

# Closure factory: create posterior_predict function with captured constants
.sdai_make_posterior_predict <- function(n_conf, n_src, n_margins,
                                         margin_suffixes, model_type) {
  force(n_conf)
  force(n_src)
  force(n_margins)
  force(margin_suffixes)
  force(model_type)

  function(i, prep, ...) {
    params <- .sdai_extract_params(i, prep, n_conf, n_margins,
                                    margin_suffixes, model_type)
    counts <- .sdai_get_counts(i, prep, n_conf, n_src)

    if (model_type == "gumbel") {
      p_tabs <- .sdai_tabs_probs_gumbel(params$thresholds)
      p_tpres <- .sdai_tpres_probs_gumbel(params$thresholds, params$dprime)
    } else {
      p_tabs <- .sdai_tabs_probs_gaussian(params$thresholds)
      p_tpres <- .sdai_tpres_probs_gaussian(params$thresholds, params$dprime,
                                             params$sdratio)
    }

    n_draws <- nrow(p_tabs)
    out <- matrix(NA_real_, nrow = n_draws, ncol = n_src * n_conf)

    for (d in seq_len(n_draws)) {
      tabs_draw <- stats::rmultinom(1, sum(counts$tabs_counts), p_tabs[d, ])
      tpres_draw <- stats::rmultinom(1, sum(counts$tpres_counts), p_tpres[d, ])
      out[d, ] <- c(tabs_draw, tpres_draw)
    }
    out
  }
}
