############################################################################# !
# ORDERED-THRESHOLD PARAMETERIZATIONS                                     ####
############################################################################# !

# Shared by every model that places K - 1 ordered thresholds on a latent
# axis (sdt_rating(), sdt_cdp(), likert_cumulative()). The builders are
# .make_thresholds() in R/distributions.R and make_thresholds() in
# inst/stan_chunks/threshold_functions.stan.

# Threshold-parameterization registry: the position defines the integer
# thresh_type code passed to Stan -- reordering entries changes the R <-> Stan
# contract (see make_thresholds in threshold_functions.stan)
.threshold_types <- c("parsimonious", "equidistant", "log_distance",
                      "log_ratio", "softmax")

# Parameter spec for one threshold parameterization. parsimonious/equidistant
# need a single spacing; log_distance/log_ratio need K-2 interval deltas, named
# in interval order (delta_i lies between thresholds i and i + 1); softmax needs
# spacing plus K-3 allocation deltas.
# Every threshold parameter is a log-scale quantity -- the Stan builders read
# exp(spacing) and exp(delta) -- so its sd default follows sdratio's rate 2
# rather than d's rate 1: at rate 2 an individual's spacing stays within a
# factor of about 2.8 of the group value at the one-sided 95% quantile (90%
# central), at rate 1 within a factor of 8.
.threshold_parameter_parts <- function(n_ratings, threshold_type) {
  parameters <- list()
  default_priors <- list()
  param_links <- list()
  delta_prior <- list(
    main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)"
  )

  if (threshold_type %in% c("equidistant", "parsimonious")) {
    parameters$spacing <- paste0(
      "Threshold spacing: controls distance between adjacent thresholds ",
      "(exp(spacing) ensures positive spacing)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"
    )
    param_links$spacing <- "identity"
  } else if (threshold_type == "log_distance") {
    for (i in seq_len(n_ratings - 2L)) {
      pname <- paste0("delta", i)
      parameters[[pname]] <- glue(
        "Log width of the interval between thresholds {i} and {i + 1}"
      )
      default_priors[[pname]] <- delta_prior
      param_links[[pname]] <- "identity"
    }
  } else if (threshold_type == "log_ratio") {
    labels <- .log_ratio_delta_labels(n_ratings)
    for (i in seq_along(labels)) {
      pname <- paste0("delta", i)
      parameters[[pname]] <- labels[i]
      default_priors[[pname]] <- delta_prior
      param_links[[pname]] <- "identity"
    }
  } else if (threshold_type == "softmax") {
    parameters$spacing <- paste0(
      "Average log spacing across threshold intervals ",
      "(exp(spacing) is the mean interval size)"
    )
    default_priors$spacing <- list(
      main = "normal(0, 0.5)", effects = "normal(0, 0.3)", sd = "exponential(2)"
    )
    param_links$spacing <- "identity"

    n_deltas <- max(0L, n_ratings - 3L)
    for (i in seq_len(n_deltas)) {
      pname <- paste0("delta", i)
      parameters[[pname]] <- glue(
        "Softmax threshold allocation parameter for interval {i}, relative ",
        "to the last interval"
      )
      default_priors[[pname]] <- delta_prior
      param_links[[pname]] <- "identity"
    }
  }

  nlist(parameters, default_priors, param_links)
}

# The log_ratio deltas are not exchangeable (see .log_ratio_widths()), so
# each label names the role its interval plays.
.log_ratio_delta_labels <- function(n_ratings) {
  n_gaps <- n_ratings - 2L
  labels <- character(n_gaps)
  if (n_ratings %% 2L == 0L) {
    m <- n_ratings %/% 2L
    labels[m] <- "Log spread: width of the interval just above the criterion"
    labels[m - 1L] <- paste0(
      "Log ratio of the interval just below the criterion to the spread above it"
    )
    first <- c(rep(m - 1L, max(0L, m - 2L)), rep(m, max(0L, n_gaps - m)))
    rest <- setdiff(seq_len(n_gaps), c(m - 1L, m))
  } else {
    g <- (n_ratings - 1L) %/% 2L
    labels[g] <- "Log spread: width of the middle category, centred on the criterion"
    labels[c(g - 1L, g + 1L)] <- glue(
      "Log ratio of interval {c(g - 1L, g + 1L)} to the middle category"
    )
    rest <- setdiff(seq_len(n_gaps), c(g - 1L, g, g + 1L))
    first <- ifelse(rest < g, g - 1L, g + 1L)
  }
  labels[rest] <- glue("Log ratio of interval {rest} to interval {first}")
  labels
}

.threshold_delta_names <- function(object) {
  grep("^delta", names(object$parameters), value = TRUE)
}

# aliases, remove when no branch reads them
.sdt_threshold_types <- .threshold_types
.sdt_threshold_parameter_parts <- .threshold_parameter_parts
.sdt_log_ratio_delta_labels <- .log_ratio_delta_labels
.sdt_threshold_delta_names <- .threshold_delta_names
