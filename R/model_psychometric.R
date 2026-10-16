############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_psychometric <- function(response = NULL, intensity = NULL,
                                n_trials = NULL, sigmoid = "normal",
                                guess = NULL, lapse = NULL, call = NULL, ...) {
  log_x <- .psychometric_sigmoids[[sigmoid]]$log_x
  parameters <- list(
    midpoint = paste0(
      "Stimulus intensity at which the sigmoid F is 0.5, halfway between ",
      "the asymptotes, in the units of the intensity",
      if (log_x) "; log link, so effects are multiplicative" else ""
    ),
    width = paste0(
      "Distance on the intensity axis from F = 0.05 to F = 0.95",
      if (log_x) ", in log units of the intensity" else "",
      "; log link"
    ),
    guess = paste0(
      "Guess rate: the lower asymptote, the response probability far below ",
      "the midpoint; logit link"
    ),
    lapse = paste0(
      "Lapse rate: the proportion of trials on which the observer does not ",
      "attend and responds at the guess rate; the upper asymptote is ",
      "1 - lapse * (1 - guess); logit link"
    )
  )
  links <- list(
    midpoint = if (log_x) "log" else "identity",
    width = "log", guess = "logit", lapse = "logit"
  )
  fixed_rates <- c(guess = !is.null(guess), lapse = !is.null(lapse))
  parameters <- parameters[setdiff(names(parameters), names(fixed_rates)[fixed_rates])]
  links <- links[names(parameters)]

  # guess and lapse both take a logit-normal fitted by least squares to the
  # logit quantiles (2.5% to 97.5%) of psignifit 4's Beta(1, 10) default,
  # (-2.97, 1.31): median 0.047, 95% below 0.39. lapse here is the share of
  # lapse trials, so the same prior on the Wichmann-Hill lapse, lapse * (1 -
  # guess), is half as wide in 2AFC. midpoint and width are placeholders for an
  # intensity range of 0 to 1; check_model() rescales them to the data.
  default_priors <- c(
    .psychometric_shape_priors(c(0, 1)),
    list(
      guess = list(main = "normal(-3, 1.3)", effects = "normal(0, 1)", sd = "exponential(1)"),
      lapse = list(main = "normal(-3, 1.3)", effects = "normal(0, 1)", sd = "exponential(1)")
    )
  )[names(parameters)]

  requirements <- glue(
    "Provide data with one row per trial or per stimulus level:", "\n\n",
    "  - Response (response): 0/1 per trial, or the number of positive ",
    "responses out of n_trials", "\n",
    "  - Stimulus intensity (intensity): the level of the manipulated stimulus ",
    "dimension", "\n",
    "  - Number of trials (n_trials, optional): trials per row for counts"
  )

  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(intensity, n_trials, sigmoid, guess, lapse),
      domain = "Psychophysics",
      task = "Detection or discrimination at varying stimulus intensity",
      name = "Psychometric function",
      citation = c(
        glue(
          "Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: ",
          "I. Fitting, sampling, and goodness of fit. Perception & ",
          "Psychophysics, 63(8), 1293-1313. https://doi.org/10.3758/BF03194544"
        ),
        glue(
          "Sch\u00fctt, H. H., Harmeling, S., Macke, J. H., & Wichmann, F. A. ",
          "(2016). Painfree and accurate Bayesian estimation of psychometric ",
          "functions for (potentially) overdispersed data. Vision Research, ",
          "122, 105-123. https://doi.org/10.1016/j.visres.2016.02.002"
        )
      ),
      version = "NA",
      requirements = requirements,
      parameters = parameters,
      links = links,
      fixed_parameters = list(),
      default_priors = default_priors,
      init_ranges = c(
        .psychometric_shape_inits(c(0, 1), log_x),
        list(guess = c(0.02, 0.1), lapse = c(0.01, 0.05))
      )[names(parameters)]
    ),
    class = c("bmmodel", "psychometric"),
    call = call
  )
  set_links(out, NULL)
}

# The midpoint and width priors are written on the links the model fixes, so
# they would no longer fit the scale they describe under another link, and the
# log-scale sigmoids need a positive midpoint
#' @exportS3Method
settable_links.psychometric <- function(model) {
  character(0)
}


#' @title Psychometric Function Model
#' @name psychometric
#' @description Estimates the psychometric function, the probability of a
#'   positive response as a function of stimulus intensity, with its midpoint,
#'   width, guess rate and lapse rate, from trial-level binary responses or
#'   from response counts per stimulus level.
#' @details `r model_docs(.model_psychometric())`
#'
#' The probability of a positive response ("yes", "correct", "longer") at
#' intensity \eqn{x} is
#' \deqn{\psi(x) = \gamma + (1 - \gamma)(1 - \lambda) F(x),}{psi(x) = guess + (1 - guess) * (1 - lapse) * F(x),}
#' where \eqn{\gamma} is `guess`, \eqn{\lambda} is `lapse` and \eqn{F} the
#' sigmoid. A lapse is a trial on which the observer does not attend and so
#' responds at the guess rate, which keeps \eqn{\psi} increasing for every
#' value of the two rates.
#'
#' @section Midpoint and width:
#' The sigmoid is placed by `midpoint`, the intensity at which \eqn{F = 0.5},
#' and scaled by `width`, the distance on the intensity axis between
#' \eqn{F = 0.05} and \eqn{F = 0.95} (Schütt et al., 2016). Both mean the same
#' thing for every sigmoid, so estimates are comparable across sigmoids.
#'
#' `midpoint` is not the threshold at a target performance level: because of
#' the guess and lapse rates, \eqn{F = 0.5} is \eqn{\psi = 0.75} in 2AFC
#' without lapses, but \eqn{\psi = 0.5} in a yes/no task. Use
#' [psychometric_threshold()] for the intensity at any performance level, with
#' its posterior, and for the slope there.
#'
#' @section Sigmoids:
#' \itemize{
#'   \item `"normal"`: cumulative normal (probit), on the intensity axis
#'   \item `"logistic"`: logistic, on the intensity axis
#'   \item `"gumbel_min"`: \eqn{1 - \exp(-\exp(z))}{1 - exp(-exp(z))}, skewed
#'     towards low intensities (psignifit's `"gumbel"`)
#'   \item `"gumbel_max"`: \eqn{\exp(-\exp(-z))}{exp(-exp(-z))}, skewed towards
#'     high intensities (psignifit's `"rgumbel"`)
#'   \item `"weibull"`: `"gumbel_min"` on log intensity, the Weibull function
#'   \item `"lognormal"`, `"loglogistic"`: `"normal"` and `"logistic"` on log
#'     intensity
#' }
#' The last three read the intensity on the log scale, so it must be positive.
#' Their `midpoint` is still reported in the units of the intensity, through a
#' log link, while their `width` is in log units.
#'
#' @section Guess and lapse rates:
#' Each rate is either estimated, with a logit link and its own formula, or
#' fixed by the design. Fix them in the call, on the probability scale, not in
#' the formula: `guess = 0.5` for 2AFC, `guess = 1 / m` for m-AFC,
#' `guess = "chance"` to read a per-row guess rate from a column (for designs
#' that mix 2AFC and 4AFC trials), `lapse = 0` for no lapses. A fixed value of
#' zero is allowed. Leave `guess` estimated for yes/no tasks, where it is the
#' false-alarm rate far below the midpoint.
#'
#' The lapse rate of Wichmann and Hill (2001), and of psignifit and Palamedes,
#' is the gap between the upper asymptote and 1, which is
#' `lapse * (1 - guess)` here: in 2AFC, half the `lapse` this model reports.
#' [psychometric_threshold()] returns both.
#'
#' @section Default priors:
#' Stimulus units are arbitrary, so the priors on `midpoint` and `width` are
#' scaled to the range of the intensities in the data (on the log scale for the
#' log-scale sigmoids), as psignifit does. With \eqn{c} the centre and \eqn{R}
#' the range, the `midpoint` intercept gets `normal(c, R / 2)`, its effects
#' `normal(0, R / 4)` and its random-effect SDs `exponential(4 / R)`; the log of
#' `width` gets `normal(log(R / 2), 1)`. The priors listed above are these for a
#' range of 0 to 1; [default_prior()] shows the ones for your data.
#'
#' @param response The name of the variable in the data with the responses:
#'   0/1 (or logical) per trial when `n_trials` is `NULL`, and otherwise the
#'   number of positive responses out of `n_trials`.
#' @param intensity The name of the variable in the data with the stimulus
#'   intensity of each row.
#' @param n_trials The name of the variable in the data with the number of
#'   trials per row, for data aggregated into counts. `NULL` (default) reads
#'   `response` as one binary response per row.
#' @param sigmoid The sigmoid \eqn{F}: one of `"normal"` (default),
#'   `"logistic"`, `"gumbel_min"`, `"gumbel_max"`, `"weibull"`, `"lognormal"`
#'   or `"loglogistic"`. See the Sigmoids section.
#' @param guess The guess rate. `NULL` (default) estimates it; a number in
#'   \[0, 1) fixes it; the name of a variable in the data fixes it per row.
#' @param lapse The lapse rate. `NULL` (default) estimates it; a number in
#'   \[0, 1) fixes it; the name of a variable in the data fixes it per row.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#'
#' @references
#' Schütt, H. H., Harmeling, S., Macke, J. H., & Wichmann, F. A. (2016).
#'   Painfree and accurate Bayesian estimation of psychometric functions for
#'   (potentially) overdispersed data. \emph{Vision Research}, \emph{122},
#'   105--123. \doi{10.1016/j.visres.2016.02.002}
#'
#' Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: I.
#'   Fitting, sampling, and goodness of fit. \emph{Perception & Psychophysics},
#'   \emph{63}(8), 1293--1313. \doi{10.3758/BF03194544}
#'
#' @seealso [dpsychometric()], [rpsychometric()], [psychometric_threshold()],
#'   and [sdt_mafc()] for m-AFC accuracy at a single stimulus level.
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # 2AFC contrast detection, 40 trials per level, 10 observers
#' dat <- expand.grid(contrast = c(0.01, 0.02, 0.04, 0.08, 0.16), id = 1:10)
#' dat$n_trials <- 40L
#' dat$n_correct <- rpsychometric(nrow(dat), dat$contrast, midpoint = 0.04,
#'                                width = 2, guess = 0.5, lapse = 0.02,
#'                                n_trials = dat$n_trials, sigmoid = "weibull")
#'
#' model <- psychometric(response = "n_correct", intensity = "contrast",
#'                       n_trials = "n_trials", sigmoid = "weibull",
#'                       guess = 0.5)
#'
#' fit <- bmm(
#'   formula = bmf(midpoint ~ 1 + (1 | id), width ~ 1, lapse ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # contrast at 75% correct, and the slope there
#' psychometric_threshold(fit, p = 0.75)
#' }
psychometric <- function(response, intensity, n_trials = NULL,
                         sigmoid = c("normal", "logistic", "gumbel_min",
                                     "gumbel_max", "weibull", "lognormal",
                                     "loglogistic"),
                         guess = NULL, lapse = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  sigmoid <- match.arg(sigmoid)
  rates <- list(guess = guess, lapse = lapse)
  for (rate in names(rates)) {
    value <- rates[[rate]]
    stopif(
      !is.null(value) &&
        !(is.character(value) && length(value) == 1) &&
        !(is.numeric(value) && length(value) == 1 && is.finite(value) &&
            value >= 0 && value < 1),
      "{rate} must be NULL (estimated), a single number in [0, 1), or the \\
      name of a column in the data"
    )
  }

  .model_psychometric(response = response, intensity = intensity,
                      n_trials = n_trials, sigmoid = sigmoid,
                      guess = guess, lapse = lapse, call = call, ...)
}


############################################################################# !
# DATA-SCALED PRIORS AND INITS                                           ####
############################################################################# !

# psignifit 4 puts its threshold prior flat on the stimulus range with fall-offs
# to half a range beyond either end, and its width prior flat up to the range
# with a fall-off to three ranges. normal(c, R / 2) has its 95% interval on
# exactly that threshold support; normal(log(R / 2), 1) on log width has its
# 97.5% quantile at 3.55 R and its 2.5% quantile at 0.07 R. A condition effect
# of up to half a range and a between-subject SD of R / 4 on average are the
# scales the other two midpoint priors express.
.psychometric_shape_priors <- function(x_range) {
  centre <- signif(mean(x_range), 3)
  span <- diff(x_range)
  list(
    midpoint = list(
      main = glue("normal({centre}, {signif(span / 2, 3)})"),
      effects = glue("normal(0, {signif(span / 4, 3)})"),
      sd = glue("exponential({signif(4 / span, 3)})")
    ),
    width = list(
      main = glue("normal({signif(log(span / 2), 3)}, 1)"),
      effects = "normal(0, 0.5)",
      sd = "exponential(2)"
    )
  )
}

# on the natural scale, which for the log-scale sigmoids' midpoint is the
# intensity itself
.psychometric_shape_inits <- function(x_range, log_x) {
  span <- diff(x_range)
  midpoint <- mean(x_range) + c(-0.25, 0.25) * span
  list(
    midpoint = if (log_x) exp(midpoint) else midpoint,
    width = c(0.25, 1) * span
  )
}

# Range of the intensities on the sigmoid's axis; NULL when the column cannot
# give one, which check_data() then reports. A single level has no range, so
# the priors are centred on it with a span of its own magnitude (or 1 at 0),
# rather than left at the 0-to-1 placeholders whatever the units
.psychometric_axis_range <- function(model, data) {
  x <- data[[model$other_vars$intensity]]
  if (!is.numeric(x)) {
    return(NULL)
  }
  if (.psychometric_sigmoids[[model$other_vars$sigmoid]]$log_x) {
    x <- log(x[x > 0])
  }
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NULL)
  }
  if (length(unique(x)) == 1) {
    return(x[1] + c(-0.5, 0.5) * max(abs(x[1]), 1))
  }
  range(x)
}

#' @export
check_model.psychometric <- function(model, data = NULL, formula = NULL) {
  x_range <- .psychometric_axis_range(model, data)
  if (!is.null(x_range)) {
    log_x <- .psychometric_sigmoids[[model$other_vars$sigmoid]]$log_x
    model$default_priors[c("midpoint", "width")] <- .psychometric_shape_priors(x_range)
    model$init_ranges[c("midpoint", "width")] <- .psychometric_shape_inits(x_range, log_x)
  }
  NextMethod("check_model")
}


############################################################################# !
# CHECK_FORMULA S3 METHODS                                               ####
############################################################################# !

# The likelihood reads the intensity, so a midpoint ~ intensity term rescales
# the intensity axis exactly as width does, and the two trade off along a
# ridge; on width it bends the sigmoid into a curve that is no longer one
#' @export
check_formula.psychometric <- function(model, data, formula) {
  intensity <- model$other_vars$intensity
  uses_intensity <- vapply(rhs_vars(formula, collapse = FALSE),
                           function(x) intensity %in% x, logical(1))
  uses_intensity <- uses_intensity[names(uses_intensity) %in% c("midpoint", "width")]
  stopif(
    any(uses_intensity),
    "The formula for {collapse_comma(names(uses_intensity)[uses_intensity])} \\
    uses the intensity variable '{intensity}'. {model$name} already places \\
    every row on the sigmoid by '{intensity}', so a term in it rescales the \\
    intensity axis, which is what width does, and the model is not \\
    identified. Drop '{intensity}' from the formula."
  )
  NextMethod("check_formula")
}


############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.psychometric <- function(model, data, formula) {
  ov <- model$other_vars
  resp_var <- model$resp_vars$response

  reserved <- intersect(.psychometric_reserved_cols, colnames(data))
  warnif(length(reserved) > 0,
         "Column(s) {collapse_comma(reserved)} in your data are reserved by \\
         {model$name} and will be overwritten")

  if (is.null(ov$n_trials)) {
    data[[resp_var]] <- .psychometric_binary_response(data, resp_var)
    data$psy_trials <- 1L
  } else {
    .validate_sdt_counts(data, resp_var, ov$n_trials)
  }

  .psychometric_check_intensity(data, ov$intensity, ov$sigmoid)

  for (rate in c("guess", "lapse")) {
    if (is.character(ov[[rate]])) {
      .psychometric_check_rate_column(data, ov[[rate]], rate)
    } else {
      # an estimated rate never reads this column, which only fills its slot
      data[[paste0("psy_", rate)]] <- ov[[rate]] %||% 0
    }
  }

  spec <- .psychometric_sigmoids[[ov$sigmoid]]
  data$psy_dist <- .sdt_dist_id(spec$dist)
  data$psy_logx <- as.integer(spec$log_x)

  NextMethod("check_data")
}

.psychometric_reserved_cols <- c("psy_trials", "psy_guess", "psy_lapse",
                                 "psy_dist", "psy_logx")

.psychometric_binary_response <- function(data, resp_var) {
  stopif(!resp_var %in% colnames(data),
         "Response variable '{resp_var}' missing in the data")
  y <- data[[resp_var]]
  if (is.logical(y)) {
    y <- as.integer(y)
  }
  stopif(!is.numeric(y) || any(!y %in% c(0, 1, NA)),
         "Response variable '{resp_var}' must be 0/1 (or logical) when \\
         n_trials is NULL, one binary response per row. For counts of \\
         positive responses, name the column with the number of trials in \\
         n_trials")
  as.integer(y)
}

.psychometric_check_intensity <- function(data, intensity, sigmoid) {
  stopif(!intensity %in% colnames(data),
         "Intensity variable '{intensity}' missing in the data")
  x <- data[[intensity]]
  stopif(!is.numeric(x),
         "Intensity variable '{intensity}' must be numeric (it is {class(x)[1]})")
  stopif(any(!is.finite(x)),
         "Intensity variable '{intensity}' must be finite, without NA")
  stopif(.psychometric_sigmoids[[sigmoid]]$log_x && any(x <= 0),
         "Intensity variable '{intensity}' must be positive for sigmoid = \\
         '{sigmoid}', which reads it on the log scale. Use sigmoid = \\
         '{.psychometric_linear_twin(sigmoid)}' for intensities already on a \\
         log scale (e.g. in dB or log units)")
  warnif(length(unique(x)) < 2,
         "Intensity variable '{intensity}' has a single value, so the data \\
         cannot separate midpoint from width; only the priors will, and the \\
         default ones are centred on that value with an arbitrary span. Set \\
         your own priors on midpoint and width")
}

.psychometric_linear_twin <- function(sigmoid) {
  c(weibull = "gumbel_min", lognormal = "normal", loglogistic = "logistic")[[sigmoid]]
}

.psychometric_check_rate_column <- function(data, col, rate) {
  stopif(!col %in% colnames(data),
         "The {rate} rate column '{col}' is missing in the data")
  vals <- data[[col]]
  stopif(!is.numeric(vals) || anyNA(vals) || any(vals < 0 | vals >= 1),
         "The {rate} rate column '{col}' must be numeric, in [0, 1), \\
         without NA; for m-AFC trials it holds 1 / m")
}


############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.psychometric <- function(model, formula) {
  ov <- model$other_vars
  rate_col <- function(rate) {
    if (is.character(ov[[rate]])) ov[[rate]] else paste0("psy_", rate)
  }
  brms::bf(glue(
    "{model$resp_vars$response} | trials({ov$n_trials %||% 'psy_trials'}) + ",
    "vreal({ov$intensity}, {rate_col('guess')}, {rate_col('lapse')}) + ",
    "vint(psy_dist, psy_logx) ~ 0"
  ))
}


############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.psychometric <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  rates <- intersect(c("guess", "lapse"), names(model$parameters))
  dpars <- c("midpoint", "width", rates)
  formula$family <- brms::custom_family(
    .psychometric_family_name(rates),
    dpars = c("mu", dpars),
    links = c("identity", unlist(model$links[dpars], use.names = FALSE)),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_psychometric,
    posterior_predict = posterior_predict_psychometric,
    posterior_epred = posterior_epred_psychometric,
    vars = c("trials[n]", "vreal1[n]", "vreal2[n]", "vreal3[n]",
             "vint1[n]", "vint2[n]", "psy_zmid", "psy_zspan")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- paste(
    read_lines2(paste0(sc_path, "/sdt_dist_funs.stan")),
    read_lines2(paste0(sc_path, "/psychometric_funs.stan")),
    sep = "\n"
  )
  z <- .psychometric_z_constants(model$other_vars$sigmoid)
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(x = z[["zmid"]], name = "psy_zmid", scode = "real psy_zmid;") +
    brms::stanvar(x = z[["zspan"]], name = "psy_zspan", scode = "real psy_zspan;")

  nlist(formula, data, stanvars)
}

# psychometric_funs.stan has one family for each combination of estimated
# rates; a fixed rate reaches it as data, in vreal2 (guess) or vreal3 (lapse)
.psychometric_family_name <- function(rates) {
  suffix <- switch(paste(rates, collapse = "+"),
    "guess+lapse" = "",
    "guess" = "_fixlapse",
    "lapse" = "_fixguess",
    "_fixboth"
  )
  paste0("psychometric", suffix)
}


############################################################################# !
# LOG_LIK & POSTERIOR_PREDICT                                            ####
############################################################################# !

.psychometric_sigmoid_name <- function(dist_id, log_x) {
  hit <- vapply(.psychometric_sigmoids, function(s) {
    .sdt_dist_id(s$dist) == dist_id && s$log_x == as.logical(log_x)
  }, logical(1))
  names(.psychometric_sigmoids)[hit]
}

# A rate is a dpar when estimated, and otherwise the data column it was fixed to
.psychometric_rate <- function(prep, rate, i = NULL) {
  if (rate %in% names(prep$dpars)) {
    return(brms::get_dpar(prep, rate, i = i))
  }
  col <- prep$data[[c(guess = "vreal2", lapse = "vreal3")[[rate]]]]
  if (is.null(i)) col else col[i]
}

.psychometric_row_log_probs <- function(i, prep) {
  .psychometric_log_probs(
    intensity = prep$data$vreal1[i],
    midpoint = brms::get_dpar(prep, "midpoint", i = i),
    width = brms::get_dpar(prep, "width", i = i),
    guess = .psychometric_rate(prep, "guess", i),
    lapse = .psychometric_rate(prep, "lapse", i),
    sigmoid = .psychometric_sigmoid_name(prep$data$vint1[i], prep$data$vint2[i])
  )
}

log_lik_psychometric <- function(i, prep) {
  lp <- .psychometric_row_log_probs(i, prep)
  y <- prep$data$Y[i]
  n <- prep$data$trials[i]
  lchoose(n, y) + times_nonzero(y, lp$one) + times_nonzero(n - y, lp$zero)
}

posterior_predict_psychometric <- function(i, prep, ...) {
  lp <- .psychometric_row_log_probs(i, prep)
  stats::rbinom(prep$ndraws, prep$data$trials[i], exp(lp$one))
}

posterior_epred_psychometric <- function(prep) {
  rate <- function(name) {
    if (name %in% names(prep$dpars)) {
      prep$dpars[[name]]
    } else {
      .epred_data(.psychometric_rate(prep, name), prep)
    }
  }
  lp <- .psychometric_log_probs(
    intensity = .epred_data(prep$data$vreal1, prep),
    midpoint = prep$dpars$midpoint,
    width = prep$dpars$width,
    guess = rate("guess"),
    lapse = rate("lapse"),
    sigmoid = .psychometric_sigmoid_name(prep$data$vint1[1], prep$data$vint2[1])
  )
  .epred_matrix(.epred_data(prep$data$trials, prep) * exp(lp$one), prep)
}


############################################################################# !
# POST-PROCESSING                                                        ####
############################################################################# !

#' Thresholds and slopes of a fitted psychometric function
#'
#' Computes, per posterior draw and condition, the stimulus intensity at which
#' a [psychometric()] fit reaches a target response probability, the slope of
#' the psychometric function there, and the guess and lapse rates, including the
#' lapse rate as Wichmann and Hill (2001) define it.
#'
#' The threshold at target probability \eqn{p} solves
#' \eqn{\gamma + (1 - \gamma)(1 - \lambda) F(x) = p}{guess + (1 - guess) * (1 - lapse) * F(x) = p}.
#' It exists only between the lower asymptote, `guess`, and the upper asymptote,
#' `1 - lapse * (1 - guess)`. A draw whose asymptotes do not enclose `p` has no
#' threshold: its threshold and slope are `NA`, the summary reports the share of
#' such draws as `undefined`, and the function warns. Asking for 0.75 in a 2AFC
#' task with `guess = 0.5` is the usual case; asking for 0.5 is not, because no
#' 2AFC observer falls below chance.
#'
#' The slope is \eqn{d\psi / dx}{dpsi / dx} at the threshold, in probability
#' per unit of intensity.
#'
#' @param fit A `bmmfit` object returned by [bmm()] with a [psychometric()]
#'   model.
#' @param p Numeric vector of target response probabilities, each strictly
#'   between 0 and 1. The default 0.75 is the conventional 2AFC threshold.
#' @param conditions A data frame of predictor values at which to evaluate the
#'   parameters. The default uses every unique combination of the
#'   population-level predictors in the data, and of the guess or lapse column
#'   when a rate is fixed per row. Group-level effects are left out.
#' @param probs Numeric vector of length 2. Lower and upper quantiles for the
#'   credible interval in the summary. Default `c(0.025, 0.975)`.
#' @param ... Further arguments passed to [brms::posterior_linpred()], such as
#'   `draw_ids`.
#'
#' @return A data frame of class `"bmm_psychometric_threshold"` with one row
#'   per draw, measure and condition, and columns `measure` (`"threshold"`,
#'   `"slope"`, `"guess"`, `"lapse"`, `"lapse_wh"`), `p` (`NA` for the three
#'   rates, which do not depend on it), `value`, `.draw` and the condition
#'   columns. `lapse_wh` is `lapse * (1 - guess)`, the gap between the upper
#'   asymptote and 1, which psignifit and Palamedes report as the lapse rate.
#'   Attribute `summary` holds the posterior mean, the credible interval and
#'   the share of `undefined` draws per measure, `p` and condition.
#'
#' @references
#' Wichmann, F. A., & Hill, N. J. (2001). The psychometric function: I.
#'   Fitting, sampling, and goodness of fit. \emph{Perception & Psychophysics},
#'   \emph{63}(8), 1293--1313. \doi{10.3758/BF03194544}
#'
#' @seealso [psychometric()]
#' @export
#' @examples
#' \dontrun{
#' dat <- expand.grid(contrast = c(0.01, 0.02, 0.04, 0.08, 0.16),
#'                    cond = c("a", "b"))
#' dat$n_trials <- 60L
#' dat$n_correct <- rpsychometric(nrow(dat), dat$contrast,
#'                                midpoint = ifelse(dat$cond == "a", 0.03, 0.05),
#'                                width = 2, guess = 0.5, lapse = 0.02,
#'                                n_trials = 60, sigmoid = "weibull")
#' fit <- bmm(
#'   bmf(midpoint ~ 0 + cond, width ~ 1, lapse ~ 1), dat,
#'   psychometric("n_correct", "contrast", "n_trials", sigmoid = "weibull",
#'                guess = 0.5),
#'   cores = 4, backend = "cmdstanr"
#' )
#' thr <- psychometric_threshold(fit, p = c(0.75, 0.9))
#' thr
#'
#' # the threshold difference between the conditions, per draw
#' wide <- subset(thr, measure == "threshold" & p == 0.75)
#' diff_ab <- wide$value[wide$cond == "b"] - wide$value[wide$cond == "a"]
#' quantile(diff_ab, c(0.025, 0.5, 0.975))
#' }
psychometric_threshold <- function(fit, p = 0.75, conditions = NULL,
                                   probs = c(0.025, 0.975), ...) {
  .sdt_check_args(fit, conditions, probs, ...)
  model <- fit$bmm$model
  stopif(!inherits(model, "psychometric"),
         "psychometric_threshold() is only available for psychometric() fits")
  stopif(!is.numeric(p) || length(p) == 0 || anyNA(p) || any(p <= 0 | p >= 1),
         "p must be target probabilities strictly between 0 and 1")

  conditions <- .psychometric_resolve_conditions(fit, conditions)
  midpoint <- .psychometric_par_draws(fit, "midpoint", conditions, ...)
  n_draws <- nrow(midpoint)
  pars <- c(
    list(midpoint = midpoint),
    lapply(stats::setNames(nm = c("width", "guess", "lapse")), function(par) {
      .psychometric_par_draws(fit, par, conditions, n_draws, ...)
    })
  )
  sigmoid <- model$other_vars$sigmoid

  draws_list <- list()
  for (c_i in seq_len(max(1L, nrow(conditions)))) {
    cond_row <- conditions[c_i, , drop = FALSE]
    par <- lapply(pars, function(m) m[, c_i])
    rates <- list(guess = par$guess, lapse = par$lapse,
                  lapse_wh = par$lapse * (1 - par$guess))
    measures <- c(
      unlist(lapply(p, function(target) {
        at <- .psychometric_threshold_at(target, par, sigmoid)
        list(data.frame(measure = "threshold", p = target, value = at$threshold),
             data.frame(measure = "slope", p = target, value = at$slope))
      }), recursive = FALSE),
      lapply(names(rates), function(nm) {
        data.frame(measure = nm, p = NA_real_, value = rates[[nm]])
      })
    )
    for (df in measures) {
      df$.draw <- seq_len(n_draws)
      draws_list[[length(draws_list) + 1L]] <- .sdt_bind_cond(df, cond_row)
    }
  }
  draws <- do.call(rbind, draws_list)

  undefined <- mean(is.na(draws$value[draws$measure == "threshold"]))
  warnif(undefined > 0,
         "{signif(100 * undefined, 2)}% of the threshold draws are undefined, \\
         because p lies outside the asymptotes of their psychometric function \\
         (see the 'undefined' column of the summary)")

  structure(
    draws,
    class = c("bmm_psychometric_threshold", "data.frame"),
    summary = .psychometric_summarise(draws, names(conditions), probs),
    probs = probs,
    sigmoid = sigmoid,
    conditions = conditions
  )
}

#' @export
print.bmm_psychometric_threshold <- function(x, ...) {
  if (is.null(attr(x, "summary"))) return(NextMethod())
  cat("Psychometric function thresholds (sigmoid = ", attr(x, "sigmoid"), ")\n",
      sep = "")
  cat("  threshold = intensity at which P(positive response) = p",
      "| slope = change in P per unit of intensity there\n")
  cat("  lapse_wh = Wichmann-Hill lapse rate, lapse * (1 - guess)\n")
  print(attr(x, "summary"), digits = 3, row.names = FALSE)
  invisible(x)
}

# Intensity at which psi reaches p, and the slope of psi there, per draw; NA
# where the asymptotes do not enclose p
.psychometric_threshold_at <- function(p, par, sigmoid) {
  gain <- (1 - par$guess) * (1 - par$lapse)
  f <- (p - par$guess) / gain
  f[f <= 0 | f >= 1] <- NA_real_
  spec <- .psychometric_sigmoids[[sigmoid]]
  dist <- .sdt_dists[[spec$dist]]
  threshold <- .psychometric_x_at(f, par$midpoint, par$width, sigmoid)
  dz_dx <- .psychometric_z_constants(sigmoid)[["zspan"]] / par$width
  if (spec$log_x) dz_dx <- dz_dx / threshold
  list(threshold = threshold, slope = gain * dist$pdf(dist$qf(f)) * dz_dx)
}

.psychometric_summarise <- function(draws, cond_cols, probs) {
  key_cols <- c("measure", "p", cond_cols)
  key <- do.call(paste, c(lapply(draws[key_cols], as.character), sep = "\r"))
  first <- !duplicated(key)
  values <- split(draws$value, factor(key, levels = key[first]))
  keys <- draws[first, key_cols, drop = FALSE]
  rows <- lapply(seq_len(nrow(keys)), function(k) {
    v <- values[[k]]
    data.frame(
      keys[k, , drop = FALSE],
      mean = mean(v, na.rm = TRUE),
      lower = unname(stats::quantile(v, probs[1L], na.rm = TRUE)),
      upper = unname(stats::quantile(v, probs[2L], na.rm = TRUE)),
      undefined = mean(is.na(v))
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# The default conditions are every combination of the population-level
# predictors and of a per-row guess or lapse column, which moves the threshold
# just as a predictor does
.psychometric_resolve_conditions <- function(fit, conditions) {
  model <- fit$bmm$model
  rate_cols <- Filter(is.character, model$other_vars[c("guess", "lapse")])
  cols <- union(unlist(.population_vars(fit$bmm$user_formula)), unlist(rate_cols))
  cols <- setdiff(intersect(cols, names(fit$data)),
                  c(model$other_vars$intensity, "Intercept"))
  if (is.null(conditions)) {
    return(.sdt_unique_subset(fit$data, cols))
  }
  # a column left out would be read from the first data row without notice
  missing_cols <- setdiff(cols, names(conditions))
  stopif(length(missing_cols) > 0,
         "conditions must give a value for every column the model reads; \\
         missing: {collapse_comma(missing_cols)}")
  as.data.frame(conditions)
}

# Natural-scale draws x conditions of one parameter: an estimated one through
# its inverse link, a fixed rate as its value in each condition, repeated over
# n_draws rows
.psychometric_par_draws <- function(fit, par, conditions, n_draws = NULL, ...) {
  model <- fit$bmm$model
  newdata <- fit$data[rep(1L, max(1L, nrow(conditions))), , drop = FALSE]
  rownames(newdata) <- NULL
  cols <- intersect(names(conditions), names(newdata))
  newdata[cols] <- conditions[cols]

  if (par %in% names(model$parameters)) {
    linpred <- brms::posterior_linpred(fit, dpar = par, newdata = newdata,
                                       re_formula = NA, allow_new_levels = TRUE, ...)
    return(link_transform(linpred, model$links[[par]], inverse = TRUE))
  }
  fixed <- model$other_vars[[par]]
  values <- if (is.character(fixed)) newdata[[fixed]] else rep(fixed, nrow(newdata))
  matrix(values, nrow = n_draws, ncol = length(values), byrow = TRUE)
}

# midpoint is the 50% point of the sigmoid, which is a threshold only in
# yes/no tasks, and the lapse rate psignifit reports is lapse * (1 - guess)
#' @export
summary_notes.psychometric <- function(model, x) {
  paste(
    "Note: midpoint is where the sigmoid F is 0.5, not a performance threshold,",
    "and lapse is the share of\n      lapse trials; the Wichmann-Hill lapse",
    "rate is lapse * (1 - guess). psychometric_threshold()\n      returns",
    "thresholds and slopes at any target probability, and both lapse rates."
  )
}
