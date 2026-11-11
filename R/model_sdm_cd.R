############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The parameters, links, priors and initial ranges are sdm()'s plus criterion,
# so that a change-detection and a continuous-reproduction fit name the same
# memory parameters; a joint fit of both tasks relies on that. There is a
# single precision state, so the knowledge-limited and knowledge-rich
# observers of the circular mixtures coincide and sdm_cd() takes no knowledge
# argument.
.model_sdm_cd <- function(response = NULL, probe = NULL,
                          links = NULL, version = "simple", call = NULL, ...) {
  de <- .model_sdm()
  spec <- .cd_add_criterion(list(
    parameters = de$parameters,
    links = de$links,
    fixed_parameters = de$fixed_parameters,
    priors = de$default_priors,
    init_ranges = de$init_ranges
  ))
  out <- structure(
    list(
      resp_vars = nlist(response, probe),
      other_vars = nlist(),
      domain = "Working memory",
      task = "Detection, recognition and confidence judgments",
      name = "Signal Discrimination Model (SDM) for Change Detection",
      citation = c(
        de$citation,
        .cd_citation
      ),
      version = version,
      # paste0() rather than glue(), which would trim the indent that keeps
      # the requirements in one list
      requirements = paste0(
        "- The response variable should be coded 0 ('same') and 1 ('change')\n",
        "  - The probe should be in radians and relative to the target: the ",
        "feature shown at test minus the feature shown at the probed location, ",
        "which is 0 on a 'same' trial"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = spec$fixed_parameters,
      default_priors = spec$priors,
      init_ranges = spec$init_ranges
    ),
    class = c("bmmodel", "change_detection", "sdm_cd", paste0("sdm_cd_", version)),
    call = call
  )
  set_links(out, links)
}

# As for sdm(), whose links none can set, so that both tasks report the
# memory parameters on the same scales
#' @exportS3Method
settable_links.sdm_cd <- function(model) {
  character(0)
}

# user facing alias

#' @title `r .model_sdm_cd()$name`
#' @description Fits the Signal Discrimination Model of [sdm()] to the "same"
#'   or "change" responses of a single-probe change-detection task.
#' @details `r model_docs(.model_sdm_cd())`
#'
#' @section The observer:
#' On each trial one item of the memory array is probed, and the participant
#' says whether the probe is the same as the item shown at that location or a
#' change. The observer retrieves a feature from the same memory as in
#' [sdm()], the SDM density with memory strength `c` and precision `kappa`
#' around the probed item's feature. It then compares the retrieved feature
#' with the probe (Lin & Oberauer, 2022, Eq. 8, with the SDM density in place
#' of their von Mises): "change" is the answer when the log-likelihood ratio of
#' a change against no change exceeds `criterion`. The "same" responses
#' therefore fall in an arc around the probe whose half-width depends on `c`,
#' `kappa` and `criterion`; unlike the von Mises models, no parameter drops out
#' of the boundary at `criterion = 0`. The probability of a "same" response is
#' the retrieval mass inside the arc. There is one precision state, so the
#' observer has nothing more to know about the trial, and the model has no
#' `knowledge` argument.
#'
#' The memory parameters have the same names, links and default priors as in
#' [sdm()], plus `criterion`. Whether a person gives the same values in both
#' tasks is an empirical question rather than something the model guarantees:
#' the precision estimated from continuous reproduction includes response
#' noise that change detection does not have, for example. Lin and Oberauer
#' (2022) did not fit the SDM to change-detection data, so this model extends
#' their decision rule rather than reproducing one of their fits.
#'
#' @section Decision criterion:
#' `criterion` is fixed to 0 by default, the unbiased observer. Give it a
#' formula, such as `criterion ~ 1`, to estimate it; positive values favour
#' "same" responses. An estimated criterion is kept inside the range of
#' evidence each trial can produce: it equals Lin and Oberauer's (2022)
#' criterion near 0 and is compressed toward the ends of that range, so both
#' responses stay possible on every trial. Without the bound, a criterion past
#' the range makes one response impossible and the likelihood drops off a
#' cliff the sampler cannot cross. Large estimated values are therefore not
#' directly comparable with Lin and Oberauer's unbounded criterion. bmm prints
#' a message saying so when the criterion is estimated.
#'
#' @section Fitting:
#' As for [sdm()], the normalising constant of the SDM density is computed
#' once for every run of trials with the same predictor values, so sorting the
#' data by the predictors speeds up the fit; bmm says so when they are not
#' sorted. Within-chain threading (the `threads` argument) is not supported.
#'
#' @param response The name of the variable in the data containing the
#'   change-detection response, coded 0 ("same") and 1 ("change"). A logical
#'   variable is accepted and coded `TRUE` = "change".
#' @param probe The name of the variable containing the probe in radians,
#'   relative to the target: the feature shown at test minus the feature shown
#'   at the probed location in the memory array, as [wrap()] returns it. It is 0
#'   on "same" trials.
#' @param version Character. The version of the model to use. Currently only
#'   "simple" is supported, as for [sdm()].
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @seealso [sdm()] for continuous reproduction, [sdm_cd_dist] for the
#'   distribution functions.
#' @keywords bmmodel
#' @references Oberauer, K. (2023). Measurement models for visual working
#'   memory - A factorial model comparison. Psychological Review, 130(3),
#'   841-852. \doi{10.1037/rev0000328}
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. Cognitive
#'   Psychology, 133, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # the probe is relative to the target: half of the trials show the target
#' # itself (a probe of 0), the others a random colour
#' design <- data.frame(probe = c(rep(0, 1000), runif(1000, -pi, pi)))
#' model <- sdm_cd(response = "change", probe = "probe")
#' dat <- simulate_data(model, pars = list(c = 1.5, kappa = 2), design = design)
#'
#' fit <- bmm(
#'   formula = bmmformula(c ~ 1, kappa ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#'
#' # criterion is estimated once it has a formula
#' criterion_formula <- bmmformula(c ~ 1, kappa ~ 1, criterion ~ 1)
#' @export
sdm_cd <- function(response, probe, version = "simple", ...) {
  call <- match.call()
  .check_constructor_dots(list(...), known = "links")
  stop_missing_args()
  version <- match.arg(version)
  .model_sdm_cd(
    response = response, probe = probe, version = version,
    call = call, ...
  )
}

############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

# Runs of rows with equal predictor values share one normalising constant, as
# in sdm(); sorted data make the runs as long as they can be
#' @export
check_data.sdm_cd <- function(model, data, formula) {
  data <- order_data_query(model, data, formula)
  attr(data, "sdm_run_metadata") <- sdm_run_metadata(data, formula, model)
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

#' @export
bmf2bf.sdm_cd <- function(model, formula) {
  brms::bf(.circmix_aterm(model$resp_vars$response, vreal = model$resp_vars$probe))
}

# Not the core's family: log Z costs up to 28 panels of 32 nodes, so it is
# computed once per run of predictor values rather than per row, which makes
# a gradient about 4 times cheaper (#136). The vectorised likelihood needs
# loop = FALSE. Threading is refused: brms slices Y and the dpars per thread
# but passes vreal1 and the run metadata whole, and sdm()'s threaded route
# does not port, because log Z sits inside P("same") and the arc's half-width
# rather than in a term of its own.
#' @export
configure_model.sdm_cd <- function(model, data, formula) {
  stopif(
    brms_slices_likelihood(),
    "sdm_cd() does not support within-chain threading yet. Fit it without \\
    the 'threads' argument (and without options(brms.threads))."
  )
  runs <- attr(data, "sdm_run_metadata") %||% sdm_run_metadata(data, formula, model)
  dpars <- c("mu", "c", "kappa", "criterion")
  bounds <- .circmix_bounds(model$links[dpars])
  formula <- bmf2bf(model, formula)
  formula$family <- brms::custom_family(
    paste0("sdm_cd_", model$version),
    dpars = dpars,
    links = unlist(model$links[dpars], use.names = FALSE),
    lb = bounds$lb, ub = bounds$ub,
    type = "int",
    vars = c("vreal1", "G_sdm_runs", "sdm_run_start", "sdm_run_count", "cd_gl_x", "cd_gl_w"),
    loop = FALSE,
    log_lik = log_lik_sdm_cd_simple,
    posterior_predict = posterior_predict_sdm_cd_simple,
    posterior_epred = posterior_epred_sdm_cd_simple
  )

  tdata <- read_lines2(system.file("stan_chunks", "sdm_simple_tdata.stan", package = "bmm"))
  stanvars <- .circmix_chunk_stanvar("circmix_functions.stan") +
    .cd_stanvars() +
    .circmix_chunk_stanvar("sdm_cd_functions.stan") +
    brms::stanvar(scode = tdata, block = "tdata") +
    brms::stanvar(x = runs$G_sdm_runs, name = "G_sdm_runs") +
    sdm_stanvar_int_array(runs$sdm_run_start, "sdm_run_start", "G_sdm_runs") +
    sdm_stanvar_int_array(runs$sdm_run_count, "sdm_run_count", "G_sdm_runs")

  nlist(formula, data, stanvars)
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

# P("same") of observation i for every posterior draw
.sdm_cd_prep_psame <- function(prep, i) {
  .sdm_cd_psame(
    .cd_prep_probe(prep, i),
    mu = brms::get_dpar(prep, "mu", i = i),
    c = brms::get_dpar(prep, "c", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    criterion = brms::get_dpar(prep, "criterion", i = i)
  )
}

log_lik_sdm_cd_simple <- function(i, prep) {
  .cd_bernoulli_ld(prep$data$Y[i], .sdm_cd_prep_psame(prep, i))
}

posterior_predict_sdm_cd_simple <- function(i, prep, ...) {
  p_same <- .sdm_cd_prep_psame(prep, i)
  stats::rbinom(length(p_same), 1, 1 - p_same)
}

# The expected response is the probability of a "change" response. Every dpar
# arrives as a draws x observations matrix (a fixed one as one value per
# draw, which recycles down the columns), so the cells are evaluated as one
# long vector.
posterior_epred_sdm_cd_simple <- function(prep) {
  dpars <- prep$dpars
  p_same <- .sdm_cd_psame(
    .epred_data(prep$data$vreal1, prep),
    mu = as.vector(dpars$mu),
    c = as.vector(dpars$c),
    kappa = as.vector(dpars$kappa),
    criterion = as.vector(dpars$criterion)
  )
  .epred_matrix(1 - p_same, prep)
}

############################################################################# !
# SIMULATE_DATA S3 METHODS                                               ####
############################################################################# !

# Draws through the observer, .sdm_cd_observe(), never through the P("same")
# of the likelihood
#' @export
simulate_data.sdm_cd <- function(model, pars, design = NULL, scale = c("link", "native"), ...) {
  vars <- model$resp_vars
  design <- .design_or_n(model, design, vars$probe, ...)
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  y <- .sdm_cd_observe(
    probe = wrap(design[[vars$probe]]),
    mu = p$mu, c = p$c, kappa = p$kappa, criterion = p$criterion
  )
  .simulate_result(design, stats::setNames(list(y), vars$response), sim)
}

############################################################################# !
# DENSITY HELPERS                                                        ####
############################################################################# !
# R twins of inst/stan_chunks/sdm_cd_functions.stan, vectorised over rows.
#
# The sdm density is far more peaked than a von Mises, so the arc integral
# cannot use the single 32-node rule of the circular mixtures (it errs by 149%
# against an adaptive reference; 128 nodes still by 46%, measured on the
# August change-detection branch). The integral is taken in log space,
# because the activation reaches exp(709) inside the sampled range; it is split
# at every periodic image of the peak, because an arc may span the whole
# circle; and at multiples of the spike width, so the peak sits on a panel
# boundary where Gauss-Legendre clusters its nodes. With all three the rule
# held to 3.8e-10 over c <= 60, kappa <= 150 and arcs from 0.05 rad to the
# circle.
#
# The normalising constant uses the same rule (decided for #136 when sdm()'s
# Chebyshev rule erred by up to 7.5 in log Z): against integrate() it errs by
# at most 1.6e-6 in log Z for log c <= 5 and kappa <= 300, at the corner cells
# (log c 5, kappa 200), and by at most 3.3e-10 for kappa <= 30. It is the
# costliest step of the likelihood, which therefore takes it once per run of
# predictor values (check_data.sdm_cd()).

.sdm_cd_log_peak <- function(c, kappa) {
  log(c) + 0.5 * log(kappa / (2 * pi))
}

.sdm_cd_spike_width <- function(c, kappa) {
  pmin(pi, sqrt(2 * pmax(.sdm_cd_log_peak(c, kappa), 1) / kappa))
}

# log of the integral of exp(c sqrt(kappa / 2 pi) exp(kappa (cos(u - mu) - 1)))
# over [lo, hi]
.sdm_cd_log_int <- function(lo, hi, mu, c, kappa, gl = .cd_gl_rule()) {
  args <- .circmix_recycle(lo = lo, hi = hi, mu = mu, c = c, kappa = kappa)
  n <- length(args$lo)
  w <- .sdm_cd_spike_width(args$c, args$kappa)
  peak <- exp(.sdm_cd_log_peak(args$c, args$kappa))
  mult <- c(-3, -2, -1, -0.5, 0, 0.5, 1, 2, 3)
  images <- outer(args$mu, rep(c(-2, 0, 2) * pi, each = length(mult)), "+") +
    outer(w, rep(mult, 3))
  cuts <- cbind(args$lo, args$hi, pmin(pmax(images, args$lo), args$hi))
  cuts <- matrix(cuts[order(row(cuts), cuts)], n, byrow = TRUE)

  panels <- matrix(-Inf, n, ncol(cuts) - 1)
  for (i in seq_len(ncol(panels))) {
    rows <- which(cuts[, i + 1] > cuts[, i])
    if (!length(rows)) next
    half <- (cuts[rows, i + 1] - cuts[rows, i]) / 2
    u <- (cuts[rows, i + 1] + cuts[rows, i]) / 2 + outer(half, gl$x) - args$mu[rows]
    terms <- log(half) + rep(log(gl$w), each = length(rows)) +
      peak[rows] * exp(args$kappa[rows] * (cos(u) - 1))
    panels[rows, i] <- matrixStats::rowLogSumExps(terms)
  }
  matrixStats::rowLogSumExps(panels)
}

# The core's bound with the observer's prior of storage at 1: the decision
# density is the sdm density itself, largest at the probe and smallest
# opposite it, so log_min and log_max (of 2 pi times the density) follow from
# the activation at 0 and pi and log Z
.sdm_cd_bounded_criterion <- function(criterion, c, kappa, log_z) {
  peak <- exp(.sdm_cd_log_peak(c, kappa))
  .circmix_cd_bound(
    criterion, rep(1, length(criterion)),
    log_min = log(2 * pi) + peak * exp(-2 * kappa) - log_z,
    log_max = log(2 * pi) + peak - log_z
  )
}

# half-width of the "same" arc: where the sdm density equals
# exp(-criterion) / (2 pi), with the criterion already bounded
.sdm_cd_crit_angle <- function(c, kappa, criterion, log_z) {
  log_peak <- .sdm_cd_log_peak(c, kappa)
  thresh <- log_z - criterion - log(2 * pi)
  log_thresh <- log(pmax(thresh, 0))
  out <- acos(pmin(pmax(1 + (log_thresh - log_peak) / kappa, -1), 1))
  out[log_thresh >= log_peak] <- 0
  out[log_thresh <= log_peak - 2 * kappa] <- pi
  out
}

# P("same") for one row per observation or posterior draw; twin of the row
# loop of sdm_cd_simple_lpmf() in sdm_cd_functions.stan, with log Z per row
# where Stan takes it once per run of rows that share c and kappa
.sdm_cd_psame <- function(probe, mu, c, kappa, criterion) {
  args <- .circmix_recycle(probe = probe, mu = mu, c = c, kappa = kappa, criterion = criterion)
  log_z <- .sdm_cd_log_int(-pi, pi, 0, args$c, args$kappa)
  bounded <- .sdm_cd_bounded_criterion(args$criterion, args$c, args$kappa, log_z)
  hw <- .sdm_cd_crit_angle(args$c, args$kappa, bounded, log_z)
  out <- as.numeric(hw >= pi)
  arc <- which(hw > 0 & hw < pi)
  if (length(arc)) {
    out[arc] <- exp(.sdm_cd_log_int(
      args$probe[arc] - hw[arc], args$probe[arc] + hw[arc], args$mu[arc],
      args$c[arc], args$kappa[arc]
    ) - log_z[arc])
  }
  # the arc integral of an unknown arc is empty, which would read as P("same") 0
  out[is.na(args$probe) | is.na(args$mu)] <- NA
  out
}

# The observer of sdm_cd, trial by trial: draw the remembered feature x from
# the sdm density around mu, then answer "change" when the log-likelihood
# ratio -log(2 pi f(x - probe)) exceeds the bounded criterion, f being the sdm
# density. It shares no code with .sdm_cd_psame(): x comes from rsdm()'s
# rejection sampler and log Z from integrate(), not from the Gauss-Legendre
# arc rule, so the two check each other; simulate_data() draws through it for
# the same reason.
.sdm_cd_observe <- function(probe, mu, c, kappa, criterion) {
  args <- .circmix_recycle(probe = probe, mu = mu, c = c, kappa = kappa, criterion = criterion)
  n <- length(args$probe)
  x <- rsdm(n, mu = args$mu, c = args$c, kappa = args$kappa)
  peak <- args$c * sqrt(args$kappa / (2 * pi))
  log_z <- .sdm_cd_observer_log_z(args$c, args$kappa)
  # log of 2 pi times the density, at the probe's distance from x and at the
  # two ends of the trial's range of evidence
  log_2pi_f <- function(d) log(2 * pi) + peak * exp(args$kappa * (cos(d) - 1)) - log_z
  llr <- -log_2pi_f(x - args$probe)
  llr_max <- -log_2pi_f(pi)
  llr_min <- -log_2pi_f(0)
  crit <- args$criterion
  bounded <- crit
  up <- crit > 0
  bounded[up] <- llr_max[up] * tanh(crit[up] / llr_max[up])
  down <- crit < 0
  bounded[down] <- -llr_min[down] * tanh(crit[down] / -llr_min[down])
  as.integer(llr > bounded)
}

# log Z by integrate() on [0, pi], split at multiples of the spike width,
# once per distinct (c, kappa)
.sdm_cd_observer_log_z <- function(c, kappa) {
  key <- unique(data.frame(c = c, kappa = kappa))
  key$log_z <- mapply(function(c, kappa) {
    peak <- c * sqrt(kappa / (2 * pi))
    width <- 1 / sqrt(kappa * max(peak, 1))
    breaks <- sort(unique(c(0, pmin(pi, width * c(1, 4, 16, 64)), pi)))
    f <- function(u) exp(peak * expm1(kappa * (cos(u) - 1)))
    total <- sum(vapply(seq_len(length(breaks) - 1), function(i) {
      stats::integrate(f, breaks[i], breaks[i + 1], rel.tol = 1e-12, abs.tol = 0,
                       subdivisions = 1000L)$value
    }, numeric(1)))
    peak + log(2 * total)
  }, key$c, key$kappa)
  key$log_z[match(paste(c, kappa), paste(key$c, key$kappa))]
}
