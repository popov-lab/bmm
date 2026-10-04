############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The parameters, links, priors and initial ranges are those of sdm() plus
# criterion, so that a change-detection and a continuous-reproduction fit name
# the same memory parameters. There is a single precision state, so the
# knowledge-limited and knowledge-rich observers of the circular mixtures
# coincide and sdm_cd() takes no knowledge argument.
.model_sdm_cd <- function(response = NULL, probe = NULL, target = NULL,
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
      resp_vars = nlist(response, probe, target),
      other_vars = nlist(),
      domain = "Visual working memory (change detection)",
      task = "Change detection",
      name = "Signal Discrimination Model (SDM) for change detection",
      citation = c(
        de$citation,
        glue(
          "Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual \\
          working memory: Applications to the change detection task. Cognitive \\
          Psychology, 133, 101463. https://doi.org/10.1016/j.cogpsych.2022.101463"
        )
      ),
      version = version,
      requirements = glue(
        "- The response variable is binary, coded 0 = 'same' and 1 = 'change' \\
        (or logical, TRUE = 'change')
        - The probe and the target are features in radians; the model uses the \\
        probe relative to the target
        - The decision rule of Lin & Oberauer (2022) is applied to the SDM \\
        density, which is not one of the models they compared"
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
  out <- set_links(out, links)
  out
}

# As in sdm(), configure_model.sdm_cd declares the family links itself and the
# log link of `c` is written into the Stan chunk, so no link can be set
#' @exportS3Method
settable_links.sdm_cd <- function(model) {
  character(0)
}

#' @title `r .model_sdm_cd()$name`
#' @name sdm_cd
#' @details The observer retrieves a feature from the SDM density of the
#'   probed item and answers "change" when the log-likelihood ratio of a
#'   change against no change exceeds `criterion` (Lin & Oberauer, 2022, Eq. 8,
#'   with the SDM density in place of their von Mises). The "same" responses
#'   then fall in an arc around the probe whose half-width depends on `c`,
#'   `kappa` and `criterion`; unlike the von Mises models, no parameter drops
#'   out of the boundary at `criterion = 0`. `c` and `kappa` mean the same as in
#'   [sdm()], so fits of the two tasks report the same parameters, but whether
#'   one participant has the same values in both tasks is an empirical
#'   question, not a property of the model. The SDM is not one of the models
#'   Lin & Oberauer (2022) fitted to change-detection data, so this model
#'   extends their decision rule rather than reproducing one of their fits.
#'   Within-chain threading is not supported yet. `r model_docs(.model_sdm_cd())`
#' @param response The name of the variable in the data containing the binary
#'   response, coded 0 = "same" and 1 = "change". A logical variable is
#'   accepted and coded as 1 for `TRUE`.
#' @param probe The name of the variable in the data containing the probed
#'   feature in radians.
#' @param target The name of the variable in the data containing the feature
#'   of the probed item at encoding in radians.
#' @param version Character. The version of the model to use. Currently only
#'   "simple" is supported.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @seealso [sdm()] for continuous reproduction, [dsdm_cd()] for the
#'   distribution functions.
#' @export
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simulate change-detection responses from the model
#' dat <- data.frame(target = 0, probe = runif(1000, -pi, pi))
#' dat$change <- rsdm_cd(1000, probe = dat$probe, c = 5, kappa = 4)
#'
#' fit <- bmm(
#'   formula = bmf(c ~ 1, kappa ~ 1),
#'   data = dat,
#'   model = sdm_cd(response = "change", probe = "probe", target = "target"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
sdm_cd <- function(response, probe, target, version = "simple", ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  .model_sdm_cd(
    response = response, probe = probe, target = target, version = version,
    call = call, ...
  )
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

# The normalising constant is computed once per run of rows sharing all
# predictor values, as for sdm()
#' @export
check_data.sdm_cd <- function(model, data, formula) {
  data <- order_data_query(model, data, formula)
  attr(data, "sdm_run_metadata") <- sdm_run_metadata(data, formula, model)
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
bmf2bf.sdm_cd <- function(model, formula) {
  brms::bf(.circmix_aterm(model$resp_vars$response, vreal = "probe_centered"))
}

#' @export
configure_model.sdm_cd <- function(model, data, formula) {
  # Under threading brms slices Y and the dpars but passes vreal1 and the run
  # metadata whole, so the likelihood would pair rows with the wrong probes
  stopif(
    brms_slices_likelihood(),
    "sdm_cd() does not support within-chain threading yet. Fit it without \\
    the 'threads' argument (and without options(brms.threads))."
  )
  family <- brms::custom_family(
    name = "sdm_cd",
    dpars = c("mu", "c", "kappa", "criterion"),
    # c has a log link that the Stan chunk applies itself, as in sdm()
    links = c("tan_half", "identity", "log", "identity"),
    lb = c(NA, NA, NA, NA),
    ub = c(NA, NA, NA, NA),
    type = "int",
    vars = c("vreal1", "G_sdm_runs", "sdm_run_start", "sdm_run_count", "cd_gl_x", "cd_gl_w"),
    loop = FALSE,
    log_lik = log_lik_sdm_cd,
    posterior_predict = posterior_predict_sdm_cd,
    posterior_epred = posterior_epred_sdm_cd
  )

  run_metadata <- attr(data, "sdm_run_metadata") %||% sdm_run_metadata(data, formula, model)
  sc_path <- system.file("stan_chunks", package = "bmm")
  gl <- .cd_gl_rule()
  stanvars <- .circmix_chunk_stanvar("cd_funs.stan") +
    .circmix_chunk_stanvar("sdm_cd_funs.stan") +
    brms::stanvar(scode = read_lines2(paste0(sc_path, "/sdm_simple_tdata.stan")), block = "tdata") +
    brms::stanvar(x = run_metadata$G_sdm_runs, name = "G_sdm_runs") +
    sdm_stanvar_int_array(run_metadata$sdm_run_start, "sdm_run_start", "G_sdm_runs") +
    sdm_stanvar_int_array(run_metadata$sdm_run_count, "sdm_run_count", "G_sdm_runs") +
    brms::stanvar(x = gl$x, name = "cd_gl_x") +
    brms::stanvar(x = gl$w, name = "cd_gl_w")

  formula <- bmf2bf(model, formula)
  formula$family <- family

  nlist(formula, data, stanvars)
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

# brms samples log(c) under an identity link; declaring the log link after the
# fit makes get_dpar() return c on its natural scale, which the R densities use
#' @export
postprocess_brm.sdm_cd <- function(model, fit, ...) {
  fit$family$link_c <- "log"
  fit$formula$family$link_c <- "log"
  fit
}

#' @export
revert_postprocess_brm.sdm_cd <- function(model, fit, ...) {
  fit$family$link_c <- "identity"
  fit$formula$family$link_c <- "identity"
  fit
}

log_lik_sdm_cd <- function(i, prep) {
  p_same <- .sdm_cd_dpar_psame(prep, i)
  .cd_bernoulli_ld(rep_len(prep$data$Y[i], length(p_same)), p_same)
}

posterior_predict_sdm_cd <- function(i, prep, ...) {
  p_same <- .sdm_cd_dpar_psame(prep, i)
  stats::rbinom(length(p_same), 1, 1 - p_same)
}

# The expected response is P("change")
posterior_epred_sdm_cd <- function(prep) {
  .epred_matrix(1 - with(prep$dpars, .sdm_cd_psame(
    .epred_data(prep$data$vreal1, prep), mu, c, kappa, criterion
  )), prep)
}

.sdm_cd_dpar_psame <- function(prep, i) {
  .sdm_cd_psame(
    .cd_prep_probe(prep, i),
    mu = brms::get_dpar(prep, "mu", i = i),
    c = brms::get_dpar(prep, "c", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    criterion = brms::get_dpar(prep, "criterion", i = i)
  )
}

############################################################################# !
# DENSITY HELPERS                                                        ####
############################################################################# !
# R twins of inst/stan_chunks/sdm_cd_funs.stan, vectorised over rows; c is on
# its natural scale here and on the log scale in Stan.
#
# The sdm density is far more peaked than a von Mises, so the arc integral
# cannot use the single 32-node rule of the circular mixtures (it errs by 149%
# against an adaptive reference; 128 nodes still by 46%, measured on the August
# change-detection branch, 2df8c295). The integral is taken in log space,
# because the activation reaches exp(709) inside the sampled range; it is split
# at every periodic image of the peak, because an arc may span the whole
# circle; and at multiples of the spike width, so the peak sits on a panel
# boundary where Gauss-Legendre clusters its nodes. With all three the rule held
# to 3.8e-10 over c <= 60, kappa <= 150 and arcs from 0.05 rad to the circle.
#
# The normalising constant uses the same rule rather than the Chebyshev
# quadrature of sdm() (sdm_simple_ldenom_chquad_adaptive). Measured against
# integrate() on an R transcription of it, that one errs by up to 5.7e-4 in
# log Z for log(c) in [-1, 5] and kappa in [0.5, 30], and by 0.012 to 7.5 where
# its node count reaches the cap of 200 (log(c) >= 4 at kappa = 100, >= 2 at
# kappa = 300). Here an error in log Z moves both the boundary and P("same").

.sdm_cd_log_peak <- function(c, kappa) {
  log(c) + 0.5 * log(kappa / (2 * pi))
}

.sdm_cd_spike_width <- function(c, kappa) {
  pmin(pi, sqrt(2 * pmax(.sdm_cd_log_peak(c, kappa), 1) / kappa))
}

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

.sdm_cd_crit_angle <- function(c, kappa, criterion, log_z) {
  log_peak <- .sdm_cd_log_peak(c, kappa)
  thresh <- log_z - criterion - log(2 * pi)
  log_thresh <- log(pmax(thresh, 0))
  out <- acos(pmin(pmax(1 + (log_thresh - log_peak) / kappa, -1), 1))
  out[log_thresh >= log_peak] <- 0
  out[log_thresh <= log_peak - 2 * kappa] <- pi
  out
}

.sdm_cd_psame <- function(probe, mu, c, kappa, criterion) {
  args <- .circmix_recycle(probe = probe, mu = mu, c = c, kappa = kappa, criterion = criterion)
  log_z <- .sdm_cd_log_int(-pi, pi, 0, args$c, args$kappa)
  hw <- .sdm_cd_crit_angle(args$c, args$kappa, args$criterion, log_z)
  out <- as.numeric(hw >= pi)
  arc <- which(hw > 0 & hw < pi)
  if (length(arc)) {
    out[arc] <- exp(.sdm_cd_log_int(
      args$probe[arc] - hw[arc], args$probe[arc] + hw[arc], args$mu[arc],
      args$c[arc], args$kappa[arc]
    ) - log_z[arc])
  }
  out
}
