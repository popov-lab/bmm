############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The retrieval side is imm's version table, so the two models share their
# parameter names, links, priors and inits; a joint fit of both tasks relies on
# that. Only criterion is added.
.model_imm_cd <- function(response = NULL, probe = NULL,
                          nt_features = NULL, nt_distances = NULL,
                          set_size = NULL, regex = FALSE, version = "full",
                          variable_precision = FALSE, knowledge = "limited",
                          links = NULL, call = NULL, ...) {
  spec <- .imm_version_spec(version)
  spec$fixed_parameters <- list(mu = 0, b = 0)
  spec <- .cd_add_criterion(spec)

  out <- structure(
    list(
      resp_vars = nlist(response, probe),
      other_vars = nlist(nt_features, nt_distances, set_size),
      domain = "Working memory",
      task = "Detection, recognition and confidence judgments",
      name = "Interference Measurement Model for Change Detection",
      citation = c(
        glue(
          "Oberauer, K., & Lin, H.-Y. (2017). An interference model \\
          of visual working memory. Psychological Review, 124(1), 21-59. \\
          https://doi.org/10.1037/rev0000044"
        ),
        .cd_citation
      ),
      version = version,
      # paste0() rather than glue(), which would trim the indent that keeps
      # the requirements in one list
      requirements = paste0(
        "- The response variable should be coded 0 ('same') and 1 ('change')\n",
        "  - The probe should be in radians and relative to the target: the ",
        "feature shown at test minus the feature shown at the probed location, ",
        "which is 0 on a 'same' trial\n",
        "  - The non-target features should be in radians and be centered ",
        "relative to the target",
        if (version != "abc") {
          paste0(
            "\n  - The non-target distances should be positive, in the unit ",
            "the spatial gradient `s` is to be read in"
          )
        }
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = spec$fixed_parameters,
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      variable_precision = variable_precision,
      knowledge = knowledge
    ),
    regex = regex,
    regex_vars = if (version == "abc") "nt_features" else c("nt_features", "nt_distances"),
    class = c(
      "bmmodel", "change_detection", "non_targets", "imm_cd",
      paste0("imm_cd_", version)
    ),
    call = call
  )
  set_links(out, links)
}

# as in imm(): b is fixed to 0 on its log link, the background the other
# activations are weighed against
#' @exportS3Method
settable_links.imm_cd <- function(model) {
  setdiff(names(model$links), "b")
}

# user facing alias

#' @title `r .model_imm_cd()$name`
#' @name imm_cd
#' @description Fits the interference measurement model of [imm()] to the
#'   "same" or "change" responses of a single-probe change-detection task.
#' @details `r model_docs(.model_imm_cd(), components = c('domain', 'task', 'name', 'citation'))`
#' #### Version: `full`
#' `r model_docs(.model_imm_cd(version = "full"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `bsc`
#' `r model_docs(.model_imm_cd(version = "bsc"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `abc`
#' `r model_docs(.model_imm_cd(version = "abc"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#'
#' @section The observer:
#' Each trial shows a memory array and then a probe at one location; the
#' response is whether the probe differs from the item that was there. The
#' observer retrieves a feature `x` for the probed location from the same
#' memory as [imm()]: every item is activated by the cue in proportion to
#' `c * exp(-s * distance) + a` (the target at distance 0), the background by
#' `b`, and `x` comes from a von Mises distribution with precision `kappa`
#' around the item that wins, or uniformly when the background wins. The
#' observer then compares `x` with the probe (Lin & Oberauer, 2022, Eq. 8) and
#' answers "change" when the log-likelihood ratio of a change against no change,
#' `log[(1 / 2pi) / (p_s vM(x | probe, kappa) + (1 - p_s) / 2pi)]`, exceeds
#' `criterion`. This puts a "same" region around the probe, and the probability
#' of a "same" response is the retrieval mass inside it.
#'
#' `p_s` is the observer's prior that `x` comes from the probed item, here the
#' probability of retrieving the target: the observer does not know which item
#' won, so a non-target retrieval looks to it like a guess. At the default
#' `criterion = 0`, `p_s` drops out of the decision.
#'
#' Non-targets enter only through retrieval, so a probe that matches a
#' non-target is harder to reject than a new feature, the more so the nearer
#' that non-target was to the target: that is the intrusion cost. To estimate
#' `a` and `s` the design needs such intrusion probes, beside probes of the
#' target itself and of new features. A trial with a single item has no
#' non-target, so it says nothing about `s` and informs `a` and `c` only
#' through their sum; the model warns when a formula asks for set-size-1
#' levels of these parameters.
#'
#' The memory parameters have the same names and meanings as in [imm()], plus
#' `criterion`, so the two models can share parameters in a joint fit of both
#' tasks. Whether a person gives the same values in both tasks is an empirical
#' question rather than something the model guarantees: the precision estimated
#' from continuous reproduction includes response noise that change detection
#' does not have, for example.
#'
#' @section Decision criterion:
#' `criterion` is fixed to 0 by default, the unbiased observer, for whom `p_s`
#' drops out of the decision. Give it a formula (`criterion ~ 1`) to estimate
#' it; positive values favour "same" responses. An estimated criterion is kept
#' inside the range of evidence each trial can produce: it equals Lin and
#' Oberauer's criterion near 0 and is compressed toward the ends of that range,
#' so both responses stay possible on every trial. Without the bound, a
#' criterion past the range makes one response impossible on some trials, and
#' the likelihood has edges the sampler cannot cross. Large estimated values are
#' therefore not directly comparable with Lin and Oberauer's (2022) unbounded
#' criterion; bmm prints a message when the criterion is estimated.
#'
#' @param response The name of the variable in the dataset containing the
#'   change-detection response, coded 0 ("same") and 1 ("change"). A logical
#'   variable is accepted and coded `TRUE` = "change".
#' @param probe The name of the variable containing the probe in radians,
#'   relative to the target: the feature shown at test minus the feature shown
#'   at the probed location in the memory array, as [wrap()] returns it. It is 0
#'   on "same" trials.
#' @param nt_features A character vector with the names of the non-target
#'   feature values. The non-target feature values should be in radians and
#'   centered relative to the target. Alternatively, if `regex = TRUE`, a
#'   regular expression can be used to match the non-target feature columns in
#'   the dataset.
#' @param nt_distances A vector of names of the columns containing the distances
#'   of the non-target items to the target item. Alternatively, if
#'   `regex = TRUE`, a regular expression can be used to match the non-target
#'   distance columns in the dataset. Only necessary for the `"bsc"` and
#'   `"full"` versions.
#' @param set_size Name of the column containing the set size variable (if
#'   set_size varies) or a numeric value for the set_size, if the set_size is
#'   fixed.
#' @param regex Logical. If TRUE, the `nt_features` and `nt_distances`
#'   arguments are interpreted as regular expressions to match the non-target
#'   feature and distance columns in the dataset.
#' @param version Which version of the IMM to use: `"full"` (default), `"bsc"`
#'   or `"abc"`. They are the versions of [imm()] and have the same parameters.
#' @param variable_precision Must be `FALSE`. Variable precision is not
#'   available for change detection yet; it is planned for bmm 1.6.0.
#' @param knowledge What the observer knows about the precision of the current
#'   trial when precision varies across trials. With `"limited"` (default), the
#'   observer knows only the distribution of precision and uses one decision
#'   boundary for all trials; with `"rich"`, the observer knows the precision of
#'   the trial and sets the boundary accordingly. Precision does not vary across
#'   trials in any version of `imm_cd()` until bmm 1.6.0 brings variable
#'   precision, so until then `knowledge` does not change its predictions. Lin
#'   and Oberauer (2022) found the limited observer described their data better
#'   under variable precision; their comparison for the interference model
#'   concerned a different knowledge, whether the target is in the focus of
#'   attention, which [imm()] does not model, and favoured the rich observer.
#'   `"rich"` warns, because it has no effect in any version yet.
#' @param links A named list of link functions for the model parameters, e.g.
#'   `list(kappa = "softplus")`. The link of `b` cannot be changed, because `b`
#'   is fixed on its log scale, nor that of `criterion`, because its fixed value
#'   and its bound are written for the identity scale.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @seealso [imm()], the same memory model for continuous reproduction;
#'   [imm_cd_dist] for the probability and random generation functions.
#' @references Oberauer, K., & Lin, H.-Y. (2017). An interference model of
#'   visual working memory. Psychological Review, 124(1), 21-59.
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. Cognitive
#'   Psychology, 133, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # three items at random distances from the target, every feature relative
#' # to the target: a third of the probes show the target (a probe of 0), a
#' # third the first non-target and a third a new colour
#' n <- 1500
#' design <- data.frame(
#'   nt1 = runif(n, -pi, pi), nt2 = runif(n, -pi, pi),
#'   dist1 = runif(n, 0.3, 2), dist2 = runif(n, 0.3, 2), set_size = 3,
#'   type = rep(c("same", "intrusion", "new"), length.out = n)
#' )
#' design$probe <- ifelse(design$type == "same", 0,
#'   ifelse(design$type == "intrusion", design$nt1, runif(n, -pi, pi))
#' )
#' model <- imm_cd(
#'   response = "change", probe = "probe",
#'   nt_features = c("nt1", "nt2"), nt_distances = c("dist1", "dist2"),
#'   set_size = "set_size"
#' )
#' dat <- simulate_data(model,
#'   pars = list(kappa = 2, c = 1.4, a = -0.7, s = 0.7), design = design
#' )
#'
#' fit <- bmm(
#'   formula = bmmformula(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#'
#' # "bsc" has no `a` and "abc" no `s`, nor distances; criterion is estimated
#' # once it has a formula
#' bsc <- imm_cd("change", "probe", c("nt1", "nt2"), c("dist1", "dist2"), "set_size",
#'   version = "bsc"
#' )
#' abc <- imm_cd("change", "probe", c("nt1", "nt2"), set_size = "set_size", version = "abc")
#' abc_formula <- bmmformula(kappa ~ 1, c ~ 1, a ~ 1, criterion ~ 1)
#' @export
imm_cd <- function(response, probe, nt_features, nt_distances, set_size,
                   regex = FALSE, version = c("full", "bsc", "abc"),
                   variable_precision = FALSE, knowledge = c("limited", "rich"),
                   links = NULL, ...) {
  call <- match.call()
  .check_constructor_dots(list(...))
  version <- match.arg(version)
  knowledge <- match.arg(knowledge)
  if (version == "abc") nt_distances <- NULL
  stop_missing_args()
  .cd_check_variable_precision(variable_precision)
  .cd_warn_inert_knowledge(knowledge, version, acts = FALSE)

  .model_imm_cd(
    response = response, probe = probe,
    nt_features = nt_features, nt_distances = nt_distances,
    set_size = set_size, regex = regex, version = version,
    variable_precision = variable_precision, knowledge = knowledge,
    links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

# check_data.non_targets() runs first and adds ss_numeric
#' @export
check_data.imm_cd <- function(model, data, formula) {
  stopif(
    anyNA(data$ss_numeric),
    "The set size variable '{model$other_vars$set_size}' contains missing \\
    values. Remove the trials without a set size before fitting."
  )
  NextMethod("check_data")
}

#' @export
check_data.imm_cd_bsc <- function(model, data, formula) {
  data <- .check_data_imm_dist(model, data, formula)
  NextMethod("check_data")
}

#' @export
check_data.imm_cd_full <- function(model, data, formula) {
  data <- .check_data_imm_dist(model, data, formula)
  NextMethod("check_data")
}

############################################################################# !
# CHECK_FORMULA METHODS                                                  ####
############################################################################# !

#' @export
check_formula.imm_cd <- function(model, data, formula) {
  .imm_warn_set_size_one(model, data, formula)
  NextMethod("check_formula")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

# the probe comes first, so the wrapper receives it as vreal1; the non-target
# features and distances follow as for imm()
#' @export
bmf2bf.imm_cd <- function(model, formula) {
  covariates <- c(model$resp_vars$probe, model$other_vars$nt_features)
  if (.imm_version_spec(model$version)$needs_distances) {
    covariates <- c(covariates, model$other_vars$nt_distances)
  }
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = "ss_numeric", vreal = covariates
  ))
}

#' @export
configure_model.imm_cd <- function(model, data, formula) {
  spec <- .imm_version_spec(model$version)
  n_nt <- attr(data, "max_set_size") - 1
  groups <- c(
    list(probe = 1, nt = n_nt),
    if (spec$needs_distances) list(dist = n_nt)
  )

  formula <- bmf2bf(model, formula)
  formula$family <- .cd_custom_family(
    model, paste0("imm_cd_", model$version), spec$weight_parameters,
    vint = TRUE, n_vreal = sum(unlist(groups)),
    log_lik = .imm_cd_log_lik(model$version),
    posterior_predict = .imm_cd_posterior_predict(model$version),
    posterior_epred = .imm_cd_posterior_epred(model$version)
  )

  nlist(
    formula, data,
    stanvars = .cd_model_stanvars(
      model, formula$family, c("imm_functions.stan", "imm_cd_functions.stan"),
      vint = "ss", vreal = groups
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

# P("same") of observation i for every posterior draw. The non-targets follow
# the probe in vreal, the distances after them.
.imm_cd_prep_psame <- function(prep, i, version) {
  dpar <- function(name) brms::get_dpar(prep, name, i = i)
  covariates <- .cd_prep_nt(prep, i)
  with_distances <- version != "abc"
  half <- if (with_distances) length(covariates) / 2 else length(covariates)
  .imm_cd_psame(
    probe = .cd_prep_probe(prep, i),
    mu = dpar("mu"),
    kappa = dpar("kappa"),
    c = dpar("c"),
    a = if (version != "bsc") dpar("a"),
    s = if (with_distances) dpar("s"),
    b = dpar("b"),
    criterion = dpar("criterion"),
    set_size = prep$data$vint1[i],
    nt = covariates[seq_len(half)],
    dist = if (with_distances) covariates[half + seq_len(half)],
    version = version
  )
}

.imm_cd_log_lik <- function(version) {
  force(version)
  function(i, prep) {
    .cd_bernoulli_ld(prep$data$Y[i], .imm_cd_prep_psame(prep, i, version))
  }
}

.imm_cd_posterior_predict <- function(version) {
  force(version)
  function(i, prep, ...) {
    p_same <- .imm_cd_prep_psame(prep, i, version)
    stats::rbinom(length(p_same), 1, 1 - p_same)
  }
}

# The expected response is the probability of a "change" response. Each
# observation has its own set size, non-targets and distances, so the
# observations are evaluated one at a time.
.imm_cd_posterior_epred <- function(version) {
  force(version)
  function(prep) {
    .epred_matrix(vapply(seq_len(prep$nobs), function(i) {
      1 - .imm_cd_prep_psame(prep, i, version)
    }, numeric(prep$ndraws)), prep)
  }
}

log_lik_imm_cd_full <- .imm_cd_log_lik("full")
log_lik_imm_cd_bsc <- .imm_cd_log_lik("bsc")
log_lik_imm_cd_abc <- .imm_cd_log_lik("abc")

posterior_predict_imm_cd_full <- .imm_cd_posterior_predict("full")
posterior_predict_imm_cd_bsc <- .imm_cd_posterior_predict("bsc")
posterior_predict_imm_cd_abc <- .imm_cd_posterior_predict("abc")

posterior_epred_imm_cd_full <- .imm_cd_posterior_epred("full")
posterior_epred_imm_cd_bsc <- .imm_cd_posterior_epred("bsc")
posterior_epred_imm_cd_abc <- .imm_cd_posterior_epred("abc")

############################################################################# !
# SIMULATE_DATA S3 METHODS                                               ####
############################################################################# !

# The responses come from .imm_cd_observe(), which computes the activations,
# draws the retrieved item and the remembered feature and applies the decision
# rule to it trial by trial, independently of the arc integral dimm_cd() and
# the likelihood evaluate. The fixed background b comes from the model as in
# the fit.
#' @export
simulate_data.imm_cd <- function(model, pars, design = NULL,
                                 scale = c("link", "native"), ...) {
  vars <- model$resp_vars
  other <- model$other_vars
  design <- .design_or_n(
    model, design, c(vars$probe, .set_size_column(model)),
    ragged = c(other$nt_features, other$nt_distances)
  )
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  set_size <- check_var_set_size(other$set_size, design)$ss_numeric
  .check_nt_na(design, other$nt_features, set_size)
  .check_nt_na(design, other$nt_distances, set_size)
  y <- .imm_cd_observe(
    probe = wrap(design[[vars$probe]]),
    mu = p$mu, nt = unname(as.matrix(design[other$nt_features])),
    dist = if (length(other$nt_distances)) unname(as.matrix(design[other$nt_distances])),
    set_size = set_size, kappa = p$kappa, c = p$c, a = p$a, s = p$s, b = p$b,
    criterion = p$criterion
  )
  .simulate_result(design, stats::setNames(list(y), vars$response), sim)
}
