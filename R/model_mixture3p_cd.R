############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The retrieval side is mixture3p's version table, so the two models share
# their parameter names, links, priors and inits; a joint fit of both tasks
# relies on that. Only criterion is added.
.model_mixture3p_cd <- function(response = NULL, probe = NULL,
                                nt_features = NULL, set_size = NULL,
                                regex = FALSE, links = NULL, version = "simple",
                                variable_precision = FALSE,
                                knowledge = "limited", call = NULL, ...) {
  spec <- .cd_add_criterion(.mixture3p_version_table[[version]])

  out <- structure(
    list(
      resp_vars = nlist(response, probe),
      other_vars = nlist(nt_features, set_size),
      domain = "Working memory",
      task = "Detection, recognition and confidence judgments",
      name = "Three-Parameter Mixture Model for Change Detection",
      citation = c(
        .mixture3p_citation,
        if (version != "simple") .mixture3p_capacity_citation,
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
        "relative to the target\n",
        "  - In the 'slot' and 'slot_averaging' versions, the set size should ",
        "vary for the capacity rule to be testable"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = c(list(mu = 0), spec$fixed_parameters),
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      variable_precision = variable_precision,
      knowledge = knowledge
    ),
    regex = regex,
    regex_vars = c("nt_features"),
    class = c(
      "bmmodel", "change_detection", "non_targets", "mixture3p_cd",
      paste0("mixture3p_cd_", version)
    ),
    call = call
  )
  set_links(out, links)
}

# as in mixture3p(), the softmax over thetat and thetant is part of the
# likelihood
#' @exportS3Method
settable_links.mixture3p_cd <- function(model) {
  names(model$links)[model$links != "softmax"]
}

# user facing alias

#' @title `r .model_mixture3p_cd()$name`
#' @name mixture3p_cd
#' @description Fits the three-parameter mixture model of [mixture3p()] to the
#'   "same" or "change" responses of a single-probe change-detection task.
#' @details `r model_docs(.model_mixture3p_cd(version = "slot"), components = c('domain', 'task', 'name', 'citation', 'requirements'))`
#' #### Version: `simple`
#' `r model_docs(.model_mixture3p_cd(version = "simple"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot`
#' `r model_docs(.model_mixture3p_cd(version = "slot"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot_averaging`
#' `r model_docs(.model_mixture3p_cd(version = "slot_averaging"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#'
#' @section The observer:
#' Each trial shows a memory array and then a probe at one location; the
#' response is whether the probe differs from the item that was there. The
#' observer retrieves a feature `x` for the probed location from the same
#' memory as [mixture3p()]: from a von Mises distribution with precision
#' `kappa` around the target, around one of the non-targets (a swap), or
#' uniformly (a guess). The observer then compares `x` with the probe (Lin &
#' Oberauer, 2022, Eq. 8) and answers "change" when the log-likelihood ratio
#' of a change against no change,
#' `log[(1 / 2pi) / (p_s vM(x | probe, kappa) + (1 - p_s) / 2pi)]`, exceeds
#' `criterion`. This puts a "same" region around the probe, and the probability
#' of a "same" response is the retrieval mass inside it.
#'
#' `p_s` is the observer's prior that `x` comes from the probed item. Under
#' `"simple"` and `"slot"` it is the probability of retrieving the target, as
#' in Lin and Oberauer's (2022) observers that do not know whether the item is
#' stored: a swap looks to the observer like a guess. Under `"slot_averaging"`
#' an item holds `floor(K / set_size)` slots or one more. An item with at least
#' one slot is known to be stored (`p_s = 1`) and its precision is that of the
#' slots it holds; an item with no slot may still receive the extra one, and
#' the observer's prior that it is stored is the probability of that. At the
#' default `criterion = 0`, `p_s` drops out of the decision.
#'
#' Non-targets enter only through retrieval, so a probe that matches a
#' non-target is harder to reject than a new feature: that is the intrusion
#' cost. To estimate the swap parameter (`thetant` or `ptarget`) the design
#' needs such intrusion probes, beside probes of the target itself and of new
#' features. A trial with a single item has no non-target, so the set-size-1
#' level of the swap parameter is not identified; the model warns when a
#' formula asks for it.
#'
#' The memory parameters have the same names and meanings as in [mixture3p()],
#' plus `criterion`, so the two models can share parameters in a joint fit of
#' both tasks. Whether a person gives the same values in both tasks is an
#' empirical question rather than something the model guarantees: the precision
#' estimated from continuous reproduction includes response noise that change
#' detection does not have, for example.
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
#' @param set_size Name of the column containing the set size variable (if
#'   set_size varies) or a numeric value for the set_size, if the set_size is
#'   fixed. In the `"slot"` and `"slot_averaging"` versions a single set size
#'   cannot test the capacity rule: `K` is then estimated from the responses at
#'   that set size, and only below it.
#' @param regex Logical. If TRUE, the `nt_features` argument is interpreted as
#'   a regular expression to match the non-target feature columns in the
#'   dataset.
#' @param version Which storage model to use: `"simple"` (default), `"slot"`
#'   or `"slot_averaging"`. They are the versions of [mixture3p()] and have the
#'   same parameters.
#' @param variable_precision Must be `FALSE`. Variable precision is not
#'   available for change detection yet; it is planned for bmm 1.6.0.
#' @param knowledge What the observer knows about the precision of the current
#'   trial when precision varies across trials. With `"limited"` (default), the
#'   observer knows only the distribution of precision and uses one decision
#'   boundary for all trials; with `"rich"`, the observer knows the precision of
#'   the trial and sets the boundary accordingly. Until bmm 1.6.0 brings
#'   variable precision, precision varies only under `version =
#'   "slot_averaging"`, for an item that holds one slot or one more, and
#'   `knowledge` changes the predictions of that version only. Lin and Oberauer
#'   (2022) found the limited observer described their data better under
#'   variable precision; for slot averaging they varied a different knowledge,
#'   whether the item is stored, which this argument does not model. `"rich"`
#'   warns in the versions where it has no effect.
#' @param links A named list of link functions for the model parameters, e.g.
#'   `list(kappa = "softplus")`. The links of `thetat` and `thetant` cannot be
#'   changed, because the softmax over them is part of the likelihood, nor that
#'   of `criterion`, because its fixed value and its bound are written for the
#'   identity scale.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @seealso [mixture3p()], the same memory model for continuous reproduction;
#'   [mixture3p_cd_dist] for the probability and random generation functions.
#' @references Bays, P. M., Catalao, R. F. G., & Husain, M. (2009). The
#'   precision of visual working memory is set by allocation of a shared
#'   resource. Journal of Vision, 9(10), Article 7.
#'
#'   Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations
#'   in visual working memory. Nature, 453(7192), 233-235.
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. Cognitive
#'   Psychology, 133, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # three items, every feature relative to the target: a third of the probes
#' # show the target (a probe of 0), a third the first non-target and a third a
#' # new colour
#' n <- 1500
#' design <- data.frame(
#'   nt1 = runif(n, -pi, pi), nt2 = runif(n, -pi, pi), set_size = 3,
#'   type = rep(c("same", "intrusion", "new"), length.out = n)
#' )
#' design$probe <- ifelse(design$type == "same", 0,
#'   ifelse(design$type == "intrusion", design$nt1, runif(n, -pi, pi))
#' )
#' model <- mixture3p_cd(
#'   response = "change", probe = "probe",
#'   nt_features = c("nt1", "nt2"), set_size = "set_size"
#' )
#' dat <- simulate_data(model,
#'   pars = list(kappa = 2, thetat = 1, thetant = -0.5), design = design
#' )
#'
#' fit <- bmm(
#'   formula = bmmformula(kappa ~ 1, thetat ~ 1, thetant ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#'
#' # the capacity versions estimate K and ptarget instead of the two weights;
#' # criterion is estimated once it has a formula
#' slot <- mixture3p_cd("change", "probe", c("nt1", "nt2"), "set_size", version = "slot")
#' slot_averaging <- mixture3p_cd(
#'   "change", "probe", c("nt1", "nt2"), "set_size",
#'   version = "slot_averaging", knowledge = "rich"
#' )
#' capacity_formula <- bmmformula(kappa ~ 1, K ~ 1, ptarget ~ 1, criterion ~ 1)
#' @export
mixture3p_cd <- function(response, probe, nt_features, set_size,
                         regex = FALSE,
                         version = c("simple", "slot", "slot_averaging"),
                         variable_precision = FALSE,
                         knowledge = c("limited", "rich"), links = NULL, ...) {
  call <- match.call()
  .check_constructor_dots(list(...))
  stop_missing_args()
  version <- match.arg(version)
  knowledge <- match.arg(knowledge)
  .cd_check_variable_precision(variable_precision)
  .cd_warn_inert_knowledge(knowledge, version, acts = version == "slot_averaging")

  .model_mixture3p_cd(
    response = response, probe = probe,
    nt_features = nt_features, set_size = set_size, regex = regex,
    version = version, variable_precision = variable_precision,
    knowledge = knowledge, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

# The steps of check_data.mixture3p(), which this class does not inherit;
# check_data.non_targets() runs first and adds ss_numeric
#' @export
check_data.mixture3p_cd <- function(model, data, formula) {
  nt_features <- model$other_vars$nt_features
  ss <- data$ss_numeric
  stopif(
    anyNA(ss),
    "The set size variable '{model$other_vars$set_size}' contains missing \\
    values. Remove the trials without a set size before fitting."
  )
  for_prediction <- isTRUE(attr(data, "prediction_newdata"))
  set_sizes <- unique(ss)
  warnif(
    !for_prediction && model$version != "simple" && length(set_sizes) == 1,
    "The set size does not vary, so the data cannot test the capacity rule of \\
    the '{model$version}' version. K is estimated from the responses at set \\
    size {set_sizes}, and reflects only its prior where it reaches {set_sizes}."
  )
  NextMethod("check_data")
}

############################################################################# !
# CHECK_FORMULA METHODS                                                  ####
############################################################################# !

#' @export
check_formula.mixture3p_cd <- function(model, data, formula) {
  .mixture3p_warn_set_size_one(model, data, formula)
  NextMethod("check_formula")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

# the probe comes first, so the wrapper receives it as vreal1 and the
# non-targets after it
#' @export
bmf2bf.mixture3p_cd <- function(model, formula) {
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = "ss_numeric",
    vreal = c(model$resp_vars$probe, model$other_vars$nt_features)
  ))
}

#' @export
configure_model.mixture3p_cd <- function(model, data, formula) {
  spec <- .mixture3p_version_table[[model$version]]
  n_nt <- attr(data, "max_set_size") - 1
  family <- paste0("mixture3p_cd_", model$version)

  # the softmax over the weights is part of the likelihood, so brms passes the
  # weights through on the scale they are sampled on
  family_model <- model
  family_model$links[family_model$links == "softmax"] <- "identity"

  formula <- bmf2bf(model, formula)
  formula$family <- .cd_custom_family(
    family_model, family, spec$weight_parameters,
    vint = TRUE, n_vreal = 1 + n_nt,
    log_lik = .mixture3p_cd_log_lik(model$version),
    posterior_predict = .mixture3p_cd_posterior_predict(model$version),
    posterior_epred = .mixture3p_cd_posterior_epred(model$version)
  )

  nlist(
    formula, data,
    stanvars = .cd_model_stanvars(
      model, formula$family,
      c("mixture3p_functions.stan", "mixture3p_cd_functions.stan"),
      vint = "ss", vreal = list(probe = 1, nt = n_nt)
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

# P("same") of observation i for every posterior draw
.mixture3p_cd_prep_psame <- function(prep, i, version) {
  dpar <- function(name) brms::get_dpar(prep, name, i = i)
  capacity <- version != "simple"
  .mixture3p_cd_psame(
    probe = .cd_prep_probe(prep, i),
    mu = dpar("mu"),
    kappa = dpar("kappa"),
    w1 = dpar(if (capacity) "K" else "thetat"),
    w2 = dpar(if (capacity) "ptarget" else "thetant"),
    criterion = dpar("criterion"),
    set_size = prep$data$vint1[i],
    nt = .cd_prep_nt(prep, i),
    version = version,
    rich = .cd_prep_rich(prep)
  )
}

.mixture3p_cd_log_lik <- function(version) {
  force(version)
  function(i, prep) {
    .cd_bernoulli_ld(prep$data$Y[i], .mixture3p_cd_prep_psame(prep, i, version))
  }
}

.mixture3p_cd_posterior_predict <- function(version) {
  force(version)
  function(i, prep, ...) {
    p_same <- .mixture3p_cd_prep_psame(prep, i, version)
    stats::rbinom(length(p_same), 1, 1 - p_same)
  }
}

# The expected response is the probability of a "change" response. Each
# observation has its own set size and non-targets, so the observations are
# evaluated one at a time.
.mixture3p_cd_posterior_epred <- function(version) {
  force(version)
  function(prep) {
    .epred_matrix(vapply(seq_len(prep$nobs), function(i) {
      1 - .mixture3p_cd_prep_psame(prep, i, version)
    }, numeric(prep$ndraws)), prep)
  }
}

log_lik_mixture3p_cd_simple <- .mixture3p_cd_log_lik("simple")
log_lik_mixture3p_cd_slot <- .mixture3p_cd_log_lik("slot")
log_lik_mixture3p_cd_slot_averaging <- .mixture3p_cd_log_lik("slot_averaging")

posterior_predict_mixture3p_cd_simple <- .mixture3p_cd_posterior_predict("simple")
posterior_predict_mixture3p_cd_slot <- .mixture3p_cd_posterior_predict("slot")
posterior_predict_mixture3p_cd_slot_averaging <-
  .mixture3p_cd_posterior_predict("slot_averaging")

posterior_epred_mixture3p_cd_simple <- .mixture3p_cd_posterior_epred("simple")
posterior_epred_mixture3p_cd_slot <- .mixture3p_cd_posterior_epred("slot")
posterior_epred_mixture3p_cd_slot_averaging <-
  .mixture3p_cd_posterior_epred("slot_averaging")

############################################################################# !
# SIMULATE_DATA S3 METHODS                                               ####
############################################################################# !

# The responses come from .mixture3p_cd_observe(), which draws the remembered
# feature and applies the decision rule to it trial by trial, independently of
# the arc integral dmixture3p_cd() and the likelihood evaluate. The simple
# version weighs its components with the likelihood's softmax over the link
# values, which drops the non-targets at set size 1 as the fit does.
#' @export
simulate_data.mixture3p_cd <- function(model, pars, design = NULL,
                                       scale = c("link", "native"), ...) {
  vars <- model$resp_vars
  nt_features <- model$other_vars$nt_features
  design <- .design_or_n(
    model, design, c(vars$probe, .set_size_column(model)),
    ragged = nt_features
  )
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  set_size <- check_var_set_size(model$other_vars$set_size, design)$ss_numeric
  .check_nt_na(design, nt_features, set_size)
  nt <- unname(as.matrix(design[nt_features]))
  capacity <- model$version != "simple"
  y <- .mixture3p_cd_observe(
    probe = wrap(design[[vars$probe]]),
    mu = p$mu, nt = nt, set_size = set_size, kappa = p$kappa,
    w1 = if (capacity) p$K else sim$link$thetat,
    w2 = if (capacity) p$ptarget else sim$link$thetant,
    criterion = p$criterion, version = model$version,
    rich = model$knowledge == "rich"
  )
  .simulate_result(design, stats::setNames(list(y), vars$response), sim)
}
