############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The retrieval side is mixture2p's version table, so the two models share
# their parameter names, priors, links and inits; a joint fit of both tasks
# relies on that.
.model_mixture2p_cd <- function(response = NULL, probe = NULL, target = NULL,
                                set_size = NULL, links = NULL, version = "simple",
                                variable_precision = FALSE, knowledge = "limited",
                                call = NULL, ...) {
  spec <- .cd_add_criterion(.mixture2p_version_table[[version]])

  out <- structure(
    list(
      resp_vars = nlist(response, probe, target),
      other_vars = nlist(set_size),
      domain = "Working memory",
      task = "Detection, recognition and confidence judgments",
      name = "Two-Parameter Mixture Model for Change Detection",
      citation = c(
        glue(
          "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
          representations in visual working memory. Nature, 453(7192), 233-235. \\
          https://doi.org/10.1038/nature06860"
        ),
        glue(
          "Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual \\
          working memory: Applications to the change detection task. Cognitive \\
          Psychology, 133, 101463. https://doi.org/10.1016/j.cogpsych.2022.101463"
        )
      ),
      version = version,
      requirements = paste0(
        "- The response variable should be coded 0 ('same') and 1 ('change')\n",
        "- The probe and the target should be in radians: the feature shown at test ",
        "and the feature shown at the probed location in the memory array\n",
        "  - The 'slot' and 'slot_averaging' versions additionally require a ",
        "set_size variable, which should vary for the capacity rule to be testable"
      ),
      parameters = spec$parameters,
      links = spec$links,
      fixed_parameters = spec$fixed_parameters,
      default_priors = spec$priors,
      init_ranges = spec$init_ranges,
      variable_precision = variable_precision,
      knowledge = knowledge
    ),
    class = c(
      "bmmodel", "change_detection", "mixture2p_cd",
      paste0("mixture2p_cd_", version)
    ),
    call = call
  )
  set_links(out, links)
}

# user facing alias

#' @title `r .model_mixture2p_cd()$name`
#' @details `r model_docs(.model_mixture2p_cd(), components = c('domain', 'task', 'name', 'citation', 'requirements'))`
#' #### Version: `simple`
#' `r model_docs(.model_mixture2p_cd(version = "simple"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot`
#' `r model_docs(.model_mixture2p_cd(version = "slot"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot_averaging`
#' `r model_docs(.model_mixture2p_cd(version = "slot_averaging"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#'
#' @section The observer:
#' On each trial one item of the memory array is probed, and the participant
#' says whether the probe is the same as the item shown at that location or a
#' change. The observer retrieves a feature from the same memory as in
#' [mixture2p()]: from memory with the probability of a memory response, and
#' as a guess otherwise. It then compares the retrieved feature with the
#' probe (Lin & Oberauer, 2022, Eq. 8). "Change" is the answer when the
#' log-likelihood ratio of a change against no change exceeds `criterion`,
#' which puts a "same" region around the probe; the probability of a "same"
#' response is the retrieval mass inside it. The observer's prior that the
#' probed item is stored is the probability of a memory response (`thetat`,
#' or `min(1, K / set_size)`).
#'
#' The memory parameters have the same names, links and default priors as in
#' [mixture2p()], plus `criterion`. Whether a person gives the same values in
#' both tasks is an empirical question rather than something the model
#' guarantees: the precision estimated from continuous reproduction includes
#' response noise that change detection does not have, for example.
#'
#' @section Decision criterion:
#' `criterion` is fixed to 0 by default, the unbiased observer, for whom the
#' prior of storage drops out of the decision. Give it a formula, such as
#' `criterion ~ 1`, to estimate it; positive values favour "same" responses.
#' An estimated criterion is kept inside the range of evidence each trial can
#' produce: it equals Lin and Oberauer's (2022) criterion near 0 and is
#' compressed toward the ends of that range, so both responses stay possible
#' on every trial. Without the bound, a criterion past the range makes one
#' response impossible and the likelihood drops off a cliff the sampler cannot
#' cross. Large estimated values are therefore not directly comparable with
#' Lin and Oberauer's unbounded criterion. bmm prints a message saying so when
#' the criterion is estimated.
#'
#' @param response The name of the variable in the data containing the
#'   change-detection response, coded 0 ("same") and 1 ("change"). A logical
#'   variable is accepted and coded `TRUE` = "change".
#' @param probe The name of the variable containing the feature shown at test,
#'   in radians.
#' @param target The name of the variable containing the feature shown at the
#'   probed location in the memory array, in radians. On "same" trials it
#'   equals the probe.
#' @param set_size Name of the variable in the data containing the set size,
#'   or a single number if the set size is constant. Required by the `"slot"`
#'   and `"slot_averaging"` versions and ignored by `"simple"`. With a single
#'   set size the capacity rule cannot be tested: `K` is then estimated from
#'   the responses at that set size, and only below it.
#' @param version Which storage model to use: `"simple"` (default), `"slot"`
#'   or `"slot_averaging"`, the versions of [mixture2p()].
#' @param variable_precision Must be `FALSE`. Variable precision is planned for
#'   change detection in bmm 1.6.0; `TRUE` stops with an error.
#' @param knowledge What the observer knows about the precision of the
#'   current trial when it varies across trials: `"limited"` (default) or
#'   `"rich"`. In this version of bmm precision varies only under
#'   `version = "slot_averaging"`, where a held item has one of two slot counts;
#'   the other versions give the same predictions with either value. The
#'   limited observer knows only the distribution of the slot counts and uses
#'   one decision boundary for all trials; the rich observer knows the slot
#'   count of the trial and sets the boundary from it. Lin and Oberauer (2022)
#'   found the limited observer fit better under variable precision; for slot
#'   averaging they varied whether the observer knows that the item is stored,
#'   which this argument does not model.
#' @param links A named list of link functions for the model parameters, e.g.
#'   `list(kappa = "softplus")`. Defaults are `"log"` for `kappa` and `K`,
#'   `"logit"` for `thetat` and `"identity"` for `criterion`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @seealso [mixture2p()] for continuous reproduction, [mixture2p_cd_dist] for
#'   the distribution functions.
#' @keywords bmmodel
#' @references Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
#'   representations in visual working memory. Nature, 453(7192), 233-235.
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. Cognitive
#'   Psychology, 133, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # half of the trials probe the target itself, the others a random colour
#' design <- data.frame(target = 0, probe = c(rep(0, 1000), runif(1000, -pi, pi)))
#' model <- mixture2p_cd(response = "change", probe = "probe", target = "target")
#' dat <- simulate_data(model, pars = list(kappa = 2, thetat = 1), design = design)
#'
#' fit <- bmm(
#'   formula = bmmformula(kappa ~ 1, thetat ~ 1),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#' @export
mixture2p_cd <- function(response, probe, target, set_size = NULL,
                         version = c("simple", "slot", "slot_averaging"),
                         variable_precision = FALSE,
                         knowledge = c("limited", "rich"), links = NULL, ...) {
  call <- match.call()
  .check_constructor_dots(list(...))
  stop_missing_args()
  version <- match.arg(version)
  knowledge <- match.arg(knowledge)
  .cd_check_variable_precision(variable_precision)
  .mixture2p_cd_check_version_args(version, set_size)

  .model_mixture2p_cd(
    response = response, probe = probe, target = target, set_size = set_size,
    version = version, variable_precision = variable_precision,
    knowledge = knowledge, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_DATA METHODS                                                     ####
############################################################################# !

# the set-size step of check_data.mixture2p(), which this class does not
# inherit
#' @export
check_data.mixture2p_cd <- function(model, data, formula) {
  if (!.mixture2p_version_table[[model$version]]$needs_set_size) {
    return(NextMethod("check_data"))
  }
  for_prediction <- isTRUE(attr(data, "prediction_newdata"))
  set_size <- model$other_vars$set_size
  stopif(
    is.character(set_size) && anyNA(data[[set_size]]),
    "The set size variable '{set_size}' contains missing values. Remove the \\
    trials without a set size before fitting."
  )
  ss <- check_var_set_size(set_size, data)
  warnif(
    !for_prediction && length(unique(ss$ss_numeric)) == 1,
    "The set size does not vary, so the data cannot test the capacity rule of \\
    the '{model$version}' version. K is estimated from the responses at set \\
    size {ss$max_set_size}, and reflects only its prior where it reaches \\
    {ss$max_set_size}."
  )
  data$ss_numeric <- ss$ss_numeric
  attr(data, "max_set_size") <- ss$max_set_size

  NextMethod("check_data")
}

#' @exportS3Method
built_data_columns.mixture2p_cd <- function(model) {
  c(
    if (.mixture2p_version_table[[model$version]]$needs_set_size) "ss_numeric",
    NextMethod("built_data_columns")
  )
}

#' @exportS3Method
revert_check_data.mixture2p_cd <- function(model, data) {
  set_size <- model$other_vars$set_size
  # brms keeps the set_size column only when a formula predicts something with
  # it; ss_numeric is check_data()'s numeric copy, kept for the vint() term
  if (is.character(set_size) && not_in(set_size, colnames(data)) &&
    "ss_numeric" %in% colnames(data)) {
    data[[set_size]] <- data$ss_numeric
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), set_size)
  }
  NextMethod("revert_check_data")
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !

#' @export
bmf2bf.mixture2p_cd <- function(model, formula) {
  needs_set_size <- .mixture2p_version_table[[model$version]]$needs_set_size
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = if (needs_set_size) "ss_numeric",
    vreal = "probe_centered"
  ))
}

#' @export
configure_model.mixture2p_cd <- function(model, data, formula) {
  spec <- .mixture2p_version_table[[model$version]]
  family <- paste0("mixture2p_cd_", model$version)

  formula <- bmf2bf(model, formula)
  formula$family <- .cd_custom_family(
    model, family,
    weight_parameters = spec$weight_parameter,
    vint = spec$needs_set_size,
    log_lik = .mixture2p_cd_log_lik(model$version),
    posterior_predict = .mixture2p_cd_posterior_predict(model$version),
    posterior_epred = .mixture2p_cd_posterior_epred(model$version)
  )

  nlist(
    formula, data,
    stanvars = .cd_model_stanvars(
      model, formula$family,
      c("mixture2p_functions.stan", "mixture2p_cd_functions.stan"),
      vint = if (spec$needs_set_size) "ss"
    )
  )
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

# P("same") of observation i for every posterior draw
.mixture2p_cd_prep_psame <- function(prep, i, version) {
  capacity <- version != "simple"
  .mixture2p_cd_psame(
    probe = .cd_prep_probe(prep, i),
    mu = brms::get_dpar(prep, "mu", i = i),
    kappa = brms::get_dpar(prep, "kappa", i = i),
    p_mem = if (!capacity) brms::get_dpar(prep, "thetat", i = i) else 1,
    criterion = brms::get_dpar(prep, "criterion", i = i),
    K = if (capacity) brms::get_dpar(prep, "K", i = i) else 1,
    set_size = if (capacity) prep$data$vint1[i] else 1,
    version = version,
    rich = .cd_prep_rich(prep)
  )
}

.mixture2p_cd_log_lik <- function(version) {
  switch(version,
    simple = log_lik_mixture2p_cd_simple,
    slot = log_lik_mixture2p_cd_slot,
    slot_averaging = log_lik_mixture2p_cd_slot_averaging
  )
}

.mixture2p_cd_posterior_predict <- function(version) {
  switch(version,
    simple = posterior_predict_mixture2p_cd_simple,
    slot = posterior_predict_mixture2p_cd_slot,
    slot_averaging = posterior_predict_mixture2p_cd_slot_averaging
  )
}

.mixture2p_cd_posterior_epred <- function(version) {
  switch(version,
    simple = posterior_epred_mixture2p_cd_simple,
    slot = posterior_epred_mixture2p_cd_slot,
    slot_averaging = posterior_epred_mixture2p_cd_slot_averaging
  )
}

.mixture2p_cd_ll <- function(version) {
  function(i, prep) {
    .cd_bernoulli_ld(prep$data$Y[i], .mixture2p_cd_prep_psame(prep, i, version))
  }
}

.mixture2p_cd_pp <- function(version) {
  function(i, prep, ...) {
    p_same <- .mixture2p_cd_prep_psame(prep, i, version)
    stats::rbinom(length(p_same), 1, 1 - p_same)
  }
}

# The expected response is the probability of a "change" response. Every dpar
# arrives as a draws x observations matrix, so the cells are evaluated as one
# long vector and the data columns are lined up with them.
.mixture2p_cd_epred <- function(version) {
  function(prep) {
    capacity <- version != "simple"
    dpars <- prep$dpars
    p_same <- .mixture2p_cd_psame(
      probe = .epred_data(prep$data$vreal1, prep),
      mu = as.vector(dpars$mu),
      kappa = as.vector(dpars$kappa),
      p_mem = if (!capacity) as.vector(dpars$thetat) else 1,
      criterion = as.vector(dpars$criterion),
      K = if (capacity) as.vector(dpars$K) else 1,
      set_size = if (capacity) .epred_data(prep$data$vint1, prep) else 1,
      version = version,
      rich = .cd_prep_rich(prep)
    )
    .epred_matrix(1 - p_same, prep)
  }
}

log_lik_mixture2p_cd_simple <- .mixture2p_cd_ll("simple")
log_lik_mixture2p_cd_slot <- .mixture2p_cd_ll("slot")
log_lik_mixture2p_cd_slot_averaging <- .mixture2p_cd_ll("slot_averaging")

posterior_predict_mixture2p_cd_simple <- .mixture2p_cd_pp("simple")
posterior_predict_mixture2p_cd_slot <- .mixture2p_cd_pp("slot")
posterior_predict_mixture2p_cd_slot_averaging <- .mixture2p_cd_pp("slot_averaging")

posterior_epred_mixture2p_cd_simple <- .mixture2p_cd_epred("simple")
posterior_epred_mixture2p_cd_slot <- .mixture2p_cd_epred("slot")
posterior_epred_mixture2p_cd_slot_averaging <- .mixture2p_cd_epred("slot_averaging")

############################################################################# !
# SIMULATE_DATA S3 METHODS                                               ####
############################################################################# !

# Draws through the observer, .mixture2p_cd_observe(), never through the
# P("same") of the likelihood
#' @export
simulate_data.mixture2p_cd <- function(model, pars, design = NULL, scale = c("link", "native"), ...) {
  capacity <- .mixture2p_version_table[[model$version]]$needs_set_size
  vars <- model$resp_vars
  design <- .design_or_n(model, design, c(
    vars$probe, vars$target, if (capacity) .set_size_column(model)
  ), ...)
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  y <- .mixture2p_cd_observe(
    probe = wrap(design[[vars$probe]] - design[[vars$target]]),
    mu = p$mu, kappa = p$kappa,
    p_mem = if (!capacity) p$thetat else 1,
    criterion = p$criterion,
    K = if (capacity) p$K else 1,
    set_size = if (capacity) {
      check_var_set_size(model$other_vars$set_size, design)$ss_numeric
    } else {
      1
    },
    version = model$version, rich = identical(model$knowledge, "rich")
  )
  .simulate_result(design, stats::setNames(list(y), vars$response), sim)
}
