############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The retrieval side is mixture2p's version table, so the two models share
# their parameter names and meanings; a joint fit of both tasks relies on that.
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
        "- The probe and the target should be in radians: the feature shown at ",
        "test, and the feature shown at the probed location in the memory array\n",
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
#' @name mixture2p_cd
#' @details `r model_docs(.model_mixture2p_cd(), components = c('domain', 'task', 'name', 'citation', 'requirements'))`
#' #### Version: `simple`
#' `r model_docs(.model_mixture2p_cd(version = "simple"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot`
#' `r model_docs(.model_mixture2p_cd(version = "slot"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `slot_averaging`
#' `r model_docs(.model_mixture2p_cd(version = "slot_averaging"), components = c('parameters', 'fixed_parameters', 'links', 'prior'))`
#'
#' @section The observer:
#' Each trial shows a memory array and then a probe at one location; the
#' response is whether the probe differs from the item that was there. The
#' observer retrieves a feature `x` for the probed location from the same
#' memory as [mixture2p()]: with the probability of a memory response (`thetat`,
#' or `min(1, K / set_size)` under `"slot"`) from a von Mises distribution
#' around the target with precision `kappa`, otherwise uniformly. The observer
#' then compares `x` with the probe (Lin & Oberauer, 2022, Eq. 8) and answers
#' "change" when the log-likelihood ratio of a change against no change,
#' `log[(1 / 2pi) / (p_s vM(x | probe, kappa) + (1 - p_s) / 2pi)]`, exceeds
#' `criterion`. Here `p_s` is the observer's prior that the probed item is in
#' memory, the probability of a memory response. This puts a "same" region
#' around the probe, and the probability of a "same" response is the retrieval
#' mass inside it.
#'
#' Under `"slot_averaging"` an item holds `floor(K / set_size)` slots or one
#' more. An item with at least one slot is known to be stored (`p_s = 1`), and
#' its precision is that of the slots it holds; an item with no slot may still
#' receive the extra one, and the observer's prior that it is stored is the
#' probability of that.
#'
#' The memory parameters have the same names and meanings as in [mixture2p()],
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
#' @param probe The name of the variable containing the feature shown at test,
#'   in radians.
#' @param target The name of the variable containing the feature shown at the
#'   probed location in the memory array, in radians. On "same" trials it
#'   equals the probe.
#' @param set_size Name of the variable in the dataset containing the set size,
#'   or a single number if the set size is constant. Required by the `"slot"`
#'   and `"slot_averaging"` versions and ignored by `"simple"`. With a single
#'   set size the capacity rule cannot be tested: `K` is then estimated from the
#'   responses at that set size, and only below it.
#' @param version Which storage model to use: `"simple"` (default), `"slot"`
#'   or `"slot_averaging"`. They are the versions of [mixture2p()].
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
#'   whether the item is stored, which this argument does not model.
#' @param links A named list of link functions for the model parameters, e.g.
#'   `list(kappa = "softplus")`. Defaults are `"log"` for `kappa` and `K`,
#'   `"logit"` for `thetat` and `"identity"` for `criterion`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @seealso [mixture2p()], the same memory model for continuous reproduction;
#'   [mixture2p_cd_dist] for the probability and random generation functions.
#' @references Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution
#'   representations in visual working memory. Nature, 453(7192), 233-235.
#'
#'   Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. Cognitive
#'   Psychology, 133, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' # half of the trials probe the target itself, the others a random colour
#' dat <- data.frame(
#'   target = 0,
#'   probe = c(rep(0, 1000), runif(1000, -pi, pi))
#' )
#' dat$change <- rmixture2p_cd(2000, dat$probe - dat$target, kappa = 6, p_mem = 0.7)
#'
#' fit <- bmm(
#'   formula = bmmformula(kappa ~ 1, thetat ~ 1),
#'   data = dat,
#'   model = mixture2p_cd(response = "change", probe = "probe", target = "target"),
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

# The set-size step of check_data.mixture2p(), which this class does not
# inherit
#' @export
check_data.mixture2p_cd <- function(model, data, formula) {
  if (!.mixture2p_version_table[[model$version]]$needs_set_size) {
    return(NextMethod("check_data"))
  }
  ss <- check_var_set_size(model$other_vars$set_size, data)
  stopif(
    anyNA(ss$ss_numeric),
    "The set size variable '{model$other_vars$set_size}' contains missing values."
  )
  # newdata for a prediction may hold one set size; the warning is about the
  # data K is estimated from
  warnif(
    !isTRUE(attr(data, "prediction_newdata")) && length(unique(ss$ss_numeric)) == 1,
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
  brms::bf(.circmix_aterm(
    model$resp_vars$response,
    vint = if (.mixture2p_version_table[[model$version]]$needs_set_size) "ss_numeric",
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

# P("same") of observation i for every posterior draw; the version decides
# which of thetat and K the family estimates
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
  force(version)
  function(i, prep) {
    .cd_bernoulli_ld(prep$data$Y[i], .mixture2p_cd_prep_psame(prep, i, version))
  }
}

.mixture2p_cd_posterior_predict <- function(version) {
  force(version)
  function(i, prep, ...) {
    p_same <- .mixture2p_cd_prep_psame(prep, i, version)
    stats::rbinom(length(p_same), 1, 1 - p_same)
  }
}

# The expected response is the probability of a "change" response. Every dpar
# arrives as a draws x observations matrix, so the cells are evaluated as one
# long vector with the data columns lined up with them.
.mixture2p_cd_posterior_epred <- function(version) {
  force(version)
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

log_lik_mixture2p_cd_simple <- .mixture2p_cd_log_lik("simple")
log_lik_mixture2p_cd_slot <- .mixture2p_cd_log_lik("slot")
log_lik_mixture2p_cd_slot_averaging <- .mixture2p_cd_log_lik("slot_averaging")

posterior_predict_mixture2p_cd_simple <- .mixture2p_cd_posterior_predict("simple")
posterior_predict_mixture2p_cd_slot <- .mixture2p_cd_posterior_predict("slot")
posterior_predict_mixture2p_cd_slot_averaging <-
  .mixture2p_cd_posterior_predict("slot_averaging")

posterior_epred_mixture2p_cd_simple <- .mixture2p_cd_posterior_epred("simple")
posterior_epred_mixture2p_cd_slot <- .mixture2p_cd_posterior_epred("slot")
posterior_epred_mixture2p_cd_slot_averaging <-
  .mixture2p_cd_posterior_epred("slot_averaging")

############################################################################# !
# SIMULATE_DATA S3 METHODS                                               ####
############################################################################# !

# The responses come from rmixture2p_cd(), which draws the remembered feature
# and applies the decision rule to it trial by trial, independently of the arc
# integral dmixture2p_cd() and the likelihood evaluate
#' @export
simulate_data.mixture2p_cd <- function(model, pars, design = NULL,
                                       scale = c("link", "native"), ...) {
  needs_set_size <- .mixture2p_version_table[[model$version]]$needs_set_size
  vars <- model$resp_vars
  design <- .design_or_n(
    model, design, c(vars$probe, vars$target, if (needs_set_size) .set_size_column(model)), ...
  )
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  args <- list(
    n = nrow(design), probe = wrap(design[[vars$probe]] - design[[vars$target]]),
    mu = p$mu, kappa = p$kappa, criterion = p$criterion,
    version = model$version, knowledge = model$knowledge
  )
  args <- if (needs_set_size) {
    c(args, list(K = p$K, set_size = check_var_set_size(model$other_vars$set_size, design)$ss_numeric))
  } else {
    c(args, list(p_mem = p$thetat))
  }
  y <- do.call(rmixture2p_cd, args)
  .simulate_result(design, stats::setNames(list(y), vars$response), sim)
}
