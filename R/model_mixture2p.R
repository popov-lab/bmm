############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_mixture2p <- function(resp_error = NULL, links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Two-parameter mixture model by Zhang and Luck (2008).",
      version = "NA",
      citation = glue(
        "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
        representations in visual working memory. Nature, 453(7192), 233-235"
      ),
      requirements = glue(
        "- The response vairable should be in radians and \\
        represent the angular error relative to the target"
      ),
      parameters = list(
        mu1 = glue(
          "Location parameter of the von Mises distribution for memory responses \\
          (in radians). Fixed internally to 0 by default."
        ),
        kappa = "Concentration parameter of the von Mises distribution",
        thetat = "Mixture weight for target responses"
      ),
      links = list(
        mu1 = "tan_half",
        kappa = "log",
        thetat = "logit"
      ),
      fixed_parameters = list(mu1 = 0, mu2 = 0, kappa2 = -100),
      default_priors = list(
        mu1 = list(main = "student_t(1, 0, 1)"),
        kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
        thetat = list(main = "logistic(0, 1)")
      ),
      void_mu = FALSE
    ),
    class = c("bmmodel", "circular", "mixture2p"),
    call = call
  )
  out$links[names(links)] <- links
  out
}

.model_mixture2p_cd <- function(response = NULL, probe = NULL, target = NULL,
                                links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(response, probe, target),
      other_vars = nlist(),
      domain = "Visual working memory",
      task = "Change detection",
      name = "Two-parameter mixture model for single-probe change detection.",
      version = "NA",
      citation = glue(
        "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
        representations in visual working memory. Nature, 453(7192), 233-235; \\
        Lin, H.Y., & Oberauer, K. (2022). An interference model for visual \\
        working memory: Applications to the change detection task. \\
        Cognitive Psychology, 133, 101463."
      ),
      requirements = glue(
        "- response: binary, coded 0 = 'same' and 1 = 'change'
        - probe: the probed feature in radians
        - target: the feature shown at the probed location in radians"
      ),
      parameters = list(
        mu = glue(
          "Bias of the retrieval distribution relative to the target \\
          (in radians). Fixed internally to 0 by default."
        ),
        kappa = "Concentration parameter of the von Mises distribution",
        thetat = "Probability of retrieving the target from memory",
        criterion = glue(
          "Decision criterion. A 'change' response is given when the \\
          log-likelihood ratio exceeds the criterion, so larger values make 'change' \\
          responses less likely. Fixed to 0 (unbiased) by default; supply a \\
          formula for criterion to estimate it."
        )
      ),
      links = list(
        mu = "tan_half",
        kappa = "log",
        thetat = "logit",
        criterion = "identity"
      ),
      fixed_parameters = list(mu = 0, criterion = 0),
      default_priors = list(
        mu = list(main = "student_t(1, 0, 1)"),
        kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
        thetat = list(main = "logistic(0, 1)"),
        criterion = list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
      ),
      void_mu = FALSE
    ),
    class = c("bmmodel", "change_detection", "mixture2p_cd"),
    call = call
  )
  out$links[names(links)] <- links
  out
}

# user facing aliases

#' @title `r .model_mixture2p()$name`
#' @name mixture2p
#' @details
#' `mixture2p()` is a generic entry point that dispatches on the response
#' arguments you supply: `resp_error` selects the continuous reproduction model
#' and `response`/`probe`/`target` select the change detection model. The two
#' tasks require different response formats, so each has its own constructor
#' with its own parameters and likelihood — `mixture2p_de()` and
#' `mixture2p_cd()`. Call those directly when you want to be explicit.
#'
#' ## Continuous reproduction
#'
#' `r model_info(.model_mixture2p())`
#'
#' ## Change detection
#'
#' `r model_info(.model_mixture2p_cd())`
#' @param resp_error The name of the variable in the provided dataset containing
#'   the response error. The response error should code the response relative to
#'   the to-be-recalled target in radians. You can transform the response error
#'   in degrees to radian using the `deg2rad` function.
#' @param response The name of the variable in the dataset containing the binary
#'   change detection response, coded 0 = "same" and 1 = "change".
#' @param probe The name of the variable containing the probed feature in
#'   radians.
#' @param target The name of the variable containing the feature that was shown
#'   at the probed location, in radians.
#' @param links A list of links for the model parameters
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @references
#' Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution representations in
#'   visual working memory. \emph{Nature}, \emph{453}(7192), 233--235.
#'
#' Lin, H. Y., & Oberauer, K. (2022). An interference model for visual working
#'   memory: Applications to the change detection task. \emph{Cognitive
#'   Psychology}, \emph{133}, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # continuous reproduction
#' dat <- data.frame(y = rmixture2p(n = 2000))
#' fit <- bmm(
#'   formula = bmmformula(kappa ~ 1, thetat ~ 1),
#'   data = dat,
#'   model = mixture2p(resp_error = "y"),
#'   cores = 4, iter = 500, backend = "cmdstanr"
#' )
#'
#' # change detection
#' dat_cd <- data.frame(target = 0, probe = runif(2000, -pi, pi))
#' dat_cd$resp <- rmixture2p_cd(2000, dat_cd$probe, kappa = 5, p_mem = 0.7)
#' fit_cd <- bmm(
#'   formula = bmmformula(kappa ~ 1, thetat ~ 1),
#'   data = dat_cd,
#'   model = mixture2p_cd(response = "resp", probe = "probe", target = "target"),
#'   cores = 4, iter = 500, backend = "cmdstanr"
#' )
#' @export
mixture2p <- function(resp_error = NULL, response = NULL, probe = NULL,
                      target = NULL, links = NULL, ...) {
  call <- match.call()
  cd_args <- c(response, probe, target)
  stopif(
    is.null(resp_error) && length(cd_args) == 0,
    "Provide either 'resp_error' for continuous reproduction, or 'response',
    'probe' and 'target' for change detection."
  )
  stopif(
    !is.null(resp_error) && length(cd_args) > 0,
    "Provide either 'resp_error' or the change detection arguments 'response',
    'probe' and 'target', not both."
  )

  out <- if (is.null(resp_error)) {
    stopif(
      length(cd_args) < 3,
      "Change detection requires all of 'response', 'probe' and 'target'."
    )
    .model_mixture2p_cd(
      response = response, probe = probe, target = target, links = links, ...
    )
  } else {
    .model_mixture2p(resp_error = resp_error, links = links, ...)
  }
  attr(out, "call") <- call
  out
}

#' @rdname mixture2p
#' @export
mixture2p_de <- function(resp_error, links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  .model_mixture2p(resp_error = resp_error, links = links, call = call, ...)
}

#' @rdname mixture2p
#' @export
mixture2p_cd <- function(response, probe, target, links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  .model_mixture2p_cd(
    response = response, probe = probe, target = target, links = links,
    call = call, ...
  )
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.mixture2p <- function(model, data, formula) {
  # construct the brmsformula
  formula <- bmf2bf(model, formula) +
    brms::lf(kappa2 ~ 1, mu2 ~ 1) +
    brms::nlf(kappa1 ~ kappa) +
    brms::nlf(theta1 ~ thetat)

  # specify the mixture family
  formula$family <- brms::mixture("von_mises", "von_mises", order = "none")

  nlist(formula, data)
}

#' @export
configure_model.mixture2p_cd <- function(model, data, formula) {
  free_crit <- as.integer(not_in("criterion", names(model$fixed_parameters)))
  gl <- .cd_gauss_legendre()

  family <- brms::custom_family(
    name = "mixture2p_cd",
    dpars = c("mu", "kappa", "thetat", "criterion"),
    links = c("tan_half", "log", "logit", "identity"),
    lb = c(NA, 0, 0, NA),
    ub = c(NA, NA, 1, NA),
    type = "int",
    vars = c("vreal1", "cd_gl_x", "cd_gl_w", "cd_free_criterion"),
    loop = FALSE,
    log_lik = log_lik_mixture2p_cd,
    posterior_predict = posterior_predict_mixture2p_cd
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/mixture2p_cd_funs.stan"))
  stanvars <- brms::stanvar(x = gl$nodes, name = "cd_gl_x") +
    brms::stanvar(x = gl$weights, name = "cd_gl_w") +
    brms::stanvar(x = free_crit, name = "cd_free_criterion", scode = "int cd_free_criterion;") +
    brms::stanvar(scode = stan_funs, block = "functions")

  formula <- bmf2bf(model, formula)
  formula$family <- family

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.mixture2p_cd <- function(model, formula = bmmformula()) {
  brms::bf(glue(
    "{model$resp_vars$response} | vreal(probe_centered) ~ {.extract_mu_rhs(formula)}"
  ))
}

log_lik_mixture2p_cd <- function(i, prep) {
  dmixture2p_cd(
    prep$data$Y[i],
    prep$data$vreal1[i],
    kappa = brms::get_dpar(prep, "kappa", i = i),
    p_mem = brms::get_dpar(prep, "thetat", i = i),
    criterion = brms::get_dpar(prep, "criterion", i = i),
    mu = brms::get_dpar(prep, "mu", i = i),
    log = TRUE
  )
}

posterior_predict_mixture2p_cd <- function(i, prep, ...) {
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  rmixture2p_cd(
    length(kappa),
    prep$data$vreal1[i],
    kappa = kappa,
    p_mem = brms::get_dpar(prep, "thetat", i = i),
    criterion = brms::get_dpar(prep, "criterion", i = i),
    mu = brms::get_dpar(prep, "mu", i = i)
  )
}
