############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_mixture2p <- function(resp_error = NULL, response = NULL,
                             probe = NULL, target = NULL, task = "de",
                             links = NULL, call = NULL, ...) {
  domain_class <- if (task == "cd") "change_detection" else "circular"

  out <- structure(
    list(
      resp_vars = if (task == "cd") nlist(response, probe, target) else nlist(resp_error),
      other_vars = nlist(),
      domain = "Visual working memory",
      task = if (task == "cd") "Change detection" else "Continuous reproduction",
      name = "Two-parameter mixture model by Zhang and Luck (2008).",
      version = "NA",
      citation = glue(
        "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
        representations in visual working memory. Nature, 453(7192), 233-235"
      ),
      parameters = list(
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
        kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
        thetat = list(main = "logistic(0, 1)")
      ),
      fixed_parameters = list(),
      void_mu = task == "cd"
    ),
    class = c("bmmodel", domain_class, "mixture2p", paste0("mixture2p_", task)),
    call = call
  )

  if (task == "cd") {
    out$citation <- glue(
      "{out$citation}; \\
      Lin, H.Y., & Oberauer, K. (2022). An interference model for visual \\
      working memory: Applications to the change detection task. \\
      Cognitive Psychology, 133, 101463."
    )
    out$requirements <- glue(
      "- response: Binary (0='same', 1='change')
      - probe: Probe color in radians
      - target: Target color in radians"
    )
    out$parameters <- c(
      list(mu = glue(
        "Location parameter (bias) of the retrieval distribution \\
        (in radians). Fixed to 0 by default."
      )),
      out$parameters
    )
    out$links <- c(list(mu = "tan_half"), out$links)
    out$fixed_parameters$mu <- 0
    out$default_priors <- c(
      list(mu = list(main = "student_t(1, 0, 1)")),
      out$default_priors
    )
    out$parameters$beta <- glue(
      "Decision criterion (log prior odds). \\
      Fixed to 0 by default for unbiased decision."
    )
    out$links$beta <- "identity"
    out$fixed_parameters$beta <- 0
    out$default_priors$beta <- list(main = "normal(0, 0.5)")
  } else {
    out$requirements <- glue(
      "- The response vairable should be in radians and \\
      represent the angular error relative to the target"
    )
    out$parameters <- c(
      list(mu1 = glue(
        "Location parameter of the von Mises distribution for memory responses \\
        (in radians). Fixed internally to 0 by default."
      )),
      out$parameters
    )
    out$links <- c(list(mu1 = "tan_half"), out$links)
    out$fixed_parameters <- list(mu1 = 0, mu2 = 0, kappa2 = -100)
    out$default_priors <- c(
      list(mu1 = list(main = "student_t(1, 0, 1)")),
      out$default_priors
    )
  }

  out$links[names(links)] <- links
  out
}

# user facing alias

#' @title `r .model_mixture2p()$name`
#' @details `r model_info(.model_mixture2p())`
#' @param resp_error The name of the variable in the provided dataset containing
#'   the response error. The response Error should code the response relative to
#'   the to-be-recalled target in radians. You can transform the response error
#'   in degrees to radian using the `deg2rad` function. Required when
#'   `task = "de"`.
#' @param response The name of the variable in the dataset containing the binary
#'   response (0 = "same", 1 = "change"). Required when `task = "cd"`.
#' @param probe The name of the variable containing the probe color in radians.
#'   Required when `task = "cd"`.
#' @param target The name of the variable containing the target color in
#'   radians. Required when `task = "cd"`.
#' @param task Character. The experimental task: `"de"` for delayed estimation
#'   (continuous reproduction) or `"cd"` for change detection. Default is `"de"`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # generate artificial data
#' dat <- data.frame(y = rmixture2p(n = 2000))
#'
#' # define formula
#' ff <- bmmformula(kappa ~ 1, thetat ~ 1)
#'
#' model <- mixture2p(resp_error = "y")
#'
#' # fit the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   iter = 500,
#'   backend = "cmdstanr"
#' )
#' @export
mixture2p <- function(resp_error = NULL, response = NULL, probe = NULL,
                      target = NULL, task = "de", ...) {
  call <- match.call()
  task <- match.arg(task, c("de", "cd"))
  if (task == "de") {
    stopif(is.null(resp_error), "Argument 'resp_error' is required for task = 'de'.")
  } else {
    warning2("`mixture2p(task = \"cd\")` is deprecated. Please use `mixture2p_cd()` instead.")
    stopif(is.null(response), "Argument 'response' is required for task = 'cd'.")
    stopif(is.null(probe), "Argument 'probe' is required for task = 'cd'.")
    stopif(is.null(target), "Argument 'target' is required for task = 'cd'.")
  }
  .model_mixture2p(resp_error = resp_error, response = response,
                    probe = probe, target = target, task = task,
                    call = call, ...)
}

#' @rdname mixture2p
#' @export
mixture2p_cd <- function(response, probe, target, ...) {
  call <- match.call()
  stop_missing_args()
  .model_mixture2p(
    response = response,
    probe = probe,
    target = target,
    task = "cd",
    call = call,
    ...
  )
}

############################################################################# !
# CONFIGURE_MODEL METHODS                                                ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.mixture2p_de <- function(model, data, formula) {
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
  mixture2p_cd <- brms::custom_family(
    name = "mixture2p_cd",
    dpars = c("mu", "kappa", "thetat", "beta"),
    links = c("tan_half", "log", "logit", "identity"),
    lb = c(NA, 0, 0, NA),
    ub = c(NA, NA, 1, NA),
    type = "int",
    loop = FALSE,
    log_lik = log_lik_mixture2p_cd,
    posterior_predict = posterior_predict_mixture2p_cd
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/mixture2p_cd_funs.stan"))
  stan_tdata <- read_lines2(paste0(sc_path, "/mixture2p_cd_tdata.stan"))
  stan_likelihood <- read_lines2(paste0(sc_path, "/mixture2p_cd_likelihood.stan"))
  stanvars <- brms::stanvar(x = data$probe_centered, name = "probe_cd") +
    brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")

  formula <- bmf2bf(model, formula)
  formula$family <- mixture2p_cd

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.mixture2p_cd <- function(model, formula = bmmformula()) {
  resp_name <- model$resp_vars$response
  mu_rhs <- .extract_mu_rhs(formula)
  brms::bf(glue("{resp_name} ~ {mu_rhs}"))
}

log_lik_mixture2p_cd <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  thetat <- brms::get_dpar(prep, "thetat", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- .extract_cd_probe(i, prep)
  y <- prep$data$Y[i]
  dmixture2p_cd(y, probe, kappa = kappa, thetat = thetat, beta = beta,
                mu = mu, log = TRUE)
}

posterior_predict_mixture2p_cd <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  thetat <- brms::get_dpar(prep, "thetat", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- .extract_cd_probe(i, prep)
  rmixture2p_cd(length(kappa), probe, kappa = kappa, thetat = thetat,
                beta = beta, mu = mu)
}
