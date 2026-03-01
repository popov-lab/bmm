############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_mixture2p <- function(resp_error = NULL, response = NULL,
                             probe = NULL, target = NULL, task = "de",
                             links = NULL, call = NULL, ...) {
  domain_class <- if (task == "cd") "change_detection" else "circular"

  if (task == "cd") {
    out <- structure(
      list(
        resp_vars = nlist(response, probe, target),
        other_vars = nlist(),
        domain = "Visual working memory",
        task = "Change detection",
        name = "Two-parameter mixture model by Zhang and Luck (2008).",
        version = "NA",
        citation = glue(
          "Zhang, W., & Luck, S. J. (2008). Discrete fixed-resolution \\
          representations in visual working memory. Nature, 453(7192), 233-235; \\
          Lin, H.Y., & Oberauer, K. (2022). An interference model for visual \\
          working memory: Applications to the change detection task. \\
          Cognitive Psychology, 133, 101463."
        ),
        requirements = glue(
          "- response: Binary (0='same', 1='change')
          - probe: Probe color in radians
          - target: Target color in radians"
        ),
        parameters = list(
          kappa = "Concentration parameter of the von Mises distribution",
          thetat = "Mixture weight for target responses",
          beta = glue(
            "Decision criterion (log prior odds). \\
            Fixed to 0 by default for unbiased decision."
          )
        ),
        links = list(
          kappa = "log",
          thetat = "identity",
          beta = "identity"
        ),
        fixed_parameters = list(beta = 0),
        default_priors = list(
          kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
          thetat = list(main = "logistic(0, 1)"),
          beta = list(main = "normal(0, 0.5)")
        ),
        void_mu = TRUE
      ),
      class = c("bmmodel", domain_class, "mixture2p", "mixture2p_cd"),
      call = call
    )
  } else {
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
          thetat = "identity"
        ),
        fixed_parameters = list(mu1 = 0, mu2 = 0, kappa2 = -100),
        default_priors = list(
          mu1 = list(main = "student_t(1, 0, 1)"),
          kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)"),
          thetat = list(main = "logistic(0, 1)")
        ),
        void_mu = FALSE
      ),
      class = c("bmmodel", domain_class, "mixture2p", "mixture2p_de"),
      call = call
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
    stopif(is.null(response), "Argument 'response' is required for task = 'cd'.")
    stopif(is.null(probe), "Argument 'probe' is required for task = 'cd'.")
    stopif(is.null(target), "Argument 'target' is required for task = 'cd'.")
  }
  .model_mixture2p(resp_error = resp_error, response = response,
                    probe = probe, target = target, task = task,
                    call = call, ...)
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
  resp_name <- model$resp_vars$response
  probe_var <- attr(data, "probe_var")

  mixture2p_cd <- brms::custom_family(
    name = "mixture2p_cd",
    dpars = c("mu", "kappa", "thetat", "beta"),
    links = c("identity", "log", "identity", "identity"),
    lb = c(NA, 0, NA, NA),
    ub = c(NA, NA, NA, NA),
    type = "int",
    loop = TRUE,
    vars = "vreal1[n]",
    log_lik = log_lik_mixture2p_cd,
    posterior_predict = posterior_predict_mixture2p_cd
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/mixture2p_cd_funs.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  formula <- bmf2bf(model, formula)
  formula$family <- mixture2p_cd

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.mixture2p_cd <- function(model, formula = bmmformula()) {
  resp_name <- model$resp_vars$response
  probe_var <- "probe_centered"
  brms_formula <- brms::bf(
    glue("{resp_name} | vreal({probe_var}) ~ 1")
  )
  components <- lapply(formula, function(x) {
    if (is_nl(x)) brms::nlf(x) else brms::lf(x)
  })
  Reduce(`+`, components, init = brms_formula)
}

log_lik_mixture2p_cd <- function(i, prep) {
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  thetat <- brms::get_dpar(prep, "thetat", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- prep$data$vreal1[i]
  y <- prep$data$Y[i]
  dmixture2p_cd(y, probe, kappa = kappa, thetat = thetat, beta = beta,
                log = TRUE)
}

posterior_predict_mixture2p_cd <- function(i, prep, ...) {
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  thetat <- brms::get_dpar(prep, "thetat", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- prep$data$vreal1[i]
  rmixture2p_cd(length(kappa), probe, kappa = kappa, thetat = thetat,
                beta = beta)
}
