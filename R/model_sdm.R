############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdm <- function(resp_error = NULL, response = NULL, probe = NULL,
                       target = NULL, links = NULL, version = "simple",
                       task = "de", call = NULL, ...) {
  domain_class <- if (task == "cd") "change_detection" else "circular"

  out <- structure(
    list(
      resp_vars = if (task == "cd") nlist(response, probe, target) else nlist(resp_error),
      other_vars = nlist(),
      domain = "Visual working memory",
      task = if (task == "cd") "Change detection" else "Continuous reproduction",
      name = "Signal Discrimination Model (SDM) by Oberauer (2023)",
      version = version,
      citation = glue(
        "Oberauer, K. (2023). Measurement models for visual working memory - \\
        A factorial model comparison. Psychological Review, 130(3), 841-852"
      ),
      parameters = list(
        c = "Memory strength parameter of the SDM distribution",
        kappa = "Precision parameter of the SDM distribution"
      ),
      links = list(c = "log", kappa = "log"),
      default_priors = list(
        kappa = list(main = "student_t(5, 1.75, 0.75)", effects = "normal(0, 1)"),
        c = list(main = "student_t(5, 2, 0.75)", effects = "normal(0, 1)")
      ),
      fixed_parameters = list(),
      void_mu = task == "cd"
    ),
    class = c("bmmodel", domain_class, "sdm", paste0("sdm_", version),
              paste0("sdm_", version, "_", task)),
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
      "- The response variable should be in radians and represent the angular \\
      error relative to the target"
    )
    out$parameters <- c(
      list(mu = glue("Location parameter of the SDM distribution (in radians; \\
                      by default fixed internally to 0)")),
      out$parameters
    )
    out$links <- c(list(mu = "tan_half"), out$links)
    out$fixed_parameters <- list(mu = 0)
    out$default_priors <- c(
      list(mu = list(main = "student_t(1, 0, 1)")),
      out$default_priors
    )
    out$init_ranges <- list(mu = c(-0.5, 0.5), kappa = c(2.5, 3.5), c = c(4, 6))
  }

  out$links[names(links)] <- links
  out
}

# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_sdm_simple(NA)$info

#' @title `r .model_sdm()$name`
#' @name sdm
#' @details see [the online article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for a detailed description of the model
#'   and how to use it. `r model_info(.model_sdm())`
#' @param resp_error The name of the variable in the dataset containing the
#'   response error. The response error should code the response relative to the
#'   to-be-recalled target in radians. You can transform the response error in
#'   degrees to radians using the `deg2rad` function.
#' @param version Character. The version of the model to use. Currently only
#'   "simple" is supported.
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
#' @export
#' @keywords bmmodel
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # simulate data from the model
#' dat <- data.frame(y = rsdm(n = 1000, c = 4, kappa = 3))
#'
#' # specify formula
#' ff <- bmf(
#'   c ~ 1,
#'   kappa ~ 1
#' )
#'
#' # specify the model
#' fit <- bmm(
#'   formula = ff,
#'   data = dat,
#'   model = sdm(resp_error = "y"),
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
sdm <- function(resp_error = NULL, response = NULL, probe = NULL,
                target = NULL, version = "simple", task = "de", ...) {
  call <- match.call()
  task <- match.arg(task, c("de", "cd"))
  if (task == "de") {
    stopif(is.null(resp_error), "Argument 'resp_error' is required for task = 'de'.")
  } else {
    warning2("`sdm(task = \"cd\")` is deprecated. Please use `sdm_cd()` instead.")
    stopif(is.null(response), "Argument 'response' is required for task = 'cd'.")
    stopif(is.null(probe), "Argument 'probe' is required for task = 'cd'.")
    stopif(is.null(target), "Argument 'target' is required for task = 'cd'.")
  }
  .model_sdm(resp_error = resp_error, response = response, probe = probe,
             target = target, version = version, task = task,
             call = call, ...)
}

#' @rdname sdm
#' @keywords deprecated
#' @export
sdmSimple <- function(resp_error, version = "simple", task = "de", ...) {
  warning2("The function `sdmSimple()` is deprecated. Please use `sdm()` instead.")
  call <- match.call()
  stop_missing_args()
  .model_sdm(resp_error = resp_error, version = version, task = task,
             call = call, ...)
}

#' @rdname sdm
#' @export
sdm_cd <- function(response, probe, target, version = "simple", ...) {
  call <- match.call()
  stop_missing_args()
  .model_sdm(
    response = response,
    probe = probe,
    target = target,
    version = version,
    task = "cd",
    call = call,
    ...
  )
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdm <- function(model, data, formula) {
  if (!inherits(model, "change_detection")) {
    data <- order_data_query(model, data, formula)
  }
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.sdm_simple_de <- function(model, data, formula) {
  # note - c has a log link, but I've coded it manually for computational efficiency
  sdm_simple <- brms::custom_family(
    name = "sdm_simple",
    dpars = c("mu", "c", "kappa"),
    links = c("tan_half", "identity", "log"),
    lb = c(NA, NA, NA),
    ub = c(NA, NA, NA),
    type = "real", loop = FALSE,
    log_lik = log_lik_sdm_simple,
    posterior_predict = posterior_predict_sdm_simple
  )

  # prepare initial stanvars to pass to brms, model formula and priors
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdm_simple_funs.stan"))
  stan_tdata <- read_lines2(paste0(sc_path, "/sdm_simple_tdata.stan"))
  stan_likelihood <- read_lines2(paste0(sc_path, "/sdm_simple_likelihood.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")

  # construct main brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  formula$family <- sdm_simple

  # set initial values to be sampled between [-1,1] to avoid extreme SDs that
  # can cause the sampler to fail
  init <- 1

  nlist(formula, data, stanvars, init)
}

#' @export
configure_model.sdm_simple_cd <- function(model, data, formula) {
  sdm_simple_cd <- brms::custom_family(
    name = "sdm_simple_cd",
    dpars = c("mu", "c", "kappa", "beta"),
    links = c("tan_half", "log", "log", "identity"),
    lb = c(NA, NA, 0, NA),
    ub = c(NA, NA, NA, NA),
    type = "int",
    loop = FALSE,
    log_lik = log_lik_sdm_simple_cd,
    posterior_predict = posterior_predict_sdm_simple_cd
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- read_lines2(paste0(sc_path, "/sdm_simple_cd_funs.stan"))
  stan_tdata <- read_lines2(paste0(sc_path, "/sdm_simple_cd_tdata.stan"))
  stan_llik <- read_lines2(paste0(sc_path, "/sdm_simple_cd_likelihood.stan"))
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(scode = stan_llik, block = "likelihood", position = "end")

  formula <- bmf2bf(model, formula)
  formula$family <- sdm_simple_cd

  nlist(formula, data, stanvars)
}

#' @export
bmf2bf.sdm_simple_cd <- function(model, formula = bmmformula()) {
  resp_name <- model$resp_vars$response
  mu_rhs <- .extract_mu_rhs(formula)
  brms::bf(glue("{resp_name} | vreal(probe_centered) ~ {mu_rhs}"))
}

log_lik_sdm_simple_cd <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  c_par <- brms::get_dpar(prep, "c", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- prep$data$vreal1[i]
  y <- prep$data$Y[i]
  dsdm_cd(y, probe, c = c_par, kappa = kappa, beta = beta, mu = mu,
          log = TRUE)
}

posterior_predict_sdm_simple_cd <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  c_par <- brms::get_dpar(prep, "c", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  beta <- brms::get_dpar(prep, "beta", i = i)
  probe <- prep$data$vreal1[i]
  rsdm_cd(length(c_par), probe, c = c_par, kappa = kappa, beta = beta,
          mu = mu)
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

#' @export
postprocess_brm.sdm_simple_de <- function(model, fit, ...) {
  # manually set link_c to "log" since I coded it manually
  fit$family$link_c <- "log"
  fit$formula$family$link_c <- "log"
  fit
}

#' @export
revert_postprocess_brm.sdm_simple_de <- function(model, fit, ...) {
  fit$family$link_c <- "identity"
  fit$formula$family$link_c <- "identity"
  fit
}

log_lik_sdm_simple <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  c <- brms::get_dpar(prep, "c", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  y <- prep$data$Y[i]
  dsdm(y, mu, c, kappa, log = T)
}

posterior_predict_sdm_simple <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  c <- brms::get_dpar(prep, "c", i = i)
  kappa <- brms::get_dpar(prep, "kappa", i = i)
  rsdm(length(mu), mu, c, kappa)
}
