############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_sdm <- function(resp_error = NULL, links = NULL, version = "simple", call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      other_vars = nlist(),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Signal Discrimination Model (SDM) by Oberauer (2023)",
      citation = glue(
        "Oberauer, K. (2023). Measurement models for visual working memory - \\
        A factorial model comparison. Psychological Review, 130(3), 841-852"
      ),
      version = version,
      requirements = glue(
        "- The response variable should be in radians and represent the angular \\
        error relative to the target"
      ),
      parameters = list(
        mu = glue("Location parameter of the SDM distribution (in radians; \\
                  by default fixed internally to 0)"),
        c = "Memory strength parameter of the SDM distribution",
        kappa = "Precision parameter of the SDM distribution"
      ),
      links = list(
        mu = "tan_half",
        c = "log",
        kappa = "log"
      ),
      fixed_parameters = list(mu = 0),
      default_priors = list(
        mu = list(main = "student_t(1, 0, 1)"),
        kappa = list(main = "student_t(5, 1.75, 0.75)", effects = "normal(0, 1)"),
        c = list(main = "student_t(5, 2, 0.75)", effects = "normal(0, 1)")
      ),
      init_ranges = list(
        mu = c(-0.5,0.5),
        kappa = c(2.5,3.5),
        c = c(4,6)
      ),
      void_mu = FALSE
    ),
    class = c("bmmodel", "circular", "sdm", paste0("sdm_", version)),
    call = call
  )
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
sdm <- function(resp_error, version = "simple", ...) {
  call <- match.call()
  stop_missing_args()
  .model_sdm(resp_error = resp_error, version = version, call = call, ...)
}

#' @rdname sdm
#' @keywords deprecated
#' @export
sdmSimple <- function(resp_error, version = "simple", ...) {
  warning2("The function `sdmSimple()` is deprecated. Please use `sdm()` instead.")
  call <- match.call()
  stop_missing_args()
  .model_sdm(resp_error = resp_error, version = version, call = call, ...)
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.sdm <- function(model, data, formula) {
  # data sorted by predictors is necessary for speedy computation of normalizing constant
  data <- order_data_query(model, data, formula)
  attr(data, "sdm_run_metadata") <- sdm_run_metadata(data, formula)
  NextMethod("check_data")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.sdm <- function(model, data, formula) {
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
  run_metadata <- attr(data, "sdm_run_metadata")
  if (is.null(run_metadata)) {
    run_metadata <- sdm_run_metadata(data, formula)
  }
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(scode = stan_tdata, block = "tdata") +
    brms::stanvar(x = run_metadata$G_sdm_runs, name = "G_sdm_runs") +
    sdm_stanvar_int_array(run_metadata$sdm_run_start, "sdm_run_start", "G_sdm_runs") +
    sdm_stanvar_int_array(run_metadata$sdm_run_count, "sdm_run_count", "G_sdm_runs") +
    brms::stanvar(scode = stan_likelihood, block = "likelihood", position = "end")

  # construct main brms formula from the bmm formula
  formula <- bmf2bf(model, formula)
  formula$family <- sdm_simple

  # set initial values to be sampled between [-1,1] to avoid extreme SDs that
  # can cause the sampler to fail
  init <- 1

  nlist(formula, data, stanvars, init)
}

############################################################################# !
# POSTPROCESS METHODS                                                    ####
############################################################################# !

#' @export
postprocess_brm.sdm <- function(model, fit, ...) {
  # manually set link_c to "log" since I coded it manually
  fit$family$link_c <- "log"
  fit$formula$family$link_c <- "log"
  fit
}

#' @export
revert_postprocess_brm.sdm <- function(model, fit, ...) {
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

sdm_run_metadata <- function(data, formula) {
  predictors <- sdm_predictor_vars(data, formula)
  if (length(predictors) == 0L) {
    return(list(
      G_sdm_runs = 1L,
      sdm_run_start = 1L,
      sdm_run_count = nrow(data)
    ))
  }

  run_id <- interaction(data[predictors], drop = TRUE, lex.order = TRUE)
  run_start <- c(1L, which(run_id[-1] != run_id[-length(run_id)]) + 1L)
  run_count <- diff(c(run_start, nrow(data) + 1L))
  list(
    G_sdm_runs = length(run_start),
    sdm_run_start = as.integer(run_start),
    sdm_run_count = as.integer(run_count)
  )
}

sdm_predictor_vars <- function(data, formula) {
  dpars <- names(formula)
  predictors <- rhs_vars(formula)
  predictors <- predictors[not_in(predictors, dpars)]
  predictors[predictors %in% colnames(data)]
}

sdm_stanvar_int_array <- function(x, name, size) {
  out <- brms::stanvar(x = as.integer(x), name = name)
  out[[name]]$scode <- paste0("array[", size, "] int ", name, ";")
  out[[name]]$pll_args <- paste("data array[] int", name)
  out
}
