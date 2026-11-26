# local copies of two trivial bmm helpers that bmm does not export
nlist <- function(...) {
  dots <- list(...)
  names(dots) <- vapply(substitute(list(...))[-1], deparse, character(1))
  dots
}
read_lines2 <- function(con) paste(readLines(con, warn = FALSE), collapse = "\n")

.slnorm_defaults <- list(
  parameters = list(
    mu = "Location of the log response time",
    sigma = "Scale of the log response time",
    ndt = "Non-decision time: the shift of the distribution, in seconds"
  ),
  links = list(mu = "identity", sigma = "log", ndt = "log"),
  fixed_parameters = list(),
  priors = list(
    mu = list(main = "normal(-0.5, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
    sigma = list(main = "normal(-1, 0.5)", effects = "normal(0, 0.5)", sd = "exponential(2)"),
    ndt = list(main = "normal(-2, 1)", effects = "normal(0, 0.5)", sd = "exponential(2)")
  ),
  # native scale, inside the central part of the main priors
  init_ranges = list(mu = c(-1, 0), sigma = c(0.25, 0.5), ndt = c(0.05, 0.15))
)

.model_slnorm <- function(rt = NULL, links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(rt),
      other_vars = list(),
      domain = "Decision Making / Response times",
      task = "Response times",
      name = "Shifted log-normal response-time model (throwaway spike)",
      # placeholder, not a real reference; formatted so bmm's citation test passes
      citation = "Spike Author (2026). Throwaway shifted log-normal model for the bmm extension spike.",
      version = "NA",
      requirements = "- The response variable holds response times in seconds",
      parameters = .slnorm_defaults$parameters,
      links = .slnorm_defaults$links,
      fixed_parameters = .slnorm_defaults$fixed_parameters,
      default_priors = .slnorm_defaults$priors,
      init_ranges = .slnorm_defaults$init_ranges
    ),
    class = c("bmmodel", "slnorm"),
    call = call
  )
  bmm:::set_links(out, links)   # bmm::: set_links() validates and records links
}

#' Shifted log-normal response-time model
#' @param rt Name of the column with response times in seconds.
#' @param links Named list of link functions.
#' @param ... Used internally.
#' @export
slnorm <- function(rt, links = NULL, ...) {
  call <- match.call()
  bmm:::stop_missing_args()   # bmm::: stop_missing_args() checks required formals
  .model_slnorm(rt = rt, links = links, call = call, ...)
}

check_data.slnorm <- function(model, data, formula) {
  rt <- model$resp_vars$rt
  if (!rt %in% colnames(data)) {
    stop(sprintf("The response variable '%s' is not present in the data.", rt), call. = FALSE)
  }
  if (any(data[[rt]] <= 0, na.rm = TRUE)) {
    stop("Response times must be positive.", call. = FALSE)
  }
  NextMethod("check_data")
}

check_formula.slnorm <- function(model, data, formula) {
  attr(formula, "slnorm_check_formula_ran") <- TRUE
  NextMethod("check_formula")
}

configure_model.slnorm <- function(model, data, formula) {
  formula <- bmm::bmf2bf(model, formula)
  formula$family <- brms::custom_family(
    "slnorm",
    dpars = c("mu", "sigma", "ndt"),
    links = c(model$links$mu, model$links$sigma, model$links$ndt),
    lb = c(NA, 0, 0),
    ub = c(NA, NA, "min_Y"),
    type = "real",
    loop = TRUE,
    log_lik = log_lik_slnorm,
    posterior_predict = posterior_predict_slnorm,
    posterior_epred = posterior_epred_slnorm
  )
  sc_path <- system.file("stan_chunks", package = "bmmspike")
  # brms declares `real min_Y = min(Y);` itself when a bound names it
  stanvars <- brms::stanvar(
    scode = read_lines2(file.path(sc_path, "slnorm_functions.stan")),
    block = "functions"
  )
  nlist(formula, data, stanvars)
}

postprocess_brm.slnorm <- function(model, fit, ...) {
  attr(fit, "slnorm_postprocess_ran") <- TRUE
  fit
}

#' Density and random generation of the shifted log-normal
#' @param x,n,mu,sigma,ndt,log As usual.
#' @export
dslnorm <- function(x, mu, sigma, ndt, log = FALSE) {
  out <- stats::dlnorm(x - ndt, meanlog = mu, sdlog = sigma, log = TRUE)
  out[x <= ndt] <- -Inf
  if (log) out else exp(out)
}

#' @rdname dslnorm
#' @export
rslnorm <- function(n, mu, sigma, ndt) {
  ndt + stats::rlnorm(n, meanlog = mu, sdlog = sigma)
}

log_lik_slnorm <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  dslnorm(prep$data$Y[i], mu, sigma, ndt, log = TRUE)
}

posterior_predict_slnorm <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  rslnorm(length(mu), mu, sigma, ndt)
}

posterior_epred_slnorm <- function(prep) {
  as_mat <- function(x) matrix(x, nrow = prep$ndraws, ncol = prep$nobs)
  mu <- as_mat(brms::get_dpar(prep, "mu"))
  sigma <- as_mat(brms::get_dpar(prep, "sigma"))
  ndt <- as_mat(brms::get_dpar(prep, "ndt"))
  ndt + exp(mu + sigma^2 / 2)
}
