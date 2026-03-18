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
        thetat = "identity"
      ),
      fixed_parameters = list(mu1 = 0, mu2 = 0, kappa2 = -100),
      default_priors = list(
        mu1 = list(main = "student_t(1, 0, 1)", effects = "normal(0,1)"),
        kappa = list(main = "normal(2, 1)", effects = "normal(0, 1)",
                     sd = "exponential(1)"),
        thetat = list(main = "logistic(0, 1)", effects = "normal(0, 0.5)",
                      sd = "exponential(1)")
      ),
      void_mu = FALSE
    ),
    class = c("bmmodel", "circular", "mixture2p"),
    call = call
  )
  out$links[names(links)] <- links
  out
}

# user facing alias

#' @title `r .model_mixture2p()$name`
#' @details `r model_info(.model_mixture2p())`
#' @param resp_error The name of the variable in the provided dataset containing
#'   the response error. The response Error should code the response relative to
#'   the to-be-recalled target in radians. You can transform the response error
#'   in degrees to radian using the `deg2rad` function.
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
mixture2p <- function(resp_error, ...) {
  call <- match.call()
  stop_missing_args()
  .model_mixture2p(resp_error = resp_error, call = call, ...)
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
