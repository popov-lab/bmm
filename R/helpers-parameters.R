#' Transform kappa of the von Mises distribution to the circular standard
#' deviation
#' @description This function transforms the precision parameter kappa of the
#'   von Mises distribution to the circular standard deviation. Adapted from
#'   Matlab code by Paul Bays (https://www.paulbays.com/code.php)
#'
#' @param K numeric. A vector of kappa values.
#' @return A vector of sd values.
#' @keywords transform
#' @export
#' @examples
#' kappas <- runif(1000, 0.01, 100)
#'
#' # calcualte SD (in radians)
#' SDs <- k2sd(kappas)
#'
#' # transform SDs from radians to degrees
#' SDs_degress <- SDs * 180 / pi
#'
#' # plot the relationship between kappa and circular SD
#' plot(kappas, SDs)
#' plot(kappas, SDs_degress)
k2sd <- function(K) {
  log_bessel_ratio <- log(besselI(K, 1, expon.scaled = T)) - log(besselI(K, 0, expon.scaled = T))
  S <- sqrt(-2 * log_bessel_ratio)
  S[K == 0] <- Inf
  S[is.infinite(K)] <- 0
  S
}


#' Convert between parametrizations of the c parameter of the SDM distribution
#'
#' @name c_parametrizations
#' @inheritParams SDMdist
#' @return A numeric vector of the same length as `c` and `kappa`.
#' @details
#' `c_bessel2sqrtexp` converts the memory strength parameter (c)
#'   from the bessel parametrization to the sqrtexp parametrization,
#'   `c_sqrtexp2bessel` converts from the sqrtexp parametrization to the
#'   bessel parametrization.
#' @keywords transform
#' @details See [the online article](https://venpopov.com/bmm/articles/bmm_sdm_simple.html) for details on the
#'   parameterization. The sqrtexp parametrization is the default in the
#'   `bmm` package.
#' @export
#'
#' @examples
#' c_bessel <- c_sqrtexp2bessel(c = 4, kappa = 3)
#' c_sqrtexp <- c_bessel2sqrtexp(c = c_bessel, kappa = 3)
#'
c_sqrtexp2bessel <- function(c, kappa) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")
  c * besselI(kappa, 0, expon.scaled = TRUE) * sqrt(2 * pi * kappa)
}

#' @rdname c_parametrizations
#' @keywords transform
#' @export
c_bessel2sqrtexp <- function(c, kappa) {
  stopif(isTRUE(any(kappa < 0)), "kappa must be non-negative")
  stopif(isTRUE(any(c < 0)), "c must be non-negative")
  c / (besselI(kappa, 0, expon.scaled = TRUE) * sqrt(2 * pi * kappa))
}

#' @title Transform values from the native to parameter scale or vice versa according to a link function
#'
#' @description
#' This function transforms a vector of values from the native scale to the parameter scale,
#' according to the specified link function. The function is mainly used internally, to ensure proper
#' initial values
#'
#' @param values A numerical vector of values the transformation should be applied to
#' @param link A character specifying the link to be applied.
#' Available options are: "identity", "log", "log1p", "logm1", "inverse", "sqrt", "logit", "probit", "tan_half", and "cloglog".
#' @param inverse A Boolean value indicating if values should be transformed from the native to
#' the parameter scale (FALSE), or from the parameter scale to the native scale (TRUE)
#'
#' @noRd
link_transform <- function(values, link, inverse = FALSE) {
  # Handle NULL or missing link as identity (no transformation)
  if (is.null(link)) link <- "identity"

  stopifnot(is.numeric(values), is.character(link), length(link) == 1L, is.logical(inverse), length(inverse) == 1L)
  if(inverse) {
    switch(
      link,
      identity = values,
      log = exp(values),
      log1p = expm1(values),
      logm1 = brms::expp1(values),
      inverse = 1 / values,
      sqrt = values^2,
      logit = plogis(values),
      probit = pnorm(values),
      tan_half = 2 * atan(values),
      loglog = exp(-exp(values)),
      cloglog = 1 - exp(-exp(values)),
      stop2("Link not recognized.")
    )
  } else {
    switch(
      link,
      identity = values,
      log = log(values),
      log1p = log1p(values),
      logm1 = brms::logm1(values),
      inverse = 1 / values,
      sqrt = sqrt(values),
      logit = qlogis(values),
      probit = qnorm(values),
      tan_half = tan(values / 2),
      loglog = log(-log(values)),
      cloglog = log(-log1p(-values)),
      stop2("Link not recognized.")
    )
  }
}


#' Get parameter classification info for a bmmfit parameter
#'
#' @description
#' Determines whether a named model parameter is a distributional (`dpar`) or
#' non-linear (`nlpar`) parameter in the underlying brms model, retrieves its
#' link function, and checks for softmax transformation.
#'
#' @param bmmfit A bmmfit object
#' @param par Character string. Parameter name (e.g., `"kappa"`, `"c"`)
#'
#' @return A list with elements:
#'   \describe{
#'     \item{`type`}{Character: `"dpar"` or `"nlpar"`}
#'     \item{`model_name`}{Character: the parameter name as specified in the bmmodel}
#'     \item{`brms_name`}{Character: the parameter name as used in brms}
#'     \item{`link`}{Character: the link function (e.g., `"log"`, `"identity"`)}
#'     \item{`softmax`}{Logical: whether the parameter uses softmax transformation}
#'   }
#'
#' @keywords internal
#' @noRd
.get_parameter_info <- function(bmmfit, par) {
  model <- bmmfit$bmm$model
  model_pars <- names(model$parameters)
  bterms <- brms::brmsterms(bmmfit$formula)

  if (!par %in% model_pars) {
    stop2(
      "Parameter '{par}' not found in model.\n",
      "Available parameters: {paste(model_pars, collapse = ', ')}"
    )
  }

  link <- model$links[[par]] %||% "identity"
  softmax <- .is_softmax_param(par, model)

  if (!is.null(bterms$dpars) && par %in% names(bterms$dpars)) {
    type <- "dpar"
  } else if (!is.null(bterms$nlpars) && par %in% names(bterms$nlpars)) {
    type <- "nlpar"
  } else {
    # parameters not in bterms (e.g. fixed or auxiliary params) default to nlpar
    type <- "nlpar"
  }

  list(
    type = type,
    model_name = par,
    brms_name = par,
    link = link,
    softmax = softmax
  )
}


#' Check if parameter uses softmax transformation
#'
#' @param par Character string. Parameter name
#' @param model A bmmodel object
#' @return Logical. TRUE if parameter uses softmax transformation
#'
#' @keywords internal
#' @noRd
.is_softmax_param <- function(par, model) {
  "mixture3p" %in% class(model) && par %in% c("thetat", "thetant")
}
