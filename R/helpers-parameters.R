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


#' Get parameter information for a bmm model
#'
#' @description Returns a data frame with information about the model
#'   parameters, including their descriptions, whether they are fixed,
#'   their link functions, and optionally their default priors.
#'
#' @param x A \code{bmmodel} object (e.g., \code{sdm(resp_error = "y")}) or a
#'   \code{bmmfit} object (a fitted model returned by \code{\link{bmm}}).
#' @param formula An optional \code{bmmformula} object. Only relevant for
#'   M3 custom models, where additional parameters are discovered from
#'   the formula. Ignored for all other models.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame of class \code{bmm_parameters} with one row per
#'   parameter and columns: \code{parameter}, \code{description},
#'   \code{fixed}, \code{value}, and \code{link}.
#'
#' @export
#' @examples
#' # For an unfitted model
#' parameters(sdm(resp_error = "y"))
#'
#' # For an M3 model
#' parameters(m3(
#'   resp_cats = c("corr", "other", "npl"),
#'   num_options = c(1, 4, 5),
#'   version = "ss"
#' ))
parameters <- function(x, ...) {
  UseMethod("parameters")
}


#' @rdname parameters
#' @export
parameters.bmmodel <- function(x, formula = NULL, ...) {
  model <- x

  if (inherits(model, "m3_custom") && !is.null(formula)) {
    user_pars <- rhs_vars(formula[is_nl(formula)])
    user_pars <- setdiff(user_pars, names(formula[is_nl(formula)]))
    user_pars <- setdiff(user_pars, names(model$parameters))
    model$parameters <- c(model$parameters, setNames(
      as.list(user_pars), user_pars
    ))
  }

  pars <- names(model$parameters)
  if (length(pars) == 0) {
    message("This model has no parameters defined.")
    return(invisible(data.frame()))
  }

  fixed <- pars %in% names(model$fixed_parameters)
  values <- rep(NA_character_, length(pars))
  values[fixed] <- as.character(model$fixed_parameters[pars[fixed]])

  links <- vapply(pars, function(p) {
    model$links[[p]] %||% "identity"
  }, character(1))

  out <- data.frame(
    parameter = pars,
    description = as.character(model$parameters),
    fixed = fixed,
    value = values,
    link = links,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  m3_note <- NULL
  if (inherits(model, "m3_custom") && is.null(formula)) {
    m3_note <- paste(
      "Note: This is a custom M3 model. Pass a formula to",
      "discover additional parameters."
    )
  }

  class(out) <- c("bmm_parameters", "data.frame")
  attr(out, "model_name") <- model$name
  attr(out, "m3_note") <- m3_note
  out
}


#' @rdname parameters
#' @export
parameters.bmmfit <- function(x, ...) {
  x <- restructure(x)
  parameters(x$bmm$model, formula = x$bmm$user_formula, ...)
}


#' @export
print.bmm_parameters <- function(x, max_desc_width = 50, ...) {
  model_name <- attr(x, "model_name")
  m3_note <- attr(x, "m3_note")

  if (!is.null(model_name) && nzchar(model_name)) {
    cat(style("purple1")("Model: "), model_name, "\n\n")
  }

  print_df <- x
  if ("description" %in% names(print_df)) {
    print_df$description <- vapply(print_df$description, function(d) {
      d <- gsub("\\s+", " ", d)
      if (nchar(d) > max_desc_width) {
        paste0(substr(d, 1, max_desc_width - 3), "...")
      } else {
        d
      }
    }, character(1))
  }

  for (col in names(print_df)) {
    if (is.character(print_df[[col]])) {
      print_df[[col]][is.na(print_df[[col]])] <- "--"
    }
  }

  if ("fixed" %in% names(print_df)) {
    print_df$fixed <- ifelse(print_df$fixed, "yes", "no")
  }

  print.data.frame(print_df, right = FALSE, row.names = FALSE)

  if (!is.null(m3_note)) {
    cat("\n", m3_note, "\n")
  }

  invisible(x)
}
