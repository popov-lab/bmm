#' Restructure Old \code{bmmfit} Objects
#'
#' Restructure old \code{bmmfit} objects to work with
#' the latest \pkg{bmm} version. This function is called
#' internally when applying post-processing methods.
#'
#' @aliases restructure
#'
#' @param x An object of class \code{bmmfit}.
#' @param ... Currently ignored.
#'
#' @return A \code{bmmfit} object compatible with the latest version
#'   of \pkg{bmm} and \pkg{brms}.
#' @keywords transform
#' @export
#' @examplesIf file.exists("bmmfit_old.rds")
#' # Load an old bmmfit object
#' old_fit <- readRDS("bmmfit_old.rds")
#' new_fit <- restructure(old_fit)
#' @importFrom utils packageVersion
restructure.bmmfit <- function(x, ...) {
  version <- x$version$bmm
  if (is.null(version)) {
    version <- as.package_version("0.2.1")
    x$version$bmm <- version
  }
  if (!inherits(x, "bmmfit")) {
    class(x) <- c("bmmfit", class(x))
  }
  current_version <- packageVersion("bmm")
  restr_version <- restructure_version.bmm(x)
  # not gated by version: development fits carry the current version number
  # and can still predate the function
  x <- add_posterior_epred(x)

  if (restr_version >= current_version) {
    x <- NextMethod("restructure")
    return(x)
  }

  if (restr_version < "0.3.0") {
    class(x) <- c("bmmfit", class(x))
    x <- add_bmm_info(x)
  }

  if (restr_version < "0.4.3") {
    info <- x$bmm$model$info
    x$bmm$model$info <- NULL
    x$bmm$model[names(info)] <- info
    class(x$bmm$model)[1] <- "bmmodel"
    cl_len <- length(class(x$bmm$model))
    class(x$bmm$model)[cl_len] <- switch(class(x$bmm$model)[cl_len],
       sdmSimple = "sdm",
       IMMfull = "imm",
       IMMabc = "imm",
       IMMbsc = "imm",
       class(x$bmm$model)[cl_len]
    )
    x <- add_links(x)
    x$bmm$user_formula <- assign_nl_attr(x$bmm$user_formula)
  }

  if (restr_version < "0.4.4") {
    x$bmm$fit_args <- NULL
  }

  if (restr_version < "0.4.5") {
    prev_model <- x$bmm$model
    prev_model_name <- class(prev_model)[length(class(prev_model))]
    new_model <- get_model(prev_model_name)
    new_model <- brms::do_call(new_model, c(prev_model$resp_vars, prev_model$other_vars))
    x$bmm$model <- check_model(new_model, x$data, x$bmm$user_formula)
  }

  if (restr_version < "0.4.8") {
    cl <- class(x$bmm$model)
    cl[1] <- "bmmodel"
    if ("sdmSimple" %in% cl) {
      cl[cl == "sdmSimple"] <- "sdm"
      cl <- c(cl, "sdm_simple")
    } else if ("IMMfull" %in% cl) {
      cl[cl == "IMMfull"] <- "imm"
      cl <- c(cl, "imm_full")
    } else if ("IMMabc" %in% cl) {
      cl[cl == "IMMabc"] <- "imm"
      cl <- c(cl, "imm_abc")
    } else if ("IMMbsc" %in% cl) {
      cl[cl == "IMMbsc"] <- "imm"
      cl <- c(cl, "imm_bsc")
    }
    if ("nontargets" %in% cl) {
      cl[cl == "nontargets"] <- "non_targets"
      cl <- c(cl, "non_targets")
    }
    class(x$bmm$model) <- cl
  }

  x$version$bmm_restructure <- current_version
  NextMethod("restructure")
}

restructure_version.bmm <- function(x) {
  out <- x$version$bmm_restructure
  if (!is.package_version(out)) {
    out <- x$version$bmm
  }
  out
}

# Custom families saved without a posterior_epred function (#475). Without
# one, brms looks up posterior_epred_<family name> in the family environment,
# which reset_env() points at the global environment; bmm defines its
# functions under those names, so the same lookup is done here in bmm.
add_posterior_epred <- function(x) {
  family <- x$formula$family
  if (!inherits(family, "customfamily") || is.function(family$posterior_epred)) {
    return(x)
  }
  fun <- get0(paste0("posterior_epred_", family$name), envir = asNamespace("bmm"),
              mode = "function", inherits = FALSE)
  if (is.null(fun)) {
    return(x)
  }
  x$formula$family$posterior_epred <- fun
  if (inherits(x$family, "customfamily")) {
    x$family$posterior_epred <- fun
  }
  x
}

# Before bmm 1.5.0, mixture2p, mixture3p and imm were brms mixture families
# with non-linear parameters, and their location was mu1. restructure() leaves
# such a fit as it was fitted, because the extractors read the fitted
# parameters (mu1, mu2, kappa2) from its model object. update() refits it under
# the custom family instead, from a model object rebuilt by the current
# constructor. Told apart by the family rather than by the version, because
# development builds before the custom families carry the same version number
is_mixture_family_fit <- function(x) {
  any(c("mixture2p", "mixture3p", "imm") %in% class(x$bmm$model)) &&
    inherits(x$formula$family, "mixfamily")
}

upgrade_mixture_family_model <- function(x) {
  model <- x$bmm$model
  name <- intersect(c("mixture2p", "mixture3p", "imm"), class(model))
  # fits from 0.4.5 to 0.4.7 have non_targets after the version class
  version <- if (name == "imm") {
    sub("^imm_", "", grep("^imm_(full|bsc|abc)$", class(model), value = TRUE)[1])
  } else {
    "simple"
  }
  new <- brms::do_call(get_model(name), c(
    model$resp_vars, model$other_vars,
    list(version = version, regex = isTRUE(attr(model, "regex")), call = attr(model, "call"))
  ))
  # the links stay the constructor's: bmm 1.4 refused other links for these
  # models, and a stored link means what it meant in the release that wrote it
  # (0.5.3 stored thetat as identity, with the logit inside the formula)
  for (a in c("regex_vars", "regex_applied", "user_vars")) {
    attr(new, a) <- attr(model, a)
  }
  x$bmm$model <- update_model_fixed_parameters(new, x$bmm$user_formula)
  x
}

# The prior stored on such a fit holds what bmm passed to brms: the user's
# priors and bmm's own defaults, written for non-linear parameters, plus rows
# for the second mixture component (mu2, kappa2) and the constants that pinned
# the set-size-1 level of mixture3p's thetant and imm's a and s. The custom
# families have no second component, and since 1.5.0 the set-size-1 level
# samples its prior with a warning instead of being pinned (D6), so those rows
# are dropped. A constant the user wrote in exactly the form bmm generated
# cannot be told apart from it and is dropped too
translate_stored_prior <- function(old_model, model, formula, prior, fitted_with) {
  bterms <- brms::brmsterms(formula)
  component <- nzchar(prior$dpar) &
    not_in(prior$dpar, c(names(bterms$dpars), names(model$deprecated_parameters)))
  pinned <- set_size_one_constants(model, prior)
  translated <- translate_nlpar_priors(model, formula, prior[!component & !pinned, ])
  name <- intersect(c("mixture2p", "mixture3p", "imm"), class(model))
  pinned_pars <- unique(prior$nlpar[pinned])
  sentences <- c(
    glue(
      "This fit was made with bmm {fitted_with}, when {name}() was built on a brms \\
      mixture family. update() samples it again under the custom family of bmm \\
      {utils::packageVersion('bmm')} and translates its priors from `nlpar =` to \\
      `dpar =`, so your call needs no change."
    ),
    translated$details,
    if (length(pinned_pars) > 0) {
      glue(
        "The constant priors that bmm {fitted_with} set on the set-size-1 level of \\
        {collapse_comma(pinned_pars)} are dropped, so bmm no longer fixes that level."
      )
    },
    if ("mu1" %in% names(old_model$parameters) && not_in("mu1", names(old_model$fixed_parameters))) {
      glue(
        "In the refitted fit, the draws of mu1 are named after mu, e.g. b_Intercept \\
        instead of b_mu1_Intercept."
      )
    }
  )
  list(prior = translated$prior, message = paste(sentences, collapse = " "))
}

# the rows bmm 1.4 generated to pin the set-size-1 level (constrain_set_size1_*)
set_size_one_constants <- function(model, prior) {
  pars <- if (inherits(model, "mixture3p")) "thetant" else if (inherits(model, "imm")) c("a", "s")
  value <- if (inherits(model, "mixture3p")) "constant(-100)" else "constant(0)"
  set_size <- model$other_vars$set_size
  if (length(pars) == 0 || !is.character(set_size)) {
    return(rep(FALSE, nrow(prior)))
  }
  prior$nlpar %in% pars & prior$coef == paste0(set_size, 1) &
    ((prior$class == "b" & prior$prior == value) |
      (prior$class == "sd" & prior$prior == "constant(1e-8)"))
}

add_links <- function(x) {
  UseMethod("add_links")
}

#' @export
add_links.bmmfit <- function(x) {
  x$bmm$model <- add_links(x$bmm$model)
  x
}

#' @export
add_links.bmmodel <- function(x) {
  # the last class of a versioned model is <model>_<version>, which has no
  # constructor of its own, so walk back to the most specific class that does
  model_name <- Find(
    function(name) exists(paste0(".model_", name), mode = "function"),
    rev(class(x))
  )
  new_model <- get_model(model_name)()
  x$links <- new_model$links
  x
}

add_bmm_info <- function(x) {
  env <- x$family$env
  stopif(
    is.null(env),
    "Unable to restructure the object for use with the latest version of bmm. Please refit."
  )
  pforms <- env$formula$pforms
  names(pforms) <- NULL
  user_formula <- brms::do_call("bmf", pforms)
  model <- env$model
  model$resp_vars <- list(resp_error = env$formula$resp)
  model$other_vars <- list()
  if (inherits(model, "sdm")) {
    model$info$parameters$mu <- glue(
      "Location parameter of the SDM distribution \\
      (in radians; by default fixed internally to 0)"
    )
  } else {
    model$info$parameters$mu1 <- glue(
      "Location parameter of the von Mises distribution for memory responses \\
       (in radians). Fixed internally to 0 by default."
    )
  }

  x$bmm <- nlist(model, user_formula)
  x
}
