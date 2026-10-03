#' References for a measurement model
#'
#' Returns the published source(s) of a model, for the reference list of a
#' paper that uses it.
#'
#' @param x A model object (for example `sdm(resp_error = "y")`) or a `bmmfit`
#'   object returned by [bmm()].
#' @param ... Currently ignored.
#'
#' @details The references come from the model constructor of the installed
#'   version of bmm, so a fit made with an older version gets the current,
#'   corrected references. When the installed version no longer has a
#'   constructor for the model, or its references are empty, the copy stored
#'   in the fit or model object is returned instead.
#'
#'   Some versions of a model add the reference that introduced them, for
#'   example `sdt_rating(version = "dpsdt")`.
#'
#'   A custom `m3()` model is cited with the paper that introduced the M3
#'   framework. The model structure itself is defined by the user, so its
#'   source cannot be known to bmm. The returned vector then has the attribute
#'   `uncited`, naming the part of the model you have to cite yourself.
#'
#'   bmm itself is cited separately, see `citation("bmm")`.
#'
#' @return A character vector with one reference in APA style per element,
#'   or `character(0)` if no reference is available.
#' @seealso [citation()]
#' @keywords extract_info
#' @examples
#' model_citation(sdm(resp_error = "y"))
#' model_citation(
#'   sdt_rating(response = paste0("r", 1:4), stimulus = "old", version = "dpsdt")
#' )
#' @export
model_citation <- function(x, ...) {
  UseMethod("model_citation")
}

#' @export
model_citation.default <- function(x, ...) {
  stop2(
    "model_citation() needs a bmmodel or bmmfit object, \\
    not an object of class '{class(x)[1]}'."
  )
}

#' @export
model_citation.bmmfit <- function(x, ...) {
  model_citation(x$bmm$model)
}

#' @export
model_citation.bmmodel <- function(x, ...) {
  refs <- as.character(current_constructor(x)$citation)
  if (!any(nzchar(refs))) {
    refs <- split_stored_citation(x$citation)
  }
  structure(refs[nzchar(refs)], uncited = uncited_part(x))
}

# models stored before bmm 1.4.0 hold all references in one string, joined by
# a newline and "- ", with line breaks inside a reference
split_stored_citation <- function(citation) {
  trimws(gsub("\\s*\n\\s*", " ", unlist(strsplit(as.character(citation), "\n\\s*-\\s+"))))
}

# the part of a model whose source bmm cannot know, because the user defines
# it; a method on the model's class, as the bmmodel class comes first in every
# model and would shadow a model_citation() method of its own. Code that needs
# this (report_methods() in #432) calls the generic directly, because `[`,
# c() and unique() drop the uncited attribute of model_citation()
uncited_part <- function(model) {
  UseMethod("uncited_part")
}

#' @export
uncited_part.default <- function(model) {
  NULL
}

#' @export
uncited_part.m3_custom <- function(model) {
  "the user-defined model structure"
}

# Sampler settings and package versions as the fit stores them. The versions
# describe the compilation: update() without recompiling keeps the stored brms
# and Stan versions, while the settings and the date describe the last run.
# brms stores no R version at all. ndraws_stored counts the post-warmup draws
# that summary() and as_draws_array() see; brms::ndraws() reads n_save instead,
# which can disagree with the stored draws
fit_settings <- function(fit) {
  sim <- if (has_draws(fit)) fit$fit@sim else list()
  version <- fit$version
  backend <- fit$backend %||% NA_character_
  list(
    backend = backend,
    algorithm = fit$algorithm %||% NA_character_,
    chains = first_or_na(sim[["chains"]]),
    iter = first_or_na(sim[["iter"]]),
    warmup = first_or_na(sim[["warmup"]]),
    thin = first_or_na(sim[["thin"]]),
    ndraws_stored = if (has_draws(fit)) {
      sum(vapply(sim$samples, function(chain) length(chain[[1]]), numeric(1)) - sim$warmup2)
    } else {
      NA_real_
    },
    versions = c(
      bmm = version_string(version$bmm),
      brms = version_string(version$brms),
      stan = stan_version(fit),
      backend = version_string(if (backend %in% c("cmdstanr", "rstan")) version[[backend]])
    ),
    date = if (has_draws(fit)) fit$fit@date else NA_character_
  )
}

# rstan fits store the version of the StanHeaders R package, whose patch level
# is not Stan's; the stanc3 version that compiled the model is recorded in the
# model's C++, and older compilers do not write it there
stan_version <- function(fit) {
  switch(fit$backend %||% "",
    cmdstanr = version_string(fit$version$cmdstan),
    rstan = stanc_version(fit$fit@stanmodel@model_cpp$model_cppcode) %||%
      version_string(fit$version$stanHeaders),
    NA_character_
  )
}

stanc_version <- function(cpp) {
  match <- regmatches(cpp, regexec("stanc_version = stanc3 v([0-9][0-9.]*)", cpp))
  if (length(match) && length(match[[1]])) match[[1]][2]
}

# R-hat and bulk and tail ESS per variable, without lp__, lprior and constant
# variables (R-hat NA). The caller decides which parameter classes to report
convergence_summary <- function(fit) {
  # without draws, brms fails further down with an error about `@` applied to
  # a number, which does not say that the fit is a mock
  stopif(!has_draws(fit), "The fit contains no posterior draws (a mock fit?).")
  algorithm <- fit$algorithm
  if (!identical(algorithm, "sampling")) {
    out <- data.frame(
      parameter = character(), rhat = numeric(), ess_bulk = numeric(), ess_tail = numeric()
    )
    return(structure(out, algorithm = algorithm))
  }
  conv <- posterior::summarise_draws(
    brms::as_draws_array(fit), "rhat", "ess_bulk", "ess_tail"
  )
  conv <- conv[!conv$variable %in% c("lp__", "lprior") & !is.na(conv$rhat), ]
  max_treedepth <- brms::control_params(fit)$max_treedepth
  treedepth <- brms::nuts_params(fit, pars = "treedepth__")$Value
  structure(
    data.frame(
      parameter = conv$variable, rhat = conv$rhat,
      ess_bulk = conv$ess_bulk, ess_tail = conv$ess_tail
    ),
    algorithm = algorithm,
    divergent = sum(brms::nuts_params(fit, pars = "divergent__")$Value),
    max_treedepth_hits = if (is.null(max_treedepth)) NA_real_ else sum(treedepth >= max_treedepth)
  )
}

# a stanfit without draws comes from chains = 0 or empty = TRUE
has_draws <- function(fit) {
  methods::is(fit$fit, "stanfit") && length(fit$fit@sim) > 0
}

first_or_na <- function(x) {
  if (length(x)) as.numeric(x[[1]]) else NA_real_
}

version_string <- function(x) {
  if (is.null(x)) NA_character_ else as.character(x)
}

# the model as the installed version of bmm builds it; NULL when no current
# constructor matches the stored model's classes and version. A stored model
# without a version predates versions, and an unversioned model stores "NA";
# both get the default version, should the model have gained versions since
current_constructor <- function(model) {
  name <- intersect(rev(class(model)), supported_models(print_call = FALSE))[1]
  if (is.na(name)) {
    return(NULL)
  }
  versions <- model_versions(name)
  if (all(model$version %in% c(NA, "NA", "")) || all(is.na(versions))) {
    return(get_model(name)())
  }
  if (model$version %in% versions) get_model(name)(version = model$version)
}
