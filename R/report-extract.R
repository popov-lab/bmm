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
#' model_citation(sdt_rating(response = "rating", stimulus = "old", version = "dpsdt"))
#' @export
model_citation <- function(x, ...) {
  UseMethod("model_citation")
}

#' @export
model_citation.default <- function(x, ...) {
  stop2("model_citation() needs a bmmodel or bmmfit object, not an object of class '{class(x)[1]}'.")
}

#' @export
model_citation.bmmfit <- function(x, ...) {
  model_citation(x$bmm$model)
}

#' @export
model_citation.bmmodel <- function(x, ...) {
  refs <- as.character(current_constructor(x)$citation)
  if (!any(nzchar(refs))) {
    refs <- as.character(x$citation)
  }
  structure(refs[nzchar(refs)], uncited = uncited_part(x))
}

# the part of a model whose source bmm cannot know, because the user defines
# it; a method on the model's class, as the bmmodel class comes first in every
# model and would shadow a model_citation() method of its own
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

# the model as the installed version of bmm builds it; NULL when no constructor
# of this version matches the stored model's classes
current_constructor <- function(model) {
  name <- intersect(rev(class(model)), supported_models(print_call = FALSE))[1]
  if (is.na(name)) {
    return(NULL)
  }
  if (isTRUE(model$version %in% model_versions(name))) {
    get_model(name)(version = model$version)
  } else {
    get_model(name)()
  }
}
