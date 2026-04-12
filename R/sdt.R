############################################################################# !
# SDT MODEL WRAPPER                                                      ####
############################################################################# !

.sdt_version_aliases <- c(
  auto = "auto",
  binary = "binary",
  rating = "rating",
  dp = "dp",
  dpsdt = "dp",
  metad = "metad",
  meta_d = "metad",
  `meta-d` = "metad",
  cdp = "cdp",
  ranking = "ranking",
  mafc = "mafc",
  m_afc = "mafc"
)

.match_sdt_version <- function(version, allow_auto = TRUE) {
  stopif(length(version) != 1 || !is.character(version),
    "`version` must be a single character string"
  )

  key <- .sdt_version_aliases[[version]]
  stopif(is.null(key),
    paste0(
      "Unknown SDT version '{version}'. Supported values are ",
      "{collapse_comma(unique(unname(.sdt_version_aliases)))}"
    )
  )
  stopif(!allow_auto && key == "auto",
    "`version = 'auto'` is only supported for `sdt()`"
  )
  key
}

.sdt_model_dispatch <- function(version) {
  switch(version,
    binary = sdt_binary,
    rating = sdt_rating,
    dp = sdt_dp,
    metad = sdt_metad,
    cdp = sdt_cdp,
    ranking = sdt_ranking,
    mafc = sdt_mafc
  )
}

.infer_sdt_model_version <- function(args) {
  if (all(c("new_response", "old_know", "old_remember") %in% names(args))) {
    return("cdp")
  }
  if ("rank" %in% names(args)) {
    return("ranking")
  }
  if ("m" %in% names(args) && "n_trials" %in% names(args) &&
      !"stimulus" %in% names(args)) {
    return("mafc")
  }
  if ("response" %in% names(args)) {
    if (length(args$response) > 1) {
      return("rating")
    }
    if (all(c("stimulus", "n_trials") %in% names(args))) {
      return("binary")
    }
  }

  stop2(
    "Could not infer the SDT model type from the supplied arguments. ",
    "Set `version` explicitly."
  )
}

#' @title Unified SDT Model Constructor
#' @description Convenience wrapper that dispatches to the specialized SDT
#'   model constructors. Existing constructors such as [sdt_binary()] and
#'   [sdt_rating()] remain fully supported.
#' @param ... Arguments forwarded to the selected SDT model constructor.
#' @param version Character. Which SDT model to construct. Use `"auto"`
#'   (default) to infer the model from the supplied arguments.
#' @return An object of class `bmmodel`.
#' @keywords bmmodel
#' @export
sdt <- function(..., version = "auto") {
  args <- list(...)
  version <- .match_sdt_version(version, allow_auto = TRUE)
  if (version == "auto") {
    version <- .infer_sdt_model_version(args)
  }
  do.call(.sdt_model_dispatch(version), args)
}
