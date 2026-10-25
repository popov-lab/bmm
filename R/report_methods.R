#' Draft a methods paragraph for a fitted model
#'
#' Drafts the text that describes a model fitted with [bmm()] in the methods
#' section of a paper: the model, the software and sampler settings, the
#' priors, and the convergence diagnostics. The draft is meant to be edited,
#' not pasted unread.
#'
#' @param fit A `bmmfit` object returned by [bmm()]. It must contain posterior
#'   draws, so mock fits and fits made with `chains = 0` or `empty = TRUE` are
#'   refused.
#' @param sections The parts of the paragraph to draft, in the order they
#'   should appear: `"model"` (the model, the versions of bmm, brms, Stan, the
#'   backend and R, and the sampler settings), `"priors"` (the sentences of
#'   [prior_info()] with `format = "text"`) and `"convergence"` (the sentence
#'   of [report_convergence()]).
#'
#' @details Every value in the text is read from the fit, and every reference
#'   comes from the model's constructor (see [model_citation()]) or from
#'   [citation()] of an installed package. Nothing else enters the text:
#'
#'   - The versions are those the fit stores. The R version is stored by fits
#'     made with bmm 1.5.0 or later; for older fits it is left out and a note
#'     asks you to add it.
#'   - The references are those of the packages installed now. A note says
#'     when an installed version differs from the version the fit was made
#'     with.
#'   - The convergence sentence gives numbers, not verdicts: the largest
#'     R-hat and the smallest bulk and tail effective sample sizes over all
#'     estimated parameters (constants, `lp__` and `lprior` excluded), and the
#'     number of divergent transitions after warmup. Where the summary of
#'     \pkg{brms} would warn (an R-hat above 1.05 among the parameters it shows,
#'     or a divergent transition), a note says so.
#'   - For algorithms other than `"sampling"`, the convergence sentence is
#'     replaced by a note, as there are no Markov chains to diagnose.
#'
#'   In the text, the numbers in square brackets refer to the reference list.
#'   A model without a reference gets `[CITATION NEEDED]` in that list.
#'
#' @return An object of class `bmm_report_methods`, a list with
#'   - `paragraph`: the drafted text, one string;
#'   - `formula`: the bmm formula with its constants, one line per parameter
#'     (empty without the `"model"` section);
#'   - `references`: the references, numbered by position;
#'   - `notes`: what you have to check or add, not meant for the paper;
#'   - `sections`: the sections drafted, in order.
#'
#'   `print()` shows all of it; `format()` returns the paragraph, the formula
#'   and the numbered references as one string to paste into a document. The
#'   references keep the Markdown emphasis marks that [citation()] writes
#'   (`_journal_`, `*volume*`): Quarto renders them as italics; in Word, remove
#'   them.
#' @seealso [report_convergence()], [prior_info()], [model_citation()]
#' @keywords extract_info
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' fit <- bmm(
#'   bmf(c ~ 0 + set_size, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' report_methods(fit)
#' report_methods(fit, sections = c("model", "convergence"))
#' cat(format(report_methods(fit)))
#' @export
report_methods <- function(fit, sections = c("model", "priors", "convergence")) {
  stopif(
    !inherits(fit, "bmmfit"),
    "report_methods() needs a bmmfit object returned by bmm(), not an object of class '{class(fit)[1]}'."
  )
  stopif(!has_draws(fit), no_draws_message())
  stopif(
    !length(sections) || anyNA(sections),
    "sections must name at least one of 'model', 'priors' or 'convergence'."
  )
  sections <- unique(match.arg(sections, several.ok = TRUE))
  fit <- restructure(fit)
  references <- report_references(fit)
  markers <- reference_markers(references)
  conv <- if ("convergence" %in% sections) convergence_summary(fit)
  texts <- vapply(sections, function(section) {
    switch(section,
      model = model_report_text(fit, markers),
      priors = prior_report_text(prior_info(fit, format = "text")),
      convergence = if (identical(attr(conv, "algorithm"), "sampling")) convergence_report_text(conv) else ""
    )
  }, character(1))
  structure(
    list(
      paragraph = paste(texts[nzchar(texts)], collapse = " "),
      formula = if ("model" %in% sections) {
        summarise_formula(fit$bmm$user_formula, fit$bmm$model)
      } else {
        character(0)
      },
      references = unlist(references, use.names = FALSE),
      notes = report_notes(fit, sections, references, markers, conv),
      sections = sections
    ),
    class = "bmm_report_methods"
  )
}

#' Convergence diagnostics of a fitted model
#'
#' Returns the convergence diagnostics of a model fitted with [bmm()], per
#' parameter or as the sentence that [report_methods()] uses.
#'
#' @inheritParams report_methods
#' @param format `"table"` for one row per parameter, `"text"` for the
#'   sentence.
#'
#' @details The table holds the rank-normalized R-hat and the bulk and tail
#'   effective sample sizes (ESS) of every estimated parameter, computed with
#'   [posterior::summarise_draws()]. Constants, `lp__`, `lprior` and the draws
#'   from the prior (`prior_*`) are left out. The column `in_summary` marks
#'   the parameters that `summary()` shows; \pkg{brms} warns when one of them
#'   has an R-hat above 1.05. The algorithm, the number of divergent
#'   transitions after warmup and the number of transitions that reached the
#'   maximum tree depth are attributes of the table, printed below it.
#'
#'   The sentence gives the largest R-hat, the smallest bulk and tail ESS and
#'   the number of divergent transitions, without verdicts. For algorithms
#'   other than `"sampling"` the table is empty and the text says that there
#'   are no diagnostics.
#'
#' @return An object of class `bmm_report_convergence`: a `data.frame` with
#'   columns `parameter`, `rhat`, `ess_bulk`, `ess_tail` and `in_summary` for
#'   `format = "table"`, a character string for `format = "text"`. Subsetting
#'   the table with `[` returns a plain `data.frame`.
#' @seealso [report_methods()]
#' @keywords extract_info
#' @examplesIf identical(Sys.getenv("BMM_EXAMPLES"), "true")
#' fit <- bmm(
#'   bmf(c ~ 0 + set_size, kappa ~ 1),
#'   data = oberauer_lin_2017,
#'   model = sdm(resp_error = "dev_rad")
#' )
#' report_convergence(fit)
#' report_convergence(fit, format = "text")
#' @export
report_convergence <- function(fit, format = c("table", "text")) {
  stopif(
    !inherits(fit, "bmmfit"),
    "report_convergence() needs a bmmfit object returned by bmm(), not an object of class '{class(fit)[1]}'."
  )
  stopif(!has_draws(fit), no_draws_message())
  format <- match.arg(format)
  conv <- convergence_summary(restructure(fit))
  if (format == "text") {
    return(structure(convergence_report_text(conv), class = c("bmm_report_convergence", "character")))
  }
  structure(conv, class = c("bmm_report_convergence", "data.frame"))
}

no_draws_message <- function() {
  "The fit contains no posterior draws: it is a mock fit, or was made with chains = 0 or empty = TRUE."
}

#' @export
print.bmm_report_methods <- function(x, ...) {
  label <- style("purple1")
  if (nzchar(x$paragraph)) {
    cat(label("Methods paragraph (draft):"), "\n", sep = "")
    cat(strwrap(x$paragraph, width = 80), sep = "\n")
  }
  if (length(x$formula)) {
    cat("\n", label("Formula:"), "\n", sep = "")
    cat(paste0("  ", strwrap(x$formula, width = 78, exdent = 4)), sep = "\n")
  }
  cat("\n", label("References:"), "\n", sep = "")
  for (reference in numbered_references(x$references)) {
    cat(strwrap(reference, width = 80, exdent = 4), sep = "\n")
  }
  cat("\n", label("Notes (for you, not for the paper):"), "\n", sep = "")
  for (note in x$notes) {
    cat(strwrap(paste("-", note), width = 80, exdent = 2), sep = "\n")
  }
  invisible(x)
}

#' @export
format.bmm_report_methods <- function(x, ...) {
  paste(
    c(
      if (nzchar(x$paragraph)) x$paragraph,
      if (length(x$formula)) paste(x$formula, collapse = "\n"),
      paste(numbered_references(x$references), collapse = "\n")
    ),
    collapse = "\n\n"
  )
}

# citation() text style breaks lines inside a reference, which a pasted
# reference must not keep
numbered_references <- function(references) {
  paste0("[", seq_along(references), "] ", gsub("\\s*\n\\s*", " ", references))
}

#' @export
print.bmm_report_convergence <- function(x, ...) {
  if (is.character(x)) {
    cat(strwrap(x, width = 80), sep = "\n")
    return(invisible(x))
  }
  label <- style("purple1")
  if (nrow(x) == 0) {
    cat(strwrap(convergence_report_text(x), width = 80), sep = "\n")
    return(invisible(x))
  }
  print.data.frame(
    data.frame(
      parameter = x$parameter, rhat = format_rhat(x$rhat),
      ess_bulk = format_count(x$ess_bulk), ess_tail = format_count(x$ess_tail),
      in_summary = x$in_summary
    ),
    right = FALSE, row.names = FALSE
  )
  footer <- c(
    "Algorithm: " = attr(x, "algorithm"),
    "Divergent transitions after warmup: " = attr(x, "divergent"),
    "Transitions at the maximum tree depth: " = attr(x, "max_treedepth_hits")
  )
  footer <- footer[!is.na(footer)]
  cat("\n")
  cat(paste0(label(names(footer)), footer), sep = "\n")
  cat("in_summary marks the parameters that brms's summary() shows.\n")
  invisible(x)
}

#' @export
"[.bmm_report_convergence" <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) class(out) <- "data.frame"
  out
}

# the references in citation order: the model's, then bmm, brms, the backend
# and R, each a character vector; a reference that cannot be read is a
# placeholder, so the numbering in the text does not shift
report_references <- function(fit) {
  model <- as.vector(model_citation(fit))
  list(
    model = if (length(model)) model else "[CITATION NEEDED]",
    bmm = package_reference("bmm"),
    brms = package_reference("brms"),
    backend = package_reference(fit$backend),
    R = package_reference("base")
  )
}

package_reference <- function(package) {
  tryCatch(
    format(utils::citation(package), style = "text"),
    error = function(e) "[CITATION NEEDED]"
  )
}

reference_markers <- function(references) {
  ends <- cumsum(lengths(references))
  starts <- ends - lengths(references) + 1
  stats::setNames(
    Map(function(from, to) paste0("[", paste(from:to, collapse = ", "), "]"), starts, ends),
    names(references)
  )
}

# some constructors store the name with the period of its reference ("... by
# Zhang and Luck (2008)."), which would end the sentence before the marker;
# the report methods of PR 478 apply the same rule
model_display_name <- function(model) {
  sub("\\.$", "", sub("^The ", "", model$name))
}

model_report_text <- function(fit, markers) {
  settings <- fit_settings(fit)
  versions <- settings$versions
  model <- fit$bmm$model
  model_version <- if (all(model$version %in% c(NA, "NA", ""))) "" else glue(", version \"{model$version}\",")
  paste(
    glue(
      "We fitted the {model_display_name(model)} {markers$model}{model_version} ",
      "with the R package bmm{version_phrase(versions[['bmm']])} {markers$bmm}."
    ),
    glue(
      "bmm translated the model into a brms model{version_phrase(versions[['brms']], 'brms version')} ",
      "{markers$brms}, which was estimated in Stan{version_phrase(versions[['stan']])} through the R ",
      "package {settings$backend}{version_phrase(versions[['backend']])} {markers$backend}, ",
      "in R{version_phrase(versions[['R']])} {markers$R}."
    ),
    sampler_report_text(fit, settings),
    "The bmm formula of the model follows this paragraph; a parameter written with = is fixed to that constant."
  )
}

version_phrase <- function(version, label = "version") {
  if (is.na(version)) "" else glue(" ({label} {version})")
}

sampler_report_text <- function(fit, settings) {
  if (!identical(settings$algorithm, "sampling")) {
    return(glue(
      "The model was run with Stan's {settings$algorithm} algorithm, ",
      "which returned {format_count(settings$ndraws_stored)} draws."
    ))
  }
  chains <- settings$chains
  glue(
    "We ran {chains} Markov chain{if (chains > 1) 's' else ''}",
    "{if (uses_nuts(fit)) ' of the No-U-Turn sampler' else ''} with ",
    "{format_count(settings$iter)} iterations{if (chains > 1) ' each' else ''}, of which the ",
    "first {format_count(settings$warmup)} were warmup, ",
    "{if (settings$thin == 1) 'without thinning' else glue('with a thinning interval of {settings$thin}')}, ",
    "which left {format_count(settings$ndraws_stored)} post-warmup draws."
  )
}

# brms samples with Stan's default, NUTS, but the stored arguments are what
# the fit can vouch for: rstan records algorithm "NUTS", cmdstanr engine "nuts"
uses_nuts <- function(fit) {
  args <- fit$fit@stan_args[[1]]
  identical(args$algorithm, "NUTS") || identical(args$engine, "nuts")
}

convergence_report_text <- function(conv) {
  algorithm <- attr(conv, "algorithm")
  if (!identical(algorithm, "sampling")) {
    return(as.character(glue(
      "The {algorithm} algorithm runs no Markov chains, so there are no convergence diagnostics to report."
    )))
  }
  paste(
    glue(
      "Across all estimated parameters, the largest rank-normalized split R-hat was {format_rhat(max(conv$rhat))}, ",
      "the smallest bulk effective sample size (ESS) was {format_count(min(conv$ess_bulk))}, ",
      "and the smallest tail ESS was {format_count(min(conv$ess_tail))}."
    ),
    divergence_report_text(attr(conv, "divergent"))
  )
}

divergence_report_text <- function(divergent) {
  if (divergent == 0) {
    return("There were no divergent transitions after warmup.")
  }
  if (divergent == 1) {
    return("There was 1 divergent transition after warmup.")
  }
  glue("There were {format_count(divergent)} divergent transitions after warmup.")
}

format_rhat <- function(x) {
  formatC(x, format = "f", digits = 3)
}

format_count <- function(x) {
  formatC(round(x), format = "f", digits = 0, big.mark = ",")
}

report_notes <- function(fit, sections, references, markers, conv) {
  uncited <- uncited_part(fit$bmm$model)
  placeholders <- vapply(references, identical, logical(1), "[CITATION NEEDED]")
  c(
    paste(
      "This text is a draft built from the fit. It does not describe the data, the task or what the",
      "model parameters mean; add them, check the text against your analysis and edit it before use."
    ),
    "Numbers in square brackets refer to the references; replace them with citations in the style of your journal.",
    glue("Please cite bmm {markers$bmm} when you use it."),
    "No installed package provides a reference for Stan itself; add one if your journal asks for it.",
    if (placeholders[["model"]]) {
      glue("bmm has no reference for this model: replace [CITATION NEEDED] {markers$model} with its source.")
    },
    if (!is.null(uncited)) {
      glue(
        "Reference {markers$model} is the source bmm gives for this model, but bmm cannot know the source ",
        "of {uncited}: add its reference yourself."
      )
    },
    vapply(setdiff(names(references)[placeholders], "model"), function(group) {
      glue(
        "The reference of {reference_package(fit, group)} could not be read because it is not ",
        "installed: replace [CITATION NEEDED] {markers[[group]]}."
      )
    }, character(1), USE.NAMES = FALSE),
    if ("model" %in% sections && is.null(fit$version$R)) {
      paste(
        "The fit does not store the R version it ran in, as it was made with an older version of bmm:",
        "add the version you fitted the model with after \"in R\"."
      )
    },
    version_mismatch_note(fit),
    if ("convergence" %in% sections) convergence_notes(conv, fit$fit@sim$chains)
  )
}

reference_package <- function(fit, group) {
  switch(group,
    backend = fit$backend,
    R = "R",
    group
  )
}

# the references come from the installed packages, so they can name other
# versions than the ones the fit was made with
version_mismatch_note <- function(fit) {
  packages <- c("bmm", "brms", fit$backend)
  installed <- vapply(packages, function(package) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) NA_character_)
  }, character(1))
  used <- vapply(packages, function(package) version_string(fit$version[[package]]), character(1))
  if (!is.null(fit$version$R)) {
    installed[["R"]] <- as.character(getRversion())
    used[["R"]] <- as.character(fit$version$R)
  }
  differ <- !is.na(installed) & !is.na(used) & installed != used
  if (!any(differ)) {
    return(NULL)
  }
  pkgs <- names(installed)[differ]
  glue(
    "The versions in the text are those the fit stores, but the references are those of the packages ",
    "installed now, which differ: ",
    "{paste0(pkgs, ' ', installed[differ], ' (fit: ', pkgs, ' ', used[differ], ')', collapse = ', ')}."
  )
}

convergence_notes <- function(conv, chains) {
  if (!identical(attr(conv, "algorithm"), "sampling")) {
    return(convergence_report_text(conv))
  }
  above <- conv$parameter[conv$in_summary & conv$rhat > 1.05]
  divergent <- attr(conv, "divergent")
  c(
    paste(
      "R-hat and ESS were computed with posterior::summarise_draws() over all estimated parameters,",
      "including group-level effects and the centred intercepts brms samples; report_convergence()",
      "lists them per parameter."
    ),
    if (chains == 1) {
      paste(
        "With a single chain, R-hat compares the two halves of that chain only;",
        "with several chains it also compares the chains."
      )
    },
    if (length(above)) {
      glue("brms's summary() warns about this fit: R-hat above 1.05 for {collapse_comma(above)}.")
    },
    if (divergent > 0) {
      glue(
        "brms's summary() warns about the {format_count(divergent)} divergent ",
        "transition{if (divergent > 1) 's' else ''} after warmup."
      )
    }
  )
}
