#' Simulate data from a bmm model
#'
#' @description
#' `simulate_data()` draws data from a model given the values of its parameters
#' and a design, in the layout [bmm()] expects for that model: the design's
#' columns with the response columns added under the names the model was
#' specified with. The result can be passed to [bmm()] as it is, which is what
#' a parameter recovery study needs.
#'
#' @param model A `bmmodel` object, such as `mixture2p(resp_error = "y")`.
#' @param pars A named list of parameter values, named as [parameter_info()]
#'   lists the model's parameters, or a data frame of them with one row or one
#'   row per row of `design`, such as the attribute `simulate_pars` of an
#'   earlier result. Each value is a single number or a vector with one element
#'   per row of `design`, element `i` going to row `i`. Parameters the model
#'   fixes may be left out, and then take their fixed value; a value given for
#'   one of them replaces it. Of alternative parameters, such as `criterion`
#'   and `bias` of [sdt_yn()], give one.
#' @param design A data frame with one row per observation of the model's data:
#'   one trial for trial-level models, one cell for the models of response
#'   counts. It holds the columns the model reads (see the section Design) and
#'   any others, such as condition labels or a subject identifier, which are
#'   returned unchanged. `NULL` for the models that read no column, and then
#'   `n_trials` in `...` gives the number of observations.
#' @param scale The scale of the values in `pars`. `"link"` (the default) is
#'   the scale [bmm()] estimates the parameters on and `summary()` reports, for
#'   example `kappa` on the log scale. `"native"` is the scale after the
#'   inverse of the model's link, the scale of [native_parameters()]: `kappa`
#'   as a precision, `thetat` of [mixture2p()] as the probability of recall.
#'   The weights of a softmax group, such as `thetat` and `thetant` of
#'   [mixture3p()], are given as the probabilities themselves. They are the
#'   weights of a display with non-targets; where the set size is 1, the
#'   non-target weight goes to the target and to guessing in proportion, as in
#'   the fit, so `thetat = 0.6` with `thetant = 0.2` gives a probability of
#'   recall of 0.75 there.
#' @param ... `n_trials`, the number of observations when `design` is `NULL`;
#'   `formula`, the [bmmformula()] with the activation functions, for the
#'   custom version of [m3()], whose parameters exist only in that formula.
#'
#' @details
#' The values are converted from the link scale with the inverse of each
#' parameter's link, by the code [native_parameters()] uses, so that a softmax
#' group and the set-size gate of the non-target weights of [mixture3p()] (no
#' non-target weight where no non-target was shown) are applied as in the fit.
#' Some parameters are defined on a transformed scale behind an identity link,
#' such as `sdratio` of [sdt_rating()], the log of the ratio of the standard
#' deviations. For those, `"native"` is the value itself; the transformation
#' to the argument of the random number generator (here `exp(sdratio)`) is
#' applied on both scales.
#'
#' The data are drawn with the model's `r*()` function, such as [rmixture2p()]
#' or [rsdt_rating()]. They are not independent of the model's density
#' functions: a recovery study that has to rule out an error shared by the two
#' needs a generator written separately. The change-detection models (`*_cd`)
#' are the exception: their data come from an observer that retrieves a feature
#' on every trial and applies the decision rule to it, which shares no code with
#' the response probabilities the likelihood integrates.
#'
#' @section Design:
#' The columns each model reads from `design`, under the names the model was
#' specified with:
#'
#' | Model | Columns |
#' |:--|:--|
#' | [sdm()], [ddm()], [cswald()] | none |
#' | [sdm_cd()] | `probe`, relative to the target |
#' | [mixture2p()] | `set_size` (if it names a column) for the versions `"slot"` and `"slot_averaging"`; none for `"simple"` |
#' | [mixture3p()] | `set_size` (if it names a column), `nt_features`, of which a row's first `set_size - 1` are read |
#' | [imm()] | as [mixture3p()], plus `nt_distances` for the versions `"full"` and `"bsc"` |
#' | [m3()] | `n_trials`, the number of trials of each row; `num_options` when it names columns; for the custom version, the columns its activation functions read besides the parameters, as numbers (code a factor as 0/1 indicators) |
#' | [mpt()] | `n_trials`; `tree_id` and `covariates` when the model has them |
#' | [ezdm()] | `n_trials` |
#' | [sdt_yn()] | `stimulus`, `n_trials` |
#' | [sdt_mafc()] | `n_trials`; `m` when it names a column |
#' | [sdt_ranking()] | `n_trials`; `m` when it names a column |
#' | [sdt_rating()], [sdt_cdp()] | `stimulus`, `n_trials` |
#'
#' The simplex groups of [mpt()] are given as probabilities that sum to 1 on
#' both scales.
#'
#' @return A data frame: `design` with the response columns added. It carries
#'   the attribute `simulate_pars`, a data frame with the value of every
#'   parameter (fixed ones included) for every row on the link scale, the
#'   values a recovery study compares the estimates with, and `simulate_call`,
#'   the model's class and the scale. For [ezdm()], rows whose simulated mean
#'   response time is not positive are dropped with a warning, because [bmm()]
#'   rejects them; their row numbers are in the attribute `dropped_rows`. A
#'   mean response time below the non-decision time is a valid draw and is
#'   kept.
#'
#' @family model information
#' @keywords distribution
#' @export
#'
#' @examples
#' # circular responses: two-parameter mixture model, 100 trials
#' data <- bmm::simulate_data(mixture2p(resp_error = "y"),
#'   pars = list(kappa = log(8), thetat = 0.5), n_trials = 100
#' )
#' head(data)
#' attr(data, "simulate_pars")[1, ]
#'
#' # the same on the native scale: kappa as a precision, thetat as P(memory)
#' data <- bmm::simulate_data(mixture2p(resp_error = "y"),
#'   pars = list(kappa = 8, thetat = 0.62), scale = "native", n_trials = 100
#' )
#'
#' # response times and decisions: one value of the drift rate per trial
#' design <- data.frame(cond = rep(c("easy", "hard"), each = 50))
#' data <- bmm::simulate_data(ddm(rt = "rt", response = "response"),
#'   pars = list(drift = ifelse(design$cond == "easy", 2, 0.5), bound = 0, ndt = -1.5),
#'   design = design
#' )
#' head(data)
#'
#' # counts: yes/no recognition, 50 trials per stimulus type
#' design <- data.frame(stimulus = c(0, 1), n_trials = 50)
#' bmm::simulate_data(sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"),
#'   pars = list(d = 1.5, criterion = 0.5), design = design
#' )
#'
#' # categorical responses: memory measurement model, 100 trials per row
#' model <- m3(
#'   resp_cats = c("corr", "other", "npl"), num_options = c(1, 4, 5),
#'   choice_rule = "softmax", version = "ss"
#' )
#' bmm::simulate_data(model, pars = list(c = 2, a = 0.5), design = data.frame(n_trials = c(100, 100)))
simulate_data <- function(model, pars, design = NULL, scale = c("link", "native"), ...) {
  dots <- list(...)
  stopif(length(dots) > 0 && (is.null(names(dots)) || !all(nzchar(names(dots)))),
         "Arguments in `...` must be named, such as `n_trials = 100`.")
  stopif(!is.null(design) && !is.null(dots[["n_trials"]]),
         "With a `design`, the number of observations is its number of rows: drop \\
         `n_trials =`. The models of response counts read the trials of each row \\
         from a column of the design (see ?simulate_data, section Design).")
  UseMethod("simulate_data")
}

#' @rdname simulate_data
#' @export
simulate_data.default <- function(model, pars, design = NULL, scale = c("link", "native"), ...) {
  name <- if (is_bmmodel(model)) model$name else class(model)[1]
  stop2("Model '{name}' has no simulate_data() method. See ?simulate_data.")
}

# the design a method simulates on: the user's data frame, checked for the
# columns the model reads, or n_trials rows without columns. `ragged` columns
# may hold NA, where a display has fewer items than the design has columns.
# `trials` and `counts` are read columns that the generators take as whole
# numbers, the numbers of trials of a row and other counts (options, set sizes)
.design_or_n <- function(model, design, columns = character(0), ..., ragged = character(0),
                         trials = character(0), counts = character(0)) {
  columns <- c(columns, trials, counts)
  if (is.null(design)) {
    stopif(length(c(columns, ragged)) > 0,
           "{model$name} reads the column(s) {collapse_comma(c(columns, ragged))}: give them \\
           in a `design` data frame. See ?simulate_data, section Design.")
    n_trials <- list(...)[["n_trials"]]
    stopif(length(n_trials) != 1 || !is.numeric(n_trials) || n_trials < 1 || n_trials %% 1 != 0,
           "Without a `design`, give the number of observations as `n_trials`, \\
           a single positive whole number.")
    return(data.frame(row.names = seq_len(n_trials)))
  }
  stopif(!is.data.frame(design) || nrow(design) == 0,
         "`design` must be a data frame with one row per observation.")
  missing <- setdiff(c(columns, ragged), names(design))
  stopif(length(missing) > 0,
         "The design lacks the column(s) {collapse_comma(missing)} that {model$name} \\
         reads. See ?simulate_data, section Design.")
  na <- columns[vapply(design[columns], anyNA, logical(1))]
  stopif(length(na) > 0, "The design's column(s) {collapse_comma(na)} contain NA.")
  .check_whole_columns(design, trials, 1, "numbers of trials")
  .check_whole_columns(design, counts, 0, "counts")
  design
}

.check_whole_columns <- function(design, columns, min, what) {
  bad <- columns[!vapply(design[columns], function(x) is.numeric(x) && all(x >= min & x %% 1 == 0),
                         logical(1))]
  stopif(length(bad) > 0,
         "The design's column(s) {collapse_comma(bad)} hold {what}: they must be \\
         whole numbers of {min} or more.")
}

# a row's first set_size - 1 non-target columns describe the items shown; the
# columns past them may be NA
.check_nt_na <- function(design, columns, set_size) {
  nt <- as.matrix(design[columns])
  shown_na <- is.na(nt) & col(nt) < set_size
  bad <- columns[colSums(shown_na) > 0]
  rows <- which(rowSums(shown_na) > 0)
  listed <- paste(utils::head(rows, 5), collapse = ", ")
  if (length(rows) > 5) listed <- glue::glue("{listed}, ... ({length(rows)} in all)")
  stopif(length(bad) > 0,
         "The design's column(s) {collapse_comma(bad)} contain NA within the set size, \\
         in row(s) {listed}: a row's first set_size - 1 non-target columns need values.")
}

#' The values of a model's parameters for every row of a design
#'
#' Fills the fixed parameters, checks names and lengths, and converts the
#' values to the native scale through `native_transform()`, the conversion
#' `native_parameters()` uses, with every parameter as a 1 x n matrix.
#'
#' @return A list: `link`, a data frame of the link-scale values, one column
#'   per parameter; `native`, a named list of native-scale vectors; `call`,
#'   the provenance attribute.
#' @keywords internal
#' @noRd
.simulate_pars <- function(model, pars, design, scale, formula = NULL) {
  scale <- match.arg(scale, c("link", "native"))
  n <- nrow(design)
  stopif(!is.list(pars) || is.null(names(pars)) || !all(nzchar(names(pars))) ||
           is.data.frame(pars) && !nrow(pars) %in% c(1, n),
         "`pars` must be a named list of parameter values, or a data frame of them \\
         with one row or one row per row of the design ({n}).")
  pars <- as.list(pars)
  model_pars <- .simulate_par_names(model, pars, formula)

  unknown <- setdiff(names(pars), model_pars)
  stopif(length(unknown) > 0,
         "Unknown parameter(s) {collapse_comma(unknown)} for {model$name}. Its \\
         parameters are {collapse_comma(model_pars)}.")
  fixed <- setdiff(intersect(model_pars, names(model$fixed_parameters)), names(pars))
  missing <- setdiff(model_pars, c(names(pars), fixed))
  stopif(length(missing) > 0, "No value for the parameter(s) {collapse_comma(missing)}.")
  na <- names(pars)[vapply(pars, anyNA, logical(1))]
  stopif(length(na) > 0, "The value(s) of {collapse_comma(na)} contain NA.")
  not_numeric <- names(pars)[!vapply(pars, is.numeric, logical(1))]
  stopif(length(not_numeric) > 0, "The value(s) of {collapse_comma(not_numeric)} must be numeric.")
  wrong <- names(pars)[!lengths(pars) %in% c(1L, n)]
  stopif(length(wrong) > 0,
         "The value(s) of {collapse_comma(wrong)} must be a single value or one per \\
         row of the design ({n}); found length {paste(lengths(pars[wrong]), collapse = ', ')}.")

  given <- lapply(pars, function(x) rep_len(as.numeric(x), n))
  if (scale == "native") given <- .simulate_native_to_link(model, given)
  link <- c(given, lapply(model$fixed_parameters[fixed], function(x) rep_len(as.numeric(x), n)))[model_pars]

  # a fixed value is on the link the model reports for the parameter as fixed
  fixed_links <- model$links_fixed[intersect(fixed, names(model$links_fixed))]
  model$links[names(fixed_links)] <- fixed_links
  linpred <- lapply(link, matrix, nrow = 1)
  native <- lapply(native_transform(model, linpred, design)[model_pars], as.vector)
  list(link = as.data.frame(link), native = native,
       call = list(model = class(model), scale = scale))
}

# the parameters a call simulates: the model's, without the alternatives the
# call does not use and without the stick-breaking components of an mpt
# simplex, plus those a custom version's formula introduces. Of a group of
# alternatives the one in `pars` is used, as the one in a formula is
.simulate_par_names <- function(model, pars, formula) {
  model_pars <- if (is.null(formula)) names(model$parameters) else parameter_info(model, formula)$parameter
  for (group in model$alternative_parameters) {
    stopif(sum(group %in% names(pars)) > 1, "Give one of {collapse_comma(group)}, not both.")
  }
  setdiff(model_pars, c(unused_alternative_parameters(model, pars), internal_parameters(model)))
}

.simulate_native_to_link <- function(model, given) {
  softmax_pars <- .np_softmax_pars(model, names(given))
  for (par in setdiff(names(given), softmax_pars)) {
    given[[par]] <- link_transform(given[[par]], model$links[[par]])
  }
  if (length(softmax_pars) == 0) {
    return(given)
  }
  # softmaxinv() with the reference category last at 0, the link of brms'
  # unpredicted mixture component, for every row at once
  reference <- 1 - Reduce(`+`, given[softmax_pars])
  stopif(any(reference <= 0),
         "The weights {collapse_comma(softmax_pars)} must sum to less than 1, to \\
         leave weight for the reference category.")
  given[softmax_pars] <- lapply(given[softmax_pars], function(p) log(p) - log(reference))
  given
}

# the design with the response columns added and the attributes of the call
.simulate_result <- function(design, responses, sim) {
  design[names(responses)] <- responses
  attr(design, "simulate_pars") <- sim$link
  attr(design, "simulate_call") <- sim$call
  design
}

# Row indices grouped by identical values of the columns of `key`, for the
# generators that take one value per argument and call (rmixture3p(), rimm(),
# rm3(), rmpt()); 17 significant digits keep distinct doubles apart
.row_groups <- function(key) {
  key <- lapply(key, function(x) {
    if (is.numeric(x)) formatC(x, digits = 17, format = "g") else as.character(x)
  })
  split(seq_along(key[[1]]), do.call(paste, c(key, sep = "\r")))
}

# one draw per row from a generator that takes one display per call: rows
# with the same values in `key` share a call, draw(i, n) for the group of row i
.simulate_by_group <- function(key, draw) {
  y <- numeric(nrow(key))
  for (idx in .row_groups(key)) y[idx] <- draw(idx[1], length(idx))
  y
}

# the set-size column a non_targets model reads, if it names one
.set_size_column <- function(model) {
  set_size <- model$other_vars$set_size
  if (is.character(set_size)) set_size
}

# the decisions models' generators return rt and response
.rt_response <- function(model, out) {
  stats::setNames(out[c("rt", "response")], c(model$resp_vars$rt, model$resp_vars$response))
}
