############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_tcc_cd <- function(response = NULL, distance = NULL, similarity = NULL,
                          units = "radians", links = NULL, call = NULL, ...) {
  tcc_model <- .model_tcc()
  sdt_model <- .model_sdt_yn()
  out <- structure(
    list(
      resp_vars = nlist(response),
      other_vars = nlist(distance, similarity, units),
      domain = "Visual working memory (change detection)",
      task = "Change detection",
      name = "Target Confusability Competition (TCC) model for single-probe change detection",
      citation = c(
        tcc_model$citation,
        glue(
          "Williams, J. R., Robinson, M. M., Schurgin, M. W., Wixted, J. T., & Brady, \\
          T. F. (2022). You cannot 'count' how many items people remember in visual \\
          working memory: The importance of signal detection-based measures for \\
          understanding change detection performance. Journal of Experimental \\
          Psychology: Human Perception and Performance, 48(12), 1390-1409. \\
          https://doi.org/10.1037/xhp0001055"
        )
      ),
      version = "NA",
      requirements = glue(
        "- The response is binary, coded 0 = 'same' and 1 = 'change' (or logical, \\
        TRUE = 'change'), one row per trial\n",
        "- `distance` holds the probe minus the target in radians (binned to the \\
        nearest step of the similarity) or, with units = 'bins', in whole steps\n",
        "- `similarity` is the same vector as for tcc(): the similarity f of a \\
        probe to the target at each distance 0, 1, ..., n / 2 steps, with f = 1 at \\
        distance 0"
      ),
      parameters = list(
        d = glue(
          "Memory strength d, as in tcc(): the mean familiarity of a probe that \\
          matches the target, in units of the noise SD; a probe at distance x has \\
          mean d * f(x)"
        ),
        criterion = glue(
          "Decision criterion, measured from the familiarity of a probe with no \\
          similarity to the target, in noise-SD units: the observer answers \\
          'same' when the probe's familiarity exceeds it"
        ),
        bias = glue(
          "Alternative to criterion: the decision criterion measured from the \\
          midpoint between 0 and d (0 = neutral, positive = more 'change' \\
          responses). Use one of the two"
        )
      ),
      links = list(d = tcc_model$links$d, criterion = "identity", bias = "identity"),
      fixed_parameters = list(),
      alternative_parameters = list(c("criterion", "bias")),
      # d is tcc()'s memory strength and the criterion is sdt_yn()'s, so each
      # takes the default priors of the model it comes from
      default_priors = c(tcc_model$default_priors, sdt_model$default_priors[c("criterion", "bias")]),
      init_ranges = c(tcc_model$init_ranges, sdt_model$init_ranges[c("criterion", "bias")])
    ),
    class = c("bmmodel", "tcc_cd"),
    call = call
  )
  set_links(out, links)
}

#' @title `r .model_tcc_cd()$name`
#' @name tcc_cd
#' @description Single-probe change detection under the target confusability
#'   competition model: the probe produces a Gaussian familiarity signal with
#'   unit variance and mean `d * f(x)`, where `f(x)` is the psychophysical
#'   similarity of a probe `x` steps away from the remembered target, and the
#'   observer answers "same" when that signal exceeds a criterion. `d` is the
#'   memory strength of [tcc()]; the criterion is that of [sdt_yn()].
#' @details `r model_docs(.model_tcc_cd())`
#'
#' ## The decision rule
#'
#'   Schurgin et al. (2020) model memory as a familiarity signal that spreads
#'   from the target to other features according to their similarity: a
#'   feature `x` steps away has mean familiarity `d * f(x)` and unit noise.
#'   In continuous report the observer picks the strongest of many such
#'   signals ([tcc()]). With a single probe there is one signal, and the
#'   observer compares it with a criterion, as in yes/no signal detection.
#'   The probability of answering "change" is
#'
#'   \deqn{P(\mathrm{change} \mid x) = \Phi(criterion - d f(x)).}{P(change | x) = Phi(criterion - d * f(x)).}
#'
#'   This is the likelihood of [sdt_yn()] with equal variance, one trial per
#'   row, and the probe's similarity in place of its 0/1 stimulus code: a
#'   probe that matches the target is a signal trial, a probe with `f(x) = 0`
#'   a noise trial. With changes of 180 degrees only, as in Williams et al.
#'   (2022), `tcc_cd()` and [sdt_yn()] fit the same model, with "same" as
#'   "old". Probes at intermediate distances add operating points along `f`.
#'   Schurgin et al. (2020) do not write this rule out; it is their
#'   familiarity model with the decision stage of yes/no recognition. The
#'   observer here knows nothing about the state of memory on a trial. For
#'   the decision rule of Lin and Oberauer (2022), which weighs the
#'   likelihood of a change, see the article on change-detection decision
#'   rules on the package website.
#'
#'   `criterion` is measured from the familiarity of a probe with no
#'   similarity to the target, `bias` from the midpoint between 0 and `d`;
#'   they are related by `bias = criterion - d / 2`. Write one of them in the
#'   formula. A formula that names neither estimates `criterion ~ 1`; one that
#'   names both is refused. The choice matters when `d` differs between rows
#'   that share a criterion term, as described for [sdt_yn()] in its section
#'   "Choosing between `criterion` and `bias`".
#'
#' ## Similarity and distance
#'
#'   `similarity` is the vector [tcc()] takes: one value per distance in
#'   steps, in the order 0, 1, ..., `n / 2`, starting with 1. Its length sets
#'   the grid: 181 values are `n = 360` steps of 1 degree. With `units =
#'   "radians"`, `distance` is the probe minus the target in radians and is
#'   assigned to the nearest step; with `units = "bins"` it is a whole number
#'   of steps. Only the size of the distance matters, not its sign.
#'
#' ## Comparing `d` with `tcc()`
#'
#'   In the model both functions fit, `d` is one parameter, and data from both
#'   tasks generated with the same `d` and `f` return the same value. On real
#'   data, expect `tcc_cd()` to report the larger `d`. [tcc()] fits the
#'   uncorrelated-noise model, which treats the 360 options of a colour wheel
#'   as 360 independent signals. Schurgin et al. (2020, Extended Data Fig. 2)
#'   report that only the correlated-noise d' is interchangeable with a d'
#'   from a same/different task, and that the uncorrelated d' is that value
#'   scaled by about 0.65. A single probe has no such count of channels, so
#'   `tcc_cd()`'s `d` corresponds to their correlated d'. A joint fit of both
#'   tasks that estimates one `d` therefore needs the scale factor as well.
#'
#' ## Set size and other predictors
#'
#'   Set size is not an argument: name its column in the formula, e.g.
#'   `bmf(d ~ 0 + set_size, criterion ~ 1)`. `d` is on the log scale, as in
#'   [tcc()], so effects on it are multiplicative and priors on its
#'   coefficients are priors on `log(d)`.
#'
#' ## Identification
#'
#'   A cell whose probes all have the same similarity gives one rate for two
#'   parameters, so `d` and the criterion trade off. Every cell in which both
#'   are estimated needs probes of at least two similarities, such as probes
#'   at the target and probes far from it. `bmm()` warns when the whole data
#'   set has only one.
#'
#' @param response The name of the variable in the data holding the response,
#'   coded 0 = "same" and 1 = "change". A logical variable is accepted, with
#'   `TRUE` = "change".
#' @param distance The name of the variable in the data holding the probe's
#'   distance from the target: the probe minus the target in radians, or with
#'   `units = "bins"` in whole steps of the similarity. 0 is a probe that
#'   matches the target.
#' @param similarity A numeric vector of `n / 2 + 1` values in \[0, 1\], for
#'   an even `n` of at least 4: the similarity of a probe to the target at
#'   each distance 0, 1, ..., `n / 2` steps, starting with 1 at distance 0.
#'   The same vector as for [tcc()].
#' @param units Either `"radians"` or `"bins"`, see `distance`.
#' @param links A named list of link functions for `d`, `criterion` or `bias`,
#'   e.g. `links = list(d = "identity")`. `d` defaults to the log link,
#'   `criterion` and `bias` to the identity link.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Lin, H.-Y., & Oberauer, K. (2022). An interference model for visual
#'   working memory: Applications to the change detection task. \emph{Cognitive
#'   Psychology}, \emph{133}, 101463. \doi{10.1016/j.cogpsych.2022.101463}
#' @seealso [tcc()] for continuous report, [sdt_yn()] for yes/no detection with
#'   two stimulus classes, [dtcc_cd()] for the distribution functions.
#' @keywords bmmodel
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' f <- exp(-(0:180) / 20)
#' f <- (f - min(f)) / (1 - min(f))
#' # half the probes match the target, the rest lie 10 to 180 degrees away
#' dat <- data.frame(dist = c(rep(0, 500), runif(500, 10, 180)) * pi / 180)
#' dat$change <- rtcc_cd(1000, dat$dist, d = 2, similarity = f, criterion = 1)
#' fit <- bmm(bmf(d ~ 1, criterion ~ 1), dat, tcc_cd("change", "dist", similarity = f),
#'            cores = 4)
tcc_cd <- function(response, distance, similarity, units = c("radians", "bins"),
                   links = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  units <- match.arg(units)
  .tcc_cd_check_similarity(similarity)
  .model_tcc_cd(response = response, distance = distance, similarity = similarity,
                units = units, links = links, call = call, ...)
}

# the vector tcc() takes, read without n_options: its length sets the grid
.tcc_cd_check_similarity <- function(similarity) {
  stopif(length(similarity) < 3,
         "`similarity` must hold at least 3 values, one for each distance from 0 to n / 2 steps.")
  .tcc_check_similarity(similarity, .tcc_cd_n_steps(similarity))
}

.tcc_cd_n_steps <- function(similarity) {
  2 * (length(similarity) - 1)
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.tcc_cd <- function(model, data, formula) {
  resp <- model$resp_vars$response
  dist <- model$other_vars$distance
  units <- model$other_vars$units
  similarity <- .tcc_similarity(model)
  stopif(not_in(resp, colnames(data)),
         "The response variable '{resp}' is not present in the data.")
  if (is.logical(data[[resp]])) {
    data[[resp]] <- as.integer(data[[resp]])
  }
  stopif(!is.numeric(data[[resp]]) || any(!data[[resp]] %in% c(0, 1, NA)),
         "The response variable '{resp}' must be coded 0 = 'same' and 1 = 'change'.")
  stopif(not_in(dist, colnames(data)),
         "The distance variable '{dist}' is not present in the data.")
  stopif(!is.numeric(data[[dist]]), "The distance variable '{dist}' must be numeric.")
  stopif(units == "radians" && isTRUE(max(abs(data[[dist]]), na.rm = TRUE) > 2 * pi),
         "The distance variable '{dist}' holds values beyond 2*pi, so it is not in \\
         radians. Convert it with deg2rad(), or set units = 'bins' if it counts steps.")
  stopif(units == "bins" && any(data[[dist]] %% 1 != 0, na.rm = TRUE),
         "With units = 'bins', the distance variable '{dist}' must hold whole numbers.")
  warnif("tcc_cd_dist" %in% colnames(data),
         "Column 'tcc_cd_dist' in your data is reserved by the TCC change-detection \\
         model and will be overwritten.")
  data$tcc_cd_dist <- abs(.tcc_bin(data[[dist]], .tcc_cd_n_steps(similarity), units))
  warnif(length(unique(stats::na.omit(similarity[data$tcc_cd_dist + 1]))) == 1,
         "Every probe in the data has the same similarity to its target, so there is \\
         only one similarity, which gives one rate for 'd' and the criterion together. \\
         Neither is identified; add probes at other distances (e.g. at the target \\
         and far from it).")
  data[[resp]] <- as.integer(data[[resp]])
  NextMethod("check_data")
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.tcc_cd <- function(model, formula) {
  brms::bf(paste0(model$resp_vars$response, " | vint(tcc_cd_dist) ~ 0"))
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.tcc_cd <- function(model, data, formula) {
  anchor <- .sdt_anchor_parameter(formula)
  formula <- bmf2bf(model, formula)
  # the dpar names must match the parameter the formula estimates, so each
  # parameterization is its own family, as for sdt_yn()
  formula$family <- brms::custom_family(
    if (anchor == "bias") "tcc_cd_bias" else "tcc_cd",
    dpars = c("mu", "d", anchor),
    links = c("identity", model$links$d, model$links[[anchor]]),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_tcc_cd,
    posterior_predict = posterior_predict_tcc_cd,
    posterior_epred = posterior_epred_tcc_cd,
    vars = c("vint1[n]", "tcc_cd_similarity")
  )
  nlist(formula, data, stanvars = .tcc_cd_stanvars(model))
}

# The yes/no SDT distribution functions and a kernel of their own beside them,
# so that the sdt_yn programs stay as they are (#566)
.tcc_cd_stanvars <- function(model) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  stan_funs <- paste(read_lines2(file.path(sc_path, "sdt_dist_funs.stan")),
                     read_lines2(file.path(sc_path, "tcc_cd_funs.stan")), sep = "\n")
  brms::stanvar(scode = stan_funs, block = "functions") +
    brms::stanvar(x = .tcc_similarity(model), name = "tcc_cd_similarity")
}

############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT AND POSTERIOR_EPRED                         ####
############################################################################# !

# log P of the response (1 = "change") at a probe of the given similarity, the
# R twin of tcc_cd_lpmf(); "same" is sdt_yn's "old". A single response is
# scored against every draw, so the arguments recycle to the longest
.tcc_cd_log_p <- function(change, strength, d, criterion) {
  eta <- d * strength - criterion
  n <- if (length(change) && length(eta)) max(length(change), length(eta)) else 0L
  eta <- rep_len(eta, n)
  out <- .sdt_log_p_old(eta, "normal")
  is_change <- which(rep_len(change, n) == 1)
  out[is_change] <- .sdt_log_p_new(eta[is_change], "normal")
  out
}

# the similarity of a probe at each distance in radians, binned as check_data()
# bins it with units = "radians"
.tcc_cd_strength <- function(distance, similarity) {
  similarity[abs(.tcc_bin(distance, .tcc_cd_n_steps(similarity), "radians")) + 1]
}

.tcc_cd_prep_strength <- function(prep, i = NULL) {
  dist <- if (is.null(i)) prep$data$vint1 else prep$data$vint1[i]
  as.numeric(prep$data$tcc_cd_similarity)[dist + 1]
}

log_lik_tcc_cd <- function(i, prep) {
  d <- brms::get_dpar(prep, "d", i = i)
  criterion <- .sdt_criterion_dpar(prep, d, 1, i = i)
  .tcc_cd_log_p(prep$data$Y[i], .tcc_cd_prep_strength(prep, i), d, criterion)
}

posterior_predict_tcc_cd <- function(i, prep, ...) {
  d <- brms::get_dpar(prep, "d", i = i)
  criterion <- .sdt_criterion_dpar(prep, d, 1, i = i)
  stats::rbinom(length(d), 1, exp(.tcc_cd_log_p(1, .tcc_cd_prep_strength(prep, i), d, criterion)))
}

# the expected response is the probability of answering "change"
posterior_epred_tcc_cd <- function(prep) {
  d <- prep$dpars$d
  criterion <- .sdt_criterion_dpar(prep, d, 1)
  strength <- .epred_data(.tcc_cd_prep_strength(prep), prep)
  .epred_matrix(exp(.tcc_cd_log_p(1, strength, d, criterion)), prep)
}

############################################################################# !
# STORED FRAME                                                           ####
############################################################################# !

# brms keeps the binned distance but not the user's column, which is rebuilt
# at the centre of its step
#' @exportS3Method
revert_check_data.tcc_cd <- function(model, data) {
  dist <- model$other_vars$distance
  if (!dist %in% colnames(data)) {
    data[[dist]] <- if (model$other_vars$units == "radians") {
      data$tcc_cd_dist * 2 * pi / .tcc_cd_n_steps(.tcc_similarity(model))
    } else {
      data$tcc_cd_dist
    }
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), dist)
  }
  data$tcc_cd_dist <- NULL
  NextMethod("revert_check_data")
}
