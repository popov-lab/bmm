# Pure computations ------------------------------------------------------------
# Everything below works on matrices of posterior draws (draws x persons) and
# never touches brms, so the reliability formulas can be tested exactly against
# conjugate normal posteriors. Whether a model was fitted to trial-level or to
# aggregated data does not matter here: only the person-level posterior enters.

# Posterior (EAP / marginal) reliability of a person-level score.
# The error variance is the mean posterior variance, not the squared mean
# posterior SD: under unequal precision across persons the latter understates
# the error and overstates the reliability.
#
# The interval reflects the uncertainty about the true-score variance V. Each
# person's sampling variance (noise) is backed out of their posterior variance
# once, under the normal approximation e2 = V * noise / (V + noise), and the
# shrinkage is recomputed in every draw, so the reliability per draw,
# mean(V_s / (V_s + noise)), stays within (0, 1) and equals the point estimate
# at the average V. Holding e2 fixed while V varies would not.
.rel_core <- function(theta, prob = 0.95) {
  e2 <- matrixStats::colVars(theta)
  V <- matrixStats::rowVars(theta)
  Vbar <- mean(V)
  noise <- .rel_noise(e2, Vbar)
  draws <- .rel_draws(V, noise)
  probs <- c((1 - prob) / 2, (1 + prob) / 2)
  bounds <- stats::quantile(draws, probs, names = FALSE)
  list(
    reliability = mean(Vbar / (Vbar + noise)),
    lower = bounds[1],
    upper = bounds[2],
    draws = draws,
    person = Vbar / (Vbar + noise),
    m = colMeans(theta),
    e2 = e2,
    noise = noise,
    V = V,
    Vbar = Vbar
  )
}

# Sampling variance of each person's score: the error variance their posterior
# variance implies once the shrinkage towards the group is undone. A posterior
# variance at or above V means the data say nothing about the person.
.rel_noise <- function(e2, Vbar) {
  e2 <- pmin(e2, Vbar * (1 - 1e-8))
  e2 * Vbar / (Vbar - e2)
}

.rel_draws <- function(V, noise) {
  vapply(V, function(v) mean(v / (v + noise)), numeric(1))
}

# Reliability after multiplying every person's trials by f, which divides
# their sampling variance by f. Persons keep their unequal trial counts; a
# pooled Spearman-Brown correction overstates the gain when precision varies.
.rel_project_personwise <- function(core, f, prob = 0.95) {
  noise <- core$noise / f
  probs <- c((1 - prob) / 2, (1 + prob) / 2)
  bounds <- stats::quantile(.rel_draws(core$V, noise), probs, names = FALSE)
  list(reliability = mean(core$Vbar / (core$Vbar + noise)), lower = bounds[1], upper = bounds[2])
}

# True-score stability (Pearson) and absolute agreement (ICC(A,1)) between two
# occasions, computed per draw across persons. Persons are the columns, in the
# same order in both matrices.
.rel_retest <- function(theta1, theta2) {
  N <- ncol(theta1)
  c1 <- theta1 - rowMeans(theta1)
  c2 <- theta2 - rowMeans(theta2)
  covariance <- rowSums(c1 * c2) / (N - 1)
  v1 <- rowSums(c1^2) / (N - 1)
  v2 <- rowSums(c2^2) / (N - 1)
  shift <- rowMeans(theta1) - rowMeans(theta2)
  list(
    stability = covariance / sqrt(v1 * v2),
    agreement = 2 * covariance / (v1 + v2 + shift^2)
  )
}

# G-study from the universe score (predictions with only the person-level
# terms). Its posterior reliability is the G coefficient of the design as
# fitted; subtracting the person x facet variances from the implied relative
# error leaves the error of the finest person-level score, which a D-study then
# divides by the new numbers of facet levels and trials.
#
# `interactions` and `mains` are lists of list(facets = <character>, var =
# <variance draws>) for the person x facet and the facet-only random effects;
# `n` holds the number of levels of each facet per person.
.rel_gstudy <- function(universe, interactions, mains, n, prob = 0.95) {
  core <- .rel_core(universe, prob)
  delta <- core$Vbar * (1 - core$reliability) / core$reliability
  residual_facets <- unique(unlist(lapply(interactions, `[[`, "facets")))
  explained <- .rel_error_sum(interactions, n, stats::median)
  g <- list(
    core = core,
    interactions = interactions,
    mains = mains,
    n = n,
    residual_facets = residual_facets,
    s2 = max(prod(n[residual_facets]) * (delta - explained), 0),
    prob = prob
  )
  c(g, .rel_dstudy(g, n, 1))
}

.rel_dstudy <- function(g, n, f) {
  n_all <- g$n
  n_all[names(n)] <- n
  residual <- g$s2 / (f * prod(n_all[g$residual_facets]))
  rel_error <- .rel_error_sum(g$interactions, n_all, identity) + residual
  abs_error <- rel_error + .rel_error_sum(g$mains, n_all, identity)
  rel_point <- .rel_error_sum(g$interactions, n_all, stats::median) + residual
  abs_point <- rel_point + .rel_error_sum(g$mains, n_all, stats::median)
  V <- g$core$V
  Vbar <- g$core$Vbar
  probs <- c((1 - g$prob) / 2, (1 + g$prob) / 2)
  rel_bounds <- stats::quantile(V / (V + rel_error), probs, names = FALSE)
  abs_bounds <- stats::quantile(V / (V + abs_error), probs, names = FALSE)
  list(
    g_relative = Vbar / (Vbar + rel_point),
    g_relative_lower = rel_bounds[1],
    g_relative_upper = rel_bounds[2],
    g_absolute = Vbar / (Vbar + abs_point),
    g_absolute_lower = abs_bounds[1],
    g_absolute_upper = abs_bounds[2]
  )
}

# Sum of variance components, each divided by the number of levels of its
# facets. `summarise` is identity for draws and median for the point estimate.
.rel_error_sum <- function(terms, n, summarise) {
  Reduce(`+`, lapply(terms, function(term) summarise(term$var) / prod(n[term$facets])), 0)
}

# Smallest trial multiplier at which project(f) reaches target; Inf if no
# number of trials can, because the remaining error does not shrink with trials.
.rel_solve_factor <- function(project, target, range = c(1e-4, 1e6)) {
  if (project(range[2]) < target) {
    return(Inf)
  }
  if (project(range[1]) >= target) {
    return(range[1])
  }
  root <- stats::uniroot(
    function(log_f) project(exp(log_f)) - target,
    log(range), tol = 1e-10
  )
  exp(root$root)
}


# Random-effect structure ------------------------------------------------------

# brms stores "" rather than an empty frame for a parameter without
# group-level terms
.rel_re_frame <- function(bterms, par) {
  re <- (bterms$dpars[[par]] %||% bterms$nlpars[[par]])$re
  if (is.data.frame(re)) re else NULL
}

# The person variable is the grouping variable named in `group`, or else the
# only plain grouping variable, or else the one plain grouping variable that
# every interaction grouping (id:session, id:stimulus) contains.
.rel_person_var <- function(bterms, pars, group) {
  labels <- unique(unlist(lapply(pars, function(par) .rel_re_frame(bterms, par)$group)))
  stopif(
    length(labels) == 0,
    "None of the parameters {collapse_comma(pars)} has group-level effects, so \\
    there is no between-person variation to separate from measurement error."
  )
  vars <- unique(unlist(strsplit(labels, ":", fixed = TRUE)))
  if (!is.null(group)) {
    stopif(
      length(group) != 1 || !group %in% vars,
      "'group' must be one of the grouping variables {collapse_comma(vars)}."
    )
    return(group)
  }
  plain <- labels[!grepl(":", labels, fixed = TRUE)]
  candidates <- plain
  if (length(plain) > 1) {
    interactions <- strsplit(labels[grepl(":", labels, fixed = TRUE)], ":", fixed = TRUE)
    if (length(interactions) > 0) {
      candidates <- intersect(plain, Reduce(intersect, interactions))
    }
  }
  stopif(
    length(candidates) != 1,
    "Cannot tell which grouping variable identifies the persons among \\
    {collapse_comma(plain)}. Name it with 'group'."
  )
  candidates
}

# Sort the random-effect terms of one parameter by how they involve the person
# variable: the person term itself, person x facet interactions (id:session),
# and facet terms without the person (session, stimulus). A G-study needs the
# facet terms to be plain random intercepts.
.rel_terms <- function(bterms, par, group) {
  re <- .rel_re_frame(bterms, par)
  if (is.null(re) || NROW(re) == 0) {
    return(NULL)
  }
  parts <- strsplit(re$group, ":", fixed = TRUE)
  with_person <- vapply(parts, function(p) group %in% p, logical(1))
  person <- vapply(parts, function(p) identical(p, group), logical(1))
  intercept_only <- vapply(re$form, function(f) {
    tt <- stats::terms(f)
    length(attr(tt, "term.labels")) == 0 && attr(tt, "intercept") == 1
  }, logical(1))
  facet_term <- function(i) list(label = re$group[i], facets = setdiff(parts[[i]], group))
  facets <- !person
  list(
    person = re[person, , drop = FALSE],
    interactions = lapply(which(with_person & !person), facet_term),
    mains = lapply(which(!with_person), facet_term),
    facet_vars = unique(unlist(lapply(parts[facets], setdiff, group))),
    simple_facets = all(intercept_only[facets] & re$gtype[facets] == "")
  )
}

# re_formula with only the person terms, each with all its coefficients:
# brms subsets coefficients too, so ~ (1 | id) would drop the slope of
# (1 + cond | id) from the universe score
.rel_person_re_formula <- function(person_terms) {
  terms <- vapply(seq_len(NROW(person_terms)), function(i) {
    glue::glue("({deparse1(person_terms$form[[i]][[2]])} | {person_terms$group[i]})")
  }, character(1))
  stats::as.formula(paste("~", paste(terms, collapse = " + ")))
}

.rel_sd_draws <- function(x, label, par, draw_ids) {
  coef <- if (par == "mu") "Intercept" else paste0(par, "_Intercept")
  name <- paste0("sd_", label, "__", coef)
  stopif(
    !name %in% brms::variables(x),
    "Cannot find the standard deviation '{name}' in the fitted model."
  )
  # draw ids index the draws with chains merged, as in posterior_linpred()
  as.numeric(brms::as_draws_matrix(x, variable = name))[draw_ids]^2
}

# Harmonic mean over persons of the number of levels of one facet a person was
# observed at, within the rows of one cell. Terms crossing several facets are
# divided by the product of these, which assumes the facets are fully crossed.
.rel_facet_levels <- function(data, rows, group, facets) {
  sub <- data[rows, c(group, facets), drop = FALSE]
  combos <- unique(sub)
  counts <- table(combos[[group]])
  counts <- counts[counts > 0]
  length(counts) / sum(1 / counts)
}


# Trial counts -----------------------------------------------------------------

# Number of trials each row of the model data contributes. Rows of trial-level
# models count once; aggregated models carry the count in the brms trials()
# addition term, so most models need no method. A model that aggregates
# without trials() needs a method returning one count per row.
person_trials <- function(model, data, formula) {
  UseMethod("person_trials")
}

#' @exportS3Method
person_trials.default <- function(model, data, formula) {
  trials <- brms::brmsterms(formula)$adforms$trials
  if (is.null(trials)) {
    return(rep(1, nrow(data)))
  }
  as.numeric(eval(trials[[2]][[2]], data))
}

#' @exportS3Method
person_trials.ezdm_4par <- function(model, data, formula) {
  as.numeric(data[[model$other_vars$n_trials]])
}


# Exported functions -----------------------------------------------------------

#' Reliability and generalizability of person-level parameters
#'
#' @description
#' Estimates how reliably a fitted model separates persons on each parameter,
#' from the posterior of the person-level parameters. Because only the
#' person-level posterior enters, the estimates are the same whether the model
#' was fitted to trial-level or to aggregated data, and they work for every
#' model in the package, including models added later.
#'
#' When a parameter has random effects for facets crossed with persons, such as
#' sessions or stimuli (`(1 | id) + (1 | id:session) + (1 | stimulus)`), the
#' function also reports the coefficients of generalizability theory.
#'
#' @param x A `bmmfit` object.
#' @param pars Character vector of parameters. Defaults to all parameters that
#'   are not fixed to a constant.
#' @param group Name of the grouping variable that identifies persons. Needed
#'   only when the formula has several plain grouping variables that could.
#' @param generalizability Logical. If `TRUE` (the default), random effects for
#'   facets crossed with persons are treated as facets of a generalizability
#'   study (see Details). If `FALSE`, every random effect counts as part of the
#'   person's score in each design cell.
#' @param newdata Optional data frame of design cells and persons to evaluate,
#'   completed as in [native_parameters()]. By default all observed cells are
#'   used. Supply it for models with continuous predictors.
#' @param contrast Name of a predictor. If given, the reliability of each
#'   person's difference between every level of this predictor and its first
#'   level is reported instead of the reliability per level.
#' @param retest Name of a predictor coding measurement occasions, entered with
#'   person-specific effects (for example `kappa ~ 0 + session + (0 + session |
#'   id)`). If given, test-retest statistics for every pair of occasions are
#'   returned instead of the reliability per occasion.
#' @param scale `"link"` (the default, also accepted as `"sampling"`) for the
#'   scale on which the parameter is sampled and its group-level effects are
#'   normal, or `"native"` for the scale of the model's parameters.
#'   Generalizability coefficients exist on the link scale only.
#' @param ndraws Number of posterior draws to use. Defaults to all.
#' @param prob Width of the credible intervals.
#' @param ... Further arguments passed to [brms::posterior_linpred()].
#'
#' @details
#' # What is estimated
#'
#' For each parameter and design cell, the person-level parameter
#' \eqn{\theta_i} has a posterior mean \eqn{m_i}, the score one would extract
#' for person \eqn{i}, and a posterior variance \eqn{e_i^2}, its squared
#' standard error of measurement. The reliability is
#' \deqn{\rho = 1 - \bar{e^2} / \bar{V},}
#' where \eqn{\bar{V}} is the average over draws of the variance of
#' \eqn{\theta_i} across persons. This is the empirical (marginal, or EAP)
#' reliability of item response theory (Bock & Mislevy, 1982; Adams, 2005), and
#' one minus the pooling factor of Gelman and Pardoe (2006). The credible
#' interval comes from evaluating the ratio in every draw.
#'
#' The estimate refers to the posterior means, which are shrunk towards the
#' group mean. When persons contribute different numbers of trials these scores
#' are more reliable than unpooled estimates would be. The error variance is
#' the average of the posterior variances, not the square of the average
#' posterior SD, which would overstate the reliability when precision differs
#' between persons. How much it differs is shown by the person-specific
#' reliabilities \eqn{1 - e_i^2 / \bar{V}} (Raudenbush & Bryk, 2002; Williams et
#' al., 2022), summarised in the `rel_person_*` columns and listed in
#' `attr(result, "persons")`.
#'
#' # Scale
#'
#' Reliability is not invariant under a nonlinear transformation, so the value
#' on the native scale (for example `kappa`) differs from the value on the log
#' scale on which `kappa` is sampled (Nakagawa & Schielzeth, 2010). The link
#' scale is the default because the group-level effects are normal and additive
#' there, and the correlations of a multivariate model are estimated there.
#'
#' # Generalizability theory
#'
#' With random facets, the parameter is decomposed into a person variance
#' \eqn{\sigma^2_p}, person-by-facet variances \eqn{\sigma^2_{pF}}, facet
#' variances \eqn{\sigma^2_F} and the estimation error \eqn{s^2} of the finest
#' person-level score. The relative G coefficient and the absolute
#' (dependability) coefficient are
#' \deqn{E\rho^2 = \sigma^2_p / (\sigma^2_p + \sum \sigma^2_{pF} / n_F +
#' s^2 / \prod n_F)}
#' \deqn{\Phi = \sigma^2_p / (\sigma^2_p + \sum \sigma^2_{pF} / n_F +
#' s^2 / \prod n_F + \sum \sigma^2_F / n_F)}
#' (Cronbach et al., 1972; Brennan, 2001), where \eqn{n_F} is the number of
#' levels of a facet per person (a harmonic mean over persons). The G
#' coefficient of the design as fitted is the posterior reliability of the
#' universe score, the prediction with only the person terms; \eqn{s^2} is
#' what remains of its error after the person-by-facet variances.
#' [reliability_for_design()] uses these components for decision studies.
#' Facet terms must be random intercepts. Facets with few levels give
#' prior-sensitive variance components (Li, 2026); with two or three sessions
#' prefer `retest`.
#'
#' # Test-retest
#'
#' With `retest`, the stability is the correlation of the persons' true scores
#' between two occasions, computed in every draw (Chen et al., 2021), and the
#' agreement is the corresponding ICC(A,1), which also penalises a change of
#' the mean (McGraw & Wong, 1996). Stability says how much the persons' true
#' scores change, not how precisely they are measured (Katahira et al., 2024).
#' The correlation one would observe between the occasions' scores is about
#' `stability * sqrt(reliability_1 * reliability_2)`, reported as
#' `expected_retest`. The model has to estimate the correlation between
#' occasions; with uncorrelated person effects (`||`) stability is biased
#' towards zero.
#'
#' @return A data frame of class `bmm_reliability` with one row per parameter
#'   and design cell (or contrast, or pair of occasions). Attributes `persons`
#'   and `components` hold the person-specific reliabilities and the variance
#'   components.
#'
#' @references
#' Adams, R. J. (2005). Reliability as a measurement design effect. *Studies in
#' Educational Evaluation, 31*(2-3), 162-172.
#'
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability in a
#' microcomputer environment. *Applied Psychological Measurement, 6*(4),
#' 431-444.
#'
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' Chen, G., Pine, D. S., Brotman, M. A., Smith, A. R., Cox, R. W., & Haller,
#' S. P. (2021). Trial and error: A hierarchical modeling approach to
#' test-retest reliability. *NeuroImage, 245*, 118647.
#'
#' Cronbach, L. J., Gleser, G. C., Nanda, H., & Rajaratnam, N. (1972). *The
#' dependability of behavioral measurements: Theory of generalizability for
#' scores and profiles*. Wiley.
#'
#' Gelman, A., & Pardoe, I. (2006). Bayesian measures of explained variance and
#' pooling in multilevel (hierarchical) models. *Technometrics, 48*(2), 241-251.
#'
#' Katahira, K., Oba, T., & Toyama, A. (2024). Does the reliability of
#' computational models truly improve with hierarchical modeling? Some
#' recommendations and considerations for the assessment of model parameter
#' reliability. *Psychonomic Bulletin & Review, 31*(6), 2465-2486.
#'
#' Li, G. (2026). Influence of uninformative prior distributions for MCMC
#' method on estimating variance components in generalizability theory.
#' *Applied Psychological Measurement, 50*(6), 325-355.
#'
#' McGraw, K. O., & Wong, S. P. (1996). Forming inferences about some intraclass
#' correlation coefficients. *Psychological Methods, 1*(1), 30-46.
#'
#' Nakagawa, S., & Schielzeth, H. (2010). Repeatability for Gaussian and
#' non-Gaussian data: A practical guide for biologists. *Biological Reviews,
#' 85*(4), 935-956.
#'
#' Raudenbush, S. W., & Bryk, A. S. (2002). *Hierarchical linear models:
#' Applications and data analysis methods* (2nd ed.). Sage.
#'
#' Williams, D. R., Martin, S. R., & Rast, P. (2022). Putting the individual
#' into reliability: Bayesian testing of homogeneous within-person variance in
#' hierarchical models. *Behavior Research Methods, 54*(3), 1272-1290.
#'
#' @seealso [reliability_for_design()], [native_parameters()]
#' @keywords extract_info
#' @export
bmm_reliability <- function(x, ...) {
  UseMethod("bmm_reliability")
}

#' @rdname bmm_reliability
#' @export
bmm_reliability.bmmfit <- function(x, pars = NULL, group = NULL,
                                   generalizability = TRUE, newdata = NULL,
                                   contrast = NULL, retest = NULL,
                                   scale = c("link", "native", "sampling"),
                                   ndraws = NULL, prob = 0.95, ...) {
  x <- restructure(x)
  scale <- match.arg(scale)
  if (scale == "sampling") {
    scale <- "link"
  }
  model <- x$bmm$model
  model_pars <- names(model$parameters)
  pars <- pars %||% setdiff(model_pars, names(model$fixed_parameters))
  unknown <- pars[not_in(pars, model_pars)]
  stopif(
    length(unknown) > 0,
    "Unknown parameter(s) {collapse_comma(unknown)}. \\
    The parameters of this model are {collapse_comma(model_pars)}."
  )
  stopif(
    !is.numeric(prob) || length(prob) != 1L || is.na(prob) || prob <= 0 || prob >= 1,
    "'prob' must be a single number between 0 and 1, not {prob}."
  )
  stopif(
    !is.null(ndraws) &&
      (!is.numeric(ndraws) || length(ndraws) != 1L || is.na(ndraws) ||
        ndraws < 2 || ndraws > brms::ndraws(x)),
    "'ndraws' must be a single number between 2 and {brms::ndraws(x)}."
  )
  stopif(
    !is.null(contrast) && !is.null(retest),
    "Supply either 'contrast' or 'retest', not both."
  )
  for (arg in c("contrast", "retest")) {
    value <- get(arg)
    stopif(
      !is.null(value) && !(is.character(value) && length(value) == 1),
      "'{arg}' must be the name of one predictor."
    )
  }
  stopif(
    !is.logical(generalizability) || length(generalizability) != 1 || is.na(generalizability),
    "'generalizability' must be TRUE or FALSE."
  )
  stopif(
    "re_formula" %in% names(list(...)),
    "bmm_reliability() chooses the group-level terms itself; remove 're_formula'."
  )
  draw_ids <- sort(if (is.null(ndraws)) seq_len(brms::ndraws(x)) else sample.int(brms::ndraws(x), ndraws))

  bterms <- brms::brmsterms(x$formula)
  group <- .rel_person_var(bterms, model_pars, group)
  settings <- nlist(
    x, model_pars, group, generalizability, newdata, contrast, retest, scale,
    draw_ids, prob, bterms, dots = list(...)
  )
  results <- lapply(pars, .rel_parameter, settings = settings)
  rows <- .rel_bind_rows(lapply(results, `[[`, "rows"))
  structure(
    rows,
    class = c("bmm_reliability", "data.frame"),
    fits = unlist(lapply(results, `[[`, "fits"), recursive = FALSE),
    persons = .rel_bind_rows(lapply(results, `[[`, "persons")),
    components = .rel_bind_rows(lapply(results, `[[`, "components")),
    group = group,
    scale = scale,
    type = if (!is.null(retest)) "retest" else if (!is.null(contrast)) "contrast" else "cell"
  )
}


# One parameter: classify its random effects, predict the person scores per
# cell, and compute the requested statistics for every cell.
.rel_parameter <- function(par, settings) {
  x <- settings$x
  group <- settings$group
  terms <- .rel_terms(settings$bterms, par, group)
  if (is.null(terms) || NROW(terms$person) == 0) {
    note <- if (par %in% names(x$bmm$model$fixed_parameters)) {
      "fixed to a constant"
    } else {
      glue::glue("no random effect over '{group}'")
    }
    return(list(rows = data.frame(parameter = par, note = note)))
  }

  g_study <- settings$generalizability && length(terms$facet_vars) > 0
  if (g_study && !terms$simple_facets) {
    message2(
      "The facet terms of '{par}' are not all random intercepts, so no \\
      generalizability study is done for it; its reliability is per cell."
    )
    g_study <- FALSE
  }
  # a variable with a person-specific slope defines cells even when it has no
  # population-level effect, as in c ~ 1 + (1 + cond | id)
  person_vars <- unlist(lapply(terms$person$form, all.vars))
  fixed_vars <- setdiff(union(.np_grid_vars(x, par, NA), person_vars), c(group, terms$facet_vars))
  clash <- intersect(fixed_vars, .rel_reserved_names)
  stopif(
    length(clash) > 0,
    "Predictor(s) {collapse_comma(clash)} share a name with a column of the \\
    output. Rename them in the data and refit, or pass 'newdata' with \\
    different names."
  )
  for (var in c(settings$contrast, settings$retest)) {
    stopif(
      !var %in% fixed_vars,
      "'{var}' is not a predictor of '{par}'. Its predictors are \\
      {collapse_comma(fixed_vars)}."
    )
  }
  # with a G-study the person score is the universe score (person terms only);
  # without one, every random effect is part of the person's score in a cell
  cell_vars <- c(fixed_vars, if (!g_study) terms$facet_vars)
  re_formula <- if (g_study) .rel_person_re_formula(terms$person) else NULL
  scores <- .rel_cell_scores(settings, par, c(cell_vars, group), re_formula)
  cells <- .rel_split(scores, cell_vars, group)
  # variance components exist on the link scale only, and the person x facet
  # intercepts cancel from a within-person difference, so contrasts of
  # universe scores get the posterior reliability alone
  g_components <- g_study && settings$scale == "link" && is.null(settings$contrast)
  if (g_components) {
    terms <- .rel_component_draws(settings, par, terms)
  }

  if (!is.null(settings$retest)) {
    return(.rel_retest_rows(par, cells, settings))
  }
  if (!is.null(settings$contrast)) {
    cells <- .rel_contrast_cells(cells, settings$contrast)
  }

  data <- x$data
  trials <- person_trials(x$bmm$model, data, x$formula)
  fits <- list()
  rows <- list()
  persons <- list()
  components <- list()
  for (k in seq_along(cells)) {
    cell <- cells[[k]]
    stopif(
      length(cell$persons) < 3,
      "A cell of '{par}' ({.rel_cell_label(cell$key)}) has {length(cell$persons)} \\
      person(s); reliability needs at least 3. A continuous predictor makes one \\
      cell per observed value: supply 'newdata' with the cells you need."
    )
    data_rows <- .rel_data_rows(data, cell$data_key %||% cell$key)
    K <- tapply(trials[data_rows], data[[group]][data_rows], sum)[as.character(cell$persons)]
    row <- data.frame(parameter = par, cell$key, n_persons = length(cell$persons))
    if (g_components && length(data_rows) > 0) {
      fit <- .rel_g_cell(settings, par, terms, cell, data, data_rows)
      row <- cbind(row, .rel_summary_row(fit$core, K, FALSE), fit$row)
      components[[k]] <- cbind(data.frame(parameter = par, cell$key), fit$components)
    } else {
      fit <- list(core = .rel_core(cell$theta, settings$prob))
      row <- cbind(row, .rel_summary_row(fit$core, K, TRUE))
    }
    fit$K <- K
    fit$key <- data.frame(parameter = par, cell$key)
    fit$facets <- if (g_study) terms$facet_vars
    fit$native_facets <- g_study && settings$scale == "native"
    fits[[k]] <- fit
    rows[[k]] <- row
    persons[[k]] <- data.frame(
      parameter = par, cell$key, person = cell$persons,
      estimate = fit$core$m, se = sqrt(fit$core$e2), reliability = fit$core$person,
      trials = as.numeric(K), row.names = NULL
    )
  }
  rows <- .rel_bind_rows(rows)
  rows$scale <- settings$scale
  if (g_study && settings$scale == "native") {
    rows$note <- "generalizability coefficients need scale = 'link'"
  }
  list(rows = rows, fits = fits, persons = .rel_bind_rows(persons), components = .rel_bind_rows(components))
}

.rel_reserved_names <- c(
  "parameter", "n_persons", "reliability", "lower", "upper", "note", "scale",
  "contrast", "person", "estimate", "se", "trials", "occasion_1", "occasion_2",
  "target", "absolute"
)

.rel_cell_label <- function(key) {
  if (ncol(key) == 0) {
    return("all observations")
  }
  paste(names(key), "=", vapply(key, as.character, character(1)), collapse = ", ")
}

# Draws x rows of the parameter for every observed combination of score_vars
.rel_cell_scores <- function(settings, par, score_vars, re_formula) {
  x <- settings$x
  newdata <- .np_newdata(x, score_vars, settings$newdata)
  predict_pars <- if (settings$scale == "native") settings$model_pars else par
  draws <- .np_draws(
    x, predict_pars, par, newdata, re_formula, settings$scale,
    settings$draw_ids, settings$dots
  )[[par]]
  list(draws = draws, grid = newdata[, score_vars, drop = FALSE])
}

# Split the score draws into one draws x persons matrix per cell, the persons
# in sorted order so that cells can be paired person by person
.rel_split <- function(scores, cell_vars, group) {
  grid <- scores$grid
  key <- if (length(cell_vars) == 0) rep("all", nrow(grid)) else interaction(grid[cell_vars], drop = TRUE, lex.order = TRUE)
  lapply(split(seq_len(nrow(grid)), key, drop = TRUE), function(cols) {
    persons <- grid[[group]][cols]
    stopif(
      anyDuplicated(persons) > 0,
      "Persons appear more than once in a cell. Supply 'newdata' with one row \\
      per person and cell."
    )
    ord <- order(as.character(persons))
    cols <- cols[ord]
    key <- grid[cols[1], cell_vars, drop = FALSE]
    row.names(key) <- NULL
    list(
      key = key,
      persons = persons[ord],
      theta = scores$draws[, cols, drop = FALSE]
    )
  }) |> unname()
}

# Rows of the model data in the cells of `key`, one cell per row of `key`
.rel_data_rows <- function(data, key) {
  vars <- intersect(names(key), names(data))
  if (length(vars) == 0) {
    return(seq_len(nrow(data)))
  }
  matches <- Reduce(`|`, lapply(seq_len(nrow(key)), function(i) {
    Reduce(`&`, lapply(vars, function(var) {
      as.character(data[[var]]) == as.character(key[[var]][i])
    }))
  }))
  which(matches)
}

# Summary columns of the posterior reliability of one cell. The signal-to-noise
# ratio is the between-person SD in units of the per-trial noise backed out of
# each person's posterior variance (Rouder & Haaf, 2019).
.rel_summary_row <- function(core, K, snr) {
  data.frame(
    reliability = core$reliability,
    lower = core$lower,
    upper = core$upper,
    sd_between = sqrt(core$Vbar),
    se_rms = sqrt(mean(core$e2)),
    rel_person_min = min(core$person),
    rel_person_median = stats::median(core$person),
    rel_person_max = max(core$person),
    trials_median = stats::median(K),
    snr = if (snr) sqrt(core$Vbar / mean(core$noise * K)) else NA_real_
  )
}

# Variance draws of the facet terms, read once per parameter
.rel_component_draws <- function(settings, par, terms) {
  as_component <- function(term) {
    c(term, list(var = .rel_sd_draws(settings$x, term$label, par, settings$draw_ids)))
  }
  for (term in terms$mains) {
    levels <- length(unique(settings$x$data[[term$facets[1]]]))
    warnif(
      length(term$facets) == 1 && levels < 5,
      "The facet '{term$label}' has only {levels} levels, so its variance \\
      component depends strongly on the prior (Li, 2026). With few sessions, \\
      consider coding them as a predictor with person-specific effects and \\
      using 'retest'."
    )
  }
  terms$interactions <- lapply(terms$interactions, as_component)
  terms$mains <- lapply(terms$mains, as_component)
  terms
}

.rel_g_cell <- function(settings, par, terms, cell, data, data_rows) {
  group <- settings$group
  interactions <- terms$interactions
  mains <- terms$mains
  n <- vapply(terms$facet_vars, function(f) {
    .rel_facet_levels(data, data_rows, group, f)
  }, numeric(1))
  g <- .rel_gstudy(cell$theta, interactions, mains, n, settings$prob)
  occasion <- .rel_occasion_reliability(settings, par, terms, cell, data_rows)
  g_row <- if (settings$scale == "link") {
    data.frame(
      g_relative = g$g_relative, g_relative_lower = g$g_relative_lower,
      g_relative_upper = g$g_relative_upper, g_absolute = g$g_absolute,
      g_absolute_lower = g$g_absolute_lower, g_absolute_upper = g$g_absolute_upper
    )
  } else {
    data.frame(g_relative = NA_real_, g_absolute = NA_real_)
  }
  g_row$reliability_occasion <- occasion
  g_row <- cbind(g_row, as.data.frame(as.list(stats::setNames(n, paste0("n_", names(n))))))
  summarise_var <- function(v) {
    q <- stats::quantile(v, c((1 - settings$prob) / 2, (1 + settings$prob) / 2), names = FALSE)
    c(stats::median(v), q)
  }
  comp <- rbind(
    c(summarise_var(g$core$V), 1),
    do.call(rbind, lapply(interactions, function(t) c(summarise_var(t$var), prod(n[t$facets])))),
    do.call(rbind, lapply(mains, function(t) c(summarise_var(t$var), prod(n[t$facets])))),
    c(g$s2, NA, NA, prod(n[g$residual_facets]))
  )
  components <- data.frame(
    component = c(
      group,
      vapply(interactions, `[[`, character(1), "label"),
      vapply(mains, `[[`, character(1), "label"),
      "error"
    ),
    type = c(
      "person", rep("person x facet", length(interactions)),
      rep("facet", length(mains)), "estimation error"
    ),
    variance = comp[, 1], lower = comp[, 2], upper = comp[, 3], divisor = comp[, 4]
  )
  list(core = g$core, g = g, row = g_row, components = components)
}

# Average reliability of the single-occasion scores: every random effect
# included, one cell per combination of facet levels
.rel_occasion_reliability <- function(settings, par, terms, cell, data_rows) {
  if (!is.null(settings$newdata)) {
    return(NA_real_)
  }
  score_vars <- c(names(cell$key), settings$group, terms$facet_vars)
  newdata <- settings$x$data[data_rows, , drop = FALSE]
  settings$newdata <- newdata[!duplicated(newdata[score_vars]), , drop = FALSE]
  scores <- .rel_cell_scores(settings, par, score_vars, NULL)
  sub_cells <- .rel_split(scores, terms$facet_vars, settings$group)
  sub_cells <- Filter(function(cell) length(cell$persons) > 2, sub_cells)
  if (length(sub_cells) == 0) {
    return(NA_real_)
  }
  mean(vapply(sub_cells, function(cell) .rel_core(cell$theta)$reliability, numeric(1)))
}

# Person-wise differences between each level of `contrast` and its first level,
# at every combination of the other cell variables
.rel_contrast_cells <- function(cells, contrast) {
  keys <- .rel_bind_rows(lapply(cells, `[[`, "key"))
  others <- setdiff(names(keys), contrast)
  strata <- if (length(others) == 0) rep("all", nrow(keys)) else interaction(keys[others], drop = TRUE)
  out <- list()
  for (stratum in split(seq_along(cells), strata, drop = TRUE)) {
    stopif(
      length(stratum) < 2,
      "'{contrast}' needs at least two levels within each cell \\
      ({.rel_cell_label(keys[stratum[1], others, drop = FALSE])})."
    )
    levels <- keys[[contrast]][stratum]
    order_levels <- if (is.factor(levels)) order(as.integer(levels)) else order(levels)
    stratum <- stratum[order_levels]
    ref <- cells[[stratum[1]]]
    for (k in stratum[-1]) {
      cell <- cells[[k]]
      common <- intersect(as.character(ref$persons), as.character(cell$persons))
      key <- cell$key[, others, drop = FALSE]
      key$contrast <- paste(cell$key[[contrast]], "-", ref$key[[contrast]])
      out[[length(out) + 1]] <- list(
        key = key,
        data_key = rbind(ref$key, cell$key),
        persons = cell$persons[as.character(cell$persons) %in% common],
        theta = cell$theta[, as.character(cell$persons) %in% common, drop = FALSE] -
          ref$theta[, as.character(ref$persons) %in% common, drop = FALSE]
      )
    }
  }
  out
}

.rel_retest_rows <- function(par, cells, settings) {
  retest <- settings$retest
  keys <- .rel_bind_rows(lapply(cells, `[[`, "key"))
  others <- setdiff(names(keys), retest)
  strata <- if (length(others) == 0) rep("all", nrow(keys)) else interaction(keys[others], drop = TRUE)
  probs <- c((1 - settings$prob) / 2, (1 + settings$prob) / 2)
  rows <- list()
  for (stratum in split(seq_along(cells), strata, drop = TRUE)) {
    stopif(
      length(stratum) < 2,
      "'{retest}' needs at least two occasions within each cell of '{par}'."
    )
    for (pair in utils::combn(stratum, 2, simplify = FALSE)) {
      a <- cells[[pair[1]]]
      b <- cells[[pair[2]]]
      common <- intersect(as.character(a$persons), as.character(b$persons))
      ta <- a$theta[, as.character(a$persons) %in% common, drop = FALSE]
      tb <- b$theta[, as.character(b$persons) %in% common, drop = FALSE]
      stopif(
        length(common) < 3,
        "Fewer than 3 persons were measured on both '{a$key[[retest]]}' and \\
        '{b$key[[retest]]}'."
      )
      rt <- .rel_retest(ta, tb)
      rel_a <- .rel_core(ta, settings$prob)$reliability
      rel_b <- .rel_core(tb, settings$prob)$reliability
      stab <- stats::quantile(rt$stability, probs, names = FALSE)
      agree <- stats::quantile(rt$agreement, probs, names = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        parameter = par, a$key[, others, drop = FALSE],
        occasion_1 = as.character(a$key[[retest]]),
        occasion_2 = as.character(b$key[[retest]]),
        n_persons = length(common),
        stability = stats::median(rt$stability), lower = stab[1], upper = stab[2],
        agreement = stats::median(rt$agreement),
        agreement_lower = agree[1], agreement_upper = agree[2],
        reliability_1 = rel_a, reliability_2 = rel_b,
        expected_retest = stats::median(rt$stability) * sqrt(rel_a * rel_b),
        scale = settings$scale
      )
    }
  }
  list(rows = .rel_bind_rows(rows))
}

.rel_key_in <- function(key, rows) {
  if (!all(names(key) %in% names(rows))) {
    return(FALSE)
  }
  any(Reduce(`&`, lapply(names(key), function(var) {
    as.character(rows[[var]]) %in% as.character(key[[var]]) & !is.na(rows[[var]])
  })))
}

# rbind() for data frames with different columns; missing columns become NA
.rel_bind_rows <- function(dfs) {
  dfs <- Filter(function(df) !is.null(df) && NROW(df) > 0, dfs)
  if (length(dfs) == 0) {
    return(NULL)
  }
  cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(df) {
    for (col in setdiff(cols, names(df))) df[[col]] <- NA
    df[cols]
  })
  out <- do.call(rbind, dfs)
  row.names(out) <- NULL
  out
}


#' Reliability of a planned design
#'
#' @description
#' Projects the reliability estimated by [bmm_reliability()] to a design with
#' more or fewer trials, or with other numbers of facet levels such as sessions
#' or stimuli (a decision study in generalizability theory), or finds the
#' number of trials that reaches a target reliability.
#'
#' @param object The result of [bmm_reliability()], per cell or per contrast.
#' @param trials Multipliers of the current number of trials per person and
#'   cell: `2` doubles every person's trials, `0.5` halves them.
#' @param ... Numbers of facet levels per person, named by facet, for example
#'   `session = c(2, 4)`. Applied to the parameters that have this facet.
#' @param target A reliability between 0 and 1. If given, the smallest trial
#'   multiplier that reaches it is returned for each combination of facet
#'   levels, instead of the reliability for each design.
#'
#' @details
#' Without facets, each person's per-trial noise is backed out of their
#' posterior variance under a normal approximation and divided by the new
#' number of trials, so persons keep their unequal trial counts. A single
#' Spearman-Brown correction of the average reliability would overstate the
#' gain whenever precision differs between persons. With facets, the error
#' variance \eqn{s^2} and the person-by-facet variances of the G-study are
#' divided by the new numbers of levels (Brennan, 2001).
#'
#' The projection assumes that information grows in proportion to trials and
#' that the between-person variance does not change with practice, which need
#' not hold for long tasks (Kucina et al., 2023). It is most trustworthy for
#' designs close to the one that was fitted. A target that no number of trials
#' can reach, because error from too few facet levels remains, is returned as
#' `Inf`.
#'
#' @return A data frame with one row per parameter, cell and design: the trial
#'   multiplier, the median number of trials per person it implies, and the
#'   projected reliability with its credible interval (the relative G
#'   coefficient where there are facets, with the absolute coefficient in
#'   `absolute`).
#'
#' @references
#' Brennan, R. L. (2001). *Generalizability theory*. Springer.
#'
#' Kucina, T., Wells, L., Lewis, I., de Salas, K., Kohl, A., Palmer, M. A.,
#' Sauer, J. D., Matzke, D., Aidman, E., & Heathcote, A. (2023). Calibration
#' of cognitive tests to address the reliability paradox for decision-conflict
#' tasks. *Nature Communications, 14*, 2234.
#'
#' @seealso [bmm_reliability()]
#' @keywords extract_info
#' @export
reliability_for_design <- function(object, trials = 1, ..., target = NULL) {
  stopif(
    !inherits(object, "bmm_reliability"),
    "'object' must be the result of bmm_reliability()."
  )
  stopif(
    identical(attr(object, "type"), "retest"),
    "reliability_for_design() needs the reliability per cell or per contrast, \\
    not test-retest statistics."
  )
  stopif(
    !is.numeric(trials) || length(trials) == 0 || any(is.na(trials) | trials <= 0),
    "'trials' must be positive multipliers of the current number of trials."
  )
  stopif(
    !is.null(target) &&
      (!is.numeric(target) || length(target) != 1 || is.na(target) || target <= 0 || target >= 1),
    "'target' must be a single number between 0 and 1."
  )
  facets <- list(...)
  stopif(
    length(facets) > 0 && (is.null(names(facets)) || any(names(facets) == "")),
    "Numbers of facet levels in '...' must be named by facet, for example session = 4."
  )
  # rows the user dropped by subsetting drop their posterior summaries too
  fits <- Filter(function(fit) .rel_key_in(fit$key, object), attr(object, "fits"))
  stopif(
    length(fits) == 0,
    "'object' has no rows with a reliability estimate to project."
  )
  known <- unique(unlist(lapply(fits, `[[`, "facets")))
  unknown <- setdiff(names(facets), known)
  stopif(
    length(unknown) > 0,
    "{collapse_comma(unknown)} {ifelse(length(unknown) == 1, 'is', 'are')} not a \\
    facet of any parameter. Facets: {if (length(known)) collapse_comma(known) else 'none'}."
  )
  .rel_bind_rows(lapply(fits, .rel_design_rows, trials = trials, facets = facets, target = target))
}

.rel_design_rows <- function(fit, trials, facets, target) {
  stopif(
    isTRUE(fit$native_facets),
    "Projections for '{fit$key$parameter}' need the variance components of \\
    its facets, which exist on the link scale only. Call bmm_reliability() \\
    with scale = \"link\"."
  )
  own <- facets[names(facets) %in% fit$facets]
  designs <- if (length(own) == 0) data.frame(row.names = 1) else expand.grid(own)
  project <- function(f, n) {
    if (is.null(fit$g)) {
      p <- .rel_project_personwise(fit$core, f)
      return(c(p$reliability, p$lower, p$upper, NA, NA, NA))
    }
    d <- .rel_dstudy(fit$g, n, f)
    c(d$g_relative, d$g_relative_lower, d$g_relative_upper,
      d$g_absolute, d$g_absolute_lower, d$g_absolute_upper)
  }
  trials_now <- stats::median(fit$K)
  rows <- lapply(seq_len(nrow(designs)), function(i) {
    n <- vapply(designs, function(levels) levels[i], numeric(1))
    with_key <- function(values) {
      out <- cbind(fit$key[rep(1, nrow(values)), , drop = FALSE], values)
      for (facet in names(n)) out[[paste0("n_", facet)]] <- n[[facet]]
      row.names(out) <- NULL
      out
    }
    if (!is.null(target)) {
      f <- .rel_solve_factor(function(f) project(f, n)[1], target)
      return(with_key(data.frame(
        target = target, trials = f, trials_per_person = f * trials_now
      )))
    }
    values <- t(vapply(trials, project, numeric(6), n = n))
    with_key(data.frame(
      trials = trials, trials_per_person = trials * trials_now,
      reliability = values[, 1], lower = values[, 2], upper = values[, 3],
      absolute = values[, 4], absolute_lower = values[, 5], absolute_upper = values[, 6]
    ))
  })
  out <- .rel_bind_rows(rows)
  if (is.null(fit$g)) {
    out <- out[not_in(names(out), c("absolute", "absolute_lower", "absolute_upper"))]
  }
  out
}


#' @export
print.bmm_reliability <- function(x, digits = 2, ...) {
  type <- attr(x, "type") %||% "cell"
  df <- structure(x, class = "data.frame")
  stats_cols <- switch(type,
    retest = c("n_persons", "stability", "lower", "upper", "agreement",
               "reliability_1", "reliability_2", "expected_retest"),
    c("n_persons", "reliability", "lower", "upper", "rel_person_min",
      "rel_person_max", "trials_median", "g_relative", "g_absolute",
      "reliability_occasion")
  )
  hidden <- c(
    "sd_between", "se_rms", "rel_person_median", "snr", "scale",
    "g_relative_lower", "g_relative_upper", "g_absolute_lower",
    "g_absolute_upper", "agreement_lower", "agreement_upper",
    grep("^n_", names(df), value = TRUE)
  )
  keys <- setdiff(names(df), c(stats_cols, hidden, "note"))
  shown <- intersect(c(keys, stats_cols, "note"), names(df))
  shown <- shown[vapply(df[shown], function(col) !all(is.na(col)), logical(1))]
  out <- df[shown]
  numeric_cols <- vapply(out, is.numeric, logical(1)) & names(out) != "n_persons"
  out[numeric_cols] <- lapply(out[numeric_cols], round, digits = digits)
  cat(glue::glue(
    "Reliability of person-level parameters (persons: '{attr(x, 'group') %||% '?'}', \\
    scale: {attr(x, 'scale') %||% '?'})"
  ), "\n\n")
  print(out, row.names = FALSE)
  notes <- switch(type,
    retest = c(
      "stability: correlation of the true scores between two occasions",
      "expected_retest: correlation expected between the two occasions' scores"
    ),
    c(
      "reliability: of the posterior-mean scores; [lower, upper] is the credible interval",
      "rel_person_min/max: range of the person-specific reliabilities (attr(x, \"persons\"))",
      if ("g_relative" %in% shown) c(
        "g_relative / g_absolute: generalizability to the universe score over the facets as fitted",
        "reliability_occasion: reliability of single-occasion scores",
        "variance components: attr(x, \"components\"); plan designs with reliability_for_design()"
      )
    )
  )
  cat("\n", paste0("* ", notes, collapse = "\n"), "\n", sep = "")
  invisible(x)
}
