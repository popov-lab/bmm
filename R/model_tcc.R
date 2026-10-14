############################################################################# !
# MODELS                                                                 ####
############################################################################# !

.model_tcc <- function(resp_error = NULL, similarity = NULL, n_options = 360,
                       units = "radians", d_max = 12, links = NULL, call = NULL, ...) {
  out <- structure(
    list(
      resp_vars = nlist(resp_error),
      # table = "r" passes the d table from R instead of building it in Stan;
      # an internal switch for table transforms (correlated or motor noise)
      other_vars = nlist(similarity, n_options, units, d_max, table = "stan"),
      domain = "Visual working memory",
      task = "Continuous reproduction",
      name = "Target Confusability Competition (TCC) model by Schurgin, Wixted & Brady (2020)",
      citation = glue(
        "Schurgin, M. W., Wixted, J. T., & Brady, T. F. (2020). Psychophysical \\
        scaling reveals a unified theory of visual memory strength. Nature Human \\
        Behaviour, 4(11), 1156-1172. https://doi.org/10.1038/s41562-020-00938-0"
      ),
      version = "NA",
      requirements = glue(
        "- The response is the error relative to the target, in radians (binned \\
        to the nearest of `n_options` options) or, with units = 'bins', in whole \\
        option steps\n",
        "- `similarity` gives the psychophysical similarity f of an option to the \\
        target at each distance 0, 1, ..., n_options / 2 options, scaled so that \\
        f = 1 at distance 0"
      ),
      parameters = list(
        d = glue(
          "Memory strength d': the mean memory-match signal of the target in units \\
          of the noise SD; an option at distance x has mean d * f(x)"
        )
      ),
      links = list(d = "identity"),
      fixed_parameters = list(),
      # provisional (#566): d of the uncorrelated model runs from about 0.6 to 2
      # over set sizes 1-8 in oberauer_lin_2017, and the authors tabulate d' up
      # to 4.5 for their correlated model, which they scale by ~0.65
      default_priors = list(
        d = list(main = "normal(2, 1.5)", effects = "normal(0, 1)", sd = "exponential(1)")
      ),
      init_ranges = list(d = c(1, 2))
    ),
    class = c("bmmodel", "tcc"),
    call = call
  )
  set_links(out, links)
}

#' @title `r .model_tcc()$name`
#' @name tcc
#' @description Continuous report as an `n_options`-alternative forced choice:
#'   every option on the response wheel produces a Gaussian memory-match
#'   signal with unit variance and mean `d * f(x)`, where `f(x)` is the
#'   psychophysical similarity of an option `x` options away from the target,
#'   and the response is the option with the strongest signal. This is the
#'   uncorrelated-noise model of Schurgin et al. (2020); their reported d' comes
#'   from a correlated-noise version and is about 0.65 times this `d`.
#' @details `r model_docs(.model_tcc())`
#'
#'   The likelihood reads the response probabilities from a table on a grid
#'   of `d` values from -2 to `d_max` in steps of 0.02, built once per chain
#'   by quadrature and interpolated with a cubic spline (error below 1e-6 in
#'   log probability). Beyond either end of the grid the table is extended
#'   linearly, which keeps the sampler moving during warmup but is not the
#'   model: raise `d_max` if the posterior of `d` reaches it.
#'
#'   The fit is fastest when rows that share a value of `d` are adjacent, so
#'   bmm offers to sort the data by the predictors of `d`, as for [sdm()].
#'   Continuous predictors of `d` are supported at the per-trial cost.
#' @param resp_error The name of the variable in the data holding the response
#'   error. With `units = "radians"` (the default) it is the error relative to
#'   the target in radians, as for the other circular models, and is binned to
#'   the nearest of the `n_options` options; with `units = "bins"` it is the
#'   offset of the chosen option from the target in whole option steps (any
#'   integer, taken modulo `n_options`).
#' @param similarity A numeric vector of `n_options / 2 + 1` values in \[0, 1\]:
#'   the similarity of an option to the target at each distance 0, 1, ...,
#'   `n_options / 2` options, starting with 1 at distance 0.
#' @param n_options The number of response options on the wheel, an even
#'   number: 360 for a colour wheel in 1-degree steps, 180 for orientation.
#' @param units Either `"radians"` or `"bins"`, see `resp_error`.
#' @param d_max The upper end of the `d` grid of the likelihood table.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @keywords bmmodel
#' @export
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' f <- exp(-(0:180) / 20)
#' f <- (f - min(f)) / (1 - min(f))
#' dat <- data.frame(err = rtcc(1000, d = 2, similarity = f))
#' fit <- bmm(bmf(d ~ 1), dat, tcc("err", similarity = f), cores = 4)
tcc <- function(resp_error, similarity, n_options = 360, units = c("radians", "bins"),
                d_max = 12, ...) {
  call <- match.call()
  stop_missing_args()
  units <- match.arg(units)
  .tcc_check_similarity(similarity, n_options)
  stopif(!is.numeric(d_max) || length(d_max) != 1 || !isTRUE(d_max > 0),
         "`d_max` must be a single positive number.")
  .model_tcc(resp_error = resp_error, similarity = similarity, n_options = n_options,
             units = units, d_max = d_max, call = call, ...)
}

.tcc_check_similarity <- function(similarity, n_options) {
  stopif(!is.numeric(n_options) || length(n_options) != 1 || !isTRUE(n_options >= 4) ||
           n_options %% 2 != 0,
         "`n_options` must be an even number of at least 4.")
  n_values <- n_options / 2 + 1
  stopif(!is.numeric(similarity), "`similarity` must be a numeric vector.")
  stopif(length(similarity) != n_values,
         "`similarity` must hold {n_values} values, one for each distance from 0 to \\
         {n_options / 2} options, but it holds {length(similarity)}.")
  stopif(anyNA(similarity) || any(similarity < 0 | similarity > 1),
         "`similarity` must lie between 0 and 1.")
  stopif(similarity[1] != 1,
         "`similarity` must be 1 at distance 0 (the target), as d is the target's signal.")
}

############################################################################# !
# CHECK_DATA S3 METHODS                                                  ####
############################################################################# !

#' @export
check_data.tcc <- function(model, data, formula) {
  resp <- model$resp_vars$resp_error
  units <- model$other_vars$units
  stopif(not_in(resp, colnames(data)),
         "The response variable '{resp}' is not present in the data.")
  stopif(units == "radians" && isTRUE(max(abs(data[[resp]]), na.rm = TRUE) > 2 * pi),
         "The response variable '{resp}' holds values beyond 2*pi, so it is not in \\
         radians. Convert it with deg2rad(), or set units = 'bins' if it counts options.")
  stopif(units == "bins" && any(data[[resp]] %% 1 != 0, na.rm = TRUE),
         "With units = 'bins', the response variable '{resp}' must hold whole numbers.")
  warnif("tcc_bin" %in% colnames(data),
         "Column 'tcc_bin' in your data is reserved by the TCC model and will be overwritten.")
  data$tcc_bin <- .tcc_bin(data[[resp]], model$other_vars$n_options, units)
  data <- order_data_query(model, data, formula)
  NextMethod("check_data")
}

# signed offset of the chosen option from the target, in (-n_options / 2, n_options / 2]
.tcc_bin <- function(x, n_options, units) {
  steps <- if (units == "radians") floor(x * n_options / (2 * pi) + 0.5) else x
  half <- n_options / 2
  as.integer(half - (half - steps) %% n_options)
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.tcc <- function(model, formula) {
  brms::bf("tcc_bin ~ 0")
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.tcc <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)
  formula$family <- brms::custom_family(
    "tcc",
    dpars = c("mu", "d"),
    links = c("identity", model$links$d),
    type = "int",
    loop = FALSE,
    vars = c("tcc_tab", "tcc_grid_lo", "tcc_grid_step")
  )
  stanvars <- .tcc_stanvars(model, model$other_vars$table)
  nlist(formula, data, stanvars)
}

# One kernel, two sources for its table: built in transformed data from the
# similarity, or computed in R and passed as data under the same name
.tcc_stanvars <- function(model, table) {
  grid <- .tcc_grid(model)
  similarity <- .tcc_similarity(model)
  sc_path <- system.file("stan_chunks", package = "bmm")
  out <- brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_funs.stan")),
                       block = "functions") +
    brms::stanvar(x = grid$lo, name = "tcc_grid_lo") +
    brms::stanvar(x = grid$step, name = "tcc_grid_step")
  if (table == "r") {
    return(out + brms::stanvar(x = .tcc_table(similarity, grid), name = "tcc_tab"))
  }
  rule <- .gl_rule(.tcc_gl_nodes)
  out +
    brms::stanvar(x = grid$rows, name = "tcc_grid_rows") +
    brms::stanvar(x = similarity, name = "tcc_similarity") +
    brms::stanvar(x = rule$nodes, name = "tcc_gl_nodes") +
    brms::stanvar(x = log(rule$weights), name = "tcc_gl_logw") +
    brms::stanvar(scode = read_lines2(file.path(sc_path, "tcc_tdata.stan")),
                  block = "tdata", pll_args = "data matrix tcc_tab")
}

############################################################################# !
# THE LIKELIHOOD TABLE                                                   ####
############################################################################# !

# The one place the model's similarity by distance is read, so that a
# parametric similarity function changes only this function
.tcc_similarity <- function(model) {
  model$other_vars$similarity
}

# The quadrature error peaks at small |d|, where the winning signal is a
# narrow bump: the probabilities sum to 1 within 6e-7 at 100 nodes, 1e-8 at
# 120 and 2e-10 at 140 (d from -2 to 12, n_options = 360). 120 nodes cost
# 1.2 times the table build and nothing per gradient. A d step of 0.02 adds
# about 1e-7 of interpolation error (TCC cost spike, 2026-10). Row 1 and the
# last row pad the grid by one step, so that every node inside it has a
# central-difference slope
.tcc_gl_nodes <- 120L

.tcc_grid <- function(model) {
  lo <- -2
  step <- 0.02
  rows <- as.integer(ceiling(round((model$other_vars$d_max - lo) / step, 6)) + 3)
  nlist(lo, step, rows)
}

# log P of responding one particular option at each distance 0, ..., n / 2
# from the target. With z_i = x - d f_i,
#   P(r) = int phi(z_r) prod_{i != r} Phi(z_i) dx = int phi(z_r) / Phi(z_r) G(x) dx,
# where G = prod_i Phi(z_i) is shared by all options, and log G is a sum over
# distances weighted by the number of options at each (1, 2, ..., 2, 1). The
# largest mean m = max(d f) kills every integrand below m - 8 through its own
# Phi factor, and phi(z_r) kills it above d f_r + 8 <= m + 8, so the window
# [m - 8, m + 8] holds all the mass whatever the sign and size of d.
# tcc_logp_by_distance() in tcc_funs.stan is the Stan twin
.tcc_logp <- function(d, similarity) {
  n_dist <- length(similarity)
  mult <- c(1, rep(2, n_dist - 2), 1)
  means <- d * similarity
  rule <- .gl_rule(.tcc_gl_nodes)
  z <- outer(max(means) + 8 * rule$nodes, means, "-")
  lphi <- stats::pnorm(z, log.p = TRUE)
  terms <- log(8 * rule$weights) + drop(lphi %*% mult) + stats::dnorm(z, log = TRUE) - lphi
  matrixStats::colLogSumExps(terms)
}

.tcc_table <- function(similarity, grid) {
  d_values <- grid$lo + (seq_len(grid$rows) - 2) * grid$step
  t(vapply(d_values, .tcc_logp, numeric(length(similarity)), similarity = similarity))
}

############################################################################# !
# STORED FRAME                                                           ####
############################################################################# !

# brms keeps the binned response but not the user's column, which is rebuilt
# at the bin centre
#' @exportS3Method
revert_check_data.tcc <- function(model, data) {
  resp <- model$resp_vars$resp_error
  if (!resp %in% colnames(data)) {
    data[[resp]] <- if (model$other_vars$units == "radians") {
      data$tcc_bin * 2 * pi / model$other_vars$n_options
    } else {
      data$tcc_bin
    }
    attr(data, "rebuilt") <- c(attr(data, "rebuilt"), resp)
  }
  data$tcc_bin <- NULL
  NextMethod("revert_check_data")
}
