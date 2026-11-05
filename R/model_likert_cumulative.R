############################################################################# !
# MODELS                                                                 ####
############################################################################# !

# The kernel's name by noise distribution: brms's cloglog link is
# P(Y <= k) = 1 - exp(-exp(thr_k - eta)), the gumbel_min cdf, so gumbel_min is
# the complementary log-log model and gumbel_max the log-log one
.likert_kernels <- c(
  normal = "probit",
  logistic = "logit, proportional odds",
  gumbel_min = "complementary log-log",
  gumbel_max = "log-log"
)

.likert_citations <- list(
  mccullagh = glue(
    "McCullagh, P. (1980). Regression models for ordinal data. Journal of the ",
    "Royal Statistical Society: Series B (Methodological), 42(2), 109-127. ",
    "https://doi.org/10.1111/j.2517-6161.1980.tb01109.x"
  ),
  burkner_vuorre = glue(
    "Bürkner, P.-C., & Vuorre, M. (2019). Ordinal regression models in ",
    "psychology: A tutorial. Advances in Methods and Practices in Psychological ",
    "Science, 2(1), 77-101. https://doi.org/10.1177/2515245918823199"
  ),
  burkner = glue(
    "Bürkner, P.-C. (2021). Bayesian item response modeling in R with brms ",
    "and Stan. Journal of Statistical Software, 100(5). ",
    "https://doi.org/10.18637/jss.v100.i05"
  ),
  johnson = glue(
    "Johnson, T. R. (2003). On the use of heterogeneous thresholds ordinal ",
    "regression models to account for individual differences in response style. ",
    "Psychometrika, 68(4), 563-583. https://doi.org/10.1007/BF02295612"
  )
)

.model_likert_cumulative <- function(response = NULL, item = NULL, n_categories = NULL,
                                     keying = NULL, dist = "normal",
                                     threshold_type = "log_distance",
                                     response_levels = NULL, links = NULL, call = NULL, ...) {
  parameters <- list(
    activation = paste0(
      "Activation: the mean of the latent response variable of a person on an ",
      "item, in units of the noise SD, measured from the centre of the threshold ",
      "set. Person and item terms both go here"
    )
  )
  default_priors <- list(
    activation = list(main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)")
  )
  param_links <- list(activation = "identity")

  threshold_parts <- .likert_threshold_parts(n_categories, threshold_type)
  parameters <- c(parameters, threshold_parts$parameters)
  default_priors <- c(default_priors, threshold_parts$default_priors)
  param_links <- c(param_links, threshold_parts$param_links)

  init_ranges <- list(activation = c(-0.5, 0.5))
  for (p in names(threshold_parts$parameters)) {
    init_ranges[[p]] <- if (p == "spacing") c(-0.3, 0.3) else c(-0.5, 0.2)
  }

  # rows 1-4 of the cell table: only the logistic parsimonious set is nobody's
  # named model, the equidistant one being the rating scale model
  if (dist == "logistic" && threshold_type == "parsimonious") {
    message2("likert_cumulative(): two-parameter logistic-quantile thresholds (bmm's \\
             parsimonious type); the equidistant type is the rating scale model")
  }

  citation <- switch(dist,
    normal = c(.likert_citations$johnson, .likert_citations$burkner_vuorre),
    logistic = c(.likert_citations$mccullagh, .likert_citations$burkner_vuorre,
                 .likert_citations$burkner),
    .likert_citations$burkner_vuorre
  )

  requirements <- glue(
    "Provide data with one row per person and item:", "\n\n",
    "  - Response (response): the category chosen, an integer from 1 to K ",
    "with n_categories = K, or an ordered factor", "\n",
    "  - Item (item): the item identifier", "\n",
    "  - Keying (keying, optional): 1 for regular and -1 for reversed items", "\n",
    "  - A person identifier for the formula, e.g. (1 | id)"
  )

  out <- structure(
    list(
      resp_vars = nlist(response, item, keying),
      other_vars = nlist(n_categories, dist, threshold_type, response_levels),
      domain = "Personality and individual differences",
      task = "Questionnaire and rating-scale responses",
      name = glue("Cumulative Response Model for Likert Items ({.likert_kernels[[dist]]})"),
      citation = citation,
      version = "NA",
      requirements = requirements,
      parameters = parameters,
      links = param_links,
      fixed_parameters = list(),
      default_priors = default_priors,
      init_ranges = init_ranges
    ),
    class = c("bmmodel", "likert_cumulative"),
    call = call
  )
  set_links(out, links)
}

# sdt_rating()'s threshold parameters; without a known K only the parts that
# do not depend on it (spacing) are declared, and check_model() rebuilds the
# model once the data give K. log_distance deltas get a K-matched prior mean:
# five intervals of median width 1 would leave a K = 7 item's outer categories
# 4% of the prior mass, while log(2 * log(K - 1) / (K - 2)) gives the same median
# span as the parsimonious set at spacing 0
.likert_threshold_parts <- function(n_categories, threshold_type) {
  if (is.null(n_categories)) {
    parts <- .threshold_parameter_parts(5L, threshold_type)
    return(lapply(parts, function(x) x[names(x) == "spacing"]))
  }
  parts <- .threshold_parameter_parts(n_categories, threshold_type)
  if (threshold_type == "log_distance") {
    centre <- round(log(2 * log(n_categories - 1) / (n_categories - 2)), 2)
    for (p in names(parts$default_priors)) {
      parts$default_priors[[p]] <- list(
        main = paste0("normal(", centre, ", 0.5)"), effects = "normal(0, 0.3)", sd = "exponential(2)"
      )
    }
  }
  parts
}

# activation is read on the natural scale by the family, so its link can change;
# spacing and the deltas are read through exp() by the threshold builders
#' @exportS3Method
settable_links.likert_cumulative <- function(model) {
  "activation"
}


#' @title Cumulative Response Model for Likert Items
#' @name likert_cumulative
#' @details `r model_docs(.model_likert_cumulative(n_categories = 5, threshold_type = "log_distance"))`
#'
#' @section The model:
#' Each response is one category out of K ordered ones. A latent response
#' variable with mean `activation` and noise of SD 1 is cut by K - 1 ordered
#' thresholds \eqn{\tau_1 < \dots < \tau_{K-1}}, and the response is the
#' interval it falls in:
#' \deqn{P(Y = k) = F(\tau_k - \eta) - F(\tau_{k-1} - \eta),}
#' with \eqn{\eta} the activation, \eqn{\tau_0 = -\infty},
#' \eqn{\tau_K = \infty}, and \eqn{F} the cdf of the noise distribution
#' (`dist`). This is the cumulative model (McCullagh, 1980; Bürkner & Vuorre,
#' 2019): `dist = "normal"` is the ordinal probit model, `"logistic"` the
#' proportional odds model, and the Gumbel distributions the complementary
#' log-log (`"gumbel_min"`, brms's `cloglog`) and log-log (`"gumbel_max"`)
#' models. Raw responses are never recoded: the categories are the ones the
#' person chose, from 1 (lowest) to K (highest), for regular and reversed
#' items alike.
#'
#' @section Thresholds:
#' The thresholds are built around 0 from the same parameters as in
#' [sdt_rating()], with the criterion fixed at 0: with an even K the middle
#' threshold sits at 0, with an odd K the middle category is centred on 0.
#' `threshold_type` takes the same five values. The default `"log_distance"`
#' estimates the log width of each of the K - 2 inner intervals
#' (`delta1` ... `delta{K-2}`); `"parsimonious"` and `"equidistant"` estimate
#' one `spacing`; `"log_ratio"` and `"softmax"` are as in [sdt_rating()].
#'
#' Thresholds are shared by all items unless the formula says otherwise.
#' `spacing ~ 0 + item` gives each item its own threshold spread: an
#' item-specific latent scale on a fixed location (a location-scale model;
#' Tutz, 2022). It is not item discrimination in Samejima's sense, which
#' multiplies the person's trait; this model fits no discrimination parameter.
#' `delta<i> ~ 0 + item` on a free-width type gives each item its own
#' threshold set, which with `dist = "logistic"` is the graded response model
#' with unit discrimination (Bürkner, 2021).
#'
#' Lubbe and Schuster's (2020) scaled threshold model with item-specific
#' thresholds, a person scaling factor on item-specific interval widths, is not
#' fitted: the free-width types carry no `spacing` to scale them. With shared
#' thresholds it is fitted, see "Response styles".
#'
#' @section Person and item terms:
#' Every person and item term goes into `activation`. With items as fixed
#' effects and persons as random effects:
#'
#' `bmf(activation ~ 0 + item + (1 | id))`
#'
#' Each item then has its own location and the person intercepts are the
#' trait. `activation ~ 1 + (1 | item) + (1 | id)` treats the items as a
#' sample. The thresholds are centred on 0, so a formula for `activation`
#' without an item term places every item at one location; [bmm()] warns about
#' this. The default prior on `activation` is `normal(0, 1)`, which with
#' `0 + item` applies to every item location.
#'
#' @section Response styles:
#' `spacing ~ 1 + (1 | id)` (with `"parsimonious"` or `"equidistant"`
#' thresholds) gives each person their own threshold dispersion: wide
#' thresholds pull responses to the middle categories, narrow ones push them
#' to the extremes. This is the proportional thresholds structure of Johnson
#' (2003) and the scaled threshold model of Lubbe and Schuster (2020) with
#' item discrimination fixed at 1, an extreme response style. Lubbe and Schuster
#' centre the threshold set at the scale middle so that a growing scale factor
#' pulls responses toward the middle categories; bmm's thresholds are centred
#' on 0 in the same way. `delta<i> ~ 1 + (1 | id)` on `"log_distance"` gives
#' each person their own interval widths, Johnson's (2003) heterogeneous
#' thresholds without his symmetry constraint, which can also express a
#' midpoint response style.
#'
#' @section Acquiescence as a shift:
#' With a `keying` column (1 for regular, -1 for reversed items), the shift
#' account of acquiescence is a formula:
#'
#' `bmf(activation ~ 0 + item + (0 + keying | id) + (1 | id))`
#'
#' The random slope on `keying` is the content trait, which moves responses up
#' on regular and down on reversed items; the random intercept moves responses
#' up on both, the acquiescence shift (Plieninger & Heck, 2018, eq. 14). This
#' needs both keyings in the data.
#'
#' @section Relation to brms::cumulative():
#' With free thresholds (`"log_distance"`, `"log_ratio"` or `"softmax"`), the
#' model is brms's `bf(y ~ 1 + (1 | item) + (1 | id), family = cumulative())`
#' with the probit, logit or cloglog link. Two things differ. brms estimates
#' the thresholds freely and the intercept against them, whereas bmm centres the
#' thresholds on 0 and moves the latent mean, so brms's thresholds equal bmm's
#' minus the intercept of `activation`. And bmm's threshold parameters are log
#' interval widths with priors on them, rather than brms's thresholds with a
#' flat or user prior.
#'
#' @section Not yet:
#' Mixtures for careless responding and acquiescence as a response process
#' (rather than a shift) are planned.
#'
#' @section Priors:
#' `activation` has `normal(0, 1)` on its intercept or item levels,
#' `normal(0, 0.5)` on effects, and `exponential(1)` on the SDs of its random
#' effects. `spacing` has `normal(0, 0.5)`, `normal(0, 0.3)` and
#' `exponential(2)`, as in [sdt_rating()]. The `"log_distance"` deltas have
#' `normal(m, 0.5)` with \eqn{m = \log(2 \log(K - 1) / (K - 2))}, rounded to two
#' decimals (0.09 at K = 4, -0.08 at K = 5, -0.33 at K = 7), and
#' `normal(0, 0.3)` and `exponential(2)` for effects and SDs. This is the one
#' default in bmm that depends on the data: it gives the thresholds the same
#' median span as the parsimonious type at `spacing = 0`, whereas a fixed
#' `normal(0, 1)` on every delta spreads the thresholds of a 7-point scale so
#' far that its outer categories receive about 4% of the prior predictive mass.
#' The `"log_ratio"` and `"softmax"` deltas keep [sdt_rating()]'s
#' `normal(0, 1)`.
#'
#' @param response The name of the response column: an integer category from 1
#'   to K, with `n_categories = K`, or an ordered factor, whose levels give K.
#'   Rows with a missing response are dropped with a warning.
#' @param item The name of the item column. Item terms go into the formula for
#'   `activation` (see "Person and item terms").
#' @param n_categories The number of response categories K, a whole number of
#'   at least 3. Required for integer responses, because bmm does not read K from
#'   the values observed; optional for an ordered factor, where it must equal the
#'   number of levels.
#' @param keying Optional name of a column holding 1 for regular and -1 for
#'   reversed items, constant within each item. It is checked and available to
#'   the formula (see "Acquiescence as a shift"); the responses are never
#'   recoded.
#' @param dist The noise distribution of the latent response variable:
#'   `"normal"` (default, cumulative probit), `"logistic"` (cumulative logit),
#'   `"gumbel_min"` (complementary log-log) or `"gumbel_max"` (log-log).
#' @param threshold_type The threshold parameterization, as in [sdt_rating()]:
#'   `"log_distance"` (default), `"parsimonious"`, `"equidistant"`,
#'   `"log_ratio"` (needs K of at least 4) or `"softmax"`.
#' @param links A named list of link functions. Only `activation` can be
#'   set; the threshold parameters are read through `exp()`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#' @references
#' Bürkner, P.-C. (2021). Bayesian item response modeling in R with brms and
#'   Stan. \emph{Journal of Statistical Software}, \emph{100}(5).
#'   \doi{10.18637/jss.v100.i05}
#'
#' Bürkner, P.-C., & Vuorre, M. (2019). Ordinal regression models in
#'   psychology: A tutorial. \emph{Advances in Methods and Practices in
#'   Psychological Science}, \emph{2}(1), 77--101.
#'   \doi{10.1177/2515245918823199}
#'
#' Johnson, T. R. (2003). On the use of heterogeneous thresholds ordinal
#'   regression models to account for individual differences in response
#'   style. \emph{Psychometrika}, \emph{68}(4), 563--583.
#'   \doi{10.1007/BF02295612}
#'
#' Lubbe, D., & Schuster, C. (2020). A scaled threshold model for measuring
#'   extreme response style. \emph{Journal of Educational and Behavioral
#'   Statistics}, \emph{45}(1), 86--107. \doi{10.3102/1076998619859541}
#'
#' McCullagh, P. (1980). Regression models for ordinal data. \emph{Journal of
#'   the Royal Statistical Society: Series B (Methodological)}, \emph{42}(2),
#'   109--127. \doi{10.1111/j.2517-6161.1980.tb01109.x}
#'
#' Plieninger, H., & Heck, D. W. (2018). A new model for acquiescence at the
#'   interface of psychometrics and cognitive psychology. \emph{Multivariate
#'   Behavioral Research}. \doi{10.1080/00273171.2018.1469966}
#'
#' Tutz, G. (2022). Ordinal regression: A review and a taxonomy of models.
#'   \emph{WIREs Computational Statistics}, \emph{14}(2), e1545.
#'   \doi{10.1002/wics.1545}
#' @keywords bmmodel
#' @export
#' @examples
#' \dontrun{
#' # 200 persons answer 10 five-point items
#' dat <- expand.grid(item = paste0("i", 1:10), id = 1:200)
#' theta <- rnorm(200)
#' beta <- seq(-1, 1, length.out = 10)
#' dat$y <- rlikert_cumulative(nrow(dat), theta[dat$id] + beta[as.integer(dat$item)],
#'                             thresholds = c(-1.5, -0.5, 0.5, 1.5))
#'
#' model <- likert_cumulative(response = "y", item = "item", n_categories = 5)
#'
#' fit <- bmm(
#'   formula = bmf(activation ~ 0 + item + (1 | id)),
#'   data = dat,
#'   model = model,
#'   cores = 4,
#'   backend = "cmdstanr"
#' )
#'
#' # person-specific threshold dispersion (extreme response style)
#' model_ers <- likert_cumulative("y", "item", n_categories = 5,
#'                                threshold_type = "parsimonious")
#' fit_ers <- bmm(bmf(activation ~ 0 + item + (1 | id), spacing ~ 1 + (1 | id)),
#'                data = dat, model = model_ers, cores = 4, backend = "cmdstanr")
#' }
likert_cumulative <- function(response, item, n_categories = NULL, keying = NULL,
                              dist = c("normal", "logistic", "gumbel_min", "gumbel_max"),
                              threshold_type = c("log_distance", "parsimonious", "equidistant",
                                                 "log_ratio", "softmax"),
                              links = NULL, ...) {
  call <- match.call()
  .check_constructor_dots(list(...), known = "response_levels")
  stop_missing_args()
  dist <- match.arg(dist)
  threshold_type <- match.arg(threshold_type)

  is_name <- function(x) is.character(x) && length(x) == 1 && !is.na(x)
  stopif(!is_name(response), "response must name one column of the data")
  stopif(!is_name(item), "item must name one column of the data")
  stopif(!is.null(keying) && !is_name(keying),
         "keying must be NULL or the name of one column of the data")
  stopif(
    !is.null(n_categories) &&
      !(is.numeric(n_categories) && length(n_categories) == 1 && is.finite(n_categories) &&
          n_categories %% 1 == 0 && n_categories >= 3),
    "n_categories must be NULL or a single whole number of at least 3. Two \\
    categories are a binary response, which a cumulative model does not need"
  )
  stopif(threshold_type == "log_ratio" && !is.null(n_categories) && n_categories < 4,
         "log_ratio thresholds require at least 4 response categories (the anchor \\
         ratio needs an interval on each side of the centre)")

  .model_likert_cumulative(response = response, item = item,
                           n_categories = if (!is.null(n_categories)) as.integer(n_categories),
                           keying = keying, dist = dist, threshold_type = threshold_type,
                           links = links, call = call, ...)
}


############################################################################# !
# CHECK_MODEL, CHECK_DATA AND CHECK_FORMULA METHODS                      ####
############################################################################# !

# An ordered factor gives K, which decides the threshold parameters, so a model
# built without n_categories is rebuilt here from the data and checked from the
# top: check_model.bmmodel() ran on the model before the rebuild. The formula
# is the user's only here; check_formula.bmmodel() has added `activation ~ 1`
# by the time the model's own check_formula() runs
#' @export
check_model.likert_cumulative <- function(model, data = NULL, formula = NULL) {
  stopif(is_bmmformula(formula) && !"activation" %in% names(formula),
         "{model$name} needs a formula for activation, which carries the person and \\
         item terms, e.g. activation ~ 0 + item + (1 | id)")
  resp_var <- model$resp_vars$response
  if (is.null(model$other_vars$n_categories) && is.data.frame(data) && resp_var %in% names(data)) {
    y <- data[[resp_var]]
    if (!is.ordered(y)) {
      .likert_response_codes(y, resp_var, NULL)
    }
    K <- nlevels(y)
    .likert_check_k(K, model$other_vars$threshold_type, resp_var)
    ov <- model$other_vars
    model <- .model_likert_cumulative(
      response = resp_var, item = model$resp_vars$item, n_categories = K,
      keying = model$resp_vars$keying, dist = ov$dist, threshold_type = ov$threshold_type,
      response_levels = levels(y), links = model$links[settable_links(model)],
      call = attr(model, "call")
    )
    return(check_model(model, data, formula))
  }
  NextMethod("check_model")
}

.likert_check_k <- function(K, threshold_type, resp_var) {
  stopif(K < 3,
         "The response variable '{resp_var}' has {K} categories; the model needs at \\
         least 3 categories. Two categories are a binary response")
  stopif(threshold_type == "log_ratio" && K < 4,
         "log_ratio thresholds require at least 4 response categories; '{resp_var}' \\
         has {K}")
}

# the category codes 1..K of the response; K is NULL until n_categories or the
# factor levels give it
.likert_response_codes <- function(y, resp_var, K) {
  if (is.factor(y)) {
    stopif(!is.ordered(y),
           "Response variable '{resp_var}' is a factor without an order. Its level \\
           order ({collapse_comma(levels(y))}) is alphabetical unless set, so make it \\
           an ordered factor with the levels from the lowest to the highest category, \\
           e.g. factor(x, levels = c(...), ordered = TRUE)")
    stopif(!is.null(K) && nlevels(y) != K,
           "n_categories = {K}, but the response variable '{resp_var}' is an ordered \\
           factor with {nlevels(y)} levels")
    return(as.integer(y))
  }
  stopif(!is.numeric(y),
         "Response variable '{resp_var}' must hold an integer category from 1 to \\
         {K %||% 'K'} or an ordered factor; it is {class(y)[1]}")
  stopif(is.null(K),
         "likert_cumulative() needs n_categories for integer responses: the number \\
         of categories is not read from the values observed. Give n_categories, or \\
         pass the response as an ordered factor")
  stopif(any(!is.na(y) & (y %% 1 != 0 | y < 1 | y > K)),
         "Response variable '{resp_var}' must hold whole numbers from 1 to {K} \\
         (n_categories = {K})")
  as.integer(y)
}

.likert_reserved_cols <- c("lc_resp", "lc_k", "lc_dist", "lc_thresh", "lc_item", "lc_key")

#' @exportS3Method
built_data_columns.likert_cumulative <- function(model) {
  c("lc_resp", "lc_k", "lc_dist", "lc_thresh", "lc_item",
    if (!is.null(model$resp_vars$keying)) "lc_key",
    NextMethod("built_data_columns"))
}

#' @export
check_data.likert_cumulative <- function(model, data, formula) {
  ov <- model$other_vars
  resp_var <- model$resp_vars$response
  item <- model$resp_vars$item
  keying <- model$resp_vars$keying

  reserved <- intersect(.likert_reserved_cols, colnames(data))
  warnif(length(reserved) > 0,
         "Column(s) {collapse_comma(reserved)} in your data are reserved by \\
         {model$name} and will be overwritten")

  stopif(!resp_var %in% colnames(data), "Response variable '{resp_var}' missing in the data")
  y <- data[[resp_var]]
  K <- ov$n_categories %||% if (is.ordered(y)) nlevels(y)
  codes <- .likert_response_codes(y, resp_var, K)
  .likert_check_k(K, ov$threshold_type, resp_var)

  missing_resp <- is.na(codes)
  warnif(any(missing_resp),
         "Dropping {sum(missing_resp)} rows with a missing response in '{resp_var}'")
  data <- data[!missing_resp, , drop = FALSE]
  codes <- codes[!missing_resp]

  unseen <- setdiff(seq_len(K), codes)
  warnif(length(unseen) > 0,
         "{if (length(unseen) == 1) 'Category' else 'Categories'} \\
         {paste(unseen, collapse = ', ')} of {K} {if (length(unseen) == 1) 'is' else 'are'} \\
         never observed in '{resp_var}'. The thresholds around \\
         {if (length(unseen) == 1) 'it' else 'them'} are then set by the prior and the \\
         neighbouring categories alone")

  stopif(!item %in% colnames(data), "Item variable '{item}' missing in the data")
  stopif(anyNA(data[[item]]), "Item variable '{item}' must not contain NA")
  if (!is.factor(data[[item]])) {
    data[[item]] <- factor(data[[item]])
  }

  if (!is.null(keying)) {
    stopif(!keying %in% colnames(data), "Keying variable '{keying}' missing in the data")
    key <- data[[keying]]
    stopif(!is.numeric(key) || anyNA(key) || any(!key %in% c(-1, 1)),
           "Keying variable '{keying}' must be 1 (regular) or -1 (reversed) for every \\
           row, without NA")
    varies <- tapply(key, data[[item]], function(x) length(unique(x)) > 1)
    varies <- names(varies)[!is.na(varies) & varies]
    stopif(length(varies) > 0,
           "Keying variable '{keying}' must be constant within each item; it varies \\
           within {collapse_comma(varies)}")
    data$lc_key <- as.integer(key)
  }

  data$lc_resp <- codes
  data$lc_k <- as.integer(K)
  data$lc_dist <- .noise_dist_id(ov$dist)
  data$lc_thresh <- match(ov$threshold_type, .threshold_types)
  # carried into the model frame so that update() can rebuild item and keying
  # when no formula keeps them
  data$lc_item <- as.integer(data[[item]])

  NextMethod("check_data")
}

#' @export
check_formula.likert_cumulative <- function(model, data, formula) {
  item <- model$resp_vars$item
  # an offset carries item locations known in advance
  has_offset <- any(grepl("offset(", deparse(formula$activation), fixed = TRUE))
  warnif(!item %in% rhs_vars(formula$activation) && !has_offset,
         "The formula for activation has no term for the item variable '{item}', so \\
         every item shares one location: thresholds and item locations are then \\
         confounded with the intercept. Add 0 + {item} or (1 | {item})")
  .likert_threshold_notes(model, formula)
  NextMethod("check_formula")
}

# rows 5-8 of the cell table: threshold terms that turn the model into one with
# its own name, or into one that is easily misread
.likert_threshold_notes <- function(model, formula) {
  item <- model$resp_vars$item
  has_group <- function(f) any(grepl("|", deparse(f), fixed = TRUE))
  type <- model$other_vars$threshold_type
  spacing <- formula$spacing
  if (!is.null(spacing) && item %in% rhs_vars(spacing)) {
    message2("spacing with an item term fits item-specific threshold spread (a \\
             location-scale model; Tutz, 2022). Discrimination in Samejima's sense \\
             multiplies the trait and is not fitted by this model")
  }
  if (!is.null(spacing) && type %in% c("parsimonious", "equidistant") && has_group(spacing)) {
    message2("spacing with a person term is the scaled threshold model with shared \\
             thresholds (Lubbe & Schuster, 2020); item discrimination is fixed at 1")
  }
  deltas <- formula[grepl("^delta[0-9]+$", names(formula))]
  free_width <- type %in% c("log_distance", "log_ratio")
  if (free_width && any(vapply(deltas, has_group, logical(1)))) {
    message2("delta with a person term gives person-specific interval widths after \\
             Johnson (2003), without his symmetry constraint")
  }
  if (model$other_vars$dist != "logistic" &&
        any(vapply(deltas, function(f) item %in% rhs_vars(f), logical(1)))) {
    message2("item-specific thresholds on a {.likert_kernels[[model$other_vars$dist]]} \\
             kernel; with dist = 'logistic' this is the graded response model")
  }
  invisible()
}

#' @exportS3Method
response_annotations.likert_cumulative <- function(model) {
  list(response = "one row per person and item, category 1..K (integer or ordered factor)")
}

#' @exportS3Method
model_column_roles.likert_cumulative <- function(model) {
  c(
    response = "the category chosen, an integer from 1 to K or an ordered factor; one row per person and item",
    item = "item identifier",
    n_categories = NA,
    keying = "1 for regular and -1 for reversed items (optional)"
  )
}


############################################################################# !
# FORMULA, FAMILY AND STAN CODE                                          ####
############################################################################# !

# brms refuses custom-family dpars that end in a digit, so each delta<i> the
# user writes is a non-linear parameter that an nlf() hands to the dpar
# delta<roman i>
.likert_threshold_dpars <- function(model) {
  pars <- setdiff(names(model$parameters), "activation")
  is_delta <- grepl("^delta[0-9]+$", pars)
  pars[is_delta] <- paste0("delta", as.character(utils::as.roman(sub("^delta", "", pars[is_delta]))))
  pars
}

#' @export
bmf2bf.likert_cumulative <- function(model, formula) {
  carried <- c("lc_item", if (!is.null(model$resp_vars$keying)) "lc_key")
  bform <- brms::bf(glue(
    "lc_resp | vint(lc_k, lc_dist, lc_thresh, {paste(carried, collapse = ', ')}) ~ 0"
  ))
  deltas <- .threshold_delta_names(model)
  dpars <- .likert_threshold_dpars(model)
  dpars <- dpars[dpars != "spacing"]
  for (j in seq_along(deltas)) {
    bform <- bform + brms::nlf(stats::as.formula(paste(dpars[j], "~", deltas[j])))
  }
  bform
}

#' @export
configure_model.likert_cumulative <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)
  threshold_dpars <- .likert_threshold_dpars(model)
  formula$family <- brms::custom_family(
    "likert_cumulative",
    dpars = c("mu", "activation", threshold_dpars),
    links = c("identity", model$links$activation, rep("identity", length(threshold_dpars))),
    type = "int",
    loop = TRUE,
    log_lik = log_lik_likert_cumulative,
    posterior_predict = posterior_predict_likert_cumulative,
    posterior_epred = posterior_epred_likert_cumulative,
    vars = c("vint1[n]", "vint2[n]", "vint3[n]")
  )

  sc_path <- system.file("stan_chunks", package = "bmm")
  chunks <- c("sdt_dist_functions.stan", "threshold_functions.stan", "ordinal_functions.stan",
              "likert_cumulative_functions.stan")
  stan_funs <- paste(
    paste(vapply(chunks, function(f) read_lines2(paste0(sc_path, "/", f)), character(1)),
          collapse = "\n"),
    .likert_cumulative_lpmf_stan(model),
    sep = "\n"
  )
  stanvars <- brms::stanvar(scode = stan_funs, block = "functions")

  nlist(formula, data, stanvars)
}

# The family's lpmf, generated per model because the number of deltas is fixed
# by K and the threshold type; K, dist and threshold type arrive as data
.likert_cumulative_lpmf_stan <- function(model) {
  dpars <- .likert_threshold_dpars(model)
  deltas <- dpars[dpars != "spacing"]
  nd <- length(deltas)
  signature <- paste(c("int y", "real mu", "real activation", if (length(dpars)) paste("real", dpars),
                       "int K", "int dist_type", "int thresh_type"), collapse = ", ")
  delta_decl <- paste0(
    "  array[", nd, "] real deltas;\n",
    if (nd) paste0(paste(glue("  deltas[{seq_len(nd)}] = {deltas};"), collapse = "\n"), "\n")
  )
  paste0(
    "real likert_cumulative_lpmf(", signature, ") {\n",
    delta_decl,
    "  return likert_cumulative_log_prob(y, activation, ",
    if ("spacing" %in% dpars) "spacing" else "0", ", deltas, K, dist_type, thresh_type);\n",
    "}\n"
  )
}

# thresholds around 0 from the threshold parameters, one row per draw; the R
# counterpart of make_thresholds(0, ...) in the family
.likert_cumulative_thresholds <- function(n_categories, threshold_type, spacing = NULL,
                                          deltas = NULL) {
  .make_thresholds(0, n_categories, threshold_type, spacing, deltas)
}


############################################################################# !
# SIMULATE_DATA                                                          ####
############################################################################# !

#' @export
simulate_data.likert_cumulative <- function(model, pars, design = NULL,
                                            scale = c("link", "native"), ...) {
  ov <- model$other_vars
  stopif(is.null(ov$n_categories),
         "simulate_data() needs the number of categories: build the model with \\
         n_categories")
  design <- .design_or_n(model, design, c(model$resp_vars$item, model$resp_vars$keying))
  sim <- .simulate_pars(model, pars, design, scale)
  p <- sim$native
  deltas <- p[grep("^delta[0-9]+$", names(p))]
  thresholds <- .likert_cumulative_thresholds(ov$n_categories, ov$threshold_type, p$spacing,
                                              if (length(deltas)) do.call(cbind, deltas))
  response <- rlikert_cumulative(nrow(design), p$activation, thresholds, dist = ov$dist)
  .simulate_result(design, stats::setNames(list(response), model$resp_vars$response), sim)
}


############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT, POSTERIOR_EPRED                            ####
############################################################################# !

# The K category log probabilities of row i (one row per draw), or of every
# row when i is NULL (draws within observations, as as.vector() of a draws x
# observations matrix)
.likert_log_probs <- function(prep, i = NULL) {
  K <- prep$data$vint1[1]
  dist <- .noise_dist_names[prep$data$vint2[1]]
  type <- .threshold_types[prep$data$vint3[1]]
  dpar <- function(name) as.vector(brms::get_dpar(prep, name, i = i))
  activation <- dpar("activation")
  n <- length(activation)
  names <- names(prep$dpars)
  spacing <- if ("spacing" %in% names) rep_len(dpar("spacing"), n)
  n_deltas <- sum(grepl("^delta", names))
  deltas <- if (n_deltas > 0) {
    delta_dpars <- paste0("delta", as.character(utils::as.roman(seq_len(n_deltas))))
    do.call(cbind, lapply(delta_dpars, function(d) rep_len(dpar(d), n)))
  }
  thresholds <- rbind(.likert_cumulative_thresholds(K, type, spacing, deltas))
  rbind(.ordinal_category_log_probs(thresholds, activation, 1, dist))
}

log_lik_likert_cumulative <- function(i, prep) {
  .likert_log_probs(prep, i)[, prep$data$Y[i]]
}

posterior_predict_likert_cumulative <- function(i, prep, ...) {
  probs <- exp(.likert_log_probs(prep, i))
  K <- ncol(probs)
  cum <- matrixStats::rowCumsums(probs)[, -K, drop = FALSE]
  as.integer(1L + rowSums(stats::runif(nrow(probs)) > cum))
}

# the expected category index, 1 to K
posterior_epred_likert_cumulative <- function(prep) {
  probs <- exp(.likert_log_probs(prep))
  .epred_matrix(as.vector(probs %*% seq_len(ncol(probs))), prep)
}
