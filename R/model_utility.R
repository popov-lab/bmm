############################################################################# !
# MODELS                                                                 ####
############################################################################# !
# A standalone multinomial choice model for utility-theoretic designs
# (random-utility / expected-utility / welfare weights). It shares the
# multinomial-logit likelihood with the m3 model but is an independent model
# family so its user surface speaks random-utility vocabulary and its pipeline
# is self-contained. The multinomial plumbing (response matrix, base formula,
# choice kernel) is a candidate for a future shared R/helpers-multinomial.R.

.model_utility <- function(resp_cats = NULL, num_options = NULL,
                           payoffs = NULL, value_cols = NULL,
                           utility_fn = "linear", weighting = "none",
                           choice_rule = "softmax", numeraire = NULL,
                           links = NULL, default_priors = NULL,
                           call = NULL, ...) {
  if (!is.null(num_options)) {
    names(num_options) <- names(num_options) %||% paste0("n_opt_", resp_cats)
  }

  token <- paste0(
    if (utility_fn == "power") "power" else "",
    if (weighting == "prelec") "prelec" else ""
  )
  if (token == "") token <- "linear"

  # numeraire fixes the scaling parameter: "b" (=> 1), c(name = value), or the
  # path default (welfare weights anchor at 1; value designs reuse the m3
  # background 0 for softmax / 0.1 for the positive simple rule)
  if (is.null(numeraire)) {
    num_name <- "b"
    num_value <- if (!is.null(payoffs)) 1 else if (choice_rule == "softmax") 0 else 0.1
  } else if (is.character(numeraire)) {
    num_name <- numeraire
    num_value <- 1
  } else {
    num_name <- names(numeraire)
    num_value <- unname(numeraire)
  }

  pars <- .utility_parameter_spec(payoffs, num_name, utility_fn, weighting)

  out <- structure(
    list(
      resp_vars = nlist(resp_cats),
      other_vars = nlist(num_options, choice_rule, payoffs, value_cols,
                         utility_fn, weighting),
      domain = "Value-based / economic decision making",
      task = "n-alternative categorical choice",
      name = "Random-Utility Choice Model",
      citation = glue(
        "Gross, J., Gotz, F., Reher, F., & Toscano, H. (2025). Nested social \\
        dilemmas. Communications Psychology, 3. McFadden, D. (1974). \\
        Conditional logit analysis of qualitative choice behavior."
      ),
      version = token,
      requirements = paste0(
        "- Provide the count columns for each response category (`resp_cats`).\n",
        "  - Supply either a `payoffs` matrix/data.frame (known coefficients) or\n",
        "    per-trial `value_cols` (observed value columns).\n",
        "  - Predict the utility parameters (welfare weights / gamma / rho / alpha)\n",
        "    at least by a fixed intercept, plus any predictors from your data.\n"
      ),
      parameters = pars$parameters,
      fixed_parameters = setNames(list(num_value), num_name),
      links = pars$links,
      default_priors = pars$default_priors,
      void_mu = FALSE
    ),
    class = c("bmmodel", "utility", paste0("utility_", token)),
    call = call
  )

  out$links[names(links)] <- links
  out$default_priors[names(default_priors)] <- default_priors
  out
}

# Parameter inventory, links, and default priors for the two entry paths. Kept
# separate from the constructor only because it is the piece a maintainer reads
# to understand the parameter set per version.
.utility_parameter_spec <- function(payoffs, num_name, utility_fn, weighting) {
  if (!is.null(payoffs)) {
    weight_pars <- setdiff(colnames(payoffs), num_name)
    parameters <- c(
      setNames(list("Numeraire weight. Fixed for scaling; anchors the utility unit."), num_name),
      setNames(
        as.list(glue("Welfare weight on payoff component '{weight_pars}'.")),
        weight_pars
      )
    )
    links <- setNames(as.list(rep("identity", length(weight_pars))), weight_pars)
    default_priors <- setNames(
      lapply(weight_pars, function(w) list(main = "normal(0, 1)", effects = "normal(0, 0.5)")),
      weight_pars
    )
    return(nlist(parameters, links, default_priors))
  }

  parameters <- setNames(
    list("Background activation. Fixed for scaling; added to each response category."),
    num_name
  )
  links <- list()
  default_priors <- list()

  parameters$gamma <- "Value sensitivity: marginal (dis)utility per unit value."
  links$gamma <- "identity"
  default_priors$gamma <- list(main = "normal(0, 1)", effects = "normal(0, 0.5)")

  if (utility_fn == "power") {
    parameters$rho <- "Power-utility curvature. Natural scale is exp(rho)."
    links$rho <- "log"
    default_priors$rho <- list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
  }
  if (weighting == "prelec") {
    parameters$alpha <- "Prelec probability-weighting curvature. Natural scale is exp(alpha)."
    links$alpha <- "log"
    default_priors$alpha <- list(main = "normal(0, 0.5)", effects = "normal(0, 0.5)")
  }

  nlist(parameters, links, default_priors)
}

#' @title `r .model_utility()$name`
#' @name utility
#'
#' @description
#' A random-utility choice model for multinomial (category-count) data. Each
#' response category carries a systematic utility built from either a known
#' `payoffs` matrix (welfare-weight designs, e.g. economic games) or per-trial
#' `value_cols` (expected-utility designs), and the category is chosen through a
#' softmax (multinomial-logit / random-utility) or simple (Luce ratio) rule.
#'
#' This model covers **category-count choice** — you observe how often each of a
#' fixed set of categories was chosen. It does **not** cover per-option
#' attribute / binary-lottery choice (prospect theory a la Nilsson et al., 2011);
#' that needs a different data interface and is out of scope.
#'
#' @param resp_cats Character vector of the data columns holding the response
#'   counts for each choice category.
#' @param num_options Either a single number / integer vector giving the number
#'   of candidates per category (constant across trials), or the names of the
#'   data columns holding per-trial candidate counts. Defaults to one candidate
#'   per category.
#' @param payoffs A numeric matrix (constant payoffs) or a data.frame whose cells
#'   are numbers (constant) or the names of per-trial data columns. Rows are named
#'   by `resp_cats`, columns by the weight parameters, and must include the
#'   numeraire column (default `"b"`). Generates the per-category activation
#'   formulas. Mutually exclusive with `value_cols`.
#' @param value_cols A named character vector mapping response categories to
#'   per-trial value columns, e.g. `c(corr = "V")`. Mutually exclusive with
#'   `payoffs`.
#' @param utility_fn Value function for the `value_cols` path: `"linear"`
#'   (`gamma * V`) or `"power"` (`gamma * V^rho`).
#' @param weighting Probability weighting for the `value_cols` path: `"none"` or
#'   `"prelec"` (adds `alpha`; requires `choice_rule = "softmax"`).
#' @param choice_rule `"softmax"` (random-utility logit) or `"simple"` (Luce ratio
#'   over raw utilities).
#' @param numeraire Fixes the scaling weight: `"b"` (fixed to 1), `c(b = value)`,
#'   or omitted for the path default.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#'
#' @details `r model_info(.model_utility(), components = c('domain', 'task', 'name', 'citation'))`
#'
#' Built on the same multinomial-logit likelihood as the [m3()] model.
#'
#' @keywords bmmodel
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # Welfare weights from a known payoff matrix (economic game)
#' model <- utility(
#'   resp_cats = c("nkeep", "ningroup", "nuniversal"),
#'   payoffs = rbind(
#'     nkeep      = c(b = 1.0, wi = 0.0, wo = 0.0),
#'     ningroup   = c(b = 0.5, wi = 1.0, wo = 0.0),
#'     nuniversal = c(b = 0.3, wi = 0.6, wo = 0.9)
#'   ),
#'   choice_rule = "simple",
#'   numeraire = c(b = 1)
#' )
#'
#' fit <- bmm(bmf(wi ~ 1, wo ~ 1), data = my_data, model = model)
#'
#' @export
utility <- function(resp_cats, num_options = NULL,
                    payoffs = NULL, value_cols = NULL,
                    utility_fn = c("linear", "power"),
                    weighting = c("none", "prelec"),
                    choice_rule = c("softmax", "simple"),
                    numeraire = NULL, ...) {
  call <- match.call()
  stop_missing_args()
  utility_fn <- match.arg(utility_fn)
  weighting <- match.arg(weighting)
  choice_rule <- match.arg(choice_rule)

  stopif(
    is.null(payoffs) == is.null(value_cols),
    "Provide exactly one of `payoffs` or `value_cols`."
  )
  if (is.null(num_options)) num_options <- rep(1L, length(resp_cats))
  stopif(
    length(num_options) != length(resp_cats),
    "`num_options` must have the same length as `resp_cats`."
  )

  num_name <- if (is.null(numeraire)) "b" else if (is.character(numeraire)) numeraire else names(numeraire)

  if (!is.null(payoffs)) {
    stopif(
      is.null(rownames(payoffs)) || is.null(colnames(payoffs)),
      "`payoffs` must have row names (the response categories) and column names \\
       (the weight parameters, including the numeraire '{num_name}')."
    )
    stopif(
      nrow(payoffs) != length(resp_cats),
      "`payoffs` must have one row per response category."
    )
    stopif(
      !num_name %in% colnames(payoffs),
      "`payoffs` must contain the numeraire column '{num_name}'."
    )
    # G5: curvature / weighting are undefined for known payoff coefficients
    stopif(
      utility_fn != "linear" || weighting != "none",
      "The `payoffs` path uses known coefficients, so `utility_fn` and \\
       `weighting` must stay at their defaults. Use `value_cols` for value \\
       curvature or probability weighting."
    )
  } else {
    stopif(
      is.null(names(value_cols)) || !all(names(value_cols) %in% resp_cats),
      "`value_cols` must be a named character vector whose names are response \\
       categories (a subset of `resp_cats`); categories without a value column \\
       carry only the background."
    )
  }
  # prelec replaces the additive count term; only defined for the softmax rule
  stopif(
    weighting == "prelec" && choice_rule != "softmax",
    "`weighting = 'prelec'` requires `choice_rule = 'softmax'`."
  )
  # G6: a mistyped numeraire name would otherwise fail deep in prior construction
  stopif(
    !is.null(numeraire) && !is.null(payoffs) && !num_name %in% colnames(payoffs),
    "The `numeraire` parameter '{num_name}' is not a column of `payoffs`."
  )

  .model_utility(
    resp_cats = resp_cats, num_options = num_options,
    payoffs = payoffs, value_cols = value_cols,
    utility_fn = utility_fn, weighting = weighting,
    choice_rule = choice_rule, numeraire = numeraire,
    call = call, ...
  )
}

############################################################################# !
# CHECK_MODEL S3 methods                                                 ####
############################################################################# !

#' @export
check_model.utility <- function(model, data = NULL, formula = NULL) {
  # G4: the Luce rule takes log() of the activation, which goes NaN if an
  # identity-linked parameter drives the activation negative during sampling
  ident <- names(model$links)[vapply(model$links, identical, logical(1), "identity")]
  ident <- setdiff(ident, names(model$fixed_parameters))
  warnif(
    model$other_vars$choice_rule == "simple" && length(ident) > 0,
    "The 'simple' (Luce) choice rule needs positive activations, but parameter(s) \\
    {collapse_comma(ident)} use the identity link and may turn the activation \\
    negative during sampling (log of a non-positive value is NaN). Add \\
    lower-bounded priors, e.g. set_prior('normal(0, 1)', nlpar = '{ident[1]}', lb = 0)."
  )

  NextMethod("check_model")
}

############################################################################# !
# CHECK_DATA S3 methods                                                  ####
############################################################################# !

#' @export
check_data.utility <- function(model, data, formula) {
  resp_cats <- model$resp_vars$resp_cats
  n_opt_vect <- model$other_vars$num_options
  col_names <- colnames(data)

  missing_cats <- setdiff(resp_cats, col_names)
  stopif(
    length(missing_cats),
    "The response count column(s) {collapse_comma(missing_cats)} are missing in the data."
  )

  # multinomial response matrix (shared plumbing -> future helpers-multinomial.R)
  resp_matrix <- as.matrix(data[resp_cats])
  resp_matrix[is.na(resp_matrix)] <- 0
  data <- data[!col_names %in% resp_cats]
  data$nTrials <- rowSums(resp_matrix)
  data$Y <- resp_matrix

  # candidate-count columns: names -> existing columns; numeric -> replicate
  if (is.character(n_opt_vect)) {
    missing_opts <- setdiff(n_opt_vect, col_names)
    stopif(
      length(missing_opts),
      "The candidate-count column(s) {collapse_comma(missing_opts)} are missing in the data."
    )
    opt_vars <- n_opt_vect
  } else {
    opt_vars <- names(n_opt_vect)
    stopif(
      any(opt_vars %in% names(data)),
      "Column name collision for candidate counts: {collapse_comma(opt_vars)}. \\
      Give explicit names to your `num_options` vector."
    )
    data[opt_vars] <- rep(n_opt_vect, each = nrow(data))
  }
  names(opt_vars) <- resp_cats

  stopif(
    any(colSums(data[opt_vars]) == 0),
    "At least one response category has zero candidates across all observations. \\
    Remove it from the model, as it is not identified."
  )

  o <- model$other_vars

  # Prelec probability weighting needs the total candidate count per trial
  if (o$weighting == "prelec") {
    stopif(
      "n_opt_total" %in% names(data),
      "Column name collision: 'n_opt_total' is reserved for the Prelec weighting."
    )
    data$n_opt_total <- rowSums(data[opt_vars])
  }

  idx_vars <- paste0("Idx_", resp_cats)
  data[idx_vars] <- as.integer(data[opt_vars] > 0)
  data[opt_vars][data[opt_vars] == 0] <- 0.0001

  # data-property validation: per-trial payoff / value columns must exist
  if (!is.null(o$payoffs)) {
    cells <- as.character(as.matrix(o$payoffs))
    payoff_cols <- unique(cells[is.na(suppressWarnings(as.numeric(cells)))])
    missing_pay <- setdiff(payoff_cols, col_names)
    stopif(
      length(missing_pay),
      "Per-trial payoff column(s) not found in the data: {collapse_comma(missing_pay)}."
    )
  }
  if (!is.null(o$value_cols)) {
    missing_val <- setdiff(unname(o$value_cols), col_names)
    stopif(
      length(missing_val),
      "Value column(s) not found in the data: {collapse_comma(missing_val)}."
    )
  }

  # G3: the value slope gamma is unidentified when value is constant
  if (!is.null(o$value_cols)) {
    for (col in unname(o$value_cols)) {
      stopif(
        length(unique(data[[col]])) < 2,
        "The value slope 'gamma' is not identifiable when value is constant: \\
        column '{col}' has a single value. Provide >= 2 distinct value levels."
      )
    }
  }

  # G1: power curvature rho additionally needs a wide value range
  if (o$utility_fn == "power") {
    for (col in unname(o$value_cols)) {
      vals <- unique(data[[col]])
      ratio <- max(vals) / min(vals)
      warnif(
        ratio < 3,
        "Power-utility 'rho' is weakly identified: value column '{col}' spans a \\
        range ratio of {round(ratio, 2)} (< 3), so rho and gamma can trade off. \\
        Use a wider range (e.g. values with max/min >= 3, such as 1, 3, 5) if rho \\
        is of interest."
      )
    }
  }

  # G2: Prelec alpha needs the choice-set proportions n_i/N to vary across trials
  if (o$weighting == "prelec") {
    counts <- as.matrix(data[opt_vars])
    props <- counts / rowSums(counts)
    prop_varies <- any(apply(props, 2, function(p) length(unique(round(p, 8))) > 1))
    stopif(
      !prop_varies,
      "Prelec 'alpha' is not identifiable: the choice-set proportions n_i/N do not \\
      vary across trials. Provide per-trial `num_options` columns whose relative \\
      sizes change (a constant set size makes alpha unidentifiable)."
    )
  }

  NextMethod("check_data")
}

############################################################################# !
# CHECK_FORMULA S3 methods                                               ####
############################################################################# !

#' @export
check_formula.utility <- function(model, data, formula) {
  formula <- construct_utility_act_funs(model) + formula
  formula <- apply_links(formula, model$links)
  formula <- assign_nl_attr(formula)

  NextMethod("check_formula")
}

# Generate the per-category activation formulas from the payoff matrix or the
# value columns. Payoff cells are emitted uniformly as `{cell}*{param}` whether
# the cell is a constant or the name of a per-trial data column. The Prelec
# weighting term enters the activation (so apply_links maps alpha -> exp(alpha))
# and replaces the log(n) count term that the softmax kernel would otherwise add.
construct_utility_act_funs <- function(model) {
  resp_cats <- model$resp_vars$resp_cats
  o <- model$other_vars
  num_name <- names(model$fixed_parameters)[1]

  if (!is.null(o$payoffs)) {
    payoff_matrix <- as.matrix(o$payoffs)
    param_cols <- colnames(payoff_matrix)
    forms <- lapply(seq_len(nrow(payoff_matrix)), function(k) {
      terms <- vapply(seq_along(param_cols), function(j) {
        cell <- as.character(payoff_matrix[k, j])
        coef <- suppressWarnings(as.numeric(cell))
        if (!is.na(coef) && coef == 0) NA_character_ else glue("{cell}*{param_cols[j]}")
      }, character(1))
      terms <- terms[!is.na(terms)]
      stats::formula(glue("{resp_cats[k]} ~ {paste(terms, collapse = ' + ')}"))
    })
  } else {
    opt_names <- if (is.character(o$num_options)) o$num_options else names(o$num_options)
    names(opt_names) <- resp_cats
    forms <- lapply(resp_cats, function(cat) {
      rhs <- num_name
      if (cat %in% names(o$value_cols)) {
        value <- if (o$utility_fn == "power") glue("{o$value_cols[[cat]]}^rho") else o$value_cols[[cat]]
        rhs <- glue("{rhs} + gamma*{value}")
      }
      if (o$weighting == "prelec") {
        rhs <- glue("{rhs} + (-((-log({opt_names[cat]} / n_opt_total))^alpha))")
      }
      stats::formula(glue("{cat} ~ {rhs}"))
    })
  }

  do.call(bmf, forms)
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.utility <- function(model, formula) {
  o <- model$other_vars
  resp_cats <- model$resp_vars$resp_cats
  opt_names <- if (is.character(o$num_options)) o$num_options else names(o$num_options)
  idx_vars <- paste0("Idx_", resp_cats)
  names(idx_vars) <- resp_cats
  names(opt_names) <- resp_cats

  cat <- resp_cats[1]
  brms_formula <- brms::bf(glue(
    "Y | trials(nTrials) ~
    {idx_vars[cat]} *", utility_choice_kernel(o$choice_rule, o$weighting, cat, opt_names),
    "+ (1 - {idx_vars[cat]}) * (-100)"
  ), nl = TRUE)

  for (cat in resp_cats[-1]) {
    brms_formula <- brms_formula + glue_nlf(
      "mu{cat} ~
      {idx_vars[cat]} *", utility_choice_kernel(o$choice_rule, o$weighting, cat, opt_names),
      "+ (1 - {idx_vars[cat]}) * (-100)"
    )
  }

  brms_formula
}

# Systematic-utility -> linear-predictor kernel. softmax adds log(n) (McFadden
# RUM), simple wraps in log (Luce ratio, needs positive activation), and prelec
# passes the activation through because the count structure already entered via
# the Prelec probability-weighting term.
utility_choice_kernel <- function(choice_rule, weighting, cat, options_vars) {
  if (weighting == "prelec") {
    return(glue("({cat})"))
  }
  switch(
    choice_rule,
    simple = glue("log({cat} * {options_vars[cat]})"),
    softmax = glue("({cat} + log({options_vars[cat]}))")
  )
}

############################################################################# !
# CONFIGURE_MODEL S3 methods                                             ####
############################################################################# !

#' @export
configure_model.utility <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- model$resp_vars$resp_cats
  formula$family$dpars <- paste0("mu", model$resp_vars$resp_cats)

  if (model$other_vars$choice_rule == "simple" && any(model$links == "identity")) {
    init <- 0
  } else {
    init <- NULL
  }

  nlist(formula, data, init)
}
