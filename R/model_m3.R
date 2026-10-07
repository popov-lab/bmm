############################################################################# !
# MODELS                                                                 ####
############################################################################# !
# see file 'R/bmm_model_mixture3p.R' for an example
# Define lookup tables for parameters, links, and default priors
.m3_version_table <- list(
  ss = list(
    parameters = list(
      c = "Context activation. Added to the item cued to be recalled, that is the correct item.",
      a = "General activation. Added to all items that were presented during the current trial."
    ),
    links = list(
      simple = list(c = "log", a = "log"),
      softmax = list(c = "identity", a = "identity"),
      gaussian = list(c = "identity", a = "identity")
    ),
    priors = list(
      simple = list(
        a = list(main = "normal(0,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        c = list(main = "normal(3,1)", effects = "normal(0,0.5)", sd = "exponential(1)")
      ),
      softmax = list(
        a = list(main = "normal(3,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        c = list(main = "normal(3,1)", effects = "normal(0,0.5)", sd = "exponential(1)")
      ),
      gaussian = list(
        a = list(main = "normal(2,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        c = list(main = "normal(2,1)", effects = "normal(0,0.5)", sd = "exponential(1)")
      )
    )
  ),
  cs = list(
    parameters = list(
      c = "Context activation. Added to the item cued to be recalled, that is the correct item.",
      a = "General activation. Added to all items that were presented during the current trial.",
      f = "Filtering. This parameter captures the extent to which distractors remained in working memory."
    ),
    links = list(
      simple = list(c = "log", a = "log", f = "logit"),
      softmax = list(c = "identity", a = "identity", f = "logit"),
      gaussian = list(c = "identity", a = "identity", f = "logit")
    ),
    priors = list(
      simple = list(
        a = list(main = "normal(0,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        c = list(main = "normal(3,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        f = list(main = "logistic(0,1)", effects = "normal(0,1)", sd = "exponential(1)")
      ),
      softmax = list(
        a = list(main = "normal(3,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        c = list(main = "normal(3,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        f = list(main = "logistic(0,1)", effects = "normal(0,1)", sd = "exponential(1)")
      ),
      gaussian = list(
        a = list(main = "normal(2,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        c = list(main = "normal(2,1)", effects = "normal(0,0.5)", sd = "exponential(1)"),
        f = list(main = "logistic(0,1)", effects = "normal(0,1)", sd = "exponential(1)")
      )
    )
  )
)

.model_m3 <- function(resp_cats = NULL, num_options = NULL,
                      choice_rule = "softmax", version = "custom", links = NULL,
                      default_priors = NULL, call = NULL, ...) {
  if(!is.null(num_options)) names(num_options) <- names(num_options) %||% paste0("n_opt_",resp_cats)
  if(!is.character(choice_rule)) choice_rule <- as.character(choice_rule)
  out <- structure(
    list(
      resp_vars = nlist(resp_cats),
      other_vars = nlist(num_options, choice_rule),
      domain = "Working Memory (categorical), Categorical Decision Making",
      task = "n-AFC retrieval",
      name = "The Multinomial / Memory Measurement Model",
      citation = glue(
        "Oberauer, K., & Lewandowsky, S. (2019). Simple measurement models \\
        for complex working-memory tasks. Psychological Review, 126(6), \\
        880-932. https://doi.org/10.1037/rev0000159"
      ),
      version = version,
      requirements = paste0(
        '- Provide names for variables specifying the number of responses in a set of response categories.\n',
        '  - Specify activation sources for each response categories\n',
        '  - Include at least an activation source "b" for all response categories\n',
        '  - Predict the specified activation at least by a fixed intercept and any additional predictors from your data\n'
      ),
      parameters = c(
        list(b = "Background activation. Added to each response category. Fixed for scaling, necessary in all models."),
        .m3_version_table[[version]][["parameters"]]
      ),
      links = .m3_version_table[[version]][["links"]][[choice_rule]],
      fixed_parameters = list(
        b = if (choice_rule == "simple") 0.1 else 0
      ),
      default_priors = .m3_version_table[[version]][["priors"]][[choice_rule]]
    ),
    class = c("bmmodel", "m3", paste0("m3_", version)),
    call = call
  )

  out <- set_links(out, links)
  out$default_priors[names(default_priors)] <- default_priors
  out
}

# the parameters of a custom m3 are the activation sources of the user's
# formula, so there is no set of names to check a link target against
# (check_model.m3_custom refuses a parameter left without a link). The ss and
# cs versions build their activation functions from the version table, so their
# parameters are known here. The custom branch is defensive rather than
# load-bearing: a custom m3 has no links at construction, so names() is already
# NULL, and check_links() never runs on one because set_links() stored no
# attribute.
#' @exportS3Method
settable_links.m3 <- function(model) {
  if (model$version == "custom") NULL else names(model$links)
}

# m3 is the one model that applies its links itself, by substituting the
# inverse link into the activation formulas (apply_links -> inv_link), so the
# links it can honour are inv_link()'s, not the ones a brms family can emit
#' @exportS3Method
settable_link_functions.m3 <- function(model) {
  eval(formals(inv_link)$link)
}


# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_M3()$info

#' @title `r .model_m3()$name`
#' @name m3
#'
#' @description
#' The Multinomial / Memory Measurement Model (M3) is a measurement model that was originally introduced
#' for working memory tasks with categorical responses. It assumes that each candidate in each response
#' category is activated by a combination of sources of activation. The probability of choosing a response
#' category is determined by the activation of the candidates. The model can be used for any n-AFC categorical
#' decision task.
#'
#' @param resp_cats The variable names that contain the number of responses for each of the
#'   response categories used for the M3.
#' @param num_options Either an integer vector of the same length as `resp_cats` if the number
#'   of candidates in the respective response categories are constant across all conditions
#'   in the experiment. Or a vector specifying the variable names that contain the number of
#'   candidates in each response category. The order of these variables should be in the
#'   same order as the names of the response categories passed to `resp_cats`. Numbers
#'   named after the response categories, e.g. `c(corr = 1, other = 4)`, are matched to
#'   the categories by name. Numbers without names, or with other names, are taken in
#'   the order of `resp_cats`, and other names become the names of the columns
#'   bmm adds to the data. Column names given category names, e.g.
#'   `c(other = "n_other", corr = "n_corr")`, are matched by name as well, in any order.
#'   Custom activation formulas can use numbers by these column names; numbers
#'   without names or named after the categories are called `n_opt_<category>`.
#' @param choice_rule The choice rule that should be used for the M3. The options are "softmax",
#'   "simple", or "gaussian". The "softmax" option implements the softmax normalization of activation into
#'   probabilities for choosing the different response categories. The "simple" option implements
#'   a simple normalization of the absolute activations over the sum of all activations. For details
#'   on the differences of these choice rules please see the appendix of Oberauer & Lewandowsky (2019)
#'   "Simple measurement models for complex working memory tasks" published in Psychological Review.
#'   The "gaussian" option adds independent standard normal noise to the activation of every
#'   candidate and chooses the candidate with the largest value (a Thurstonian rule); "softmax" is
#'   the same with Gumbel noise. Activations under "gaussian" are on a smaller scale, and the
#'   difference is not one constant factor: `a` shrinks more than `c`, so a comparison of `c`
#'   with `a`, and effects of conditions that change the number of candidates (such as set
#'   size), can differ between the two rules. Fits with "gaussian" take much longer than with
#'   "softmax", because each probability is a numerical integral.
#' @param version Character. The version of the M3 model to use. Can be one of
#'  `ss`, `cs`, or `custom`. The default is `custom`.
#' @param ... used internally for testing, ignore it
#' @return An object of class `bmmodel`
#'
#' @details `r model_docs(.model_m3(), components =c('domain', 'task', 'name', 'citation'))`
#' #### Version: `ss`
#' `r model_docs(.model_m3(version = "ss"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `cs`
#' `r model_docs(.model_m3(version = "cs"), components =c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Version: `custom`
#' `r model_docs(.model_m3(version = "custom"), components = c('requirements', 'parameters', 'fixed_parameters', 'links', 'prior'))`
#' #### Missing values and reserved names
#' A missing response count (`NA`) is counted as 0. If the category has options in that
#' row, `bmm()` warns and says how many counts were replaced; if it has none (`num_options`
#' is 0 in that row), the 0 is true and there is no warning, as for `dist` in
#' [oberauer_lewandowsky_2019_e1]. A missing value in a column named in `num_options` is
#' an error: enter 0 where there were no options. Data columns named `Y`, `nTrials` or
#' `Idx_<category>` are refused because `bmm()` creates columns with these names; response
#' categories may be called `Y` or `nTrials`.
#'
#' @keywords bmmodel
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' data <- oberauer_lewandowsky_2019_e1
#'
#' # initiate the model object
#' m3_model <- m3(
#'   resp_cats = c("corr", "other", "dist", "npl"),
#'   num_options = c("n_corr", "n_other", "n_dist", "n_npl"),
#'   choice_rule = "simple"
#' )
#'
#' # specify the model formula including the activation formulas for each response category
#' m3_formula <- bmf(
#'   corr ~ b + a + c,
#'   other ~ b + a,
#'   dist ~ b + d,
#'   npl ~ b,
#'   c ~ 1 + cond + (1 + cond | ID),
#'   a ~ 1 + cond + (1 + cond | ID),
#'   d ~ 1 + (1 | ID)
#' )
#'
#' # specify links for the model parameters
#' m3_model$links <- list(
#'   c = "log",
#'   a = "log",
#'   d = "log"
#' )
#'
#' # check if the default priors are applied correctly
#' default_prior(m3_formula, data = data, model = m3_model)
#'
#' # fit the model
#' m3_fit <- bmm(
#'   formula = m3_formula,
#'   data = data,
#'   model = m3_model,
#'   cores = 4
#' )
#'
#' # print summary of the model
#' summary(m3_fit)
#'
#' @export
m3 <- function(resp_cats, num_options, choice_rule = "softmax",
               version = c("custom", "ss", "cs"), ...) {
  call <- match.call()
  stop_missing_args()
  version <- match.arg(version)
  stopif(
    !tolower(choice_rule) %in% c("softmax", "simple", "gaussian"),
    'Unsupported choice rule "{choice_rule}". Must be one of "simple", "softmax" or "gaussian"'
  )
  stopif(
    length(num_options) != length(resp_cats),
    "The option variables should have the same length as the response variables."
  )
  stopif(
    is.character(num_options) && any(num_options %in% resp_cats),
    "The number of options cannot be read from a response category column: \\
    {collapse_comma(intersect(num_options, resp_cats))}"
  )
  stopif(
    is.numeric(num_options) && anyNA(num_options),
    "`num_options` cannot contain missing values."
  )
  opt_names <- names(num_options)
  stopif(
    !is.null(opt_names) &&
      (anyNA(opt_names) || any(opt_names == "") || anyDuplicated(opt_names) > 0),
    "Name either all elements of `num_options` or none, and use each name only once."
  )
  stopif(
    any(opt_names %in% resp_cats) && !setequal(opt_names, resp_cats),
    "If `num_options` is named after the response categories, it needs one element for each of \\
    {collapse_comma(resp_cats)}"
  )

  .model_m3(
    resp_cats = resp_cats, num_options = num_options,
    choice_rule = tolower(choice_rule), version = version, call = call, ...
  )
}

############################################################################# !
# CHECK_Model S3 methods                                                 ####
############################################################################# !

#' @export
check_model.m3_custom <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    user_pars <- setdiff(
      m3_activation_symbols(model, formula),
      c(colnames(data), built_data_columns(model))
    )
    # a symbol without its own formula is more often a typo or a missing column
    # than a new parameter, and as a parameter it would be fitted silently
    no_formula <- setdiff(user_pars, names(formula))
    stopif(
      length(no_formula) > 0,
      "{collapse_comma(no_formula)} in your activation formula(s) is neither a \\
      data column nor a model parameter. Give each new parameter its own \\
      formula (e.g. {no_formula[1]} ~ 1), or add the column to the data (`Y` is \\
      reserved and cannot be a data column)."
    )
    model$parameters <- c(model$parameters, setNames(user_pars, user_pars))
  }

  missing_links <- setdiff(names(model$parameters), names(model$links))
  missing_links <- setdiff(missing_links, names(model$fixed_parameters))
  stopif(
    length(missing_links) > 0,
    "Please provide link functions for all model parameters via the `link` argument of `m3()` \\
     to ensure proper identification of your model.
     The following parameters are missing link functions: {paste0(missing_links, ' ', collapse = '')}"
  )

  # add default priors if missing
  missing_priors <- setdiff(names(model$parameters), names(model$default_priors))
  missing_priors <- setdiff(missing_priors, names(model$fixed_parameters))
  warnif(
    length(missing_priors) > 0 && getOption("bmm.default_priors"),
    "Default priors for each parameter will be specified internally based on the provided link function.
    Please check if the used priors are reasonable for your application"
  )
  additional_priors <- lapply(missing_priors, function(m) {
    if (model$other_vars$choice_rule == "simple") {
      switch(model$links[[m]],
             log = list(main = "normal(1, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             softplus = list(main = "normal(2, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             identity = list(main = "normal(10, 4)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             logit = list(main = "logistic(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             stop2("Invalid link function provided! Please use one of the following link functions: identity, log, softplus, logit")
      )
    } else {
      switch(model$links[[m]],
             log = list(main = "normal(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             softplus = list(main = "normal(1, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             identity = if (model$other_vars$choice_rule == "gaussian") {
               list(main = "normal(2, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)")
             } else {
               list(main = "normal(3, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)")
             },
             logit = list(main = "logistic(0, 1)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
             stop2("Invalid link function provided! Please use one of the following link functions: identity, log, softplus, logit")
      )
    }
  })
  model$default_priors <- c(model$default_priors, setNames(additional_priors, missing_priors))

  NextMethod("check_model")
}

# Candidates for the parameters a custom m3 adds: symbols in the activation
# and non-linear formulas that the model does not define already. brms fits
# every activation as non-linear, also one that is_nl() calls linear because
# no other formula parameter appears in it
m3_activation_symbols <- function(model, formula) {
  symbols <- union(
    rhs_vars(formula[is_nl(formula)]),
    rhs_vars(formula[intersect(model$resp_vars$resp_cats, names(formula))])
  )
  setdiff(symbols, c(names(formula[is_nl(formula)]), names(model$parameters)))
}

############################################################################# !
# CHECK_data S3 methods                                                  ####
############################################################################# !

# Counts or column names named after the response categories are labels: they
# are matched to the categories by name, and counts are stored under the same
# internal column names as unnamed counts. Used as column names they multiplied each
# category's activation by itself (#449). Fits from before the fix still carry
# those names in their stored model, so every reader goes through this helper
# rather than the constructor renaming them once
m3_num_options <- function(model) {
  num_options <- model$other_vars$num_options
  resp_cats <- model$resp_vars$resp_cats
  if (!setequal(names(num_options), resp_cats)) {
    return(num_options)
  }
  num_options <- num_options[resp_cats]
  if (is.numeric(num_options)) names(num_options) <- paste0("n_opt_", resp_cats)
  num_options
}

# Y is left out: as a matrix column it breaks the Stan code as a predictor
#' @exportS3Method
built_data_columns.m3 <- function(model) {
  num_options <- m3_num_options(model)
  c(
    if (is.numeric(num_options)) names(num_options),
    "nTrials", paste0("Idx_", model$resp_vars$resp_cats),
    NextMethod("built_data_columns")
  )
}

#' @export
check_data.m3 <- function(model, data, formula) {
  resp_name <- model$resp_vars$resp_cats
  n_opt_vect <- m3_num_options(model)
  col_names <- colnames(data)

  missing_variables <- setdiff(resp_name, col_names)
  stopif(length(missing_variables), "The response variable(s) {paste0(missing_variables, collapse = ', ')} missing in the data")

  # Y and nTrials may name a category, whose column is consumed first; brms
  # refuses `_` in category names, so Idx_<category> cannot be one
  reserved_cols <- intersect(
    c(setdiff(c("Y", "nTrials"), resp_name), paste0("Idx_", resp_name)),
    col_names
  )
  stopif(
    length(reserved_cols) > 0,
    "The data column(s) {collapse_comma(reserved_cols)} would be overwritten by \\
    the response matrix, trial counts and option indicators that bmm builds. \\
    Please rename them."
  )

  # Transfer all of the response variables to a matrix and name it 'Y'
  resp_matrix <- as.matrix(data[resp_name])
  missing_counts <- is.na(resp_matrix)
  resp_matrix[missing_counts] <- 0
  data <- data[!col_names %in% resp_name]
  data$nTrials <- rowSums(resp_matrix)
  data$Y <- resp_matrix

  if (is.character(n_opt_vect)) {
    # n_opt_vect is the *name* of the column in the data
    missing_options <- setdiff(n_opt_vect, col_names)
    stopif(length(missing_options), "The variable(s) {paste0(missing_options, collapse = ', ')} missing in the data")
    opt_vars <- n_opt_vect
    na_counts <- colSums(is.na(data[opt_vars]))
    na_counts <- na_counts[na_counts > 0]
    stopif(
      length(na_counts) > 0,
      "The option count column(s) contain missing values: \\
      {paste0(names(na_counts), ' (', na_counts, ' NA)', collapse = ', ')}. \\
      Give the number of response options for every row, and 0 where the category had none."
    )
  } else if (is.numeric(n_opt_vect)) {
    # n_opt_vect is the *number* of options for each response variable
    opt_vars <- names(n_opt_vect)
    # the counts become data columns under these names, and the activation
    # formulas refer to them by name, so any name already in use is taken to
    # mean the existing column or parameter instead of the count
    taken <- c(
      col_names, "Y", "nTrials", paste0("Idx_", resp_name),
      names(formula), names(model$parameters), names(model$fixed_parameters)
    )
    clashes <- intersect(opt_vars, taken)
    stopif(
      length(clashes) > 0,
      "The column name(s) {collapse_comma(clashes)} that `num_options` would be stored under are \\
      already taken by a data column, a model parameter, or a column bmm creates (`Y`, `nTrials`, \\
      `Idx_<category>`). Pass the numbers unnamed, name them after the response categories, or \\
      choose names that are not taken."
    )
    data[opt_vars] <- rep(n_opt_vect, each = nrow(data))
  } else {
    stop2("The number of options should be a string or a numeric vector.")
  }

  stopif(
    any(colSums(data[opt_vars]) == 0),
    "At least one of the specified number of candidates in the response categories is zero for all oberservations.
    Please remove this category from the model, as it is not identified."
  )

  # create index variables for any number of Option being zero in one row
  n_opt_idx_vars <- paste0("Idx_", resp_name)
  data[n_opt_idx_vars] <- as.integer(data[opt_vars] > 0)
  data[opt_vars][data[opt_vars] == 0] <- 0.0001

  # NA is how a category without options is usually recorded, and there it is
  # the true count; only where the category had options is a count lost
  n_missing <- sum(missing_counts & as.matrix(data[n_opt_idx_vars]) == 1)
  warnif(
    n_missing > 0,
    "The response category columns contain {n_missing} missing value(s) in rows \\
    where the category has response options. They are counted as 0 responses."
  )

  NextMethod("check_data")
}

############################################################################# !
# CHECK_Formula S3 methods                                               ####
############################################################################# !

#' @export
check_formula.m3 <- function(model, data, formula) {
  if (model$version != "custom") {
    formula <- construct_m3_act_funs(model, warnings = FALSE) + formula
  }

  formula <- apply_links(formula, model$links)
  formula <- assign_nl_attr(formula)

  NextMethod("check_formula")
}

#' @export
check_formula.m3_custom <- function(model, data, formula) {
  resp_cats <- model$resp_vars$resp_cats
  # test if activation functions for all categories are provided
  missing_act_funs <- !resp_cats %in% names(formula)
  stopif(
    any(missing_act_funs),
    "You did not provide activation functions for all response categories.
    Please provide activation functions for the following response categories in your bmmformula:
    {resp_cats[missing_act_funs]}"
  )

  # test if all activation functions contain background noise "b"
  act_funs <- formula[resp_cats]
  form_miss_b <- vapply(act_funs, function(f) !("b" %in% rhs_vars(f)), logical(1))
  stopif(
    any(form_miss_b),
    "Some of your activation functions do not contain the background noise parameter \"b\".
    The following activation functions need a background noise parameter:
    {resp_cats[form_miss_b]}"
  )

  NextMethod("check_formula")
}

############################################################################# !
# Convert bmmformula to brmsformla methods                               ####
############################################################################# !
#' @export
bmf2bf.m3 <- function(model, formula) {
  num_options <- m3_num_options(model)
  options_vars <- if (is.character(num_options)) num_options else names(num_options)
  resp_cats <- model$resp_vars$resp_cats
  n_opt_idx_vars <- paste0("Idx_", resp_cats)
  names(n_opt_idx_vars) <- resp_cats
  names(options_vars) <- resp_cats

  # set the base brmsformula based
  cat <- resp_cats[1]
  brms_formula <- brms::bf(glue(
    "Y | trials(nTrials) ~
    {n_opt_idx_vars[cat]} *", glue_choice_rule_functions(model$other_vars$choice_rule, cat, options_vars, resp_cats),
    "+ (1 - {n_opt_idx_vars[cat]}) * (-100)"
  ), nl = TRUE)

  # for each dependent parameter, check if it is used as a non-linear predictor of
  # another parameter and add the corresponding brms function
  for (cat in resp_cats[-1]) {
    brms_formula <- brms_formula + glue_nlf(
      "mu{cat} ~
      {n_opt_idx_vars[cat]} *", glue_choice_rule_functions(model$other_vars$choice_rule, cat, options_vars, resp_cats),
      "+ (1 - {n_opt_idx_vars[cat]}) * (-100)"
    )
  }

  brms_formula
}

#' @title glue the activation functions for the different choice rules
#'
#' @param choice_rule The choice rule that should be used for the M3: "softmax", "simple" or "gaussian"
#' @param cat The name of the response category for which the activation function should be generated
#' @param options_vars The variable names that contain the number of candidates in each response category
#' @param resp_cats The names of all response categories, in model order
#' @noRd
glue_choice_rule_functions <- function(choice_rule, cat, options_vars, resp_cats) {
  switch(
    choice_rule,
    simple = glue("log({cat} * {options_vars[cat]})"),
    softmax = glue("({cat} + log({options_vars[cat]}))"),
    # the softmax of log-probabilities that sum to one returns them unchanged, so
    # the Gaussian rule keeps the multinomial family and puts log P(cat) here
    gaussian = glue(
      "m3_gauss_logp({match(cat, resp_cats)}, {paste(resp_cats, collapse = ', ')}, \\
      {paste(options_vars[resp_cats], collapse = ', ')})"
    )
  )
}

# A non-linear formula can name only data columns and parameters, so the
# quadrature table and the number of categories are written into a generated
# wrapper with the same name and arguments as the R companion m3_gauss_logp()
m3_gaussian_stanvars <- function(model) {
  K <- length(model$resp_vars$resp_cats)
  gh <- .m3_gauss_rule()
  literal <- function(x) paste0("[", paste(formatC(x, digits = 17, format = "e"), collapse = ", "), "]'")
  wrapper <- glue(
    "real m3_gauss_logp(int k, {paste0('real A', 1:K, collapse = ', ')}, \\
    {paste0('real n', 1:K, collapse = ', ')}) {{
      return m3_gauss_logp_vec(k, [{paste0('A', 1:K, collapse = ', ')}]', \\
    [{paste0('n', 1:K, collapse = ', ')}]',
        {literal(gh$nodes)},
        {literal(gh$log_w)},
        {literal(gh$log_Phi)});
    }}"
  )
  sc_path <- system.file("stan_chunks", package = "bmm")
  brms::stanvar(
    scode = paste(read_lines2(paste0(sc_path, "/m3_gaussian_funs.stan")), wrapper, sep = "\n"),
    block = "functions"
  )
}

#' @title Category log-probability under the Gaussian choice rule of `m3()`
#' @description R companion to the Stan function `m3_gauss_logp` that
#'   `m3(choice_rule = "gaussian")` places in each category's activation
#'   formula. `brms` evaluates the non-linear formula in R for `log_lik()`,
#'   `posterior_predict()` and `posterior_epred()`, looking the function up on
#'   the search path; it is exported for that reason and is not meant to be
#'   called directly.
#' @param k Integer index of the response category.
#' @param ... The K category activations followed by the K option counts, as
#'   numbers or draws-by-observation matrices (as supplied by brms).
#' @return The log probability of category `k`, with the shape of the first
#'   activation.
#' @keywords internal
#' @export
m3_gauss_logp <- function(k, ...) {
  args <- list(...)
  K <- length(args) %/% 2
  len <- max(lengths(args))
  out <- .m3_gauss_logp_r(
    k,
    lapply(args[seq_len(K)], function(x) rep_len(as.vector(x), len)),
    lapply(args[K + seq_len(K)], function(x) rep_len(as.vector(x), len))
  )
  dim(out) <- dim(args[[1]])
  out
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !
# Each model should have a corresponding configure_model.* function. See
# ?configure_model for more information.

#' @export
configure_model.m3 <- function(model, data, formula) {
  # construct brms formula from the bmm formula
  formula <- bmf2bf(model, formula)

  # construct the family
  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- model$resp_vars$resp_cats
  formula$family$dpars <- paste0("mu", model$resp_vars$resp_cats)

  if (model$other_vars$choice_rule == "gaussian") {
    return(nlist(formula, data, stanvars = m3_gaussian_stanvars(model)))
  }
  nlist(formula, data)
}

#' @export
create_initfun.m3 <- function(model, data, formula, prior = NULL, ...) {
  # the "simple" choice rule with an identity link samples stably only from zero
  if (model$other_vars$choice_rule == "simple" && any(model$links == "identity")) {
    return(0)
  }
  NextMethod()
}


#' @title Get Activation Functions for different M3 versions
#'
#' @description
#' This function generates the activation functions for different versions of the Memory
#' Measurement Model (m3) implemented in the `bmm` package. If no `bmmodel` object is
#' passed then it will print the available model versions.
#'
#' @param model A bmmodel object that specifies the M3 model for which the
#'  activation functions should be generated. If no model is passed the available
#'  M3 versions will be printed to the console.
#' @param warnings Logical flag to indicate if information about the generated model formulas
#'  should be printed when the function is called.
#'
#' @return A bmmformula object with the activation functions for the m3 version specified in
#'  the model object. The activation functions use the names of the response categories
#'  specified in the model object.
#'
#' @examples
#' model <- m3(
#'  resp_cats = c("correct","other", "npl"),
#'  num_options = c(1, 4, 5),
#'  version = "ss"
#' )
#'
#' construct_m3_act_funs(model, warnings = FALSE)
#' @keywords transform
#' @export
construct_m3_act_funs <- function(model = NULL, warnings = TRUE) {
  if (is.null(model)) {
    message2(
      'Available m3 versions with pre-defined activation functions are:
          - "ss" for simple span tasks: 3 response categories (correct, other, npl)
          - "cs" for complex span tasks. 5 response categories (correct, dist_context, other, dist_other, npl)'
    )
    return(invisible())
  }

  stopif(
    !inherits(model, "m3") || !model$version %in% c("ss", "cs"),
    'Activation functions can only be generated for "m3" models "ss" and "cs"'
  )

  resp_cats <- model$resp_vars$resp_cats
  if (model$version == "ss") {
    warnif(
      warnings,
      '\nThe "ss" version of the m3 requires that response categories are ordered as follows:
      1) correct: correct responses
      2) other: other list responses
      3) npl: not presented lures'
    )

    act_funs <- bmf(
      formula(glue("{resp_cats[1]} ~ b + a + c")),
      formula(glue("{resp_cats[2]} ~ b + a")),
      formula(glue("{resp_cats[3]} ~ b"))
    )
  } else if (model$version == "cs") {
    warnif(
      warnings,
      "\nThe \"cs\" version of the m3 requires that response categories are ordered as follows:
      1) correct: correct responses
      2) dist_context: distractor responses close in context to the correct item
      3) other: other list responses
      4) dist_other: all distractor responses not close in context to the correct item
      5) npl: not presented lures"
    )

    act_funs <- bmf(
      formula(glue("{resp_cats[1]} ~ b + a + c")),
      formula(glue("{resp_cats[2]} ~ b + f * a + f * c")),
      formula(glue("{resp_cats[3]} ~ b + a")),
      formula(glue("{resp_cats[4]} ~ b + f * a")),
      formula(glue("{resp_cats[5]} ~ b"))
    )
  }

  act_funs
}


