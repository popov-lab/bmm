############################################################################# !
# MODELS                                                                 ####
############################################################################# !
# see file 'R/bmm_model_mixture3p.R' for an example

#' @title Create a tree for a Multinomial Processing Tree (MPT) model
#'
#' @description Specifies a single tree of a multinomial processing tree model
#'   by naming the tree and providing one branch-probability expression per
#'   response category. Trees created with `mpt_tree()` are combined into a
#'   full model specification with [mpt()].
#'
#' @param name Character. Label for the tree. For models with multiple trees,
#'   the label must match the values of the condition variable in the data that
#'   identify observations from this tree.
#' @param branches A named list of character strings. Each element gives the
#'   branch probability expression for one response category, and the element
#'   names are the response categories. The expressions can use latent
#'   parameters (e.g., `"D + (1 - D) * g"`), numeric constants, and declared
#'   covariates (see [mpt()]).
#'
#' @details Numeric fractions such as `1/4` are folded into decimal literals
#'   (`0.25`) when the tree is created. Stan compiles a bare integer fraction
#'   as integer division (`1/4 == 0`), which would silently corrupt the
#'   likelihood.
#'
#' @return An object of class `mpt_tree`
#'
#' @keywords transform
#'
#' @examples
#' tree_old <- mpt_tree(
#'   name = "old",
#'   branches = list(
#'     old = "D + (1 - D) * g",
#'     new = "(1 - D) * (1 - g)"
#'   )
#' )
#' tree_old
#' @export
mpt_tree <- function(name, branches) {
  stop_missing_args()
  stopif(
    !is.character(name) || length(name) != 1L || !nzchar(name),
    "The tree name must be a single non-empty character string."
  )
  stopif(
    !is.list(branches) || length(branches) == 0L,
    "The branches must be a named list with one branch probability expression \\
    per response category."
  )
  resp_cats <- names(branches)
  stopif(
    is.null(resp_cats) || any(!nzchar(resp_cats)),
    "Each branch probability expression must be named after its response category."
  )
  stopif(
    anyDuplicated(resp_cats) > 0,
    "Response category names must be unique within a tree. Branch lines that \\
    terminate in the same response category must be summed into one expression."
  )
  branches[] <- lapply(branches, .mpt_canonical_expr, tree_name = name)
  structure(nlist(name, branches), class = "mpt_tree")
}

#' @export
print.mpt_tree <- function(x, ...) {
  branch_lines <- glue("  P({names(x$branches)}) = {unlist(x$branches)}")
  cat(glue("MPT tree '{x$name}':"), branch_lines, sep = "\n")
  invisible(x)
}

.model_mpt <- function(trees = NULL, condition = NULL, covariates = NULL,
                       simplex = NULL, links = "logit", default_priors = NULL,
                       call = NULL, ...) {
  trees <- .mpt_as_tree_list(trees)
  if (length(trees)) names(trees) <- vapply(trees, `[[`, character(1), "name")
  covariates <- covariates %||% character(0)
  simplex <- .mpt_as_simplex_list(simplex)
  resp_cats <- if (length(trees)) names(trees[[1]]$branches) else character(0)
  parameters <- setdiff(
    unique(unlist(lapply(trees, .mpt_tree_vars))), covariates
  )
  simplex_pars <- unlist(simplex)
  raw_pars <- unlist(lapply(simplex, function(grp) paste0(grp[-length(grp)], "raw")))
  standard_pars <- setdiff(parameters, simplex_pars)

  # matched priors: effects use the same family as the intercept at half the
  # scale, so both links imply comparable regularization on their latent scale
  latent_prior <- switch(links,
    logit = list(main = "logistic(0, 1)", effects = "logistic(0, 0.5)"),
    probit = list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
  )

  parameter_info <- c(
    named_value_list(
      standard_pars,
      glue("Latent probability parameter from the tree branch expressions. \\
           Fitted on the {links} scale.")
    ),
    named_value_list(
      simplex_pars,
      paste0(
        "Probability parameter constrained within a simplex group via ",
        "stick-breaking. The last parameter of each group is derived as 1 ",
        "minus the sum of the other group members."
      )
    ),
    named_value_list(
      raw_pars,
      glue("Unconstrained stick-breaking component of a simplex parameter. \\
           Fitted on the {links} scale.")
    )
  )

  link_info <- c(
    named_value_list(standard_pars, links),
    named_value_list(simplex_pars, "identity"),
    named_value_list(raw_pars, "identity")
  )

  prior_info <- c(
    named_value_list(standard_pars, latent_prior),
    named_value_list(raw_pars, latent_prior)
  )

  out <- structure(
    list(
      resp_vars = nlist(resp_cats),
      other_vars = nlist(condition, covariates, link = links),
      trees = trees,
      simplex = simplex,
      domain = "Categorical decision making, memory, and reasoning",
      task = "Any task with categorical responses generated by a discrete processing-tree structure",
      name = "Multinomial Processing Tree (MPT) models",
      citation = glue(
        "Batchelder, W. H., & Riefer, D. M. (1999). Theoretical and empirical \\
        review of multinomial process tree modeling. Psychonomic Bulletin & \\
        Review, 6(1), 57-86. https://doi.org/10.3758/BF03210812"
      ),
      version = "",
      requirements = paste0(
        "- One tree per experimental condition created with mpt_tree(); ",
        "all trees share the same response categories\n",
        "  - The data contain one column with aggregated response counts per ",
        "response category, named after the branch names\n",
        "  - For multi-tree models, a condition column whose values match the tree names\n",
        "  - Data columns used inside branch expressions must be declared via the covariates argument\n"
      ),
      parameters = parameter_info,
      fixed_parameters = list(),
      links = link_info,
      default_priors = prior_info,
      void_mu = FALSE
    ),
    class = c("bmmodel", "mpt"),
    call = call
  )
  out$default_priors[names(default_priors)] <- default_priors
  out
}

# user facing alias
# information in the title and details sections will be filled in
# automatically based on the information in the .model_mpt()

#' @title `r .model_mpt()$name`
#' @name mpt
#'
#' @description
#' Multinomial Processing Tree (MPT) models measure the probabilities of
#' discrete latent cognitive states from categorical responses. A model is
#' specified as a set of trees (one per experimental condition) created with
#' [mpt_tree()]. Branch expressions that terminate in the same response
#' category are summed within each tree, and parameters shared across trees
#' are equated by giving them the same name. All latent probability
#' parameters are estimated on an unconstrained latent scale (logit or
#' probit); predictor formulas supplied via [bmmformula()] apply on that
#' latent scale.
#'
#' @param trees A single `mpt_tree` object or a list of `mpt_tree` objects.
#'   All trees must share the same set of response categories.
#' @param condition Character. Name of the data column identifying which tree
#'   an observation belongs to; its values must match the tree names. Can be
#'   omitted for single-tree models.
#' @param covariates Character vector. Names of data columns that appear in
#'   branch expressions but are not latent parameters, for example
#'   design-fixed guessing rates. Covariates pass into the model formulas
#'   unchanged and have no naming restrictions.
#' @param simplex A character vector, or a list of character vectors, naming
#'   groups of parameters that are jointly constrained to sum to 1. Each group
#'   is reparameterized via stick-breaking: the last parameter of each group
#'   is derived as 1 minus the sum of the other members, and each other member
#'   `p` gets an unconstrained component `praw` that predictor formulas for
#'   `p` are applied to.
#' @param links Character. The link function for all latent probability
#'   parameters: `"logit"` (default) or `"probit"`.
#' @param ... used internally for testing, ignore it
#'
#' @details `r model_info(.model_mpt(), components = c('domain', 'task', 'name', 'citation'))`
#'
#'   Parameters that receive a non-linear predictor formula (a formula whose
#'   right-hand side references other formula parameters, e.g.
#'   `D ~ Dmax * (1 - exp(-rate * ptime))`) are not transformed by the link
#'   function: the user-supplied expression must keep the parameter within
#'   (0, 1). Sub-parameters of such formulas (e.g., `Dmax` and `rate`) are
#'   estimated on the identity scale with `normal(0, 1)` default priors.
#'
#'   Parameter and response category names must start with a letter and may
#'   contain only letters and digits, because brms does not allow underscores
#'   or dots in non-linear parameter names.
#'
#' @return An object of class `bmmodel`
#'
#' @keywords bmmodel
#'
#' @examplesIf isTRUE(Sys.getenv("BMM_EXAMPLES"))
#' # two-high-threshold (2HTM) model of recognition memory
#' tree_old <- mpt_tree("old", list(
#'   old = "D + (1 - D) * g",
#'   new = "(1 - D) * (1 - g)"
#' ))
#' tree_new <- mpt_tree("new", list(
#'   old = "(1 - D) * g",
#'   new = "D + (1 - D) * (1 - g)"
#' ))
#'
#' model <- mpt(
#'   trees = list(tree_old, tree_new),
#'   condition = "item_type"
#' )
#'
#' # simulate data for 20 participants with D = 0.7, g = 0.5
#' data <- data.frame(
#'   id = rep(1:20, each = 2),
#'   item_type = rep(c("old", "new"), 20)
#' )
#' counts <- t(sapply(data$item_type, function(cond) {
#'   rmpt(
#'     n = 1, size = 50, pars = c(D = 0.7, g = 0.5),
#'     mpt_model = model, tree = cond, unpack = TRUE
#'   )
#' }))
#' data <- cbind(data, counts)
#'
#' # predict both parameters by a fixed and random intercept
#' formula <- bmf(
#'   D ~ 1 + (1 | id),
#'   g ~ 1 + (1 | id)
#' )
#'
#' fit <- bmm(
#'   formula = formula,
#'   data = data,
#'   model = model,
#'   cores = 4
#' )
#'
#' summary(fit)
#'
#' @export
mpt <- function(trees, condition = NULL, covariates = NULL, simplex = NULL,
                links = "logit", ...) {
  call <- match.call()
  stop_missing_args()
  links <- match.arg(links, c("logit", "probit"))
  trees <- .mpt_as_tree_list(trees)
  covariates <- covariates %||% character(0)
  simplex <- .mpt_as_simplex_list(simplex)

  stopif(
    length(trees) == 0L || !all(vapply(trees, inherits, logical(1), "mpt_tree")),
    "The trees argument must be an mpt_tree object or a list of mpt_tree \\
    objects created with mpt_tree()."
  )
  tree_names <- vapply(trees, `[[`, character(1), "name")
  stopif(
    anyDuplicated(tree_names) > 0,
    "Tree names must be unique. Duplicated: {collapse_comma(unique(tree_names[duplicated(tree_names)]))}"
  )
  bad_tree_names <- tree_names[!grepl("^[A-Za-z][A-Za-z0-9_]*$", tree_names)]
  stopif(
    length(bad_tree_names) > 0,
    "Tree names must start with a letter and contain only letters, digits, or \\
    underscores. Please rename: {collapse_comma(bad_tree_names)}"
  )

  resp_cats <- names(trees[[1]]$branches)
  cats_match <- vapply(
    trees, function(tree) setequal(names(tree$branches), resp_cats), logical(1)
  )
  stopif(
    !all(cats_match),
    "All trees must have the same response categories.
    Tree '{trees[[1]]$name}' has: {collapse_comma(resp_cats)}
    Tree '{tree_names[!cats_match][1]}' has: {collapse_comma(names(trees[!cats_match][[1]]$branches))}"
  )
  trees <- lapply(trees, function(tree) {
    tree$branches <- tree$branches[resp_cats]
    tree
  })
  names(trees) <- tree_names

  stopif(
    length(trees) > 1L && is.null(condition),
    "Models with multiple trees require the condition argument: the name of \\
    the data column whose values identify the tree each observation belongs to."
  )
  stopif(
    !is.null(condition) && (!is.character(condition) || length(condition) != 1L),
    "The condition argument must be a single character string naming a data column."
  )
  stopif(
    length(covariates) > 0L && !is.character(covariates),
    "The covariates argument must be a character vector of data column names."
  )

  parameters <- setdiff(unique(unlist(lapply(trees, .mpt_tree_vars))), covariates)
  stopif(
    length(parameters) == 0L,
    "The tree branch expressions contain no latent parameters. Symbols listed \\
    in the covariates argument are not treated as parameters."
  )

  .mpt_check_names(parameters, "parameter")
  .mpt_check_names(resp_cats, "response category")

  par_cat_overlap <- intersect(parameters, resp_cats)
  stopif(
    length(par_cat_overlap) > 0,
    "Names cannot be used for both a parameter and a response category: \\
    {collapse_comma(par_cat_overlap)}"
  )
  cov_collisions <- intersect(covariates, c(resp_cats, parameters))
  stopif(
    length(cov_collisions) > 0,
    "Covariate names collide with parameter or response category names: \\
    {collapse_comma(cov_collisions)}"
  )
  reserved <- intersect(c(parameters, resp_cats, covariates), c("Y", "nTrials"))
  stopif(
    length(reserved) > 0,
    "The names 'Y' and 'nTrials' are reserved for the response matrix and \\
    trial counts. Please rename: {collapse_comma(reserved)}"
  )

  simplex_pars <- unlist(simplex)
  stopif(
    anyDuplicated(simplex_pars) > 0,
    "Parameters cannot appear in more than one simplex group: \\
    {collapse_comma(unique(simplex_pars[duplicated(simplex_pars)]))}"
  )
  unknown_simplex <- setdiff(simplex_pars, parameters)
  stopif(
    length(unknown_simplex) > 0,
    "Simplex groups can only contain latent parameters from the tree \\
    expressions. Unknown: {collapse_comma(unknown_simplex)}"
  )
  short_groups <- lengths(simplex) < 2L
  stopif(
    any(short_groups),
    "Each simplex group needs at least two parameters."
  )
  raw_pars <- unlist(lapply(simplex, function(grp) paste0(grp[-length(grp)], "raw")))
  raw_collisions <- intersect(raw_pars, c(parameters, resp_cats, covariates))
  stopif(
    length(raw_collisions) > 0,
    "The stick-breaking components of simplex parameters are named by \\
    appending 'raw' to the parameter name, but these names are already in \\
    use: {collapse_comma(raw_collisions)}"
  )

  .mpt_validate_tree_sums(trees, parameters, covariates, simplex)

  .model_mpt(
    trees = trees, condition = condition, covariates = covariates,
    simplex = simplex, links = links, call = call, ...
  )
}

############################################################################# !
# CHECK_MODEL S3 methods                                                 ####
############################################################################# !

#' @export
check_model.mpt <- function(model, data = NULL, formula = NULL) {
  if (!is.null(formula)) {
    resp_cats <- model$resp_vars$resp_cats
    user_cat_formulas <- intersect(resp_cats, names(formula))
    stopif(
      length(user_cat_formulas) > 0,
      "The response category probabilities are fully determined by the tree \\
      branch expressions and cannot be predicted directly. Please remove the \\
      formula(s) for: {collapse_comma(user_cat_formulas)}"
    )
    model <- .mpt_apply_nl_bypass(model, formula, data)
    model <- .mpt_fixed_to_latent_scale(model)
  }
  NextMethod("check_model")
}

# parameters with a user-supplied non-linear formula own their (0,1)
# constraint, so the automatic link transformation must be switched off for
# them; sub-parameters referenced by such formulas become model parameters
.mpt_apply_nl_bypass <- function(model, formula, data) {
  nl_pars <- names(formula)[is_nl(formula)]
  bypass_pars <- intersect(nl_pars, names(model$parameters))
  simplex_bypass <- intersect(bypass_pars, unlist(model$simplex))
  stopif(
    length(simplex_bypass) > 0,
    "Non-linear predictor formulas are not supported for simplex parameters: \\
    {collapse_comma(simplex_bypass)}"
  )
  for (par in bypass_pars) {
    model$links[[par]] <- "identity"
    model$default_priors[[par]] <- NULL
  }

  sub_pars <- rhs_vars(formula[is_nl(formula)])
  sub_pars <- setdiff(sub_pars, nl_pars)
  sub_pars <- setdiff(sub_pars, names(model$parameters))
  sub_pars <- setdiff(sub_pars, colnames(data))
  sub_pars <- setdiff(
    sub_pars,
    c(model$other_vars$covariates, model$other_vars$condition)
  )
  if (length(sub_pars) > 0) {
    .mpt_check_names(sub_pars, "parameter")
    message2(
      "The parameter(s) {collapse_comma(sub_pars)} from your non-linear \\
      formulas are estimated on the identity scale with normal(0, 1) default \\
      priors. Apply any required transformation inside your formula and \\
      adjust the priors to the scale of your predictors."
    )
    model$parameters[sub_pars] <- named_value_list(
      sub_pars, "User-defined sub-parameter of a non-linear parameter formula."
    )
    model$links[sub_pars] <- named_value_list(sub_pars, "identity")
    model$default_priors[sub_pars] <- named_value_list(
      sub_pars, list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
    )
  }
  model
}

# fixed values from the bmmformula (e.g. g = 0.5) are probabilities, but the
# constant() prior applies to the latent intercept, so the values must be
# mapped through the link function
.mpt_fixed_to_latent_scale <- function(model) {
  if (isTRUE(attr(model, "mpt_fixed_on_latent_scale"))) {
    return(model)
  }
  fixed_pars <- intersect(names(model$fixed_parameters), names(model$parameters))
  fixed_simplex <- intersect(fixed_pars, unlist(model$simplex))
  stopif(
    length(fixed_simplex) > 0,
    "Fixing simplex parameters to constants is not supported: \\
    {collapse_comma(fixed_simplex)}"
  )
  for (par in fixed_pars) {
    link <- model$links[[par]]
    if (link == "identity") next
    value <- model$fixed_parameters[[par]]
    stopif(
      !is.numeric(value) || value <= 0 || value >= 1,
      "The fixed value for parameter '{par}' must be a probability strictly \\
      between 0 and 1. Provided: {value}"
    )
    latent <- switch(link, logit = stats::qlogis(value), probit = stats::qnorm(value))
    message2(
      "The fixed value {value} for parameter '{par}' is applied on the \\
      {link} scale as {round(latent, 4)}."
    )
    model$fixed_parameters[[par]] <- latent
  }
  attr(model, "mpt_fixed_on_latent_scale") <- TRUE
  model
}

############################################################################# !
# CHECK_data S3 methods                                                  ####
############################################################################# !

#' @export
check_data.mpt <- function(model, data, formula) {
  resp_cats <- model$resp_vars$resp_cats
  col_names <- colnames(data)

  missing_cats <- setdiff(resp_cats, col_names)
  stopif(
    length(missing_cats) > 0,
    "The data must contain one column of response counts per response category,
    named after the branch names of the trees.
    Expected columns: {collapse_comma(resp_cats)}
    Missing columns: {collapse_comma(missing_cats)}"
  )

  resp_matrix <- as.matrix(data[resp_cats])
  stopif(
    !is.numeric(resp_matrix) || any(resp_matrix < 0, na.rm = TRUE),
    "The response category columns must contain non-negative response counts."
  )
  resp_matrix[is.na(resp_matrix)] <- 0
  data <- data[!col_names %in% resp_cats]
  data$nTrials <- rowSums(resp_matrix)
  data$Y <- resp_matrix

  condition <- model$other_vars$condition
  if (!is.null(condition)) {
    stopif(
      !condition %in% colnames(data),
      "The condition variable '{condition}' is not present in the data."
    )
    cond_values <- as.character(data[[condition]])
    tree_names <- names(model$trees)
    unmatched <- setdiff(unique(cond_values), tree_names)
    stopif(
      length(unmatched) > 0,
      "All values of the condition variable '{condition}' must match a tree name.
      Tree names: {collapse_comma(tree_names)}
      Unmatched values: {collapse_comma(unmatched)}"
    )
    unused_trees <- setdiff(tree_names, unique(cond_values))
    warnif(
      length(unused_trees) > 0,
      "The data contain no observations for tree(s): {collapse_comma(unused_trees)}"
    )
    idx_vars <- paste0("Idx_", tree_names)
    idx_collisions <- intersect(idx_vars, colnames(data))
    stopif(
      length(idx_collisions) > 0,
      "The data contain column(s) {collapse_comma(idx_collisions)}, which are \\
      reserved for the generated tree indicator variables. Please rename them."
    )
    for (i in seq_along(tree_names)) {
      data[[idx_vars[i]]] <- as.integer(cond_values == tree_names[i])
    }
  }

  covariates <- model$other_vars$covariates
  missing_covariates <- setdiff(covariates, colnames(data))
  stopif(
    length(missing_covariates) > 0,
    "The declared covariates {collapse_comma(missing_covariates)} are missing \\
    from the data."
  )
  non_numeric <- covariates[
    !vapply(covariates, function(v) is.numeric(data[[v]]), logical(1))
  ]
  stopif(
    length(non_numeric) > 0,
    "The declared covariates must be numeric data columns. Not numeric: \\
    {collapse_comma(non_numeric)}"
  )

  .mpt_validate_covariate_sums(model, data)

  NextMethod("check_data")
}

# the construction-time branch-sum validation uses synthetic covariate values;
# once the data are known, the sum-to-1 property is re-checked row by row with
# the observed covariate values to catch data-preparation errors
.mpt_validate_covariate_sums <- function(model, data, tolerance = 1e-6) {
  covariates <- model$other_vars$covariates
  if (length(covariates) == 0L) {
    return(invisible(NULL))
  }
  parameters <- setdiff(
    unique(unlist(lapply(model$trees, .mpt_tree_vars))), covariates
  )
  par_vals <- setNames(rep(0.5, length(parameters)), parameters)
  for (grp in model$simplex) {
    par_vals[grp] <- 1 / length(grp)
  }

  condition <- model$other_vars$condition
  for (tree in model$trees) {
    rows <- if (is.null(condition)) {
      seq_len(nrow(data))
    } else {
      which(data[[paste0("Idx_", tree$name)]] == 1L)
    }
    if (length(rows) == 0L) next
    env <- c(as.list(par_vals), as.list(data[rows, covariates, drop = FALSE]))
    total <- Reduce(`+`, lapply(tree$branches, function(branch) {
      eval(str2lang(branch), envir = env)
    }))
    total <- rep(total, length.out = length(rows))
    deviates <- is.na(total) | abs(total - 1) > tolerance
    warnif(
      any(deviates),
      "With the covariate values in the data, the branch probabilities of \\
      tree '{tree$name}' do not sum to 1 for {sum(deviates)} row(s) \\
      (first: row {rows[deviates][1]}, sum = {signif(total[deviates][1], 6)}).
      Please check the covariate column(s): {collapse_comma(covariates)}"
    )
  }
  invisible(NULL)
}

############################################################################# !
# CHECK_Formula S3 methods                                               ####
############################################################################# !

#' @export
check_formula.mpt <- function(model, data, formula) {
  formula <- .mpt_move_simplex_formulas(model, formula)

  generated <- .mpt_category_formulas(model)
  simplex_formulas <- .mpt_simplex_formulas(model)
  if (!is.null(simplex_formulas)) {
    generated <- generated + simplex_formulas
  }
  formula <- generated + formula

  formula <- apply_links(formula, model$links)
  formula <- assign_nl_attr(formula)
  # category probabilities must always be non-linear formulas, even if a
  # branch expression happens to contain no latent parameters
  for (resp_cat in model$resp_vars$resp_cats) {
    attr(formula[[resp_cat]], "nl") <- TRUE
  }

  NextMethod("check_formula")
}

# predictor formulas for the free parameters of a simplex group apply to
# their unconstrained stick-breaking components; the parameters themselves
# receive generated stick-breaking formulas
.mpt_move_simplex_formulas <- function(model, formula) {
  for (grp in model$simplex) {
    free_pars <- grp[-length(grp)]
    derived_par <- grp[length(grp)]
    for (par in free_pars) {
      raw_par <- paste0(par, "raw")
      par_rhs <- .mpt_rhs_chr(formula[[par]])
      raw_rhs <- .mpt_rhs_chr(formula[[raw_par]])
      stopif(
        par_rhs != "1" && raw_rhs != "1" && par_rhs != raw_rhs,
        "Conflicting predictor formulas for the simplex parameter '{par}' and \\
        its stick-breaking component '{raw_par}'. Specify predictors for \\
        '{par}' only."
      )
      if (par_rhs != "1") {
        formula[raw_par] <- list(stats::formula(glue("{raw_par} ~ {par_rhs}")))
      }
    }
    derived_rhs <- .mpt_rhs_chr(formula[[derived_par]])
    stopif(
      derived_rhs != "1",
      "The parameter '{derived_par}' is derived as 1 minus the sum of \\
      {collapse_comma(free_pars)} and cannot have its own predictors. Specify \\
      predictors for the other parameters of the simplex group instead."
    )
    formula <- formula[setdiff(names(formula), grp)]
  }
  formula
}

.mpt_rhs_chr <- function(pform) {
  if (!is_formula(pform)) {
    return(NA_character_)
  }
  paste(deparse(pform[[3]]), collapse = " ")
}

.mpt_category_formulas <- function(model) {
  trees <- model$trees
  use_indicators <- !is.null(model$other_vars$condition)
  category_formulas <- lapply(model$resp_vars$resp_cats, function(resp_cat) {
    branch_exprs <- vapply(trees, function(tree) tree$branches[[resp_cat]], character(1))
    rhs <- if (use_indicators) {
      paste(glue("Idx_{names(trees)} * ({branch_exprs})"), collapse = " + ")
    } else {
      glue("({branch_exprs})")
    }
    stats::formula(glue("{resp_cat} ~ {rhs}"))
  })
  do.call(bmf, category_formulas)
}

.mpt_simplex_formulas <- function(model) {
  if (length(model$simplex) == 0L) {
    return(NULL)
  }
  ilink <- switch(model$other_vars$link, logit = "inv_logit", probit = "Phi")
  simplex_formulas <- list()
  for (grp in model$simplex) {
    n_grp <- length(grp)
    free_pars <- grp[-n_grp]
    sticks <- glue("{ilink}({free_pars}raw)")
    for (k in seq_len(n_grp - 1L)) {
      rhs <- if (k == 1L) {
        sticks[1]
      } else {
        remaining <- paste(glue("(1 - {sticks[seq_len(k - 1L)]})"), collapse = " * ")
        glue("{remaining} * {sticks[k]}")
      }
      simplex_formulas <- c(
        simplex_formulas, list(stats::formula(glue("{grp[k]} ~ {rhs}")))
      )
    }
    simplex_formulas <- c(
      simplex_formulas,
      list(stats::formula(glue("{grp[n_grp]} ~ 1 - ({paste(free_pars, collapse = ' + ')})")))
    )
  }
  do.call(bmf, simplex_formulas)
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.mpt <- function(model, formula) {
  resp_cats <- model$resp_vars$resp_cats

  brms_formula <- brms::bf(
    glue("Y | trials(nTrials) ~ log({resp_cats[1]})"),
    nl = TRUE
  )
  for (resp_cat in resp_cats[-1]) {
    brms_formula <- brms_formula + glue_nlf("mu{resp_cat} ~ log({resp_cat})")
  }

  brms_formula
}

############################################################################# !
# CONFIGURE_MODEL S3 METHODS                                             ####
############################################################################# !

#' @export
configure_model.mpt <- function(model, data, formula) {
  formula <- bmf2bf(model, formula)

  formula$family <- brms::multinomial(refcat = NA)
  formula$family$cats <- model$resp_vars$resp_cats
  formula$family$dpars <- paste0("mu", model$resp_vars$resp_cats)

  nlist(formula, data)
}

############################################################################# !
# HELPERS                                                                ####
############################################################################# !

.mpt_as_tree_list <- function(trees) {
  if (is.null(trees)) {
    return(list())
  }
  if (inherits(trees, "mpt_tree")) {
    return(list(trees))
  }
  trees
}

.mpt_as_simplex_list <- function(simplex) {
  if (is.null(simplex)) {
    return(list())
  }
  if (!is.list(simplex)) {
    return(list(simplex))
  }
  simplex
}

.mpt_canonical_expr <- function(expr, tree_name) {
  stopif(
    !is.character(expr) || length(expr) != 1L || !nzchar(expr),
    "Branch probability expressions must be single character strings \\
    (tree '{tree_name}')."
  )
  parsed <- try(str2lang(expr), silent = TRUE)
  stopif(
    is_try_error(parsed),
    "Cannot parse the branch expression '{expr}' in tree '{tree_name}'."
  )
  paste(deparse(.mpt_fold_numeric_division(parsed)), collapse = " ")
}

# Stan compiles a bare numeric fraction like 1/4 as integer division (= 0),
# so constant divisions are folded into decimal literals before emission
.mpt_fold_numeric_division <- function(expr) {
  if (!is.call(expr)) {
    return(expr)
  }
  for (i in seq_along(expr)[-1]) {
    expr[[i]] <- .mpt_fold_numeric_division(expr[[i]])
  }
  if (identical(expr[[1]], quote(`(`)) && is.numeric(expr[[2]])) {
    return(expr[[2]])
  }
  if (identical(expr[[1]], quote(`/`)) && length(expr) == 3L &&
        is.numeric(expr[[2]]) && is.numeric(expr[[3]])) {
    return(expr[[2]] / expr[[3]])
  }
  expr
}

.mpt_tree_vars <- function(tree) {
  unique(unlist(lapply(tree$branches, function(branch) all.vars(str2lang(branch)))))
}

.mpt_check_names <- function(names, what) {
  bad_names <- names[!grepl("^[A-Za-z][A-Za-z0-9]*$", names)]
  stopif(
    length(bad_names) > 0,
    "MPT {what} names must start with a letter and contain only letters and \\
    digits, because brms does not allow underscores or dots in non-linear \\
    parameter names. Please rename (e.g., 'd_A' -> 'dA'): {collapse_comma(bad_names)}"
  )
}

# branches of each tree must sum to 1 for any parameter values; evaluating at
# several distinct test points catches swapped-complement errors that a single
# symmetric point (e.g. all 0.5) would miss
.mpt_validate_tree_sums <- function(trees, parameters, covariates, simplex,
                                    tolerance = 1e-6) {
  symbols <- c(parameters, covariates)
  if (length(symbols) == 0L) {
    return(invisible(NULL))
  }
  test_vals <- c(0.137, 0.421, 0.683, 0.852)
  for (tree in trees) {
    for (shift in seq_along(test_vals)) {
      vals <- setNames(
        test_vals[(seq_along(symbols) + shift - 2L) %% length(test_vals) + 1L],
        symbols
      )
      for (grp in simplex) {
        vals[grp] <- vals[grp] / sum(vals[grp])
      }
      total <- sum(vapply(
        tree$branches,
        function(branch) eval(str2lang(branch), envir = as.list(vals)),
        numeric(1)
      ))
      if (abs(total - 1) > tolerance) {
        warning2(
          "The branch probabilities of tree '{tree$name}' sum to \\
          {signif(total, 6)} instead of 1 when evaluated at numeric test \\
          values. Please check the branch expressions."
        )
        break
      }
    }
  }
  invisible(NULL)
}

named_value_list <- function(names, value) {
  if (length(names) == 0L) {
    return(list())
  }
  setNames(rep(list(value), length(names)), names)
}

############################################################################# !
# IMPORT WRAPPERS                                                        ####
############################################################################# !

#' @title Create an MPT model from an MPTinR-style model string
#'
#' @description Parses a model definition in the format used by MPTinR model
#'   files and returns the corresponding [mpt()] model object. Each line is
#'   one branch probability expression; blank lines separate trees. The
#'   response category of each line is given either by an inline comment
#'   (`expression # category`) or by the `categories` argument. Lines that
#'   terminate in the same response category are summed.
#'
#' @param text Character. The model definition, either as a single string or
#'   as a vector of lines.
#' @param tree_names Character vector with one name per tree block, in the
#'   order the blocks appear. For multi-tree models the names must match the
#'   values of the condition variable in the data.
#' @param categories Character vector with the response category of each line
#'   within a tree, used for lines without an inline `# category` comment.
#'   Must have as many entries as the tree with the most lines.
#' @param condition Character. Name of the data column identifying the tree
#'   an observation belongs to (see [mpt()]).
#' @param covariates,simplex,links passed to [mpt()].
#'
#' @return An object of class `bmmodel` (see [mpt()])
#'
#' @keywords transform
#'
#' @examples
#' model_2htm <- "
#' D + (1 - D) * g        # old
#' (1 - D) * (1 - g)      # new
#'
#' (1 - D) * g            # old
#' D + (1 - D) * (1 - g)  # new
#' "
#' model <- mpt_from_string(
#'   model_2htm,
#'   tree_names = c("old", "new"),
#'   condition = "item_type"
#' )
#' model
#' @export
mpt_from_string <- function(text, tree_names, categories = NULL,
                            condition = NULL, covariates = NULL,
                            simplex = NULL, links = "logit") {
  stop_missing_args()
  lines <- trimws(unlist(strsplit(text, "\n")))
  blocks <- split(lines[nzchar(lines)], cumsum(!nzchar(lines))[nzchar(lines)])
  stopif(
    length(blocks) != length(tree_names),
    "Found {length(blocks)} tree block(s) separated by blank lines, but \\
    {length(tree_names)} tree_names were supplied."
  )

  trees <- mapply(function(block, tree_name) {
    parts <- strsplit(block, "#", fixed = TRUE)
    exprs <- trimws(vapply(parts, `[[`, character(1), 1))
    line_cats <- vapply(parts, function(part) {
      if (length(part) > 1) trimws(paste(part[-1], collapse = "#")) else ""
    }, character(1))
    if (any(!nzchar(line_cats))) {
      stopif(
        is.null(categories) || length(categories) < length(block),
        "Tree '{tree_name}' has lines without an inline '# category' comment.
        Please provide the response category of each line via the categories \\
        argument."
      )
      line_cats[!nzchar(line_cats)] <- categories[which(!nzchar(line_cats))]
    }
    mpt_tree(tree_name, .mpt_sum_branch_lines(exprs, line_cats))
  }, blocks, tree_names, SIMPLIFY = FALSE)

  mpt(
    trees = trees, condition = condition, covariates = covariates,
    simplex = simplex, links = links
  )
}

#' @title Create an MPT model from an EQN model file
#'
#' @description Reads a model definition in the EQN format used by multiTree,
#'   MPTinR, and TreeBUGS, and returns the corresponding [mpt()] model object.
#'   Each line of an EQN file contains three whitespace-separated fields --
#'   tree, response category, branch probability expression -- and branch
#'   lines that terminate in the same response category are summed.
#'
#'   In bmm, every observation (row) belongs to one tree, and all trees share
#'   the same response categories. Classical EQN files often label the
#'   response categories per tree instead (e.g., `hit`/`miss` in the old-item
#'   tree and `fa`/`cr` in the new-item tree of a recognition model). Use the
#'   `categories` argument to map such tree-specific labels onto the shared
#'   response categories (e.g., the response options "yes" and "no").
#'
#'   Because brms does not allow underscores or dots in non-linear parameter
#'   names, parameter and response category names are sanitized: a leading
#'   `<tree>_` prefix is removed from category names, and all remaining
#'   underscores and dots are stripped. The complete renaming map is reported
#'   and stored in the `mpt_renaming` attribute of the returned model.
#'
#' @param file Character. Path to the EQN file. A first line that does not
#'   contain three fields (the line-count header of classical EQN files) is
#'   skipped.
#' @param restrictions A named list or vector mapping parameter names to
#'   numeric constants. The constants are substituted into the branch
#'   expressions as decimal literals before the parameters are identified,
#'   which mirrors the `restrictions` mechanism of MPTinR/TreeBUGS. Equality
#'   restrictions are expressed by using the same parameter name in the EQN
#'   file itself.
#' @param categories A named character vector mapping response category names
#'   used in the EQN file onto the shared response categories of the model
#'   (see Details). Categories not listed keep their (sanitized) name.
#' @param condition Character. Name of the data column identifying the tree
#'   an observation belongs to (see [mpt()]). Defaults to the tree names used
#'   in the EQN file.
#' @param covariates,simplex,links passed to [mpt()]. Covariate names are
#'   excluded from the renaming.
#'
#' @return An object of class `bmmodel` (see [mpt()]), with an `mpt_renaming`
#'   attribute holding the applied renaming map.
#'
#' @keywords transform
#'
#' @examples
#' eqn_file <- tempfile(fileext = ".eqn")
#' writeLines(c(
#'   "6",
#'   "old  hit   D_o",
#'   "old  hit   (1-D_o)*g",
#'   "old  miss  (1-D_o)*(1-g)",
#'   "new  fa    (1-D_n)*g",
#'   "new  cr    D_n",
#'   "new  cr    (1-D_n)*(1-g)"
#' ), eqn_file)
#' model <- mpt_from_eqn(
#'   eqn_file,
#'   categories = c(hit = "yes", fa = "yes", miss = "no", cr = "no"),
#'   condition = "item_type"
#' )
#' attr(model, "mpt_renaming")
#' @export
mpt_from_eqn <- function(file, restrictions = NULL, categories = NULL,
                         condition = NULL, covariates = NULL, simplex = NULL,
                         links = "logit") {
  stop_missing_args()
  covariates <- covariates %||% character(0)
  lines <- trimws(readLines(file))
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  fields <- strsplit(lines, "[[:space:]]+")
  if (length(fields) > 0 && length(fields[[1]]) < 3) {
    fields <- fields[-1]
  }
  stopif(
    length(fields) == 0 || !all(lengths(fields) >= 3),
    "Each line of an EQN file must contain three whitespace-separated \\
    fields: tree, response category, and branch probability expression."
  )
  eqn <- data.frame(
    tree = vapply(fields, `[[`, character(1), 1),
    category = vapply(fields, `[[`, character(1), 2),
    expr = vapply(fields, function(f) paste(f[-(1:2)], collapse = ""), character(1))
  )

  eqn$expr <- .mpt_substitute_constants(eqn$expr, restrictions)

  symbols <- unique(unlist(lapply(eqn$expr, function(e) all.vars(str2lang(e)))))
  par_renaming <- .mpt_sanitize_names(setdiff(symbols, covariates))
  for (old in names(par_renaming)) {
    eqn$expr <- gsub(
      paste0("\\b", .mpt_escape_regex(old), "\\b"), par_renaming[[old]], eqn$expr
    )
  }

  category_clean <- eqn$category
  user_mapped <- category_clean %in% names(categories)
  category_clean[user_mapped] <- unlist(categories[category_clean[user_mapped]])
  has_tree_prefix <- !user_mapped &
    startsWith(category_clean, paste0(eqn$tree, "_"))
  category_clean[has_tree_prefix] <- substring(
    category_clean[has_tree_prefix], nchar(eqn$tree[has_tree_prefix]) + 2
  )
  cat_renaming <- .mpt_sanitize_names(unique(category_clean))
  renamed_cats <- unlist(cat_renaming[category_clean], use.names = FALSE)
  cat_renaming <- setNames(renamed_cats, eqn$category)
  cat_renaming <- cat_renaming[!duplicated(names(cat_renaming))]
  eqn$category <- unname(renamed_cats)

  renaming <- c(
    par_renaming[names(par_renaming) != unlist(par_renaming)],
    cat_renaming[names(cat_renaming) != unlist(cat_renaming)]
  )
  if (length(renaming) > 0) {
    message2(
      "Parameter and category names were sanitized for brms (underscores and \\
      dots removed):
      {paste(names(renaming), unlist(renaming), sep = ' -> ', collapse = '; ')}"
    )
  }

  tree_categories <- tapply(eqn$category, eqn$tree, unique)
  stopif(
    length(unique(lapply(tree_categories, sort))) > 1,
    "In bmm all trees must share the same response categories, but the trees \\
    in this EQN file use different category labels:
    {paste(names(tree_categories), vapply(tree_categories, collapse_comma, ''),
           sep = ': ', collapse = '\n    ')}
    Map the tree-specific labels onto shared response categories with the \\
    categories argument, e.g. categories = c(hit = 'yes', fa = 'yes', ...)."
  )

  trees <- lapply(unique(eqn$tree), function(tree_name) {
    tree_lines <- eqn[eqn$tree == tree_name, ]
    mpt_tree(tree_name, .mpt_sum_branch_lines(tree_lines$expr, tree_lines$category))
  })

  model <- mpt(
    trees = trees, condition = condition, covariates = covariates,
    simplex = simplex, links = links
  )
  attr(model, "mpt_renaming") <- renaming
  model
}

.mpt_sum_branch_lines <- function(exprs, categories) {
  summed <- tapply(exprs, categories, function(branch_lines) {
    if (length(branch_lines) == 1) {
      branch_lines
    } else {
      paste0("(", branch_lines, ")", collapse = " + ")
    }
  })
  as.list(summed)[unique(categories)]
}

.mpt_substitute_constants <- function(exprs, restrictions) {
  if (is.null(restrictions)) {
    return(exprs)
  }
  stopif(
    is.null(names(restrictions)) || !all(nzchar(names(restrictions))) ||
      !all(vapply(restrictions, is.numeric, logical(1))),
    "The restrictions argument must be a named list or vector mapping \\
    parameter names to numeric constants. Express equality restrictions by \\
    using the same parameter name in the model definition."
  )
  for (par in names(restrictions)) {
    exprs <- gsub(
      paste0("\\b", .mpt_escape_regex(par), "\\b"),
      format(restrictions[[par]], scientific = FALSE),
      exprs
    )
  }
  exprs
}

# maps each name to a version without underscores and dots (the brms nlpar
# constraint); names that are already valid map to themselves
.mpt_sanitize_names <- function(names) {
  sanitized <- gsub("[._]", "", names)
  clashes <- sanitized[duplicated(sanitized)]
  stopif(
    length(clashes) > 0,
    "Removing underscores and dots produces duplicated names: \\
    {collapse_comma(unique(clashes))}. Please rename them in the model file."
  )
  setNames(as.list(sanitized), names)
}

.mpt_escape_regex <- function(x) {
  gsub("([^A-Za-z0-9_])", "\\\\\\1", x)
}
