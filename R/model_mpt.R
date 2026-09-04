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
#'   the label identifies the tree in the data: the column named by the
#'   `tree_id` argument of [mpt()] holds one such label per observation.
#' @param branches A named list of character strings. Each element gives the
#'   branch probability expression for one response category, and the element
#'   names are the response categories. The expressions can use latent
#'   parameters (e.g., `"D + (1 - D) * g"`), numeric constants, and declared
#'   covariates (see [mpt()]).
#' @param impossible Character vector. Response categories that cannot occur in
#'   this tree, for example distractor responses when no distractors were
#'   presented. Such categories get no branch expression, and the remaining
#'   branches must sum to 1 on their own.
#'
#' @details Numeric fractions such as `1/4` are folded into decimal literals
#'   (`0.25`) when the tree is created. Stan compiles a bare integer fraction
#'   as integer division (`1/4 == 0`), which would silently corrupt the
#'   likelihood.
#'
#'   Categories listed in `impossible` keep their column in the data and their
#'   place in the response matrix, but their probability is set to zero for the
#'   observations belonging to this tree. Because the model is estimated on the
#'   log scale, "zero" is implemented as a large negative linear predictor
#'   rather than a literal zero; no probability has to be taken away from the
#'   remaining branches to compensate.
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
#'
#' # a three-alternative recognition task in which lures are shown on some
#' # trials only: `lure` responses cannot occur on the others
#' tree_lures <- mpt_tree(
#'   name = "lures",
#'   branches = list(
#'     target = "D * b + (1 - D) * (1/3)",
#'     lure   = "D * (1 - b) + (1 - D) * (1/3)",
#'     new    = "(1 - D) * (1/3)"
#'   )
#' )
#'
#' # the branches of the lure-free tree sum to 1 without the `lure` category;
#' # no probability is reserved for it
#' tree_nolures <- mpt_tree(
#'   name = "nolures",
#'   branches = list(
#'     target = "D + (1 - D) * (1/2)",
#'     new    = "(1 - D) * (1/2)"
#'   ),
#'   impossible = "lure"
#' )
#' tree_nolures
#'
#' # a `tree` column in the data holds "lures" or "nolures" per observation;
#' # effects of experimental conditions go in the parameter formulas instead,
#' # e.g. bmf(D ~ 0 + cond)
#' model <- mpt(list(tree_lures, tree_nolures), tree_id = "tree")
#' @export
mpt_tree <- function(name, branches, impossible = NULL) {
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
  impossible <- impossible %||% character(0)
  stopif(
    !is.character(impossible),
    "The impossible argument must be a character vector of response category names."
  )
  stopif(
    anyDuplicated(impossible) > 0,
    "Impossible response categories must be unique within a tree. Duplicated: \\
    {collapse_comma(unique(impossible[duplicated(impossible)]))}"
  )
  impossible_with_branch <- intersect(impossible, resp_cats)
  stopif(
    length(impossible_with_branch) > 0,
    "A response category cannot be both impossible and have a branch \\
    expression in tree '{name}': {collapse_comma(impossible_with_branch)}"
  )
  branches[] <- lapply(resp_cats, function(resp_cat) {
    expr <- branches[[resp_cat]]
    stopif(
      !is.character(expr) || length(expr) != 1L || !nzchar(expr),
      "Branch probability expressions must be single character strings \\
      (tree '{name}', category '{resp_cat}')."
    )
    parsed <- try(str2lang(expr), silent = TRUE)
    stopif(
      is_try_error(parsed),
      "Cannot parse the branch expression '{expr}' in tree '{name}'."
    )
    parsed <- .mpt_fold_numeric_division(parsed)
    scientific <- .mpt_scientific_constants(parsed)
    stopif(
      length(scientific) > 0,
      "The numeric constant(s) {collapse_comma(scientific)} in tree '{name}' \\
      are too extreme to be written into the generated Stan code (brms emits \\
      them in scientific notation, which breaks the Stan syntax). Please \\
      provide such values as a data column declared in the covariates \\
      argument, or use a larger constant."
    )
    parsed
  })
  structure(nlist(name, branches, impossible), class = "mpt_tree")
}

#' @export
print.mpt_tree <- function(x, ...) {
  branch_lines <- glue(
    "  P({names(x$branches)}) = {vapply(x$branches, deparse1, character(1))}"
  )
  if (length(x$impossible) > 0) {
    branch_lines <- c(
      branch_lines, glue("  P({x$impossible}) = 0 (structurally impossible)")
    )
  }
  cat(glue("MPT tree '{x$name}':"), branch_lines, sep = "\n")
  invisible(x)
}

.model_mpt <- function(trees = NULL, tree_id = NULL, covariates = NULL,
                       simplex = NULL, restrictions = NULL, links = "logit",
                       default_priors = NULL, call = NULL, ...) {
  trees <- .mpt_as_tree_list(trees)
  if (length(trees)) names(trees) <- vapply(trees, `[[`, character(1), "name")
  covariates <- covariates %||% character(0)
  simplex <- .mpt_as_simplex_list(simplex)
  resp_cats <- if (length(trees)) {
    c(names(trees[[1]]$branches), trees[[1]]$impossible)
  } else {
    character(0)
  }
  parameters <- setdiff(
    unique(unlist(lapply(trees, .mpt_expr_vars))), covariates
  )
  simplex_pars <- unlist(simplex)
  simplex_raw <- unlist(lapply(simplex, function(grp) {
    free_pars <- grp[-length(grp)]
    setNames(paste0(free_pars, "raw"), free_pars)
  }))
  raw_pars <- unname(simplex_raw)
  standard_pars <- setdiff(parameters, simplex_pars)

  # generated data columns: one 0/1 indicator per tree, and one per response
  # category that is impossible in some trees but not others
  tree_indicators <- if (is.null(tree_id) || length(trees) == 0L) {
    NULL
  } else {
    setNames(paste0("Idx_", names(trees)), names(trees))
  }
  guarded_cats <- unique(unlist(lapply(trees, `[[`, "impossible")))
  possible_indicators <- if (length(guarded_cats) == 0L) {
    NULL
  } else {
    setNames(paste0("Poss_", guarded_cats), guarded_cats)
  }

  # matched priors: effects use the same family as the intercept at half the
  # scale, so both links imply comparable regularization on their latent scale
  latent_prior <- switch(links,
    logit = list(main = "logistic(0, 1)", effects = "logistic(0, 0.5)"),
    probit = list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
  )

  parameter_info <- c(
    .mpt_named_list(
      standard_pars,
      glue("Latent probability parameter from the tree branch expressions. \\
           Fitted on the {links} scale.")
    ),
    .mpt_named_list(
      simplex_pars,
      paste0(
        "Probability parameter constrained within a simplex group via ",
        "stick-breaking. The last parameter of each group is derived as 1 ",
        "minus the sum of the other group members."
      )
    ),
    .mpt_named_list(
      raw_pars,
      glue("Unconstrained stick-breaking component of a simplex parameter. \\
           Fitted on the {links} scale.")
    )
  )

  link_info <- c(
    .mpt_named_list(standard_pars, links),
    .mpt_named_list(simplex_pars, "identity"),
    .mpt_named_list(raw_pars, "identity")
  )

  prior_info <- c(
    .mpt_named_list(standard_pars, latent_prior),
    .mpt_named_list(raw_pars, latent_prior)
  )

  out <- structure(
    list(
      resp_vars = nlist(resp_cats),
      other_vars = nlist(tree_id, covariates),
      trees = trees,
      simplex = simplex,
      restrictions = restrictions,
      link = links,
      indicators = list(tree = tree_indicators, possible = possible_indicators),
      simplex_raw = simplex_raw,
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
        "- One tree per distinct branch structure, created with mpt_tree(); ",
        "all trees share the same response categories\n",
        "  - The data contain one column with aggregated response counts per ",
        "response category, named after the branch names\n",
        "  - For multi-tree models, a column whose values name the tree each ",
        "observation belongs to, declared via the tree_id argument\n",
        "  - Data columns used inside branch expressions must be declared via the covariates argument\n"
      ),
      parameters = parameter_info,
      fixed_parameters = list(),
      links = link_info,
      default_priors = prior_info
    ),
    class = c("bmmodel", "mpt"),
    call = call
  )
  out$default_priors[names(default_priors)] <- default_priors
  out
}

# restrictions are already substituted into the stored trees, so they are
# not passed again
.mpt_constructor_args <- function(model) {
  list(
    trees = unname(model$trees),
    tree_id = model$other_vars$tree_id,
    covariates = model$other_vars$covariates,
    simplex = model$simplex,
    links = model$link
  )
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
#' specified as a set of trees (one per distinct branch structure) created
#' with [mpt_tree()]. Branch expressions that terminate in the same response
#' category are summed within each tree, and parameters shared across trees
#' are equated by giving them the same name. All latent probability
#' parameters are estimated on an unconstrained latent scale (logit or
#' probit); predictor formulas supplied via [bmmformula()] apply on that
#' latent scale.
#'
#' @param trees A single `mpt_tree` object or a list of `mpt_tree` objects.
#'   All trees must share the same set of response categories, counting those
#'   declared impossible.
#' @param tree_id Character. Name of the data column whose values name the tree
#'   each observation belongs to; the values must match the tree names. Can be
#'   omitted for single-tree models. Effects of experimental conditions belong
#'   in the parameter formulas, not here — see Details.
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
#' @param restrictions Parameter restrictions in the string syntax of MPTinR
#'   and TreeBUGS, e.g. `c("Dn = Do", "g = 0.5")`, or as a named list,
#'   `list(Dn = "Do", g = 0.5)`. A restriction either equates a parameter with
#'   another one (chains such as `"G1 = G2 = G3"` map all earlier names onto
#'   the last) or fixes it to a numeric constant (`"g = 0.5"`, `"g = 1/4"`).
#'   Restrictions are substituted into the branch expressions before the
#'   parameters are identified, so a restricted parameter is not part of the
#'   model. Order constraints (`"Do > Dn"`) are not supported here; see
#'   Details.
#' @param links Character. The link function for all latent probability
#'   parameters: `"logit"` (default) or `"probit"`.
#' @param ... used internally for testing, ignore it
#'
#' @details `r model_info(.model_mpt(), components = c('domain', 'task', 'name', 'citation'))`
#'
#'   A separate tree is needed only when the *branch expressions* differ — that
#'   is, when the trial determines which latent processes apply or which
#'   response categories are reachable (old versus new probes, a response set
#'   with or without distractors, inclusion versus exclusion instructions).
#'   Experimental factors whose *effect on the parameters* you want to estimate
#'   belong in the parameter formulas, exactly as in every other `bmm` model:
#'   `bmf(D ~ 0 + condition)`. The two roles are independent, and one data
#'   column can serve both when a factor happens to change the structure as
#'   well. When several levels of a factor share a branch structure, add a
#'   column naming the tree and keep the factor for the formulas.
#'
#'   Parameters that receive a non-linear predictor formula (a formula whose
#'   right-hand side references other formula parameters, e.g.
#'   `D ~ Dmax * (1 - exp(-rate * ptime))`) are not transformed by the link
#'   function: the user-supplied expression must keep the parameter within
#'   (0, 1). Sub-parameters of such formulas (e.g., `Dmax` and `rate`) are
#'   estimated on the identity scale with `normal(0, 1)` default priors.
#'
#'   A parameter can also be fixed to a probability at fit time, in the
#'   formula: `bmf(D ~ 1, g = 0.5)`. Unlike a restriction, the parameter stays
#'   part of the model and can be freed again by giving it a formula. The
#'   value is mapped to the latent scale when the constant prior is built, so
#'   [default_prior()] shows `constant(0)` for a guessing rate of 0.5 under
#'   the logit link.
#'
#'   Order constraints between parameters (`Do > Dn`) are expressed by
#'   reparameterizing the larger parameter in the model formula, e.g.
#'   `bmf(Do ~ Dn + (1 - Dn) * inv_logit(phi), Dn ~ 1, phi ~ 1)`; the section
#'   "Ordered parameter constraints" of the MPT article walks through the
#'   recipe.
#'
#'   Parameter and response category names must start with a letter and may
#'   contain only letters and digits, because brms does not allow underscores
#'   or dots in non-linear parameter names.
#'
#'   All latent parameters are non-linear parameters of the underlying brms
#'   model, and brms has no separate `Intercept` prior class for those: their
#'   intercept prior is stored as `class = "b", coef = "Intercept"`. Overriding
#'   the default prior of an intercept therefore requires `coef = "Intercept"`;
#'   a prior given as `class = "b", nlpar = "D"` reaches the remaining
#'   coefficients only. Call [default_prior()] to see which rows a given
#'   formula produces.
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
#'   tree_id = "item_type"
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
mpt <- function(trees, tree_id = NULL, covariates = NULL, simplex = NULL,
                restrictions = NULL, links = "logit", ...) {
  call <- match.call()
  stop_missing_args()
  links <- match.arg(links, c("logit", "probit"))
  trees <- .mpt_as_tree_list(trees)
  covariates <- covariates %||% character(0)
  simplex <- .mpt_as_simplex_list(simplex)
  restrictions <- .mpt_parse_restrictions(restrictions)

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

  tree_cats <- lapply(trees, function(tree) c(names(tree$branches), tree$impossible))
  resp_cats <- tree_cats[[1]]
  cats_match <- vapply(tree_cats, setequal, logical(1), resp_cats)
  stopif(
    !all(cats_match),
    "All trees must have the same response categories, counting those declared
    impossible.
    Tree '{trees[[1]]$name}' has: {collapse_comma(resp_cats)}
    Tree '{tree_names[!cats_match][1]}' has: {collapse_comma(tree_cats[!cats_match][[1]])}"
  )
  always_impossible <- Reduce(intersect, lapply(trees, `[[`, "impossible"))
  stopif(
    length(always_impossible) > 0,
    "The response category(ies) {collapse_comma(always_impossible)} are \\
    impossible in every tree and are therefore not identified. Please remove \\
    them from the model."
  )
  trees <- lapply(trees, function(tree) {
    tree$branches <- tree$branches[intersect(resp_cats, names(tree$branches))]
    tree
  })
  names(trees) <- tree_names

  stopif(
    length(trees) > 1L && is.null(tree_id),
    "Models with multiple trees require the tree_id argument: the name of the \\
    data column whose values identify the tree each observation belongs to."
  )
  stopif(
    !is.null(tree_id) && (!is.character(tree_id) || length(tree_id) != 1L),
    "The tree_id argument must be a single character string naming a data column."
  )
  stopif(
    length(covariates) > 0L && !is.character(covariates),
    "The covariates argument must be a character vector of data column names."
  )

  trees <- .mpt_restrict_trees(trees, restrictions, covariates)

  parameters <- setdiff(unique(unlist(lapply(trees, .mpt_expr_vars))), covariates)
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
    trees = trees, tree_id = tree_id, covariates = covariates,
    simplex = simplex, restrictions = restrictions, links = links,
    call = call, ...
  )
}

# restrictions are checked against the symbols of the unrestricted trees and
# then substituted, so every later step sees the restricted model only
.mpt_restrict_trees <- function(trees, restrictions, covariates) {
  if (length(restrictions) == 0L) {
    return(trees)
  }
  restricted <- names(restrictions)
  stopif(
    anyDuplicated(restricted) > 0,
    "Parameters cannot be restricted more than once: \\
    {collapse_comma(unique(restricted[duplicated(restricted)]))}"
  )
  restricted_covariates <- intersect(restricted, covariates)
  stopif(
    length(restricted_covariates) > 0,
    "Covariates are data columns and cannot be restricted: \\
    {collapse_comma(restricted_covariates)}"
  )
  parameters <- setdiff(unique(unlist(lapply(trees, .mpt_expr_vars))), covariates)
  restrictions <- .mpt_resolve_restrictions(restrictions)
  targets <- unique(unlist(lapply(restrictions, all.vars)))
  unknown <- setdiff(c(restricted, targets), parameters)
  stopif(
    length(unknown) > 0,
    "Restrictions refer to parameters that do not appear in the tree branch \\
    expressions: {collapse_comma(unknown)}"
  )
  circular <- intersect(targets, restricted)
  stopif(
    length(circular) > 0,
    "The restrictions on {collapse_comma(circular)} are circular."
  )
  scientific <- unlist(lapply(restrictions, .mpt_scientific_constants))
  stopif(
    length(scientific) > 0,
    "The restriction constant(s) {collapse_comma(scientific)} are too extreme \\
    to be written into the generated Stan code (brms emits them in scientific \\
    notation, which breaks the Stan syntax). Please provide such values as a \\
    data column declared in the covariates argument, or use a larger constant."
  )
  lapply(trees, .mpt_apply_restrictions, restrictions)
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
    .mpt_check_fixed_values(model)
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
    c(model$other_vars$covariates, model$other_vars$tree_id)
  )
  if (length(sub_pars) > 0) {
    .mpt_check_names(sub_pars, "parameter")
    message2(
      "The parameter(s) {collapse_comma(sub_pars)} from your non-linear \\
      formulas are estimated on the identity scale with normal(0, 1) default \\
      priors. Apply any required transformation inside your formula and \\
      adjust the priors to the scale of your predictors."
    )
    model$parameters[sub_pars] <- .mpt_named_list(
      sub_pars, "User-defined sub-parameter of a non-linear parameter formula."
    )
    model$links[sub_pars] <- .mpt_named_list(sub_pars, "identity")
    model$default_priors[sub_pars] <- .mpt_named_list(
      sub_pars, list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
    )
  }
  model
}

# fixed values from the bmmformula (e.g. g = 0.5) stay on the probability
# scale in the model object; configure_prior.mpt() maps them to the latent
# scale when it builds the constant priors
.mpt_check_fixed_values <- function(model) {
  fixed_simplex <- intersect(names(model$fixed_parameters), unlist(model$simplex))
  stopif(
    length(fixed_simplex) > 0,
    "Fixing simplex parameters to constants is not supported: \\
    {collapse_comma(fixed_simplex)}"
  )
  for (par in .mpt_latent_fixed_pars(model)) {
    value <- model$fixed_parameters[[par]]
    stopif(
      !is.numeric(value) || value <= 0 || value >= 1,
      "The fixed value for parameter '{par}' must be a probability strictly \\
      between 0 and 1. Provided: {value}"
    )
  }
  invisible(NULL)
}

# fixed latent probability parameters; sub-parameters of non-linear formulas
# have an identity link and take their fixed value as is
.mpt_latent_fixed_pars <- function(model) {
  fixed_pars <- intersect(names(model$fixed_parameters), names(model$parameters))
  fixed_pars[vapply(fixed_pars, function(par) {
    model$links[[par]] != "identity"
  }, logical(1))]
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

  tree_id <- model$other_vars$tree_id
  if (!is.null(tree_id)) {
    stopif(
      !tree_id %in% colnames(data),
      "The tree identifier column '{tree_id}' is not present in the data."
    )
    tree_values <- as.character(data[[tree_id]])
    tree_names <- names(model$trees)
    unmatched <- setdiff(unique(tree_values), tree_names)
    stopif(
      length(unmatched) > 0,
      "All values of the tree identifier column '{tree_id}' must match a tree name.
      Tree names: {collapse_comma(tree_names)}
      Unmatched values: {collapse_comma(unmatched)}
      Values of an experimental factor that share a branch structure belong to
      the same tree; add a column that names the tree for each observation and
      keep the factor for the parameter formulas."
    )
    unused_trees <- setdiff(tree_names, unique(tree_values))
    warnif(
      length(unused_trees) > 0,
      "The data contain no observations for tree(s): {collapse_comma(unused_trees)}"
    )
    # brms needs the indicators as data columns, so they are handed on through
    # the data rather than through the attribute bridge
    idx_vars <- model$indicators$tree
    idx_collisions <- intersect(idx_vars, colnames(data))
    stopif(
      length(idx_collisions) > 0,
      "The data contain column(s) {collapse_comma(idx_collisions)}, which are \\
      reserved for the generated tree indicator variables. Please rename them."
    )
    for (tree_name in tree_names) {
      data[[idx_vars[[tree_name]]]] <- as.integer(tree_values == tree_name)
    }
  }

  data <- .mpt_possibility_indicators(model, data)

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

# structurally impossible categories are switched off per row rather than per
# tree, because a category can be impossible in some trees and not in others
.mpt_possibility_indicators <- function(model, data) {
  poss_vars <- model$indicators$possible
  if (is.null(poss_vars)) {
    return(data)
  }
  poss_collisions <- intersect(poss_vars, colnames(data))
  stopif(
    length(poss_collisions) > 0,
    "The data contain column(s) {collapse_comma(poss_collisions)}, which are \\
    reserved for the generated indicators of structurally impossible response \\
    categories. Please rename them."
  )
  # a category impossible in every tree is rejected by mpt(), so any model that
  # reaches here has several trees and therefore tree indicator columns
  for (resp_cat in names(poss_vars)) {
    impossible_rows <- logical(nrow(data))
    for (tree in model$trees) {
      if (!resp_cat %in% tree$impossible) next
      impossible_rows <- impossible_rows |
        data[[model$indicators$tree[[tree$name]]]] == 1L
    }
    observed <- which(impossible_rows & data$Y[, resp_cat] > 0)
    stopif(
      length(observed) > 0,
      "The response category '{resp_cat}' is declared impossible for \\
      {length(observed)} observation(s), but the data record responses in it \\
      (first: row {observed[1]}, count {data$Y[observed[1], resp_cat]}). \\
      Please check the impossible argument of mpt_tree() and the response counts."
    )
    data[[poss_vars[[resp_cat]]]] <- as.integer(!impossible_rows)
  }
  data
}

# the construction-time branch-sum validation uses synthetic covariate values;
# once the data are known, the sum-to-1 property is re-checked row by row with
# the observed covariate values to catch data-preparation errors
.mpt_validate_covariate_sums <- function(model, data, tolerance = 1e-6) {
  covariates <- model$other_vars$covariates
  if (length(covariates) == 0L) {
    return(invisible(NULL))
  }
  parameters <- names(model$parameters)
  par_vals <- setNames(rep(0.5, length(parameters)), parameters)
  for (grp in model$simplex) {
    par_vals[grp] <- 1 / length(grp)
  }

  idx_vars <- model$indicators$tree
  for (tree in model$trees) {
    rows <- if (is.null(idx_vars)) {
      seq_len(nrow(data))
    } else {
      which(data[[idx_vars[[tree$name]]]] == 1L)
    }
    if (length(rows) == 0L) next
    env <- c(as.list(par_vals), as.list(data[rows, covariates, drop = FALSE]))
    total <- Reduce(`+`, lapply(tree$branches, eval, envir = env))
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
      raw_par <- model$simplex_raw[[par]]
      par_rhs <- .mpt_rhs_chr(formula[[par]])
      raw_rhs <- .mpt_rhs_chr(formula[[raw_par]])
      stopif(
        par_rhs != "1" && raw_rhs != "1" && par_rhs != raw_rhs,
        "Conflicting predictor formulas for the simplex parameter '{par}' and \\
        its stick-breaking component '{raw_par}'. Specify predictors for \\
        '{par}' only."
      )
      if (par_rhs != "1") {
        formula[raw_par] <- list(
          stats::as.formula(call("~", as.name(raw_par), formula[[par]][[3]]))
        )
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

# a missing formula is treated as the intercept-only default the pipeline
# would add, so direct calls before add_missing_parameters() cannot produce NA
.mpt_rhs_chr <- function(pform) {
  if (!is_formula(pform)) {
    return("1")
  }
  paste(deparse(pform[[3]]), collapse = " ")
}

.mpt_category_formulas <- function(model) {
  trees <- model$trees
  idx_vars <- model$indicators$tree
  category_formulas <- lapply(model$resp_vars$resp_cats, function(resp_cat) {
    # a placeholder keeps log() defined for trees where the category is
    # impossible; bmf2bf.mpt() overrides the linear predictor for those rows
    terms <- lapply(names(trees), function(tree_name) {
      branch <- call("(", trees[[tree_name]]$branches[[resp_cat]] %||% 1)
      if (is.null(idx_vars)) {
        branch
      } else {
        call("*", as.name(idx_vars[[tree_name]]), branch)
      }
    })
    .mpt_formula(resp_cat, .mpt_reduce_calls("+", terms))
  })
  do.call(bmf, category_formulas)
}

.mpt_simplex_formulas <- function(model) {
  if (length(model$simplex) == 0L) {
    return(NULL)
  }
  simplex_formulas <- list()
  for (grp in model$simplex) {
    n_grp <- length(grp)
    free_pars <- grp[-n_grp]
    sticks <- lapply(model$simplex_raw[free_pars], inv_link, link = model$link)
    for (k in seq_len(n_grp - 1L)) {
      remaining <- lapply(sticks[seq_len(k - 1L)], function(stick) {
        call("(", call("-", 1, stick))
      })
      simplex_formulas <- c(
        simplex_formulas,
        list(.mpt_formula(grp[k], .mpt_reduce_calls("*", c(remaining, sticks[k]))))
      )
    }
    free_sum <- .mpt_reduce_calls("+", lapply(free_pars, as.name))
    simplex_formulas <- c(
      simplex_formulas,
      list(.mpt_formula(grp[n_grp], call("-", 1, call("(", free_sum))))
    )
  }
  do.call(bmf, simplex_formulas)
}

.mpt_formula <- function(lhs, rhs) {
  stats::as.formula(call("~", as.name(lhs), rhs))
}

.mpt_reduce_calls <- function(op, calls) {
  Reduce(function(lhs, rhs) call(op, lhs, rhs), calls)
}

############################################################################# !
# Convert bmmformula to brmsformula methods                              ####
############################################################################# !

#' @export
bmf2bf.mpt <- function(model, formula) {
  resp_cats <- model$resp_vars$resp_cats
  poss_vars <- model$indicators$possible

  # -100 on the log scale is an effectively zero probability that stays finite;
  # the multinomial family renormalizes, so nothing has to compensate for it
  linpreds <- vapply(resp_cats, function(resp_cat) {
    if (!resp_cat %in% names(poss_vars)) {
      return(glue("log({resp_cat})"))
    }
    poss <- poss_vars[[resp_cat]]
    glue("{poss} * log({resp_cat}) + (1 - {poss}) * (-100)")
  }, character(1))

  brms_formula <- brms::bf(
    glue("Y | trials(nTrials) ~ {linpreds[1]}"),
    nl = TRUE
  )
  for (i in seq_along(resp_cats)[-1]) {
    brms_formula <- brms_formula + glue_nlf("mu{resp_cats[i]} ~ {linpreds[i]}")
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
# CONFIGURE_PRIOR S3 METHODS                                             ####
############################################################################# !

# the constant() prior applies to the latent intercept, so the fixed
# probability is mapped through the link; combined last by
# configure_prior.bmmodel, it overrides the probability-scale constant from
# fixed_pars_priors()
#' @export
configure_prior.mpt <- function(model, data, formula, user_prior, ...) {
  fixed_pars <- .mpt_latent_fixed_pars(model)
  if (length(fixed_pars) == 0L) {
    return(brms::empty_prior())
  }
  latent <- vapply(fixed_pars, function(par) {
    link_transform(model$fixed_parameters[[par]], model$links[[par]])
  }, numeric(1))
  brms::set_prior(
    glue("constant({latent})"),
    class = "b", coef = "Intercept", nlpar = fixed_pars
  )
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

.mpt_check_names <- function(names, what) {
  bad_names <- names[!grepl("^[A-Za-z][A-Za-z0-9]*$", names)]
  stopif(
    length(bad_names) > 0,
    "MPT {what} names must start with a letter and contain only letters and \\
    digits, because brms does not allow underscores or dots in non-linear \\
    parameter names. Please rename (e.g., 'd_A' -> 'dA'): {collapse_comma(bad_names)}"
  )
}

.mpt_named_list <- function(names, value) {
  if (length(names) == 0L) {
    return(list())
  }
  setNames(rep(list(value), length(names)), names)
}
