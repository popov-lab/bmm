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
#'   parameters (e.g., `"D + (1 - D) * g"`) and numeric constants.
#'
#' @details Constant arithmetic such as `1/4` or `1/(2*2)` is folded into a
#'   decimal literal (`0.25`) when the tree is created. Stan compiles a bare
#'   integer fraction as integer division (`1/4 == 0`), which would silently
#'   corrupt the likelihood.
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
      them in scientific notation, which breaks the Stan syntax). Please use \\
      a larger constant."
    )
    parsed
  })
  structure(nlist(name, branches), class = "mpt_tree")
}

#' @export
print.mpt_tree <- function(x, ...) {
  branch_lines <- glue(
    "  P({names(x$branches)}) = {vapply(x$branches, deparse1, character(1))}"
  )
  cat(glue("MPT tree '{x$name}':"), branch_lines, sep = "\n")
  invisible(x)
}

.model_mpt <- function(trees = NULL, tree_id = NULL, simplex = NULL,
                       restrictions = NULL, unrestricted_trees = NULL,
                       links = "logit",
                       default_priors = NULL, call = NULL, ...) {
  trees <- .mpt_as_tree_list(trees)
  if (length(trees)) names(trees) <- vapply(trees, `[[`, character(1), "name")
  simplex <- .mpt_as_simplex_list(simplex)
  resp_cats <- if (length(trees)) {
    names(trees[[1]]$branches)
  } else {
    character(0)
  }
  parameters <- unique(unlist(lapply(trees, .mpt_expr_vars)))
  simplex_pars <- unlist(simplex)
  simplex_raw <- unlist(lapply(simplex, function(grp) {
    free_pars <- grp[-length(grp)]
    setNames(paste0(free_pars, "raw"), free_pars)
  }))
  raw_pars <- unname(simplex_raw)
  standard_pars <- setdiff(parameters, simplex_pars)

  # generated data columns: one 0/1 indicator per tree
  tree_indicators <- if (is.null(tree_id) || length(trees) == 0L) {
    NULL
  } else {
    setNames(paste0("Idx_", names(trees)), names(trees))
  }

  latent_prior <- .mpt_latent_prior(links)

  parameter_info <- c(
    # the scale is left out because the link of each parameter can be changed
    # after construction; print() lists the links
    .mpt_named_list(
      standard_pars,
      "Latent probability parameter from the tree branch expressions."
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
      glue("Unconstrained stick-breaking component of a simplex parameter, \\
           transformed by the inverse {links} into the stick proportion of the \\
           remaining probability.")
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
      other_vars = list(
        tree_id = tree_id,
        trees = trees,
        unrestricted_trees = unrestricted_trees,
        simplex = simplex,
        restrictions = restrictions,
        link = links,
        indicators = list(tree = tree_indicators),
        simplex_raw = simplex_raw
      ),
      domain = "Categorical decision making, memory, and reasoning",
      task = "Any task with categorical responses generated by a discrete processing-tree structure",
      name = "Multinomial Processing Tree (MPT) models",
      citation = glue(
        "Batchelder, W. H., & Riefer, D. M. (1999). Theoretical and empirical \\
        review of multinomial process tree modeling. Psychonomic Bulletin & \\
        Review, 6(1), 57-86. https://doi.org/10.3758/BF03210812"
      ),
      version = "NA",
      requirements = paste0(
        "- One tree per distinct branch structure, created with mpt_tree(); ",
        "all trees share the same response categories\n",
        "  - The data contain one column with aggregated response counts per ",
        "response category, named after the branch names\n",
        "  - For multi-tree models, a column whose values name the tree each ",
        "observation belongs to, declared via the tree_id argument\n"
      ),
      parameters = parameter_info,
      links = link_info,
      fixed_parameters = list(),
      default_priors = prior_info
    ),
    class = c("bmmodel", "mpt"),
    call = call
  )
  out$default_priors[names(default_priors)] <- default_priors
  set_links(out, NULL)
}

# matched priors: effects get the intercept's own prior, because a narrower
# effects prior shrinks condition contrasts under `1 + cond` coding (up to 0.04
# in probability for realistic 2HTM designs at half the scale). The sd rate of
# the logit link is the one of the mixture weights; a probit SD is about 1/1.7
# of a logit SD, so the probit rate is about 1.7, rounded to 2
.mpt_latent_prior <- function(link) {
  switch(link,
    logit = list(
      main = "logistic(0, 1)", effects = "logistic(0, 1)", sd = "exponential(1)"
    ),
    probit = list(
      main = "normal(0, 1)", effects = "normal(0, 1)", sd = "exponential(2)"
    )
  )
}

# the link of a latent probability parameter can be switched between logit
# and probit; simplex members and their stick-breaking components carry
# identity links because the generated simplex formulas apply the model link
#' @exportS3Method
settable_links.mpt <- function(model) {
  defaults <- attr(model, "links_default") %||% model$links
  names(defaults)[defaults != "identity"]
}

#' @exportS3Method
settable_link_functions.mpt <- function(model) {
  c("logit", "probit")
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
#'   All trees must share the same set of response categories.
#' @param tree_id Character. Name of the data column whose values name the tree
#'   each observation belongs to; the values must match the tree names. Can be
#'   omitted for single-tree models. Effects of experimental conditions belong
#'   in the parameter formulas, not here — see Details.
#' @param simplex A character vector, or a list of character vectors, naming
#'   groups of parameters that are jointly constrained to sum to 1. Each group
#'   is reparameterized via stick-breaking: the last parameter of each group
#'   is derived as 1 minus the sum of the other members, and each other member
#'   `p` gets an unconstrained component `praw` that predictor formulas for
#'   `p` are applied to. The listing order sets the default prior and the
#'   meaning of the coefficients; see Details.
#' @param restrictions Parameter restrictions in the string syntax of MPTinR
#'   and TreeBUGS, e.g. `c("Dn = Do", "g = 0.5")`, or as a named list,
#'   `list(Dn = "Do", g = 0.5)`. A restriction either equates a parameter with
#'   another one (chains such as `"G1 = G2 = G3"` map all earlier names onto
#'   the last) or fixes it to a numeric constant (`"g = 0.5"`, `"g = 1/4"`).
#'   Constants must lie strictly between 0 and 1. The unnamed list of MPTinR
#'   and TreeBUGS, `list("Dn = Do", "g = 0.5")`, works as well. Restriction
#'   files are not read: in MPTinR and TreeBUGS a character vector names a
#'   restrictions file, in `bmm` it holds the restrictions. To use such a file,
#'   pass `readLines(path)`; blank lines and lines starting with `#` are
#'   skipped, and a `#` after a restriction starts a comment, as in MPTinR.
#'   Restrictions are substituted into the branch expressions before the
#'   parameters are identified, so a restricted parameter is not part of the
#'   model and cannot be a simplex member. Order constraints (`"Do > Dn"`) are
#'   not supported here; see Details.
#' @param links Character. The link function for all latent probability
#'   parameters: `"logit"` (default) or `"probit"`.
#' @param ... used internally for testing, ignore it
#'
#' @details `r model_docs(.model_mpt(), components = c('domain', 'task', 'name', 'citation'))`
#'
#'   The data hold aggregated counts: one row per participant and tree (and
#'   per condition, if conditions are crossed with the trees), with one count
#'   column per response category. The count columns must be named exactly
#'   like the branch names given to [mpt_tree()] (`old` and `new` in the
#'   example below); `bmm()` finds them by these names, so no response
#'   variable is passed. The `tree_id` column names the tree of each row, and
#'   the counts of a row must come from trials of that tree.
#'
#'   A separate tree is needed only when the *branch expressions* differ — that
#'   is, when the trial determines which latent processes apply (old versus
#'   new probes, inclusion versus exclusion instructions).
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
#'   estimated on the identity scale with `normal(0, 1)` default priors. Their
#'   random-effect SDs keep brms's default prior, `student_t(3, 0, 2.5)`,
#'   because no single rate fits every scale such a formula can give them.
#'
#'   Under the logit link, intercepts and regression coefficients get a
#'   `logistic(0, 1)` prior and random-effect SDs get `exponential(1)`. Under
#'   the probit link, intercepts and coefficients get `normal(0, 1)` and SDs
#'   get `exponential(2)`. Group-level correlation matrices get `lkj(2)`.
#'   Coefficients of the sub-parameters above get `normal(0, 0.5)`. A parameter
#'   switched to the other link after construction (`model$links$D <-
#'   "probit"`) takes that link's priors unless its default prior was
#'   replaced. [default_prior()] lists the rows a given formula produces.
#'
#'   In a simplex group with intercept-only formulas, the default prior makes
#'   each stick uniform on (0, 1) under both links. The prior means are therefore 0.50/0.25/0.25 for a
#'   group of three and 0.50/0.25/0.125/0.125 for a group of four: the member
#'   listed first always has prior mean 0.5, so the listing order matters. A
#'   symmetric (uniform Dirichlet) prior would put 1/K on each of K members. A
#'   predictor on member `p` acts on the logit (or probit) of `p`'s share of
#'   the probability left over by the members listed before it. A predictor
#'   on the first member therefore scales all later members by the same
#'   factor, leaving their ratios unchanged. Read such effects on the probability scale with
#'   [native_parameters()].
#'
#'   A parameter can also be fixed to a probability at fit time, in the
#'   formula: `bmf(D ~ 1, g = 0.5)`. Unlike a restriction, the parameter stays
#'   part of the model and can be freed again by giving it a formula. The
#'   value is mapped to the latent scale when the constant prior is built, so
#'   [default_prior()] shows `constant(0)` for a guessing rate of 0.5 under
#'   the logit link. For the same reason, `summary()` lists a fixed `g = 0.5`
#'   as `g_Intercept 0.00` under "Constant Parameters"; [prior_info()] and
#'   [parameter_info()] show the probability, 0.5.
#'
#'   Printing the model ends with an identifiability check for intercept-only
#'   formulas. It first compares the number of free parameters with the
#'   degrees of freedom (response categories minus one, summed over trees).
#'   It then computes the rank of the Jacobian of all category probabilities
#'   with respect to the free parameters, from exact derivatives at five
#'   interior test values; parameters fixed in the formula enter at their
#'   values. A simplex group counts with one free parameter fewer than its
#'   members, and the printout names the group by its members. A rank below
#'   the number of free parameters means that some combination of the listed
#'   parameters cannot be estimated from the data,
#'   even when the count passes, and its posterior follows the prior; a
#'   parameter whose derivatives are zero up to rounding at the test values
#'   appears not to affect any category probability and is named separately.
#'   Fix parameters in the formula or equate them in the branch expressions
#'   until the rank is full. `bmm()` warns about a rank deficit when no
#'   formula uses a data column as a predictor. When one does, a parameter
#'   that differs between conditions can identify the model across them, so
#'   `bmm()` only says that the model is not identified within one design
#'   cell and that the predictors were not checked. The
#'   check is local: it holds at the test values, not at the boundaries of
#'   the parameter space. A branch expression with a function that
#'   [stats::D()] cannot differentiate leaves the rank uncomputed, and the
#'   printout says so.
#'
#'   `summary()` reports intercepts and regression coefficients on the latent
#'   (logit or probit) scale. [native_parameters()] returns the posterior
#'   draws of every parameter on the probability scale for each observed
#'   combination of the predictors, e.g.
#'   `native_parameters(fit, re_formula = NA, summary = TRUE)`; a difference
#'   between conditions is the difference of these draws.
#'
#'   Person effects written as `(1 | id)` in each parameter's formula are
#'   independent across parameters. A shared label, `(1 |p| id)` in every
#'   formula, estimates their correlations (prior `lkj(2)`), which
#'   corresponds to the latent-trait MPT model.
#'
#'   For posterior predictive checks, [brms::posterior_predict()] returns
#'   simulated counts (draws x rows x categories) and [brms::posterior_epred()]
#'   the expected counts; the MPT article computes the T1 statistic and its
#'   posterior predictive p-value from them.
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
#' @examples
#' \dontrun{
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
#' # simulate data for 20 participants with D = 0.7, g = 0.5:
#' # one row per participant x tree; count columns named after the branches
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
#' # predict both parameters by a fixed intercept and correlated random
#' # intercepts (the shared label |p| estimates the correlation)
#' formula <- bmf(
#'   D ~ 1 + (1 |p| id),
#'   g ~ 1 + (1 |p| id)
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
#' # population-level D and g on the probability scale
#' native_parameters(fit, re_formula = NA, summary = TRUE)
#' }
#'
#' @export
mpt <- function(trees, tree_id = NULL, simplex = NULL, restrictions = NULL,
                links = "logit", ...) {
  call <- match.call()
  stop_missing_args()
  links <- match.arg(links, c("logit", "probit"))
  trees <- .mpt_as_tree_list(trees)
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

  tree_cats <- lapply(trees, function(tree) names(tree$branches))
  resp_cats <- tree_cats[[1]]
  cats_match <- vapply(tree_cats, setequal, logical(1), resp_cats)
  stopif(
    !all(cats_match),
    "All trees must have the same response categories.
    Tree '{trees[[1]]$name}' has: {collapse_comma(resp_cats)}
    Tree '{tree_names[!cats_match][1]}' has: {collapse_comma(tree_cats[!cats_match][[1]])}"
  )
  trees <- lapply(trees, function(tree) {
    tree$branches <- tree$branches[resp_cats]
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

  # plot() labels the restricted edges from the trees as written
  unrestricted_trees <- if (length(restrictions)) trees
  trees <- .mpt_restrict_trees(trees, restrictions)

  # a zero probability makes log(p) undefined in Stan; a dedicated argument
  # for impossible categories replaces this guard in a later release. It runs
  # after the restrictions, which can fold a branch to 0 (a = 0.5 in "1 - 2 * a").
  zero_branches <- unlist(lapply(trees, function(tree) {
    is_zero <- vapply(tree$branches, function(b) is.numeric(b) && b == 0, logical(1))
    glue("'{names(tree$branches)[is_zero]}' in tree '{tree$name}'")
  }))
  stopif(
    length(zero_branches) > 0,
    "Branch probabilities that are the constant 0: \\
    {paste(zero_branches, collapse = ', ')}. A zero probability makes the \\
    likelihood undefined. Response categories that a tree cannot produce will \\
    be declared with mpt_tree(impossible = ) in a later release; this version \\
    does not support them."
  )

  parameters <- unique(unlist(lapply(trees, .mpt_expr_vars)))
  stopif(
    length(parameters) == 0L,
    "The tree branch expressions contain no latent parameters."
  )

  .mpt_check_names(parameters, "parameter")
  .mpt_check_names(resp_cats, "response category")

  par_cat_overlap <- intersect(parameters, resp_cats)
  stopif(
    length(par_cat_overlap) > 0,
    "Names cannot be used for both a parameter and a response category: \\
    {collapse_comma(par_cat_overlap)}"
  )
  reserved <- intersect(c(parameters, resp_cats), c("Y", "nTrials"))
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
  restricted_simplex <- intersect(simplex_pars, names(restrictions))
  stopif(
    length(restricted_simplex) > 0,
    "The simplex parameter(s) {collapse_comma(restricted_simplex)} are also \\
    restricted, and the restriction removed them from the trees. Restrict \\
    simplex members only through the group itself."
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
  raw_collisions <- intersect(raw_pars, c(parameters, resp_cats))
  stopif(
    length(raw_collisions) > 0,
    "The stick-breaking components of simplex parameters are named by \\
    appending 'raw' to the parameter name, but these names are already in \\
    use: {collapse_comma(raw_collisions)}"
  )

  deviations <- .mpt_tree_sum_deviations(trees, parameters, simplex)
  for (tree_name in names(deviations)[!is.na(deviations)]) {
    stop2(
      "The branch probabilities of tree '{tree_name}' sum to \\
      {signif(deviations[[tree_name]], 6)} instead of 1 when evaluated at \\
      numeric test values. Please check the branch expressions. To equate \\
      parameters or fix one to a constant, use the restrictions argument \\
      (e.g. restrictions = 'Dn = Do'), not the branch expressions or the formula. \\
      Parameters that must sum to 1 across branches (e.g. guessing over three \\
      options) belong in the simplex argument."
    )
  }

  .model_mpt(
    trees = trees, tree_id = tree_id, simplex = simplex,
    restrictions = restrictions, unrestricted_trees = unrestricted_trees,
    links = links, call = call, ...
  )
}

# restrictions are checked against the symbols of the unrestricted trees and
# then substituted, so every later step sees the restricted model only
.mpt_restrict_trees <- function(trees, restrictions) {
  if (length(restrictions) == 0L) {
    return(trees)
  }
  restricted <- names(restrictions)
  stopif(
    anyDuplicated(restricted) > 0,
    "Parameters cannot be restricted more than once: \\
    {collapse_comma(unique(restricted[duplicated(restricted)]))}"
  )
  parameters <- unique(unlist(lapply(trees, .mpt_expr_vars)))
  restrictions <- .mpt_resolve_restrictions(restrictions)
  targets <- unique(unlist(lapply(restrictions, all.vars)))
  unknown <- setdiff(c(restricted, targets), parameters)
  stopif(
    "FE" %in% setdiff(targets, parameters),
    "FE is TreeBUGS syntax for a parameter without person effects. In bmm, \\
    leave the random effects out of that parameter's formula instead \\
    (e.g. g ~ 1 instead of g ~ 1 + (1 | id))."
  )
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
    notation, which breaks the Stan syntax). Please use a larger constant."
  )
  lapply(trees, .mpt_apply_restrictions, restrictions)
}

#' @export
print_model_details.mpt <- function(model, ...) {
  for (tree in model$other_vars$trees) {
    print(tree)
  }
  if (length(model$other_vars$restrictions) > 0) {
    cat(
      "Restrictions:",
      paste(
        names(model$other_vars$restrictions),
        vapply(model$other_vars$restrictions, deparse1, character(1)),
        sep = " = ", collapse = ", "
      ),
      "\n"
    )
  }
  # the links field lists the sticks as identity because their generated
  # formulas apply the model link; the display says which one
  for (grp in model$other_vars$simplex) {
    cat(glue(
      "Simplex:    {collapse_comma(grp)} via stick-breaking; the components \\
      {collapse_comma(model$other_vars$simplex_raw[grp[-length(grp)]])} pass \\
      through {inv_link('x', model$other_vars$link)[[1]]}()"
    ), "\n")
  }
  # the classical parameters-versus-categories bound, followed by the Jacobian
  # rank, which also catches redundant parameters when the count passes; each
  # tree contributes categories minus one. A simplex group counts through its
  # stick-breaking components, one fewer than its members.
  n_free <- length(setdiff(
    names(attr(model, "links_default")) %||% names(model$parameters),
    c(names(model$fixed_parameters), unlist(model$other_vars$simplex))
  ))
  df <- sum(lengths(lapply(model$other_vars$trees, `[[`, "branches")) - 1L)
  cat(glue(
    "Identifiability (intercept-only formulas): {n_free} free parameter(s), \\
    {df} degrees of freedom (response categories minus 1, summed over trees)"
  ), "\n")
  if (n_free > df) {
    cat(
      "  More free parameters than degrees of freedom: the model is not",
      "identified without further constraints.\n"
    )
  }
  identifiability <- .mpt_identifiability(model)
  rank_text <- if (!is.null(identifiability$error)) {
    glue(
      "Jacobian rank not computed: stats::D() cannot differentiate the branch \\
      expressions ({identifiability$error})."
    )
  } else if (identifiability$rank < identifiability$n_free) {
    .mpt_rank_deficit_text(identifiability)
  } else {
    glue(
      "Jacobian rank {identifiability$rank} of {identifiability$n_free} at \\
      interior test values: locally identified."
    )
  }
  cat(strwrap(rank_text, indent = 2, exdent = 4), sep = "\n")
  invisible(NULL)
}

# rank of the category probabilities in the free tree parameters; a parameter
# fixed in the formula enters at its value. A simplex group is free through its
# stick-breaking components, which users never write, so the text names the
# group's members instead
.mpt_identifiability <- function(model) {
  trees <- model$other_vars$trees
  simplex <- model$other_vars$simplex
  sticks <- model$other_vars$simplex_raw
  parameters <- setdiff(
    unique(unlist(lapply(trees, .mpt_expr_vars))), unlist(simplex)
  )
  fixed <- model$fixed_parameters[
    intersect(names(model$fixed_parameters), parameters)
  ]
  free <- c(
    setdiff(parameters, names(fixed)),
    setdiff(unname(sticks), names(model$fixed_parameters))
  )
  labels <- setNames(paste0("'", free, "'"), free)
  for (grp in simplex) {
    labels[intersect(free, sticks[grp])] <- glue(
      "the simplex group {collapse_comma(grp)}"
    )
  }
  c(.mpt_jacobian_rank(trees, free, fixed, simplex, sticks), list(labels = labels))
}

.mpt_rank_deficit_text <- function(identifiability) {
  free <- identifiability$free
  labels <- identifiability$labels
  absent <- identifiability$absent
  entangled <- setdiff(identifiability$involved, absent)
  n_combinations <- identifiability$n_free - identifiability$rank -
    length(absent)
  paste(c(
    glue(
      "The model is not identified: the Jacobian of the category \\
      probabilities has rank {identifiability$rank} for \\
      {identifiability$n_free} free parameters at interior test values."
    ),
    if (n_combinations > 0 && length(entangled) > 0) {
      glue(
        "{n_combinations} combination(s) of \\
        {.mpt_parameter_set(entangled, free, labels)} cannot be estimated \\
        from the data."
      )
    },
    if (length(absent) > 0) {
      glue(
        "The derivative with respect to {.mpt_labels(absent, labels)} is zero up \\
        to rounding at the test values, so these parameter(s) appear not to \\
        affect any category probability."
      )
    },
    glue(
      "Fix parameters in the formula (bmf(name = value)) or equate them in \\
      the branch expressions."
    )
  ), collapse = " ")
}

# a list longer than half the free parameters reads better as its complement
.mpt_parameter_set <- function(pars, free, labels) {
  if (length(pars) <= length(free) / 2) {
    return(.mpt_labels(pars, labels))
  }
  rest <- setdiff(free, pars)
  if (length(rest) == 0L) {
    return("all free parameters")
  }
  glue("all free parameters except {.mpt_labels(rest, labels)}")
}

.mpt_labels <- function(pars, labels) {
  paste(unique(labels[pars]), collapse = ", ")
}

############################################################################# !
# CHECK_MODEL S3 methods                                                 ####
############################################################################# !

#' @export
check_model.mpt <- function(model, data = NULL, formula = NULL) {
  model <- .mpt_undo_formula_state(model)
  model <- .mpt_match_link_priors(model)
  if (!is.null(formula)) {
    resp_cats <- model$resp_vars$resp_cats
    user_cat_formulas <- intersect(resp_cats, names(formula))
    stopif(
      length(user_cat_formulas) > 0,
      "The response category probabilities are fully determined by the tree \\
      branch expressions and cannot be predicted directly. Please remove the \\
      formula(s) for: {collapse_comma(user_cat_formulas)}"
    )
    nl_pars <- intersect(names(formula)[is_nl(formula)], names(model$parameters))
    nl_simplex <- intersect(nl_pars, unlist(model$other_vars$simplex))
    stopif(
      length(nl_simplex) > 0,
      "Non-linear predictor formulas are not supported for simplex parameters: \\
      {collapse_comma(nl_simplex)}"
    )
    sub_pars <- .mpt_nl_subparameters(model, formula, data)
    no_formula <- setdiff(sub_pars, names(formula))
    stopif(
      length(no_formula) > 0,
      "{collapse_comma(no_formula)} in your non-linear formula(s) is neither a \\
      data column nor a model parameter. Give each new parameter its own \\
      formula (e.g. {no_formula[1]} ~ 1), or add the column to the data."
    )
    .mpt_check_names(sub_pars, "parameter")
    if (length(sub_pars) > 0) {
      message2(
        "The parameter(s) {collapse_comma(sub_pars)} from your non-linear \\
        formulas are estimated on the identity scale with normal(0, 1) default \\
        priors, and their random-effect SDs keep brms's student_t(3, 0, 2.5) \\
        default. Apply any required transformation inside your formula and \\
        adjust the priors to the scale of your predictors."
      )
    }
    model <- .mpt_bypass_links(model, nl_pars)
    model <- .mpt_add_subparameters(model, sub_pars)
    # these links are the pipeline's own doing, so a later check_links() must
    # not read them as user changes
    attr(model, "links_checked") <- model$links

    # fixed values from the bmmformula (e.g. g = 0.5) stay on the probability
    # scale in the model object; configure_prior.mpt() maps them to the
    # latent scale when it builds the constant priors
    fixed_simplex <- intersect(
      names(model$fixed_parameters),
      c(unlist(model$other_vars$simplex), model$other_vars$simplex_raw)
    )
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

    # population-level predictors can identify a parameter across design
    # cells that a single cell leaves open, so the per-cell rank is the rank
    # of the fitted model only when no formula has a data predictor; with
    # predictors, a per-cell deficit is announced but not warned about
    with_predictors <- .mpt_predictor_formulas(formula, names(model$parameters))
    identifiability <- .mpt_identifiability(model)
    if (is.null(identifiability$error) &&
          identifiability$rank < identifiability$n_free) {
      if (length(with_predictors) == 0L) {
        warning2(
          "{.mpt_rank_deficit_text(identifiability)} Along the non-identified \\
          direction(s), the posterior follows the prior."
        )
      } else {
        message2(
          "The model is not identified within one design cell (Jacobian \\
          rank {identifiability$rank} for {identifiability$n_free} free \\
          parameters; print(model) names the parameters involved). Whether \\
          the predictors on {collapse_comma(with_predictors)} identify it \\
          across cells is not checked."
        )
      }
    }
  }
  NextMethod("check_model")
}

# formulas whose population-level terms use a data column; random-effect terms
# carry a bar, and parameters inside a non-linear formula are not data. A
# formula terms() cannot read counts as having predictors, which turns the
# warning into the milder message rather than risking a false one
.mpt_predictor_formulas <- function(formula, parameters) {
  formula <- formula[!is_constant(formula)]
  has_predictors <- vapply(formula, function(par_formula) {
    formula_terms <- try(stats::terms(par_formula), silent = TRUE)
    if (is_try_error(formula_terms)) {
      return(TRUE)
    }
    labels <- attr(formula_terms, "term.labels")
    population_vars <- unlist(lapply(
      labels[!grepl("|", labels, fixed = TRUE)],
      function(label) all.vars(str2lang(label))
    ))
    length(setdiff(population_vars, c(names(formula), parameters))) > 0
  }, logical(1))
  names(formula)[has_predictors]
}

# a parameter whose link was switched after construction gets the matched
# prior of its new link, unless the user replaced its default prior
.mpt_match_link_priors <- function(model) {
  model_prior <- .mpt_latent_prior(model$other_vars$link)
  for (par in settable_links(model)) {
    link <- model$links[[par]]
    if (identical(link, model$other_vars$link)) next
    if (!identical(model$default_priors[[par]], model_prior)) next
    model$default_priors[[par]] <- .mpt_latent_prior(link)
  }
  model
}

# symbols of non-linear formulas that are neither parameters, formula
# parameters nor data columns are sub-parameters the user introduces
.mpt_nl_subparameters <- function(model, formula, data) {
  nl_formulas <- formula[is_nl(formula)]
  setdiff(
    rhs_vars(nl_formulas),
    c(
      names(nl_formulas), names(model$parameters), colnames(data),
      model$other_vars$tree_id
    )
  )
}

# a model checked before (fit$bmm$model, re-checked by update()) carries the
# bypassed links and the sub-parameters of the formula it was checked with;
# the new formula may make a parameter linear again or drop a sub-parameter
.mpt_undo_formula_state <- function(model) {
  sub_pars <- setdiff(names(model$parameters), names(attr(model, "links_default")))
  model$parameters[sub_pars] <- NULL
  model$links[sub_pars] <- NULL
  model$default_priors[sub_pars] <- NULL
  # back to the model's latent prior, which check_model.mpt() then matches to
  # a switched link
  bypassed <- attr(model, "mpt_bypassed_links")
  model$links[names(bypassed)] <- bypassed
  model$default_priors[names(bypassed)] <- list(.mpt_latent_prior(model$other_vars$link))
  attr(model, "mpt_bypassed_links") <- NULL
  attr(model, "links_checked") <- model$links
  model
}

# parameters with a user-supplied non-linear formula own their (0,1)
# constraint, so the automatic link transformation is switched off for them;
# their links are kept for a later check with another formula
.mpt_bypass_links <- function(model, pars) {
  attr(model, "mpt_bypassed_links") <- model$links[pars]
  for (par in pars) {
    model$links[[par]] <- "identity"
    model$default_priors[[par]] <- NULL
  }
  model
}

# no sd default: a sub-parameter lives on whatever scale the user's formula
# gives it, so no rate fits all of them and brms's default applies
.mpt_add_subparameters <- function(model, sub_pars) {
  if (length(sub_pars) == 0L) {
    return(model)
  }
  model$parameters[sub_pars] <- .mpt_named_list(
    sub_pars, "User-defined sub-parameter of a non-linear parameter formula."
  )
  model$links[sub_pars] <- .mpt_named_list(sub_pars, "identity")
  model$default_priors[sub_pars] <- .mpt_named_list(
    sub_pars, list(main = "normal(0, 1)", effects = "normal(0, 0.5)")
  )
  model
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

  reserved_cols <- intersect(c("Y", "nTrials"), col_names)
  stopif(
    length(reserved_cols) > 0,
    "The data column(s) {collapse_comma(reserved_cols)} would be overwritten by \\
    the response matrix and trial counts that bmm builds. Please rename them."
  )

  resp_matrix <- as.matrix(data[resp_cats])
  stopif(
    !is.numeric(resp_matrix) || any(resp_matrix < 0, na.rm = TRUE),
    "The response category columns must contain non-negative response counts."
  )
  warnif(
    anyNA(resp_matrix),
    "The response count columns contain {sum(is.na(resp_matrix))} missing \\
    value(s), which are counted as 0 responses."
  )
  resp_matrix[is.na(resp_matrix)] <- 0
  data <- data[!col_names %in% resp_cats]
  # an integer column keeps conditional_effects() grids at a valid trial count
  data$nTrials <- as.integer(rowSums(resp_matrix))
  data$Y <- resp_matrix

  tree_id <- model$other_vars$tree_id
  if (!is.null(tree_id)) {
    stopif(
      !tree_id %in% colnames(data),
      "The tree identifier column '{tree_id}' is not present in the data."
    )
    tree_values <- as.character(data[[tree_id]])
    tree_names <- names(model$other_vars$trees)
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
    idx_vars <- model$other_vars$indicators$tree
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

  raw_collisions <- intersect(model$other_vars$simplex_raw, colnames(data))
  stopif(
    length(raw_collisions) > 0,
    "The data contain column(s) {collapse_comma(raw_collisions)}, which are \\
    the names of the stick-breaking components of the simplex parameters. \\
    Please rename them."
  )

  NextMethod("check_data")
}

############################################################################# !
# CHECK_Formula S3 methods                                               ####
############################################################################# !

#' @export
check_formula.mpt <- function(model, data, formula) {
  for (grp in model$other_vars$simplex) {
    free_pars <- grp[-length(grp)]
    derived_par <- grp[length(grp)]
    conflicts <- free_pars[vapply(free_pars, function(par) {
      .mpt_formulas_conflict(formula[[par]], formula[[model$other_vars$simplex_raw[[par]]]])
    }, logical(1))]
    stopif(
      length(conflicts) > 0,
      "Conflicting predictor formulas for the simplex parameter(s) \\
      {collapse_comma(conflicts)} and the stick-breaking component(s) \\
      {collapse_comma(model$other_vars$simplex_raw[conflicts])}. Specify predictors for \\
      the simplex parameters only."
    )
    stopif(
      !.mpt_intercept_only(formula[[derived_par]]),
      "The parameter '{derived_par}' is derived as 1 minus the sum of \\
      {collapse_comma(free_pars)} and cannot have its own predictors. Specify \\
      predictors for the other parameters of the simplex group instead."
    )
  }
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
# receive generated stick-breaking formulas. check_formula.bmmodel() has added
# `par ~ 1` for every parameter by the time these run, so a user's explicit
# `par ~ 1` counts as "no predictors" and the raw component's formula stands.
.mpt_move_simplex_formulas <- function(model, formula) {
  for (grp in model$other_vars$simplex) {
    for (par in grp[-length(grp)]) {
      if (.mpt_intercept_only(formula[[par]])) next
      raw_par <- model$other_vars$simplex_raw[[par]]
      formula[raw_par] <- list(.mpt_formula(raw_par, formula[[par]][[3]]))
    }
    formula <- formula[setdiff(names(formula), grp)]
  }
  formula
}

.mpt_intercept_only <- function(pform) {
  identical(pform[[3]], 1)
}

.mpt_formulas_conflict <- function(par_form, raw_form) {
  !.mpt_intercept_only(par_form) && !.mpt_intercept_only(raw_form) &&
    !identical(par_form[[3]], raw_form[[3]])
}

.mpt_category_formulas <- function(model) {
  trees <- model$other_vars$trees
  idx_vars <- model$other_vars$indicators$tree
  category_formulas <- lapply(model$resp_vars$resp_cats, function(resp_cat) {
    terms <- lapply(names(trees), function(tree_name) {
      branch <- call("(", trees[[tree_name]]$branches[[resp_cat]])
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
  if (length(model$other_vars$simplex) == 0L) {
    return(NULL)
  }
  simplex_formulas <- list()
  for (grp in model$other_vars$simplex) {
    n_grp <- length(grp)
    free_pars <- grp[-n_grp]
    sticks <- lapply(model$other_vars$simplex_raw[free_pars], inv_link, link = model$other_vars$link)
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
  linpreds <- glue("log({resp_cats})")

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

# no expected_response_defined() method: brms's multinomial posterior_epred() gives
# the expected count per category, trials times the tree's branch probabilities
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
