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
#'   values of the tree identifier column in the data.
#' @param categories Character vector with the response category of each line
#'   within a tree, used for lines without an inline `# category` comment.
#'   Must have as many entries as the tree with the most lines.
#' @param tree_id Character. Name of the data column whose values name the tree
#'   each observation belongs to (see [mpt()]).
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
#'   tree_id = "item_type"
#' )
#' model
#' @export
mpt_from_string <- function(text, tree_names, categories = NULL,
                            tree_id = NULL, covariates = NULL,
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

  model <- mpt(
    trees = trees, tree_id = tree_id, covariates = covariates,
    simplex = simplex, links = links
  )
  attr(model, "call") <- match.call()
  model
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
#' @param tree_id Character. Name of the data column whose values name the tree
#'   each observation belongs to (see [mpt()]); the tree names come from the
#'   EQN file.
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
#'   tree_id = "item_type"
#' )
#' attr(model, "mpt_renaming")
#' @export
mpt_from_eqn <- function(file, restrictions = NULL, categories = NULL,
                         tree_id = NULL, covariates = NULL, simplex = NULL,
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
  # longest names first: a dot creates a word boundary, so the pattern for
  # 'd_A' would otherwise also match the prefix of a parameter named 'd_A.x'
  for (old in names(par_renaming)[order(-nchar(names(par_renaming)))]) {
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
    trees = trees, tree_id = tree_id, covariates = covariates,
    simplex = simplex, links = links
  )
  attr(model, "mpt_renaming") <- renaming
  attr(model, "call") <- match.call()
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
  for (par in names(restrictions)[order(-nchar(names(restrictions)))]) {
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
