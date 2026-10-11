#' Reshape wide MPT frequency tables to the long format of `bmm()`
#'
#' MPTinR and TreeBUGS take one row per participant and one column per
#' response category of the EQN file. [bmm()] takes one row per participant
#' and tree, a column naming the tree, and one count column per *shared*
#' response category (see [mpt()]). `mpt_long_data()` reshapes the first
#' layout into the second for a given model.
#'
#' @param wide A data frame or matrix with one row per participant. The count
#'   columns are read as described under `columns`; with `columns = "names"`,
#'   all other columns except `id` are carried to the long data.
#' @param model An `mpt` model from [mpt()], [mpt_from_eqn()] or
#'   [mpt_from_string()]. For a model from [mpt_from_eqn()] the `mpt_renaming`
#'   attribute maps the labels of the EQN file onto the response categories
#'   of the model. For any other model, the columns of `wide` are named
#'   `<tree>.<category>`, e.g. `old.hit`, with the tree names and response
#'   categories of the model.
#' @param id Character. The name of the participant column in `wide`, or
#'   `NULL` (the default) to number the rows of `wide`, as `id = seq_len(nrow(wide))`
#'   in a column named `id`. The participants must be unique.
#' @param columns Character. How the count columns of `wide` are identified:
#'   `"names"` (the default) by their names, as in TreeBUGS data; `"mptinr"`
#'   by their position, as in MPTinR data. See Details.
#'
#' @details
#'   **Named columns (`columns = "names"`).** For a model from
#'   [mpt_from_eqn()], a column is matched to a category by the EQN label
#'   (`hit`), or, when trees share labels, by the tree name and label
#'   (`old.hit`, with the bmm name of the tree); the qualified name wins
#'   where both exist. A label used in several trees must be qualified.
#'   Labels that the `categories` argument of [mpt_from_eqn()] merged onto one
#'   response category are summed. All columns of `wide` that are neither
#'   count columns nor `id` are carried unchanged to every long row, e.g.
#'   between-participant covariates.
#'
#'   **Positional columns (`columns = "mptinr"`).** MPTinR data for an EQN
#'   model have no column names. Every column of `wide` except `id` is read as
#'   a count, in the order of `check.mpt()` in MPTinR (version 1.14.1),
#'   elements `eqn.order.trees` and `eqn.order.categories`: trees sorted by
#'   their label in the EQN file and, within each tree, the category labels of
#'   the file sorted. No other column can be carried; add covariates
#'   afterwards, e.g. with [merge()] on the participant column. The help page
#'   of `check.mpt()` states only that these elements give the order; the
#'   sorting rule is that of its internal `.get.EQN.model.order()`, which
#'   sorts the columns that `read.table()` returns for the file: numerically
#'   if all tree labels (all category labels) are numbers, else as character
#'   strings in the collation of the R session. `mpt_long_data()` applies the
#'   same rule, so run it in the session where MPTinR's order was determined
#'   if the labels mix upper and lower case. This mode needs a model from
#'   [mpt_from_eqn()]: MPTinR's own model format orders the columns by the
#'   lines of the model text, and [mpt()] reorders the branches of every tree
#'   to the categories of the first, so the line order cannot be recovered
#'   from the model. For other models, name the columns `<tree>.<category>`
#'   and use `columns = "names"`.
#'
#'   **Impossible categories.** A response category that a tree cannot
#'   produce (`impossible` in [mpt_tree()], [mpt_from_string()] and
#'   [mpt_from_eqn()]) has no column in `wide` and gets 0 in the long rows of
#'   that tree, which [bmm()] accepts.
#'
#'   **Counts.** The count columns must be numeric, without missing values,
#'   and hold non-negative whole numbers; the long data hold them as integers.
#'
#'   **Checks.** The function cannot tell a correct reshape from one that
#'   swaps columns of equal meaning. Run [bmm_data_check()] on the result, and
#'   compare one participant's long rows with the wide row (see the Examples
#'   of [mpt_from_eqn()]).
#'
#' @return A data frame with one row per participant and tree, participants
#'   first. Its columns are the participant column (`id`), the tree column
#'   named by the `tree_id` of the model (absent in a model without it), the
#'   carried columns, and one integer count column per response category of
#'   the model.
#'
#' @seealso [mpt_from_eqn()] and [mpt_from_string()] for importing the model
#'   itself; [mpt()] for the data layout [bmm()] expects; [bmm_data_check()]
#'   to check the result before fitting.
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
#'   restrictions = "D_n = D_o",
#'   categories = c(hit = "yes", fa = "yes", miss = "no", cr = "no"),
#'   tree_id = "item_type"
#' )
#'
#' # TreeBUGS data: one named column per category label of the EQN file.
#' # Columns that are no counts, here a between-participant factor, are
#' # carried to the long data
#' wide <- data.frame(
#'   subject = c("s1", "s2"), group = c("a", "b"),
#'   hit = c(30, 25), miss = c(10, 15), fa = c(8, 12), cr = c(32, 28)
#' )
#' long <- mpt_long_data(wide, model, id = "subject")
#' long
#' bmm_data_check(bmf(Do ~ 1, g ~ 1), long, model)
#'
#' # MPTinR data: unnamed columns, sorted by tree, then by label (cr, fa for
#' # tree new, then hit, miss for tree old)
#' wide_mptinr <- rbind(c(32, 8, 30, 10), c(28, 12, 25, 15))
#' identical(
#'   mpt_long_data(wide_mptinr, model, columns = "mptinr"),
#'   mpt_long_data(wide[c("cr", "fa", "hit", "miss")], model)
#' )
#'
#' # trees with response categories that only one tree can produce
#' model_pc <- mpt(
#'   list(
#'     mpt_tree("pairs", list(E1 = "c * r", E2 = "1 - c * r"),
#'              impossible = "F1"),
#'     mpt_tree("singles", list(F1 = "u", E2 = "1 - u"),
#'              impossible = "E1")
#'   ),
#'   tree_id = "item_type"
#' )
#' wide_pc <- data.frame(
#'   pairs.E1 = c(20, 22), pairs.E2 = c(20, 18),
#'   singles.F1 = c(25, 24), singles.E2 = c(15, 16)
#' )
#' mpt_long_data(wide_pc, model_pc)
#' @export
mpt_long_data <- function(wide, model, id = NULL,
                          columns = c("names", "mptinr")) {
  stop_missing_args()
  columns <- match.arg(columns)
  stopif(
    !inherits(model, "mpt"),
    "The model must be an mpt model, created with mpt(), mpt_from_eqn() or \\
    mpt_from_string()."
  )
  stopif(
    !is.data.frame(wide) && !is.matrix(wide),
    "The wide data must be a data frame or a matrix."
  )
  stopif(
    !is.null(id) && !(is.character(id) && length(id) == 1 && !is.na(id)),
    "The id argument must be NULL or the name of a single column of the wide \\
    data."
  )
  wide <- as.data.frame(wide, stringsAsFactors = FALSE)
  stopif(nrow(wide) == 0, "The wide data contain no rows.")
  stopif(
    !is.null(id) && !id %in% names(wide),
    "The id column '{id}' is not in the wide data. Columns: \\
    {collapse_comma(names(wide))}"
  )

  participants <- if (is.null(id)) seq_len(nrow(wide)) else wide[[id]]
  stopif(
    anyNA(participants) || anyDuplicated(participants) > 0,
    "The participant column '{id %||% 'id'}' must identify each row of the \\
    wide data: it has missing or repeated values."
  )
  id_name <- id %||% "id"

  layout <- .mpt_long_layout(model, columns)
  non_id <- setdiff(names(wide), id)
  count_cols <- if (columns == "names") {
    .mpt_named_count_columns(layout, non_id)
  } else {
    .mpt_positional_count_columns(layout, non_id)
  }
  counts <- wide[count_cols]
  .mpt_check_wide_counts(counts)

  carried <- wide[setdiff(non_id, count_cols)]
  tree_id <- model$other_vars$tree_id
  resp_cats <- model$resp_vars$resp_cats
  taken <- c(id_name, tree_id, resp_cats)
  collisions <- intersect(names(carried), taken)
  stopif(
    length(collisions) > 0,
    "The wide data hold column(s) {collapse_comma(collisions)} that the long \\
    data need for the participant, the tree or a response category. Please \\
    rename them. A participant column is passed as id = 'name'."
  )

  tree_names <- names(model$other_vars$trees)
  blocks <- lapply(tree_names, function(tree) {
    in_tree <- layout$tree == tree
    tree_counts <- matrix(0, nrow(wide), length(resp_cats),
                          dimnames = list(NULL, resp_cats))
    summed <- t(rowsum(
      t(as.matrix(counts[in_tree])), layout$category[in_tree], reorder = FALSE
    ))
    tree_counts[, colnames(summed)] <- summed
    storage.mode(tree_counts) <- "integer"
    data.frame(
      c(
        setNames(list(participants), id_name),
        if (!is.null(tree_id)) setNames(list(rep(tree, nrow(wide))), tree_id),
        as.list(carried),
        as.list(as.data.frame(tree_counts))
      ),
      check.names = FALSE
    )
  })
  long <- do.call(rbind, blocks)
  long <- long[order(
    rep(seq_len(nrow(wide)), times = length(tree_names)),
    rep(seq_along(tree_names), each = nrow(wide))
  ), ]
  rownames(long) <- NULL
  long
}

# one row per count column of the wide data: the tree, the response category
# of the model, and for EQN models the label of the file. For columns =
# "mptinr" the rows are in the order of MPTinR's columns.
.mpt_long_layout <- function(model, columns) {
  renaming <- attr(model, "mpt_renaming")
  if (is.null(renaming)) {
    stopif(
      columns == "mptinr",
      "columns = 'mptinr' needs a model from mpt_from_eqn(), because the \\
      column order of MPTinR follows the labels of the EQN file. For other \\
      models, name the columns '<tree>.<category>' and use columns = 'names'."
    )
    trees <- model$other_vars$trees
    return(data.frame(
      tree = rep(names(trees), lengths(lapply(trees, `[[`, "branches"))),
      category = unlist(lapply(trees, function(tree) names(tree$branches)),
                        use.names = FALSE),
      label = NA_character_
    ))
  }
  rows <- renaming[renaming$kind == "category", ]
  layout <- data.frame(
    tree = rows$tree, category = rows$bmm_name, label = rows$file_name
  )
  if (columns == "mptinr") {
    layout <- layout[.mpt_eqn_column_order(layout, renaming), ]
  }
  layout
}

# MPTinR sorts the tree labels and, within a tree, the category labels, each
# after read.table() has typed the file's column: numbers sort numerically,
# anything else as character strings
.mpt_eqn_column_order <- function(layout, renaming) {
  trees <- renaming[renaming$kind == "tree", ]
  tree_key <- utils::type.convert(trees$file_name, as.is = TRUE)
  order(
    tree_key[match(layout$tree, trees$bmm_name)],
    utils::type.convert(layout$label, as.is = TRUE)
  )
}

.mpt_named_count_columns <- function(layout, available) {
  qualified <- paste(layout$tree, ifelse(
    is.na(layout$label), layout$category, layout$label
  ), sep = ".")
  shared_label <- !is.na(layout$label) & !duplicated(layout$label) &
    !duplicated(layout$label, fromLast = TRUE)
  count_cols <- ifelse(
    qualified %in% available, qualified,
    ifelse(shared_label & layout$label %in% available, layout$label, NA)
  )
  missing_cols <- is.na(count_cols)
  stopif(
    any(missing_cols),
    "The wide data lack count columns for the following response categories \\
    (expected column names in parentheses): \\
    {paste0(layout$category[missing_cols], ' in tree ', layout$tree[missing_cols],
            ' (', qualified[missing_cols],
            ifelse(shared_label[missing_cols],
                   paste0(' or ', layout$label[missing_cols]), ''), ')',
            collapse = ', ')}"
  )
  count_cols
}

.mpt_positional_count_columns <- function(layout, available) {
  stopif(
    length(available) != nrow(layout),
    "With columns = 'mptinr', every column of the wide data except the \\
    participant column is read as a count, in the order of MPTinR. The model \\
    expects {nrow(layout)} columns ({collapse_comma(layout$category)} in the \\
    trees {collapse_comma(unique(layout$tree))}), but the wide data have \\
    {length(available)}."
  )
  available
}

.mpt_check_wide_counts <- function(counts) {
  non_numeric <- names(counts)[!vapply(counts, is.numeric, logical(1))]
  stopif(
    length(non_numeric) > 0,
    "The count columns must be numeric. Not numeric: \\
    {collapse_comma(non_numeric)}"
  )
  values <- as.matrix(counts)
  affected <- function(problem) collapse_comma(colnames(values)[colSums(problem) > 0])
  stopif(
    anyNA(values),
    "The count columns must not contain missing values. Missing values in: \\
    {affected(is.na(values))}"
  )
  not_counts <- values < 0 | values != round(values) | is.infinite(values)
  stopif(
    any(not_counts),
    "The count columns must hold non-negative whole numbers. Other values in: \\
    {affected(not_counts)}"
  )
  invisible(NULL)
}
