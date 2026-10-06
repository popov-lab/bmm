############################################################################# !
# IMPORT WRAPPERS                                                        ####
############################################################################# !

#' @title Create an MPT model from an MPTinR-style model string
#'
#' @description Reads a model definition in the "easy" format of MPTinR
#'   model files and returns the corresponding [mpt()] model object; EQN
#'   files (multiTree, MPTinR, TreeBUGS) are read by [mpt_from_eqn()]. Each
#'   line is one branch probability expression; blank lines separate trees.
#'   The response category of each line is given either by an inline comment
#'   (`expression # category`) or by the `categories` argument. Lines that
#'   terminate in the same response category are summed.
#'
#' @param text Character. The model definition, either as a single string or
#'   as a vector of lines (e.g., from [readLines()]).
#' @param tree_names Character vector with one name per tree block. The names
#'   are matched to the blocks by position: the first name goes to the first
#'   block, and so on. bmm cannot check this order, so names given in a
#'   different order than the blocks import a different model without an
#'   error. For multi-tree models the names must match the values of the tree
#'   identifier column in the data.
#' @param categories The response categories of lines without an inline
#'   `# category` comment, matched to the lines of a tree by position: the
#'   i-th entry is the category of the i-th line. Either a list by tree,
#'   `list(<tree> = c("<category>", ...))`, with one vector per tree that
#'   has such lines, or one unnamed character vector shared by all trees.
#'   A tree's vector in the list must not be longer than the tree. A shared
#'   vector must not be longer than the longest tree; a shorter tree takes
#'   its leading entries, which suits trees that lack the last categories
#'   (see `impossible`). If the trees use different categories, as in the
#'   pair-clustering model, use the list. Either form needs an entry for
#'   every line without an inline comment; an inline `# category` comment
#'   takes precedence over the entry for its line.
#' @param tree_id Character. Name of the data column whose values name the tree
#'   each observation belongs to (see [mpt()]). Add this column to the data
#'   yourself: for each row it holds the name, as given in `tree_names`, of
#'   the tree the counts come from.
#' @param covariates,simplex,restrictions,links passed to [mpt()].
#'   `print(model)` lists the symbols that became parameters.
#' @param impossible A named list with the response categories that cannot
#'   occur in a tree, one element per such tree:
#'   `list(<tree> = c("<category>", ...))` (see the `impossible` argument of
#'   [mpt_tree()]). Trees and categories are named as in the returned model,
#'   i.e. by `tree_names` and by the category labels of the lines. Trees that
#'   are not listed have no impossible categories.
#'
#' @details **Comments.** A line whose first non-blank character is `#` is a
#'   comment and is skipped. It does not separate trees; only blank lines do.
#'   On a line with an expression, the text after `#` is the response category
#'   of that line. This differs from MPTinR, which ignores everything after a
#'   `#` and treats a line that starts with `#` as a blank line, so that there
#'   comment lines can separate trees. An MPTinR model file whose trees are
#'   separated only by comment lines, or whose inline comments are not
#'   category labels, has to be edited before it is imported.
#'
#'   **Trees with different response categories.** In bmm all trees share the
#'   same response categories. If a category cannot occur in some trees (e.g.,
#'   the pair and singleton trees of the pair-clustering model), declare it
#'   with the `impossible` argument. If the trees only label the same
#'   responses differently, give those lines the same category label.
#'
#'   **Names.** Parameter names are used as written in the model text and are
#'   case-sensitive, so the formulas passed to [bmm()] name them the same way
#'   (`bmf(D ~ 1)`; `bmf()` is short for [bmmformula()]). They must be valid
#'   bmm names (letters and digits, see [mpt()]); unlike [mpt_from_eqn()],
#'   this function does not rename them. The import fails when the branches
#'   of a tree do not sum to 1 at numeric test values.
#'
#'   **Data for an imported model.** The Details of [mpt_from_eqn()] describe
#'   the layout [bmm()] expects (one row per participant and tree, a tree
#'   column, one count column per response category, 0 for impossible
#'   categories) and the checks to run before fitting; they apply to models
#'   from this function as well, with the tree names of `tree_names` and the
#'   category labels of the lines. The data MPTinR fits to a model file in
#'   this format are wide and positional: one row per participant (or one
#'   vector for aggregated counts) and one unnamed column per line of the
#'   model file, in file order, because each line is one response category.
#'   (For EQN files MPTinR sorts the columns instead; see [mpt_from_eqn()].)
#'   Name the columns by tree and category first, then stack one block per
#'   tree (see Examples); when a tree lacks some categories, start each tree
#'   from zeros, as in the Examples of [mpt_from_eqn()].
#'
#'   **Provenance.** The returned model carries an `mpt_source` attribute,
#'   `list(text = )`, with the model definition as one string.
#'
#' @return An object of class `bmmodel` (see [mpt()]), with an `mpt_source`
#'   attribute holding the imported text.
#'
#' @seealso [mpt_from_eqn()] for EQN files (multiTree, MPTinR, TreeBUGS);
#'   [mpt()] and [mpt_tree()] for building the trees in R, including the data
#'   layout `bmm()` expects.
#'
#' @keywords transform
#'
#' @examples
#' model_2htm <- "
#' # two-high-threshold model
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
#'
#' # fix the guessing rate as in MPTinR / TreeBUGS
#' mpt_from_string(
#'   model_2htm,
#'   tree_names = c("old", "new"),
#'   tree_id = "item_type",
#'   restrictions = "g = 0.5"
#' )
#'
#' # pair clustering, written as in MPTinR without category labels: pairs and
#' # singletons have different response categories
#' model_pc <- "
#' c * r
#' (1 - c) * u * u
#' 2 * (1 - c) * u * (1 - u)
#' c * (1 - r) + (1 - c) * (1 - u) * (1 - u)
#'
#' u
#' 1 - u
#' "
#' mpt_from_string(
#'   model_pc,
#'   tree_names = c("pairs", "singles"),
#'   categories = list(
#'     pairs = c("E1", "E2", "E3", "E4"),
#'     singles = c("F1", "F2")
#'   ),
#'   tree_id = "item_type",
#'   impossible = list(pairs = c("F1", "F2"), singles = c("E1", "E2", "E3", "E4"))
#' )
#'
#' # an MPTinR model file read with readLines() keeps the blank lines that
#' # separate the trees. A restrictions file is not read by bmm: drop its
#' # comments and blank lines and pass the remaining lines.
#' model_file <- tempfile(fileext = ".model")
#' writeLines(model_2htm, model_file)
#' restr_file <- tempfile(fileext = ".restr")
#' writeLines(c("# fix the guessing rate", "g = 0.5"), restr_file)
#' restr <- trimws(sub("#.*", "", readLines(restr_file)))
#' model <- mpt_from_string(
#'   readLines(model_file),
#'   tree_names = c("old", "new"),
#'   tree_id = "item_type",
#'   restrictions = restr[nzchar(restr)]
#' )
#'
#' # MPTinR data: one unnamed column per line of the model file, in file
#' # order. Name the columns by tree and category, then stack one block per
#' # tree to get bmm's one row per participant and tree.
#' wide <- rbind(c(30, 10, 8, 32), c(25, 15, 12, 28))
#' colnames(wide) <- c("old.old", "old.new", "new.old", "new.new")
#' long <- do.call(rbind, lapply(c("old", "new"), function(tree) {
#'   counts <- wide[, paste(tree, c("old", "new"), sep = "."), drop = FALSE]
#'   colnames(counts) <- c("old", "new")
#'   data.frame(id = seq_len(nrow(wide)), item_type = tree, counts)
#' }))
#' long
#' # the row sums are the trials per row (40 here), but they cannot show
#' # swapped columns or tree labels: compare one participant with the wide
#' # row, then run the data checks (bmm_data_check() prints
#' # "Hard checks: passed" when the columns, tree labels and counts fit)
#' rowSums(long[c("old", "new")])
#' long[long$id == 1, ]
#' wide[1, ]
#' bmm_data_check(bmf(D ~ 1), long, model)
#' @export
mpt_from_string <- function(text, tree_names, categories = NULL,
                            tree_id = NULL, covariates = NULL,
                            simplex = NULL, restrictions = NULL,
                            links = "logit", impossible = NULL) {
  stop_missing_args()
  # a vector of lines keeps its "" elements, which separate the trees
  text <- paste(text, collapse = "\n")
  lines <- trimws(strsplit(text, "\n", fixed = TRUE)[[1]])
  lines <- lines[!startsWith(lines, "#")]
  blocks <- split(lines[nzchar(lines)], cumsum(!nzchar(lines))[nzchar(lines)])
  stopif(
    length(blocks) != length(tree_names),
    "Found {length(blocks)} tree block(s) separated by blank lines, but \\
    {length(tree_names)} tree_names were supplied."
  )

  tree_lines <- setNames(lapply(blocks, function(block) {
    parts <- strsplit(block, "#", fixed = TRUE)
    data.frame(
      expr = trimws(vapply(parts, `[[`, character(1), 1)),
      category = vapply(parts, function(part) {
        if (length(part) > 1) trimws(paste(part[-1], collapse = "#")) else ""
      }, character(1))
    )
  }), tree_names)
  tree_sizes <- vapply(tree_lines, nrow, integer(1))
  unlabelled <- lapply(tree_lines, function(block_lines) {
    which(!nzchar(block_lines$category))
  })

  by_tree <- is.list(categories)
  stopif(
    by_tree && !is_namedlist(categories),
    "A categories list must be a list by tree, named by tree_names, e.g. \\
    categories = list({tree_names[1]} = c('a', 'b'))."
  )
  stopif(
    !by_tree && !is.null(names(categories)),
    "The categories vector of mpt_from_string() labels the lines of a tree \\
    by position, so it takes no names. Give it as an unnamed vector shared by \\
    all trees, or as a list by tree, e.g. \\
    categories = list({tree_names[1]} = c('a', 'b'))."
  )
  stopif(
    by_tree && !all(names(categories) %in% tree_names),
    "The categories list names trees that are not in the model: \\
    {collapse_comma(setdiff(names(categories), tree_names))}. The trees are: \\
    {collapse_comma(tree_names)}."
  )
  too_long <- if (by_tree) {
    names(categories)[lengths(categories) > tree_sizes[names(categories)]]
  } else {
    character(0)
  }
  stopif(
    length(too_long) > 0,
    "The categories list gives more entries than a tree has lines: \\
    {paste0(too_long, ' (', tree_sizes[too_long], ' lines, ',
            lengths(categories[too_long]), ' entries)', collapse = ', ')}."
  )
  stopif(
    !by_tree && length(categories) > max(tree_sizes),
    "The categories vector has {length(categories)} entries, but the longest \\
    tree has {max(tree_sizes)} lines."
  )

  branches <- Map(function(block_lines, tree_name, missing_cats) {
    tree_categories <- if (by_tree) {
      categories[[tree_name]]
    } else {
      categories
    }
    stopif(
      length(tree_categories) < max(missing_cats, 0),
      "Tree '{tree_name}' has lines without an inline '# category' label at \\
      positions {paste(missing_cats, collapse = ', ')}. Entry i of categories \\
      labels line i, so the tree needs {max(missing_cats)} entries, but \\
      categories gives {length(tree_categories)}."
    )
    block_lines$category[missing_cats] <- tree_categories[missing_cats]
    .mpt_sum_branch_lines(block_lines$expr, block_lines$category)
  }, tree_lines, tree_names, unlabelled)

  # one vector over trees of different lengths gives the shorter trees the
  # leading labels of the longer ones, which is rarely meant
  unlabelled_counts <- lengths(unlabelled)[lengths(unlabelled) > 0]
  .mpt_check_import_trees(
    branches, impossible,
    relabel_hint = paste(
      "If the trees label the same responses differently (e.g., hit/miss and",
      "fa/cr for yes/no), give those lines the same category label."
    ),
    lead_hint = if (!by_tree && length(unique(unlabelled_counts)) > 1) {
      glue(
        "The trees have different numbers of lines without an inline label \\
        ({paste0(names(unlabelled_counts), ': ', unlabelled_counts,
                 collapse = ', ')}), \\
        but one categories vector labels them all by position. If the trees \\
        use different categories, give categories as a list by tree: \\
        categories = list({paste0(names(unlabelled_counts), ' = c(...)',
                                  collapse = ', ')})."
      )
    }
  )
  model <- .mpt_model_from_branches(
    branches, impossible, tree_id = tree_id, covariates = covariates,
    simplex = simplex, restrictions = restrictions, links = links
  )
  attr(model, "mpt_source") <- list(text = text)
  attr(model, "call") <- match.call()
  model
}


#' @title Create an MPT model from an EQN model file
#'
#' @description Reads an EQN model file, the format used by multiTree,
#'   MPTinR, and TreeBUGS, and returns the corresponding [mpt()] model object;
#'   model files in MPTinR's own "easy" format are read by
#'   [mpt_from_string()]. Each line of an EQN file contains three
#'   whitespace-separated fields -- tree, response category, branch
#'   probability expression -- and branch lines that terminate in the same
#'   response category of a tree are summed.
#'
#' @param file Character. Path to the EQN file. How its first line is read is
#'   described in Details.
#' @param restrictions Parameter restrictions as in [mpt()], written with the
#'   parameter names used in the EQN file, e.g. `c("D_n = D_o", "g = 0.5")`
#'   or `list(D_n = "D_o", g = 0.5)`. They mirror the `restrictions` argument
#'   of MPTinR and TreeBUGS and are renamed together with the parameters.
#'   Restriction files are not read: pass their lines (the Examples of
#'   [mpt_from_string()] read one with [readLines()]).
#' @param categories A named character vector mapping response category names
#'   used in the EQN file onto the shared response categories of the model
#'   (see Details), e.g. `c(hit = "yes", fa = "yes", miss = "no", cr = "no")`.
#'   Categories not listed keep their (sanitized) name; labels that only
#'   carry their tree name as a prefix (`WC_Correct`, `LC_Correct`) need no
#'   mapping, because the prefix is removed. Labels that are not valid bmm
#'   names, such as numeric labels, must be mapped. Mapping two categories of
#'   the same tree onto one name merges them: their branch lines are summed,
#'   and a message says so.
#' @param tree_id Character. Name of the data column whose values name the tree
#'   each observation belongs to (see [mpt()]). Add this column to the data
#'   yourself: for each row it holds the name of the tree the counts come
#'   from, as printed by the model (the file's label, or its replacement
#'   from `tree_names`).
#' @param covariates,simplex,links passed to [mpt()]. Covariate names are
#'   written as in the EQN file and are excluded from the renaming; simplex
#'   groups use the bmm names. `print(model)` lists the symbols that became
#'   parameters.
#' @param impossible A named list, `list(<tree> = c("<category>", ...))`, of
#'   the response categories that cannot occur in a tree (see [mpt_tree()]),
#'   written with the bmm names of the trees and categories, i.e. the
#'   `bmm_name` column of the renaming map; `print(model)` lists the shared
#'   categories under "Response". Trees not listed have no impossible
#'   categories. When the trees do not share their categories and
#'   `impossible` does not cover the difference, the error names the
#'   categories each tree lacks.
#' @param tree_names A named character vector mapping tree labels of the EQN
#'   file onto tree names in bmm, e.g. `c("1" = "old", "2" = "new")`. bmm tree
#'   names must start with a letter and contain only letters, digits, or
#'   underscores, so files with numeric tree labels need this mapping. Labels
#'   not listed keep their name.
#'
#' @details **File format.** Blank lines and lines that start with `#` are
#'   skipped. If the first remaining line holds a single integer and nothing
#'   else, it is read as the equation count of classical EQN files and
#'   skipped, with a warning when it differs from the number of equation
#'   lines that follow. Any other first line is read as an equation. Note
#'   that TreeBUGS always ignores the first line of an EQN file, whatever it
#'   holds, and so does MPTinR with `model.type = "eqn"` (its default for
#'   `.eqn` files; with `"eqn2"` it reads the first line as the equation
#'   count). bmm therefore reads a file without a count line in full, while
#'   those programs drop its first equation; a title line produces an error
#'   that says to delete it or start it with `#`. A file that starts with the
#'   equation count, or with a `#` line, is read the same way by all three.
#'   A line with fewer than three fields, or an expression that cannot be
#'   parsed, is an error that names the line; a tree whose branches do not
#'   sum to 1 at numeric test values is an error that names the tree.
#'
#'   **Response categories.** In bmm, every observation (row) belongs to one
#'   tree, and all trees share the same response categories. Classical EQN
#'   files often label the response categories per tree instead (e.g.,
#'   `hit`/`miss` in the old-item tree and `fa`/`cr` in the new-item tree of
#'   a recognition model). Use the `categories` argument to map such
#'   tree-specific labels onto the shared response categories (e.g., the
#'   response options "yes" and "no"). If a category cannot occur in some
#'   trees (e.g., the pair and singleton trees of the pair-clustering model),
#'   declare it with the `impossible` argument instead.
#'
#'   **Names.** Because brms does not allow underscores or dots in non-linear
#'   parameter names, parameter and response category names are sanitized: a
#'   leading `<tree>_` prefix is removed from category names that are not
#'   mapped with `categories`, and all remaining underscores and dots are
#'   stripped. Names that would clash after this are an error; names are
#'   case-sensitive, so `P_G1` and `P_g1` stay distinct (`PG1`, `Pg1`). Tree
#'   labels are kept, unless `tree_names` maps them. The formulas passed to
#'   [bmm()] use the bmm names (`bmf(Do ~ 1)` for a file parameter `D_o`;
#'   `bmf()` is short for [bmmformula()]); the `parameter` rows of the
#'   renaming map translate them, and `print(model)` lists them. Of the
#'   arguments that name parameters, only `restrictions` and `covariates`
#'   use the file's names.
#'
#'   **Renaming map.** The `mpt_renaming` attribute of the returned model is
#'   a data frame with one row per parameter, per tree, and per category of
#'   each tree, renamed or not, and the columns `kind` (`"parameter"`,
#'   `"tree"`, or `"category"`), `tree` (the bmm tree name for categories,
#'   `NA` otherwise), `file_name`, `bmm_name`, and `renamed` (`TRUE` where
#'   the two names differ). The names that changed are also reported in a
#'   message.
#'
#'   **Data for an imported model.** [bmm()] expects one row per participant
#'   and tree: a column named by `tree_id` holding the bmm tree names, and
#'   one count column per shared response category, named by `bmm_name` in
#'   the renaming map (see [mpt()]). A category that is `impossible` in a
#'   tree gets 0 (or `NA`) in the rows of that tree. TreeBUGS data are wide:
#'   one row per participant, one column per category label of the EQN file,
#'   and no participant column, so add an `id` column first. The `category`
#'   rows of the renaming map hold every tree with its labels and bmm names,
#'   so the wide data can be reshaped by name: for each tree, start from a
#'   block of zeros over all shared categories and fill in the columns of
#'   that tree (see Examples). MPTinR data for an EQN model are unnamed, and
#'   their columns do not follow the file: the trees and, within each tree,
#'   the category labels are sorted (alphabetically, or numerically when they
#'   are numbers); `MPTinR::check.mpt("model.eqn")` returns the order in its
#'   elements `eqn.order.trees` and `eqn.order.categories`. (For MPTinR's own
#'   model format the columns follow the lines instead; see
#'   [mpt_from_string()].) Name the columns by tree and label,
#'   `paste(tree, label, sep = ".")` with the bmm tree names, because trees
#'   may share labels, and select them with the same key (see Examples).
#'
#'   If `categories` merged several labels of a tree onto one name, sum their
#'   columns in the wide data into the first label's column and drop the map
#'   rows of the other labels (`map[!duplicated(map[c("tree", "bmm_name")]), ]`)
#'   before reshaping; otherwise the recipe fails with "duplicate subscripts".
#'
#'   bmm cannot tell a reshape that swaps tree labels or category columns
#'   from a correct one: such data fit without an error. The row sums,
#'   `rowSums(long[cats])`, are the trials per row; they cannot show swapped
#'   columns, and they show swapped tree labels only when the trees differ
#'   in their number of trials. Comparing the response proportions per tree
#'   with the design (old items must get more "yes" responses than new
#'   items) catches swapped tree labels and response columns swapped in
#'   every tree, but not always a column slipped or a map entry wrong in one
#'   tree, nor a swap between two trees of the same kind; so also compare
#'   one participant's long rows with the wide row. [bmm_data_check()] runs
#'   the remaining checks (columns, tree labels, counts in impossible cells)
#'   without a fit and prints `Hard checks: passed` when they hold.
#'
#'   **Provenance.** The returned model carries an `mpt_source` attribute,
#'   `list(file = , md5 = , lines = )`: the full path of the EQN file, its MD5
#'   checksum, and the lines as read, so that a saved fit records which file
#'   it was built from.
#'
#' @return An object of class `bmmodel` (see [mpt()]), with the attributes
#'   `mpt_renaming` and `mpt_source` (see Details).
#'
#' @seealso [mpt_from_string()] for model text in MPTinR's own format;
#'   [mpt()] and [mpt_tree()] for building the trees in R, including the data
#'   layout `bmm()` expects.
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
#' # restrictions use the file's names; the formulas use the bmm names
#' model <- mpt_from_eqn(
#'   eqn_file,
#'   restrictions = "D_n = D_o",
#'   categories = c(hit = "yes", fa = "yes", miss = "no", cr = "no"),
#'   tree_id = "item_type"
#' )
#' attr(model, "mpt_renaming")
#'
#' # TreeBUGS data: one row per participant, columns named by the category
#' # labels of the EQN file, no participant column. bmm needs one row per
#' # participant and tree, columns named by bmm_name: reshape by name through
#' # the map, starting each tree from zeros so that categories the tree
#' # cannot produce get 0
#' wide <- data.frame(
#'   hit = c(30, 25), miss = c(10, 15), fa = c(8, 12), cr = c(32, 28)
#' )
#' wide$id <- seq_len(nrow(wide))
#' map <- attr(model, "mpt_renaming")
#' map <- map[map$kind == "category", ]
#' cats <- unique(map$bmm_name)
#' long <- do.call(rbind, lapply(split(map, map$tree), function(tree_map) {
#'   counts <- as.data.frame(matrix(0, nrow(wide), length(cats)))
#'   names(counts) <- cats
#'   counts[tree_map$bmm_name] <- wide[tree_map$file_name]
#'   data.frame(id = wide$id, item_type = tree_map$tree[1], counts)
#' }))
#' long
#' # the row sums are the trials per row (40 here). The proportions catch
#' # swapped tree labels and swapped response columns (old items must get
#' # more "yes" than new items), but not always a column slipped in one
#' # tree: also compare one participant's long rows with the wide row
#' rowSums(long[cats])
#' tapply(long$yes / rowSums(long[cats]), long$item_type, mean)
#' long[long$id == 1, ]
#' wide[1, ]
#' bmm_data_check(bmf(Do ~ 1, g ~ 1), long, model)
#'
#' # MPTinR data for an EQN model are unnamed and sorted by tree, then by
#' # label: MPTinR::check.mpt(eqn_file)$eqn.order.categories lists the
#' # labels (here cr, fa for tree new, then hit, miss for tree old). Name the
#' # columns by tree and label, so that a label used in several trees stays
#' # apart, and select with the same key
#' wide_mptinr <- rbind(c(32, 8, 30, 10), c(28, 12, 25, 15))
#' colnames(wide_mptinr) <- c("new.cr", "new.fa", "old.hit", "old.miss")
#' long2 <- do.call(rbind, lapply(split(map, map$tree), function(tree_map) {
#'   counts <- as.data.frame(matrix(0, nrow(wide_mptinr), length(cats)))
#'   names(counts) <- cats
#'   key <- paste(tree_map$tree, tree_map$file_name, sep = ".")
#'   counts[tree_map$bmm_name] <- wide_mptinr[, key]
#'   data.frame(id = seq_len(nrow(wide_mptinr)), item_type = tree_map$tree[1],
#'              counts)
#' }))
#' identical(long2, long)
#' @export
mpt_from_eqn <- function(file, restrictions = NULL, categories = NULL,
                         tree_id = NULL, covariates = NULL, simplex = NULL,
                         links = "logit", impossible = NULL,
                         tree_names = NULL) {
  stop_missing_args()
  stopif(
    !is.null(categories) && !.mpt_is_named_map(categories),
    "The categories argument must be a named character vector without \\
    missing values that maps the category names of the EQN file onto \\
    response categories, e.g. \\
    categories = c(hit = 'yes', fa = 'yes', miss = 'no', cr = 'no')."
  )
  stopif(
    !is.null(tree_names) && !.mpt_is_named_map(tree_names),
    "The tree_names argument must be a named character vector without \\
    missing values that maps the tree labels of the EQN file onto tree \\
    names, e.g. tree_names = c('1' = 'old', '2' = 'new')."
  )
  eqn_file <- .mpt_read_eqn(file)
  tree_map <- .mpt_map_eqn_trees(eqn_file$eqn, tree_names, eqn_file$title_line)
  parameters <- .mpt_rename_eqn_parameters(
    eqn_file$exprs, restrictions, covariates %||% character(0)
  )
  renamed_cats <- .mpt_rename_eqn_categories(
    eqn_file$eqn, categories, tree_map, eqn_file$title_line
  )
  eqn <- renamed_cats$eqn
  eqn$expr <- parameters$exprs

  renaming <- rbind(
    .mpt_renaming_rows(
      "parameter", NA_character_, names(parameters$renaming),
      unlist(parameters$renaming)
    ),
    .mpt_renaming_rows("tree", NA_character_, names(tree_map), tree_map),
    .mpt_renaming_rows(
      "category", renamed_cats$rows$tree, renamed_cats$rows$category,
      renamed_cats$rows$bmm_category
    ),
    make.row.names = FALSE
  )
  if (any(renaming$renamed)) {
    message2(
      "Names that differ from the EQN file (stored in \\
      attr(model, 'mpt_renaming')):
      {.mpt_format_renaming(renaming[renaming$renamed, ])}"
    )
  }

  branches <- lapply(
    split(eqn, factor(eqn$tree, levels = unique(eqn$tree))),
    function(tree_lines) {
      .mpt_sum_branch_lines(tree_lines$expr, tree_lines$bmm_category)
    }
  )
  .mpt_check_import_trees(
    branches, impossible,
    relabel_hint = paste(
      "If the trees label the same responses differently (e.g., hit/miss and",
      "fa/cr for yes/no), map the labels onto shared response categories with",
      "the categories argument, e.g. categories = c(hit = 'yes', fa = 'yes',",
      "miss = 'no', cr = 'no')."
    ),
    # a title forms a tree of its own, which lacks every other category
    lead_hint = if (!anyNA(eqn_file$title_line)) {
      .mpt_title_sentence(eqn_file$title_line)
    }
  )
  model <- .mpt_model_from_branches(
    branches, impossible, tree_id = tree_id, covariates = covariates,
    simplex = simplex, restrictions = parameters$restrictions, links = links
  )
  attr(model, "mpt_renaming") <- renaming
  attr(model, "mpt_source") <- list(
    file = normalizePath(file),
    md5 = unname(tools::md5sum(file)),
    lines = eqn_file$raw_lines
  )
  attr(model, "call") <- match.call()
  model
}

# The helpers below are called only from mpt_from_eqn(), so their errors can
# speak about the EQN file and the importer's arguments.

# reads the equation lines of an EQN file, keeping their line numbers for
# the error messages
.mpt_read_eqn <- function(file) {
  raw_lines <- readLines(file, warn = FALSE)
  lines <- trimws(raw_lines)
  line_no <- which(nzchar(lines) & !startsWith(lines, "#"))
  # classical EQN files start with the number of equations; MPTinR and
  # TreeBUGS drop the first line whatever it holds, bmm only a lone integer
  if (length(line_no) > 0 && grepl("^[0-9]+$", lines[line_no[1]])) {
    warnif(
      as.numeric(lines[line_no[1]]) != length(line_no) - 1,
      "The first line of the EQN file announces {lines[line_no[1]]} \\
      equations, but {length(line_no) - 1} equation lines follow. Please \\
      check that the file is complete."
    )
    line_no <- line_no[-1]
  }
  stopif(
    length(line_no) == 0,
    "The EQN file '{file}' contains no equation lines."
  )
  fields <- strsplit(lines[line_no], "[[:space:]]+")
  # a title read as an equation forms a tree of its own
  title_line <- if (length(fields) > 1 &&
                      !fields[[1]][1] %in% vapply(fields[-1], `[`, "", 1)) {
    line_no[1]
  } else {
    NA_integer_
  }
  short <- lengths(fields) < 3
  stopif(
    any(short),
    "Each line of an EQN file must contain three whitespace-separated \\
    fields: tree, response category, and branch probability expression. \\
    These lines do not:
    {paste0('  line ', line_no[short], ': ',
            sQuote(lines[line_no[short]], FALSE), collapse = '\\n')}\\
    {.mpt_title_hint(line_no[short], title_line)}"
  )
  eqn <- data.frame(
    tree = vapply(fields, `[[`, character(1), 1),
    category = vapply(fields, `[[`, character(1), 2),
    expr = vapply(fields, function(f) paste(f[-(1:2)], collapse = " "), character(1)),
    line = line_no
  )
  exprs <- mapply(
    .mpt_parse_eqn_expression, eqn$expr, line_no, lines[line_no], title_line,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )
  nlist(raw_lines, eqn, exprs, title_line)
}

# TreeBUGS and MPTinR skip the first line of an EQN file, so it may hold a
# title, which bmm reads as an equation; errors on that line suggest removing it
.mpt_title_hint <- function(bad_lines, title_line) {
  # title_line is NA when line 1 does not look like a title; bad_lines are
  # line numbers and never NA
  if (title_line %in% bad_lines) {
    paste0("\n", .mpt_title_sentence(title_line))
  } else {
    ""
  }
}

.mpt_title_sentence <- function(title_line) {
  paste0("If line ", title_line, " is a title, delete it or start it with #.")
}

.mpt_parse_eqn_expression <- function(expr, line_no, line, title_line) {
  tryCatch(str2lang(expr), error = function(e) {
    # R's parser reports '<text>:1:3: unexpected symbol' plus an excerpt
    reason <- sub(
      "^<text>:[0-9:]+ *", "", strsplit(conditionMessage(e), "\n")[[1]][1]
    )
    stop2(
      "Cannot parse the branch probability expression on line {line_no} of \\
      the EQN file ({reason}): '{line}'{.mpt_title_hint(line_no, title_line)}"
    )
  })
}

# maps the file's tree labels onto bmm tree names
.mpt_map_eqn_trees <- function(eqn, tree_names, title_line) {
  stopif(
    !all(names(tree_names) %in% eqn$tree),
    "The tree_names argument names trees that are not in the EQN file: \\
    {collapse_comma(setdiff(names(tree_names), eqn$tree))}. The file uses: \\
    {collapse_comma(unique(eqn$tree))}."
  )
  tree_map <- setNames(unique(eqn$tree), unique(eqn$tree))
  tree_map[names(tree_names)] <- unlist(tree_names)
  stopif(
    anyDuplicated(tree_map) > 0,
    "Each tree needs its own name, but tree_names gives several trees the \\
    same name: {collapse_comma(unique(tree_map[duplicated(tree_map)]))}."
  )
  # checked here rather than left to mpt(), whose message cannot point to
  # tree_names
  bad_tree_names <- tree_map[!grepl("^[A-Za-z][A-Za-z0-9_]*$", tree_map)]
  bad_lines <- eqn$line[match(names(bad_tree_names), eqn$tree)]
  stopif(
    length(bad_tree_names) > 0,
    "Tree names must start with a letter and contain only letters, digits, \\
    or underscores. Please rename \\
    {paste0(sQuote(bad_tree_names, FALSE), ' (first on line ', bad_lines, ')',
            collapse = ', ')} \\
    with the tree_names argument, e.g. \\
    tree_names = c('1' = 'old', '2' = 'new').\\
    {.mpt_title_hint(bad_lines, title_line)}"
  )
  tree_map
}

# renames the parameters in the branch expressions and in the restrictions,
# which are written with the file's names
.mpt_rename_eqn_parameters <- function(exprs, restrictions, covariates) {
  file_symbols <- unique(unlist(lapply(exprs, all.vars)))
  absent <- setdiff(covariates, file_symbols)
  stopif(
    length(absent) > 0,
    "The covariates argument names symbols that are not in the EQN file: \\
    {collapse_comma(absent)}. Write covariates as the file does, before \\
    underscores and dots are removed.\\
    {if (any(gsub('[._]', '', file_symbols) %in% absent)) {
      paste0(' The file has ', collapse_comma(
        file_symbols[gsub('[._]', '', file_symbols) %in% absent]
      ), '.')
    } else {
      ''
    }}"
  )
  restrictions <- .mpt_parse_restrictions(restrictions)
  unknown <- setdiff(
    c(names(restrictions), unlist(lapply(restrictions, all.vars))),
    c(file_symbols, covariates, "FE")
  )
  stopif(
    length(unknown) > 0,
    "The restrictions refer to names that are not in the EQN file: \\
    {collapse_comma(unknown)}. Write them as the file does, before \\
    underscores and dots are removed."
  )
  symbols <- unique(c(
    file_symbols, names(restrictions), unlist(lapply(restrictions, all.vars))
  ))
  renaming <- .mpt_sanitized_names(setdiff(symbols, covariates))
  .mpt_check_sanitized_names(renaming)
  symbol_map <- lapply(renaming, as.name)
  restrictions <- lapply(restrictions, function(value) {
    if (is.language(value)) .mpt_substitute_symbols(value, symbol_map) else value
  })
  names(restrictions) <- vapply(names(restrictions), function(par) {
    renaming[[par]] %||% par
  }, character(1))
  list(
    exprs = vapply(exprs, function(expr) {
      deparse1(.mpt_substitute_symbols(expr, symbol_map))
    }, character(1)),
    restrictions = restrictions,
    renaming = renaming
  )
}

# adds the bmm name of each line's category (bmm_category) and the bmm tree
# names (eqn), plus one row per category of a tree (rows); categories maps
# first, then the '<tree>_' prefix strip and the removal of underscores and
# dots apply to the unmapped labels
.mpt_rename_eqn_categories <- function(eqn, categories, tree_map, title_line) {
  stopif(
    !all(names(categories) %in% eqn$category),
    "The categories argument names categories that are not in the EQN file: \\
    {collapse_comma(setdiff(names(categories), eqn$category))}. The file \\
    uses: {collapse_comma(unique(eqn$category))}."
  )
  category_clean <- eqn$category
  user_mapped <- category_clean %in% names(categories)
  category_clean[user_mapped] <- unlist(categories[category_clean[user_mapped]])
  has_tree_prefix <- !user_mapped &
    startsWith(category_clean, paste0(eqn$tree, "_"))
  category_clean[has_tree_prefix] <- substring(
    category_clean[has_tree_prefix], nchar(eqn$tree[has_tree_prefix]) + 2
  )
  unmapped <- unique(data.frame(
    tree = eqn$tree, original = eqn$category, clean = category_clean
  )[!user_mapped, ])
  clash <- unmapped[duplicated(unmapped[c("tree", "clean")]), ]
  stopif(
    nrow(clash) > 0,
    "In tree '{clash$tree[1]}', the categories \\
    {collapse_comma(unmapped$original[unmapped$tree == clash$tree[1] &
                                      unmapped$clean == clash$clean[1]])} \\
    would all be named '{clash$clean[1]}' once the '<tree>_' prefix is \\
    removed. Rename them in the EQN file or map them with the categories \\
    argument."
  )
  cat_renaming <- .mpt_sanitized_names(unique(category_clean))
  .mpt_check_sanitized_names(cat_renaming)
  eqn$bmm_category <- unlist(cat_renaming[category_clean], use.names = FALSE)
  eqn$prefix <- ifelse(has_tree_prefix, paste0(eqn$tree, "_"), "")

  # checked here rather than left to mpt(), whose message cannot point to
  # categories; numeric labels are common in classical EQN files
  bad <- !grepl("^[A-Za-z][A-Za-z0-9]*$", eqn$bmm_category)
  bad_cats <- unique(eqn[bad, c("category", "bmm_category", "line")])
  bad_cats <- bad_cats[!duplicated(bad_cats$category), ]
  stopif(
    any(bad),
    "Response category names must start with a letter and contain only \\
    letters and digits. Please give \\
    {paste0(sQuote(bad_cats$category, FALSE),
            ifelse(bad_cats$category == bad_cats$bmm_category, ' (',
                   paste0(' (becomes ', sQuote(bad_cats$bmm_category, FALSE),
                          ', ')),
            'first on line ', bad_cats$line, ')', collapse = ', ')} \\
    valid names with the categories argument, e.g. \\
    categories = c('1' = 'hit', '2' = 'miss').\\
    {.mpt_title_hint(bad_cats$line, title_line)}"
  )
  eqn$tree <- unname(tree_map[eqn$tree])

  # only the categories map can merge categories of a tree here: merges by
  # the prefix strip or by sanitizing alone were errors above
  cat_rows <- unique(eqn[c("tree", "category", "bmm_category", "prefix")])
  merge_key <- cat_rows[c("tree", "bmm_category")]
  merged <- cat_rows[
    duplicated(merge_key) | duplicated(merge_key, fromLast = TRUE),
  ]
  if (nrow(merged) > 0) {
    message2(
      "The categories argument maps several categories of a tree onto one \\
      response category; their branch lines are summed:
      {.mpt_format_merges(merged)}"
    )
  }
  list(eqn = eqn, rows = cat_rows)
}

.mpt_is_named_map <- function(x) {
  !is.null(names(x)) && !anyNA(names(x)) && all(nzchar(names(x))) &&
    !anyNA(x)
}

# the importers check impossible and the category sets themselves, so that
# their errors can point to the importers' arguments, which mpt() does not have
.mpt_check_import_trees <- function(branches, impossible, relabel_hint,
                                    lead_hint = NULL) {
  stopif(
    length(impossible) > 0 && !is_namedlist(impossible),
    "The impossible argument must be a named list with the categories that \\
    cannot occur in each tree, e.g. impossible = list(tree1 = c('a', 'b'))."
  )
  unknown_trees <- setdiff(names(impossible), names(branches))
  stopif(
    length(unknown_trees) > 0,
    "The impossible argument names trees that are not in the model: \\
    {collapse_comma(unknown_trees)}. The trees are: \\
    {collapse_comma(names(branches))}."
  )
  model_cats <- unique(unlist(lapply(branches, names)))
  unknown_cats <- setdiff(unlist(impossible), model_cats)
  stopif(
    length(unknown_cats) > 0,
    "The impossible argument names categories that are not in the model: \\
    {collapse_comma(unknown_cats)}. The categories are: \\
    {collapse_comma(model_cats)}."
  )
  lacking <- Map(function(tree_branches, tree_name) {
    setdiff(model_cats, c(names(tree_branches), impossible[[tree_name]]))
  }, branches, names(branches))
  lacking <- lacking[lengths(lacking) > 0]
  stopif(
    length(lacking) > 0,
    "In bmm all trees share the same response categories, but some trees \\
    have no branch for categories that other trees use:
    {paste0('  ', names(lacking), ': ', vapply(lacking, collapse_comma, ''),
            collapse = '\\n')}
    {paste(Filter(nzchar, c(
      lead_hint,
      paste0(
        'If these categories cannot occur in those trees, declare them with ',
        'the impossible argument: impossible = list(',
        paste0(names(lacking), ' = c(', vapply(lacking, collapse_comma, ''),
               ')', collapse = ', '),
        ').'
      ),
      relabel_hint
    )), collapse = '\\n')}"
  )
  invisible(NULL)
}

.mpt_model_from_branches <- function(branches, impossible, ...) {
  trees <- Map(function(tree_name, tree_branches) {
    mpt_tree(tree_name, tree_branches, impossible[[tree_name]])
  }, names(branches), branches)
  mpt(trees = trees, ...)
}

.mpt_renaming_rows <- function(kind, tree, file_name, bmm_name) {
  data.frame(
    kind = rep(kind, length(file_name)),
    tree = rep(tree, length.out = length(file_name)),
    file_name = as.character(file_name),
    bmm_name = as.character(bmm_name),
    renamed = as.character(file_name) != as.character(bmm_name)
  )
}

.mpt_format_renaming <- function(renaming) {
  single <- renaming[renaming$kind != "category", ]
  cats <- renaming[renaming$kind == "category", ]
  paste(c(
    paste0(
      "  ", format(single$kind), "  ", single$file_name, " -> ", single$bmm_name,
      recycle0 = TRUE
    ),
    vapply(unique(cats$tree), function(tree) {
      in_tree <- cats[cats$tree == tree, ]
      paste0(
        "  categories in tree ", tree, ": ",
        paste(
          in_tree$file_name, in_tree$bmm_name, sep = " -> ", collapse = ", "
        )
      )
    }, character(1))
  ), collapse = "\n")
}

.mpt_format_merges <- function(merged) {
  groups <- unique(merged[c("tree", "bmm_category")])
  paste(vapply(seq_len(nrow(groups)), function(i) {
    members <- merged[
      merged$tree == groups$tree[i] &
        merged$bmm_category == groups$bmm_category[i],
    ]
    stripped <- members[nzchar(members$prefix), ]
    paste0(
      "  tree ", groups$tree[i], ": ", collapse_comma(members$category),
      " -> ", groups$bmm_category[i],
      if (nrow(stripped) > 0) {
        paste0(
          " (after removing the prefix ",
          paste(sQuote(stripped$prefix, FALSE), "from",
                sQuote(stripped$category, FALSE), collapse = ", "),
          ")"
        )
      }
    )
  }, character(1)), collapse = "\n")
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

# maps each name to a version without underscores and dots (the brms nlpar
# constraint); names that are already valid map to themselves
.mpt_sanitized_names <- function(names) {
  setNames(as.list(gsub("[._]", "", names)), names)
}

.mpt_check_sanitized_names <- function(renaming) {
  sanitized <- unlist(renaming)
  clashes <- sanitized[duplicated(sanitized)]
  stopif(
    length(clashes) > 0,
    "Removing underscores and dots produces duplicated names: \\
    {collapse_comma(unique(clashes))}. Please rename them in the model file."
  )
  invisible(NULL)
}
