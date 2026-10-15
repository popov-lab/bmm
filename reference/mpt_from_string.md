# Create an MPT model from an MPTinR-style model string

Reads a model definition in the "easy" format of MPTinR model files and
returns the corresponding
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md) model object; EQN
files (multiTree, MPTinR, TreeBUGS) are read by
[`mpt_from_eqn()`](https://venpopov.com/bmm/reference/mpt_from_eqn.md).
Each line is one branch probability expression; blank lines separate
trees. The response category of each line is given either by an inline
comment (`expression # category`) or by the `categories` argument. Lines
that terminate in the same response category are summed.

## Usage

``` r
mpt_from_string(
  text,
  restrictions = NULL,
  tree_names,
  categories = NULL,
  tree_id = NULL,
  covariates = NULL,
  simplex = NULL,
  links = "logit",
  impossible = NULL
)
```

## Arguments

- text:

  Character. The model definition, either as a single string or as a
  vector of lines (e.g., from
  [`readLines()`](https://rdrr.io/r/base/readLines.html)).

- restrictions:

  Parameter restrictions as in
  [`mpt()`](https://venpopov.com/bmm/reference/mpt.md), written with the
  parameter names of the model text, e.g. `"g = 0.5"` or
  `list(g = 0.5)`. They mirror the `restrictions` argument of MPTinR and
  TreeBUGS. Restriction files are not read: pass their lines (see
  Examples).

- tree_names:

  Character vector with one name per tree block. The names are matched
  to the blocks by position: the first name goes to the first block, and
  so on. bmm cannot check this order, so names given in a different
  order than the blocks import a different model without an error. For
  multi-tree models the names must match the values of the tree
  identifier column in the data.

- categories:

  The response categories of lines without an inline `# category`
  comment, matched to the lines of a tree by position: the i-th entry is
  the category of the i-th line. Either a list by tree,
  `list(<tree> = c("<category>", ...))`, with one vector per tree that
  has such lines, or one unnamed character vector shared by all trees. A
  tree's vector in the list must not be longer than the tree. A shared
  vector must not be longer than the longest tree; a shorter tree takes
  its leading entries, which suits trees that lack the last categories
  (see `impossible`). If the trees use different categories, as in the
  pair-clustering model, use the list. Either form needs an entry for
  every line without an inline comment; an inline `# category` comment
  takes precedence over the entry for its line.

- tree_id:

  Character. Name of the data column whose values name the tree each
  observation belongs to (see
  [`mpt()`](https://venpopov.com/bmm/reference/mpt.md)). Add this column
  to the data yourself: for each row it holds the name, as given in
  `tree_names`, of the tree the counts come from.

- covariates, simplex, links:

  passed to [`mpt()`](https://venpopov.com/bmm/reference/mpt.md).
  `print(model)` lists the symbols that became parameters.

- impossible:

  A named list with the response categories that cannot occur in a tree,
  one element per such tree: `list(<tree> = c("<category>", ...))` (see
  the `impossible` argument of
  [`mpt_tree()`](https://venpopov.com/bmm/reference/mpt_tree.md)). Trees
  and categories are named as in the returned model, i.e. by
  `tree_names` and by the category labels of the lines. Trees that are
  not listed have no impossible categories.

## Value

An object of class `bmmodel` (see
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md)), with an
`mpt_source` attribute holding the imported text.

## Details

**Comments.** A line whose first non-blank character is `#` is a comment
and is skipped. It does not separate trees; only blank lines do. On a
line with an expression, the text after `#` is the response category of
that line. This differs from MPTinR, which ignores everything after a
`#` and treats a line that starts with `#` as a blank line, so that
there comment lines can separate trees. An MPTinR model file whose trees
are separated only by comment lines, or whose inline comments are not
category labels, has to be edited before it is imported.

**Trees with different response categories.** In bmm all trees share the
same response categories. If a category cannot occur in some trees
(e.g., the pair and singleton trees of the pair-clustering model),
declare it with the `impossible` argument. If the trees only label the
same responses differently, give those lines the same category label.

**Names.** Parameter names are used as written in the model text and are
case-sensitive, so the formulas passed to
[`bmm()`](https://venpopov.com/bmm/reference/bmm.md) name them the same
way (`bmf(D ~ 1)`;
[`bmf()`](https://venpopov.com/bmm/reference/bmmformula.md) is short for
[`bmmformula()`](https://venpopov.com/bmm/reference/bmmformula.md)).
They must be valid bmm names (letters and digits, see
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md)); unlike
[`mpt_from_eqn()`](https://venpopov.com/bmm/reference/mpt_from_eqn.md),
this function does not rename them. The import fails when the branches
of a tree do not sum to 1 at numeric test values.

**Data for an imported model.** The Details of
[`mpt_from_eqn()`](https://venpopov.com/bmm/reference/mpt_from_eqn.md)
describe the layout [`bmm()`](https://venpopov.com/bmm/reference/bmm.md)
expects (one row per participant and tree, a tree column, one count
column per response category, 0 for impossible categories) and the
checks to run before fitting; they apply to models from this function as
well, with the tree names of `tree_names` and the category labels of the
lines. The data MPTinR fits to a model file in this format are wide and
positional: one row per participant (or one vector for aggregated
counts) and one unnamed column per line of the model file, in file
order, because each line is one response category. (For EQN files MPTinR
sorts the columns instead; see
[`mpt_from_eqn()`](https://venpopov.com/bmm/reference/mpt_from_eqn.md).)
Name the columns `<tree>.<category>`;
[`mpt_data_reshape()`](https://venpopov.com/bmm/reference/mpt_data_reshape.md)
then stacks the trees and gives 0 to the categories a tree cannot
produce (see Examples).

**Provenance.** The returned model carries an `mpt_source` attribute,
`list(text = )`, with the model definition as one string.

## See also

[`mpt_from_eqn()`](https://venpopov.com/bmm/reference/mpt_from_eqn.md)
for EQN files (multiTree, MPTinR, TreeBUGS);
[`mpt()`](https://venpopov.com/bmm/reference/mpt.md) and
[`mpt_tree()`](https://venpopov.com/bmm/reference/mpt_tree.md) for
building the trees in R, including the data layout
[`bmm()`](https://venpopov.com/bmm/reference/bmm.md) expects.

## Examples

``` r
model_2htm <- "
# two-high-threshold model
D + (1 - D) * g        # old
(1 - D) * (1 - g)      # new

(1 - D) * g            # old
D + (1 - D) * (1 - g)  # new
"
model <- mpt_from_string(
  model_2htm,
  tree_names = c("target", "lure"),
  tree_id = "item_type"
)
model
#> mpt_from_string(text = model_2htm,                 tree_names = c("target", "lure"),                 tree_id = "item_type") 
#> Response:   resp_cats = old, new (counts per response category) 
#> Parameters: D, g 
#> Links:      D = logit; g = logit 
#> MPT tree 'target':
#>   P(old) = D + (1 - D) * g
#>   P(new) = (1 - D) * (1 - g)
#> MPT tree 'lure':
#>   P(old) = (1 - D) * g
#>   P(new) = D + (1 - D) * (1 - g)
#> Identifiability (intercept-only formulas): 2 free parameter(s), 2 degrees of freedom (response categories minus 1, summed over trees) 
#>   Jacobian rank 2 of 2 at interior test values: locally identified.
#> Use parameter_info() for parameter descriptions, links and fixed values.

# fix the guessing rate as in MPTinR / TreeBUGS
mpt_from_string(
  model_2htm,
  tree_names = c("target", "lure"),
  tree_id = "item_type",
  restrictions = "g = 0.5"
)
#> mpt_from_string(text = model_2htm,                 restrictions = "g = 0.5",                 tree_names = c("target", "lure"),                 tree_id = "item_type") 
#> Response:   resp_cats = old, new (counts per response category) 
#> Parameters: D 
#> Links:      D = logit 
#> MPT tree 'target':
#>   P(old) = D + (1 - D) * 0.5
#>   P(new) = (1 - D) * 0.5
#> MPT tree 'lure':
#>   P(old) = (1 - D) * 0.5
#>   P(new) = D + (1 - D) * 0.5
#> Restrictions: g = 0.5 
#> Identifiability (intercept-only formulas): 1 free parameter(s), 2 degrees of freedom (response categories minus 1, summed over trees) 
#>   Jacobian rank 1 of 1 at interior test values: locally identified.
#> Use parameter_info() for parameter descriptions, links and fixed values.

# pair clustering, written as in MPTinR without category labels: pairs and
# singletons have different response categories
model_pc <- "
c * r
(1 - c) * u * u
2 * (1 - c) * u * (1 - u)
c * (1 - r) + (1 - c) * (1 - u) * (1 - u)

u
1 - u
"
mpt_from_string(
  model_pc,
  tree_names = c("pairs", "singles"),
  categories = list(
    pairs = c("E1", "E2", "E3", "E4"),
    singles = c("F1", "F2")
  ),
  tree_id = "item_type",
  impossible = list(pairs = c("F1", "F2"), singles = c("E1", "E2", "E3", "E4"))
)
#> mpt_from_string(text = model_pc,                 tree_names = c("pairs", "singles"),                 categories = list(pairs = c("E1", "E2", "E3", "E4"), singles = c("F1", "F2")),                 tree_id = "item_type",                 impossible = list(pairs = c("F1", "F2"), singles = c("E1", "E2", "E3", "E4"))) 
#> Response:   resp_cats = E1, E2, E3, E4, F1, F2 (counts per response category) 
#> Parameters: c, r, u 
#> Links:      c = logit; r = logit; u = logit 
#> MPT tree 'pairs':
#>   P(E1) = c * r
#>   P(E2) = (1 - c) * u * u
#>   P(E3) = 2 * (1 - c) * u * (1 - u)
#>   P(E4) = c * (1 - r) + (1 - c) * (1 - u) * (1 - u)
#>   P(F1) = 0 (structurally impossible)
#>   P(F2) = 0 (structurally impossible)
#> MPT tree 'singles':
#>   P(F1) = u
#>   P(F2) = 1 - u
#>   P(E1) = 0 (structurally impossible)
#>   P(E2) = 0 (structurally impossible)
#>   P(E3) = 0 (structurally impossible)
#>   P(E4) = 0 (structurally impossible)
#> Identifiability (intercept-only formulas): 3 free parameter(s), 4 degrees of freedom (response categories minus 1, summed over trees) 
#>   Jacobian rank 3 of 3 at interior test values: locally identified.
#> Use parameter_info() for parameter descriptions, links and fixed values.

# an MPTinR model file read with readLines() keeps the blank lines that
# separate the trees. A restrictions file is not read by bmm: drop its
# comments and blank lines and pass the remaining lines.
model_file <- tempfile(fileext = ".model")
writeLines(model_2htm, model_file)
restr_file <- tempfile(fileext = ".restr")
writeLines(c("# fix the guessing rate", "g = 0.5"), restr_file)
restr <- trimws(sub("#.*", "", readLines(restr_file)))
model <- mpt_from_string(
  readLines(model_file),
  tree_names = c("target", "lure"),
  tree_id = "item_type",
  restrictions = restr[nzchar(restr)]
)

# MPTinR data: one unnamed column per line of the model file, in file
# order. Name the columns by tree and category, then let mpt_data_reshape()
# stack the trees to get bmm's one row per participant and tree.
wide <- rbind(c(30, 10, 8, 32), c(25, 15, 12, 28))
colnames(wide) <- c("target.old", "target.new", "lure.old", "lure.new")
long <- mpt_data_reshape(wide, model)
long
#>   id item_type old new
#> 1  1    target  30  10
#> 2  1      lure   8  32
#> 3  2    target  25  15
#> 4  2      lure  12  28
# the row sums are the trials per row (40 here), but they cannot show
# swapped columns or tree labels: compare one participant with the wide
# row, then run the data checks (bmm_data_check() prints
# "Hard checks: passed" when the columns, tree labels and counts fit)
rowSums(long[c("old", "new")])
#> [1] 40 40 40 40
long[long$id == 1, ]
#>   id item_type old new
#> 1  1    target  30  10
#> 2  1      lure   8  32
wide[1, ]
#> target.old target.new   lure.old   lure.new 
#>         30         10          8         32 
bmm_data_check(bmf(D ~ 1), long, model)
#> Pre-fit data check
#> Model: mpt_from_string(text = readLines(model_file), restrictions = restr[nzchar(restr)], tree_names = c("target", "lure"), tree_id = "item_type") 
#> Data:  4 observations ('long') 
#> 
#> Response variables:
#>   old  numeric with 4 unique values: 8, 12, 25, 30 (expected: counts per
#>        response category)
#>   new  numeric with 4 unique values: 10, 15, 28, 32 (expected: counts per
#>        response category)
#> 
#> Findings: No issues detected.
#> 
#> Hard checks: passed
```
