# Create an MPT model from an EQN model file

Reads an EQN model file, the format used by multiTree, MPTinR, and
TreeBUGS, and returns the corresponding
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md) model
object; model files in MPTinR's own "easy" format are read by
[`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md).
Each line of an EQN file contains three whitespace-separated fields –
tree, response category, branch probability expression – and branch
lines that terminate in the same response category of a tree are summed.

## Usage

``` r
mpt_from_eqn(
  file,
  restrictions = NULL,
  categories = NULL,
  tree_id = NULL,
  covariates = NULL,
  simplex = NULL,
  links = "logit",
  impossible = NULL,
  tree_names = NULL
)
```

## Arguments

- file:

  Character. Path to the EQN file. How its first line is read is
  described in Details.

- restrictions:

  Parameter restrictions as in
  [`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md),
  written with the parameter names used in the EQN file, e.g.
  `c("D_n = D_o", "g = 0.5")` or `list(D_n = "D_o", g = 0.5)`. They
  mirror the `restrictions` argument of MPTinR and TreeBUGS and are
  renamed together with the parameters. Restriction files are not read:
  pass their lines (the Examples of
  [`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md)
  read one with [`readLines()`](https://rdrr.io/r/base/readLines.html)).

- categories:

  A named character vector mapping response category names used in the
  EQN file onto the shared response categories of the model (see
  Details), e.g. `c(hit = "yes", fa = "yes", miss = "no", cr = "no")`.
  Categories not listed keep their name. Mapped values and unlisted
  names alike are sanitized (underscores and dots removed, see Details),
  so `c(hit = "y_es")` gives `yes`. Labels that only carry their tree
  name as a prefix (`WC_Correct`, `LC_Correct`) need no mapping, because
  the prefix is removed. Labels that are not valid bmm names, such as
  numeric labels, must be mapped. Mapping two categories of the same
  tree onto one name merges them: their branch lines are summed, and a
  message says so.

- tree_id:

  Character. Name of the data column whose values name the tree each
  observation belongs to (see
  [`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md)). Add
  this column to the data yourself: for each row it holds the name of
  the tree the counts come from, as printed by the model (the file's
  label, or its replacement from `tree_names`).

- covariates, simplex, links:

  passed to
  [`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md).
  Covariate names are written as in the EQN file and are excluded from
  the renaming; simplex groups use the bmm names. `print(model)` lists
  the symbols that became parameters.

- impossible:

  A named list, `list(<tree> = c("<category>", ...))`, of the response
  categories that cannot occur in a tree (see
  [`mpt_tree()`](https://popov-lab.github.io/bmm/dev/reference/mpt_tree.md)),
  written with the bmm names of the trees and categories, i.e. the
  `bmm_name` column of the renaming map; `print(model)` lists the shared
  categories under "Response". Trees not listed have no impossible
  categories. When the trees do not share their categories and
  `impossible` does not cover the difference, the error names the
  categories each tree lacks.

- tree_names:

  A named character vector mapping tree labels of the EQN file onto tree
  names in bmm, e.g. `c("1" = "target", "2" = "lure")`. bmm tree names
  must start with a letter and contain only letters, digits, or
  underscores, so files with numeric tree labels need this mapping.
  Labels not listed keep their name.

## Value

An object of class `bmmodel` (see
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md)), with
the attributes `mpt_renaming` and `mpt_source` (see Details).

## Details

**File format.** Blank lines and lines that start with `#` are skipped.
If the first remaining line holds a single integer and nothing else, it
is read as the equation count of classical EQN files and skipped, with a
warning when it differs from the number of equation lines that follow.
Any other first line is read as an equation. Note that TreeBUGS always
ignores the first line of an EQN file, whatever it holds, and so does
MPTinR with `model.type = "eqn"` (its default for `.eqn` files; with
`"eqn2"` it reads the first line as the equation count). bmm therefore
reads a file without a count line in full, while those programs drop its
first equation; a title line produces an error that says to delete it or
start it with `#`. A file that starts with the equation count, or with
`#` lines followed by equations, is read the same way by all three. A
line with fewer than three fields, or an expression that cannot be
parsed, is an error that names the line; a tree whose branches do not
sum to 1 at numeric test values is an error that names the tree.

**Response categories.** In bmm, every observation (row) belongs to one
tree, and all trees share the same response categories. Classical EQN
files often label the response categories per tree instead (e.g.,
`hit`/`miss` in the old-item tree and `fa`/`cr` in the new-item tree of
a recognition model). Use the `categories` argument to map such
tree-specific labels onto the shared response categories (e.g., the
response options "yes" and "no"). If a category cannot occur in some
trees (e.g., the pair and singleton trees of the pair-clustering model),
declare it with the `impossible` argument instead.

**Names.** Because brms does not allow underscores or dots in non-linear
parameter names, parameter and response category names are sanitized: a
leading `<tree>_` prefix is removed from category names that are not
mapped with `categories`, and all remaining underscores and dots are
stripped. Two labels of the file, or two different values within one
tree, that become the same name are an error; mapped values that become
the same name in different trees are that one shared category. Names are
case-sensitive, so `P_G1` and `P_g1` stay distinct (`PG1`, `Pg1`). Tree
labels are kept, unless `tree_names` maps them. The formulas passed to
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) use the
bmm names (`bmf(Do ~ 1)` for a file parameter `D_o`;
[`bmf()`](https://popov-lab.github.io/bmm/dev/reference/bmmformula.md)
is short for
[`bmmformula()`](https://popov-lab.github.io/bmm/dev/reference/bmmformula.md));
the `parameter` rows of the renaming map translate them, and
`print(model)` lists them. Of the arguments that name parameters, only
`restrictions` and `covariates` use the file's names.

**Renaming map.** The `mpt_renaming` attribute of the returned model is
a data frame with one row per parameter, per tree, and per category of
each tree, renamed or not, and the columns `kind` (`"parameter"`,
`"tree"`, or `"category"`), `tree` (the bmm tree name for categories,
`NA` otherwise), `file_name`, `bmm_name`, and `renamed` (`TRUE` where
the two names differ). The names that changed are also reported in a
message.

**Data for an imported model.**
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) expects
one row per participant and tree: a column named by `tree_id` holding
the bmm tree names, and one count column per shared response category,
named by `bmm_name` in the renaming map (see
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md)). A
category that is `impossible` in a tree gets 0 (or `NA`) in the rows of
that tree. TreeBUGS data are wide: one row per participant, one column
per category label of the EQN file. MPTinR data for an EQN model are
unnamed, and their columns do not follow the file: the trees and, within
each tree, the category labels are sorted (alphabetically, or
numerically when all labels in the file are numbers);
`MPTinR::check.mpt("model.eqn")` returns the order in its elements
`eqn.order.trees` and `eqn.order.categories`. (For MPTinR's own model
format the columns follow the lines instead; see
[`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md).)
[`mpt_data_reshape()`](https://popov-lab.github.io/bmm/dev/reference/mpt_data_reshape.md)
reshapes both layouts with the renaming map, sums the labels that
`categories` merged onto one response category, and gives 0 to the
categories a tree cannot produce (see Examples).

bmm cannot tell a reshape that swaps tree labels or category columns
from a correct one: such data fit without an error. The row sums of the
count columns are the trials per row; they cannot show swapped columns,
and they show swapped tree labels only when the trees differ in their
number of trials. Comparing the response proportions per tree with the
design (targets must get more "yes" responses than lures) catches
swapped tree labels and response columns swapped in every tree, but not
always a column slipped or a map entry wrong in one tree, nor a swap
between two trees of the same kind; so also compare one participant's
long rows with the wide row.
[`bmm_data_check()`](https://popov-lab.github.io/bmm/dev/reference/bmm_data_check.md)
runs the remaining checks (columns, tree labels, counts in impossible
cells) without a fit and prints `Hard checks: passed` when they hold.

**Provenance.** The returned model carries an `mpt_source` attribute,
`list(file = , md5 = , lines = )`: the full path of the EQN file, its
MD5 checksum, and the lines as read, so that a saved fit records which
file it was built from.

## See also

[`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md)
for model text in MPTinR's own format;
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md) and
[`mpt_tree()`](https://popov-lab.github.io/bmm/dev/reference/mpt_tree.md)
for building the trees in R, including the data layout
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) expects.

## Examples

``` r
eqn_file <- tempfile(fileext = ".eqn")
writeLines(c(
  "6",
  "target  hit   D_o",
  "target  hit   (1-D_o)*g",
  "target  miss  (1-D_o)*(1-g)",
  "lure    fa    (1-D_n)*g",
  "lure    cr    D_n",
  "lure    cr    (1-D_n)*(1-g)"
), eqn_file)
# restrictions use the file's names; the formulas use the bmm names
model <- mpt_from_eqn(
  eqn_file,
  restrictions = "D_n = D_o",
  categories = c(hit = "yes", fa = "yes", miss = "no", cr = "no"),
  tree_id = "item_type"
)
#> Names that differ from the EQN file (stored in attr(model, 'mpt_renaming')):
#>   parameter  D_o -> Do
#>   parameter  D_n -> Dn
#>   categories in tree target: hit -> yes, miss -> no
#>   categories in tree lure: fa -> yes, cr -> no
attr(model, "mpt_renaming")
#>        kind   tree file_name bmm_name renamed
#> 1 parameter   <NA>       D_o       Do    TRUE
#> 2 parameter   <NA>         g        g   FALSE
#> 3 parameter   <NA>       D_n       Dn    TRUE
#> 4      tree   <NA>    target   target   FALSE
#> 5      tree   <NA>      lure     lure   FALSE
#> 6  category target       hit      yes    TRUE
#> 7  category target      miss       no    TRUE
#> 8  category   lure        fa      yes    TRUE
#> 9  category   lure        cr       no    TRUE

# TreeBUGS data: one row per participant, columns named by the category
# labels of the EQN file, no participant column. mpt_data_reshape() reshapes
# them to bmm's one row per participant and tree, through the map; the
# categories a tree cannot produce get 0
wide <- data.frame(
  hit = c(30, 25), miss = c(10, 15), fa = c(8, 12), cr = c(32, 28)
)
long <- mpt_data_reshape(wide, model)
long
#>   id item_type yes no
#> 1  1    target  30 10
#> 2  1      lure   8 32
#> 3  2    target  25 15
#> 4  2      lure  12 28
# the row sums are the trials per row (40 here). The proportions catch
# swapped tree labels and swapped response columns (targets must get
# more "yes" than lures), but not always a column slipped in one
# tree: also compare one participant's long rows with the wide row
rowSums(long[c("yes", "no")])
#> [1] 40 40 40 40
tapply(long$yes / rowSums(long[c("yes", "no")]), long$item_type, mean)
#>   lure target 
#> 0.2500 0.6875 
long[long$id == 1, ]
#>   id item_type yes no
#> 1  1    target  30 10
#> 2  1      lure   8 32
wide[1, ]
#>   hit miss fa cr
#> 1  30   10  8 32
bmm_data_check(bmf(Do ~ 1, g ~ 1), long, model)
#> Pre-fit data check
#> Model: mpt_from_eqn(file = eqn_file, restrictions = "D_n = D_o", categories = c(hit = "yes", fa = "yes", miss = "no", cr = "no"), tree_id = "item_type") 
#> Data:  4 observations ('long') 
#> 
#> Response variables:
#>   yes  numeric with 4 unique values: 8, 12, 25, 30 (expected: counts per
#>        response category)
#>   no   numeric with 4 unique values: 10, 15, 28, 32 (expected: counts per
#>        response category)
#> 
#> Findings: No issues detected.
#> 
#> Hard checks: passed

# MPTinR data for an EQN model are unnamed and sorted by tree, then by
# label: MPTinR::check.mpt(eqn_file)$eqn.order.categories lists the
# labels (here cr, fa for tree new, then hit, miss for tree old)
wide_mptinr <- rbind(c(32, 8, 30, 10), c(28, 12, 25, 15))
identical(mpt_data_reshape(wide_mptinr, model, columns = "mptinr"), long)
#> [1] TRUE
```
