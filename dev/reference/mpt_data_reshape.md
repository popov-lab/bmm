# Reshape wide MPT frequency tables to the long format of `bmm()`

MPTinR and TreeBUGS take one row per participant and one column per
response category of the EQN file.
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) takes
one row per participant and tree, a column naming the tree, and one
count column per *shared* response category (see
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md)).
`mpt_data_reshape()` reshapes the first layout into the second for a
given model.

## Usage

``` r
mpt_data_reshape(data, model, id = NULL, columns = c("names", "mptinr"))
```

## Arguments

- data:

  A data frame or matrix with one row per participant. The count columns
  are read as described under `columns`; with `columns = "names"`, all
  other columns except `id` are carried to the long data.

- model:

  An `mpt` model from
  [`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md),
  [`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md)
  or
  [`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md).
  For a model from
  [`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md)
  the `mpt_renaming` attribute maps the labels of the EQN file onto the
  response categories of the model. For any other model, the columns of
  `data` are named `<tree>.<category>`, e.g. `target.old`, with the tree
  names and response categories of the model.

- id:

  Character. The name of the participant column in `data`, or `NULL`
  (the default) to number the rows of `data`, as
  `id = seq_len(nrow(data))` in a column named `id`. The participants
  must be unique.

- columns:

  Character. How the count columns of `data` are identified: `"names"`
  (the default) by their names, as in TreeBUGS data; `"mptinr"` by their
  position, as in MPTinR data. See Details.

## Value

A data frame with one row per participant and tree, participants first.
Its columns are the participant column (`id`), the tree column named by
the `tree_id` of the model (absent in a model without it), the carried
columns, and one integer count column per response category of the
model.

## Details

**Named columns (`columns = "names"`).** For a model from
[`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md),
a column is matched to a category by the EQN label (`hit`), or, when
trees share labels, by the tree name and label (`target.hit`, with the
bmm name of the tree); the qualified name wins where both exist. A label
used in several trees must be qualified. Labels that the `categories`
argument of
[`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md)
merged onto one response category are summed. All columns of `data` that
are neither count columns nor `id` are carried unchanged to every long
row, e.g. between-participant covariates.

**Positional columns (`columns = "mptinr"`).** MPTinR data for an EQN
model have no column names. Every column of `data` except `id` is read
as a count, in the order of `check.mpt()` in MPTinR (version 1.14.1),
elements `eqn.order.trees` and `eqn.order.categories`: trees sorted by
their label in the EQN file and, within each tree, the category labels
of the file sorted. No other column can be carried; add covariates
afterwards, e.g. with [`merge()`](https://rdrr.io/r/base/merge.html) on
the participant column. The help page of `check.mpt()` states only that
these elements give the order; the sorting rule is that of its internal
`.get.EQN.model.order()`, which sorts the columns that
[`read.table()`](https://rdrr.io/r/utils/read.table.html) returns for
the file: numerically if all tree labels (all category labels) are
numbers, else as character strings in the collation of the R session.
`mpt_data_reshape()` applies the same rule, so run it in the session
where MPTinR's order was determined if the labels mix upper and lower
case. This mode needs a model from
[`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md):
MPTinR's own model format orders the columns by the lines of the model
text, and
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md) reorders
the branches of every tree to the categories of the first, so the line
order cannot be recovered from the model. For other models, name the
columns `<tree>.<category>` and use `columns = "names"`.

**Impossible categories.** A response category that a tree cannot
produce (`impossible` in
[`mpt_tree()`](https://popov-lab.github.io/bmm/dev/reference/mpt_tree.md),
[`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md)
and
[`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md))
has no column in `data` and gets 0 in the long rows of that tree, which
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) accepts.

**Counts.** The count columns must be numeric, without missing values,
and hold non-negative whole numbers; the long data hold them as
integers.

**Checks.** The function cannot tell a correct reshape from one that
swaps columns of equal meaning. Run
[`bmm_data_check()`](https://popov-lab.github.io/bmm/dev/reference/bmm_data_check.md)
on the result, and compare one participant's long rows with the wide row
(see the Examples of
[`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md)).

## See also

[`mpt_from_eqn()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_eqn.md)
and
[`mpt_from_string()`](https://popov-lab.github.io/bmm/dev/reference/mpt_from_string.md)
for importing the model itself;
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md) for the
data layout
[`bmm()`](https://popov-lab.github.io/bmm/dev/reference/bmm.md) expects;
[`bmm_data_check()`](https://popov-lab.github.io/bmm/dev/reference/bmm_data_check.md)
to check the result before fitting.

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

# TreeBUGS data: one named column per category label of the EQN file.
# Columns that are no counts, here a between-participant factor, are
# carried to the long data
wide <- data.frame(
  subject = c("s1", "s2"), group = c("a", "b"),
  hit = c(30, 25), miss = c(10, 15), fa = c(8, 12), cr = c(32, 28)
)
long <- mpt_data_reshape(wide, model, id = "subject")
long
#>   subject item_type group yes no
#> 1      s1    target     a  30 10
#> 2      s1      lure     a   8 32
#> 3      s2    target     b  25 15
#> 4      s2      lure     b  12 28
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

# MPTinR data: unnamed columns, sorted by tree, then by label (cr, fa for
# tree lure, then hit, miss for tree target)
wide_mptinr <- rbind(c(32, 8, 30, 10), c(28, 12, 25, 15))
identical(
  mpt_data_reshape(wide_mptinr, model, columns = "mptinr"),
  mpt_data_reshape(wide[c("cr", "fa", "hit", "miss")], model)
)
#> [1] TRUE

# trees with response categories that only one tree can produce
model_pc <- mpt(
  list(
    mpt_tree("pairs", list(E1 = "c * r", E2 = "1 - c * r"),
             impossible = "F1"),
    mpt_tree("singles", list(F1 = "u", E2 = "1 - u"),
             impossible = "E1")
  ),
  tree_id = "item_type"
)
wide_pc <- data.frame(
  pairs.E1 = c(20, 22), pairs.E2 = c(20, 18),
  singles.F1 = c(25, 24), singles.E2 = c(15, 16)
)
mpt_data_reshape(wide_pc, model_pc)
#>   id item_type E1 E2 F1
#> 1  1     pairs 20 20  0
#> 2  1   singles  0 15 25
#> 3  2     pairs 22 18  0
#> 4  2   singles  0 16 24
```
