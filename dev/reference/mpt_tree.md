# Create a tree for a Multinomial Processing Tree (MPT) model

Specifies a single tree of a multinomial processing tree model by naming
the tree and providing one branch-probability expression per response
category. Trees created with `mpt_tree()` are combined into a full model
specification with
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md).

## Usage

``` r
mpt_tree(name, branches, impossible = NULL)
```

## Arguments

- name:

  Character. Label for the tree. For models with multiple trees, the
  label identifies the tree in the data: the column named by the
  `tree_id` argument of
  [`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md) holds
  one such label per observation.

- branches:

  A named list of character strings. Each element gives the branch
  probability expression for one response category, and the element
  names are the response categories. The expressions can use latent
  parameters (e.g., `"D + (1 - D) * g"`), numeric constants, and
  declared covariates (see
  [`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md)).

- impossible:

  Character vector. Response categories that cannot occur in this tree,
  for example distractor responses when no distractors were presented.
  Such categories get no branch expression, and the remaining branches
  must sum to 1 on their own.

## Value

An object of class `mpt_tree`

## Details

Constant arithmetic such as `1/4` or `1/(2*2)` is folded into a decimal
literal (`0.25`) when the tree is created. Stan compiles a bare integer
fraction as integer division (`1/4 == 0`), which would silently corrupt
the likelihood.

Categories listed in `impossible` keep their column in the data and
their place in the response matrix, but their probability is set to zero
for the observations belonging to this tree. Because the model is
estimated on the log scale, "zero" is implemented as a large negative
linear predictor rather than a literal zero; no probability has to be
taken away from the remaining branches to compensate.

## Examples

``` r
tree_old <- mpt_tree(
  name = "old",
  branches = list(
    old = "D + (1 - D) * g",
    new = "(1 - D) * (1 - g)"
  )
)
tree_old
#> MPT tree 'old':
#>   P(old) = D + (1 - D) * g
#>   P(new) = (1 - D) * (1 - g)

# a three-alternative recognition task in which lures are shown on some
# trials only: `lure` responses cannot occur on the others
tree_lures <- mpt_tree(
  name = "lures",
  branches = list(
    target = "D * b + (1 - D) * (1/3)",
    lure   = "D * (1 - b) + (1 - D) * (1/3)",
    new    = "(1 - D) * (1/3)"
  )
)

# the branches of the lure-free tree sum to 1 without the `lure` category;
# no probability is reserved for it
tree_nolures <- mpt_tree(
  name = "nolures",
  branches = list(
    target = "D + (1 - D) * (1/2)",
    new    = "(1 - D) * (1/2)"
  ),
  impossible = "lure"
)
tree_nolures
#> MPT tree 'nolures':
#>   P(target) = D + (1 - D) * 0.5
#>   P(new) = (1 - D) * 0.5
#>   P(lure) = 0 (structurally impossible)

# a `tree` column in the data holds "lures" or "nolures" per observation;
# effects of experimental conditions go in the parameter formulas instead,
# e.g. bmf(D ~ 0 + cond)
model <- mpt(list(tree_lures, tree_nolures), tree_id = "tree")
```
