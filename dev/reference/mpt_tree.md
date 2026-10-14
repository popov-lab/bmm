# Create a tree for a Multinomial Processing Tree (MPT) model

Specifies a single tree of a multinomial processing tree model by naming
the tree and providing one branch-probability expression per response
category. Trees created with `mpt_tree()` are combined into a full model
specification with
[`mpt()`](https://popov-lab.github.io/bmm/dev/reference/mpt.md).

## Usage

``` r
mpt_tree(name, branches)
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
  parameters (e.g., `"D + (1 - D) * g"`) and numeric constants.

## Value

An object of class `mpt_tree`

## Details

Constant arithmetic such as `1/4` or `1/(2*2)` is folded into a decimal
literal (`0.25`) when the tree is created. Stan compiles a bare integer
fraction as integer division (`1/4 == 0`), which would silently corrupt
the likelihood.

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
```
