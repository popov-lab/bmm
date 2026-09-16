# Create a file with a template for adding a new model (for developers)

Create a file with a template for adding a new model (for developers)

## Usage

``` r
use_model_template(
  model_name,
  custom_family = FALSE,
  stanvar_blocks = c("data", "tdata", "parameters", "tparameters", "model", "likelihood",
    "genquant", "functions"),
  open_files = TRUE,
  testing = FALSE
)
```

## Arguments

- model_name:

  A string with the name of the model. The file will be named
  `model_model_name.R` and all necessary functions will be created with
  the appropriate names and structure. The file will be saved in the
  `R/` directory

- custom_family:

  Logical; Do you plan to define a brms::custom_family()? If TRUE the
  function will add a section for the custom family, placeholders for
  the stan_vars and corresponding empty .stan files in
  `inst/stan_chunks/`, that you can fill For an example, see the sdm
  model in `/R/model_sdm.R`. If FALSE (default) the function will not
  add the custom family section nor stan files.

- stanvar_blocks:

  A character vector with the names of the blocks that will be added to
  the custom family section. See
  [`brms::stanvar()`](https://paulbuerkner.com/brms/reference/stanvar.html)
  for more details. The default lists all the possible blocks, but it is
  unlikely that you will need all of them. You can specify a vector of
  only those that you need. The function will add a section for each
  block in the list

- open_files:

  Logical; If TRUE (default), the function will open the template files
  that were created in RStudio

- testing:

  Logical; If TRUE, the function will return the file content but will
  not save the file. If FALSE (default), the function will save the file

## Value

If `testing` is TRUE, the function will return the file content as a
string. If `testing` is FALSE, the function will return NULL

## Details

If you get a warning during check() about non-ASCII characters, this is
often due to the citation field. You can find what the problem is by
running

    remotes::install_github("eddelbuettel/dang")
    dang::checkPackageAsciiCode(dir = ".")

usually rewriting the numbers (issue, page numbers) manually fixes it

## Examples

``` r
if (FALSE) { # isTRUE(Sys.getenv("BMM_EXAMPLES"))
library(usethis)


# create a new model file with a brms::custom_family, three .stan files in
# inst/stan_chunks/ and open the files
use_model_template("abc",
  custom_family = TRUE,
  stanvar_blocks = c("functions", "likelihood", "tdata")
)
}
```
