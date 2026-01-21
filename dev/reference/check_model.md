# Generic S3 method for checking if the model is supported and model preprocessing

In addition for validating the model, specific methods might add
information to the model object based on the provided data and formula

## Usage

``` r
check_model(model, data = NULL, formula = NULL)
```

## Arguments

- model:

  the model argument supplied by the user

- data:

  the data argument supplied by the user

- formula:

  the formula argument supplied by the user

## Value

An object of type 'bmmodel'
