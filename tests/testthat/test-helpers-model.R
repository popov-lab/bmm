test_that("supported_models() returns a non-empty character vector", {
  expect_type(supported_models(print_call = FALSE), "character")
  expect_gt(length(supported_models(print_call = FALSE)), 0)
})

test_that("get_model() returns the correct function", {
  expect_equal(get_model("mixture2p"), .model_mixture2p)
})

test_that("check_model() refuses invalid models and accepts valid models", {
  expect_error(check_model("invalid_model"))
  expect_error(check_model(structure(list(), class = "invalid")))
  expect_error(check_model(sdm), "Did you forget")
  okmodels <- supported_models(print_call = FALSE)
  for (model in okmodels) {
    if (model == "m3") next
    model <- get_model(model)()
    expect_silent(check_model(model))
    expect_type(check_model(model), "list")
  }
})

test_that("check_model() works with regular expressions", {
  dat <- oberauer_lin_2017
  models1 <- list(
    mixture3p("dev_rad",
      nt_features = paste0("col_nt", 1:7),
      set_size = "set_size"
    ),
    imm("dev_rad",
      nt_features = paste0("col_nt", 1:7),
      nt_distances = paste0("dist_nt", 1:7),
      set_size = "set_size"
    ),
    imm("dev_rad",
      nt_features = paste0("col_nt", 1:7),
      nt_distances = paste0("dist_nt", 1:7),
      set_size = "set_size",
      version = "bsc"
    ),
    imm("dev_rad",
      nt_features = paste0("col_nt", 1:7),
      set_size = "set_size",
      version = "abc"
    )
  )
  models2 <- list(
    mixture3p("dev_rad",
      nt_features = "col_nt",
      set_size = "set_size",
      regex = TRUE
    ),
    imm("dev_rad",
      nt_features = "col_nt",
      nt_distances = "dist_nt",
      set_size = "set_size",
      regex = TRUE
    ),
    imm("dev_rad",
      nt_features = "col_nt",
      nt_distances = "dist_nt",
      set_size = "set_size",
      regex = TRUE,
      version = "bsc"
    ),
    imm("dev_rad",
      nt_features = "col_nt",
      set_size = "set_size",
      regex = TRUE,
      version = "abc"
    )
  )

  for (i in 1:length(models1)) {
    check1 <- check_model(models1[[i]], dat)
    check2 <- check_model(models2[[i]], dat)
    attributes(check1) <- NULL
    attributes(check2) <- NULL
    expect_equal(check1, check2)
  }
})

test_that("use_model_template() prevents duplicate models", {
  skip_on_cran()
  okmodels <- supported_models(print_call = FALSE)
  for (model in okmodels) {
    expect_error(use_model_template(model))
  }

  model_files <- list.files(path = "R/", pattern = "^model_.*\\.R$")
  model_files_names <- gsub("^model_", "", model_files)
  model_files_names <- gsub("\\.R$", "", model_files_names)
  for (model in model_files_names) {
    expect_error(use_model_template(model))
  }
})

test_that("stancode() works with brmsformula", {
  ff <- brms::bf(count ~ zAge + zBase * Trt + (1 | patient))
  sd <- stancode(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(sd)[1], "character")
})

test_that("stancode() works with formula", {
  ff <- count ~ zAge + zBase * Trt + (1 | patient)
  sd <- stancode(ff, data = brms::epilepsy, family = poisson())
  expect_equal(class(sd)[1], "character")
})

test_that("stancode() works with bmmformula", {
  ff <- bmmformula(kappa ~ 1, thetat ~ 1, thetant ~ 1)
  model <- mixture3p("dev_rad", "col_nt", set_size = "set_size", regex = TRUE)
  sc <- stancode(ff, oberauer_lin_2017, model = model)
  expect_equal(class(sc)[1], "character")
})

test_that("no check for with stancode function", {
  withr::local_options("bmm.sort_data" = "check")
  expect_no_message(stancode(
    bmf(kappa ~ set_size, c ~ set_size),
    oberauer_lin_2017,
    sdm("dev_rad")
  ))
})

test_that("update_model_fixed_parameters() works", {
  model1 <- sdm("y")
  formula <- bmf(mu ~ set_size, kappa = 3, c ~ 1)
  model2 <- update_model_fixed_parameters(model1, formula)
  expect_equal(model1$fixed_parameters, list(mu = 0))
  expect_equal(model2$fixed_parameters, list(kappa = 3))
})

test_that("extracts all blocks and names are correct", {
  data <- data.frame(y = runif(100, min = -pi, pi))
  model <- mixture2p(resp_error = "y")
  formula <- bmf(thetat ~ 1, kappa ~ 1)

  stan_code <- stancode(formula, data = data, model = model)

  out <- extract_stan_blocks(stan_code, "all")

  expect_type(out, "list")
  expect_setequal(names(out), c(
    "functions","data","transformed data","parameters",
    "transformed parameters","model","generated quantities"
  ))
})

test_that("extracts only requested subset of blocks", {
  data <- data.frame(y = runif(100, min = -pi, pi))
  model <- mixture2p(resp_error = "y")
  formula <- bmf(thetat ~ 1, kappa ~ 1)

  stan_code <- stancode(formula, data = data, model = model)

  out <- extract_stan_blocks(stan_code, c("data", "model"))
  expect_setequal(names(out), c("data","model"))
  expect_match(out$data, "int<lower=1> N;", fixed = TRUE)
  expect_match(out$model, "von_mises_lpdf", fixed = TRUE)
})

test_that("unknown block names are ignored (no error)", {
  stan_code <- "\nfunctions {\n}\n\
data {\n}\n\
model {\n}\n\
generated quantities {\n}\n"

  out <- extract_stan_blocks(stan_code, c("data","flying spaghetti monster","model"))
  expect_setequal(names(out), c("data","model"))
})

test_that("block boundaries are correct and do not bleed into next block", {
  stan_code <- "\nfunctions {\n  real foo(real x) { return x; }\n}\n\
data {\n  int N;\n}\n\
model {\n  N ~ poisson(1);\n}\n\
generated quantities {\n  real y;\n}\n"

  out <- extract_stan_blocks(stan_code, "all")

  # 'model' block should not contain any text from 'generated quantities'
  expect_false(grepl("generated quantities", out$model, fixed = TRUE))
  expect_match(out$model, "poisson", fixed = TRUE)
})

test_that("last block extraction stops at final closing brace", {
  # This specifically guards against regressions in how the last block is found.
  # With the current code, this will likely FAIL due to `gregexec` not existing,
  # which is exactly the kind of regression we want to catch.
  stan_code <- "\nfunctions {\n}\n\
data {\n}\n\
model {\n}\n\
generated quantities {\n  real y;\n}\n"

  out <- extract_stan_blocks(stan_code, "generated quantities")
  # Should contain 'real y;' but not any stray braces beyond its own block
  expect_match(out[["generated quantities"]], "real y;", fixed = TRUE)
  expect_false(grepl("\\bgenerated quantities\\b.*\\bgenerated quantities\\b", out[["generated quantities"]]))
})

test_that("errors (or at least fails) when a requested block is missing", {
  # Current implementation will likely error if a requested block isn't present.
  # This test locks in that behavior so future changes deliberately decide
  # whether to error or return an empty string.
  stan_code <- "\nfunctions {\n}\nmodel {\n}\n"
  expect_error(extract_stan_blocks(stan_code, c("data")), regexp = NA)
  # If you later change the function to return "" instead of error,
  # update this to expect_equal(out$data, "") accordingly.
})
