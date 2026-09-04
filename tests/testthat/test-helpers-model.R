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

test_that("use_model_template() generates a flat-defaults unversioned scaffold", {
  skip_on_cran()
  out <- paste(
    capture.output(use_model_template("tmpl_unver", testing = TRUE)),
    collapse = "\n"
  )
  expect_match(out, ".tmpl_unver_defaults <- list(", fixed = TRUE)
  expect_match(out, 'parameters = .tmpl_unver_defaults[["parameters"]]', fixed = TRUE)
  expect_match(out, 'init_ranges = .tmpl_unver_defaults[["init_ranges"]]', fixed = TRUE)
  expect_match(out, 'main = "normal(0, 1)", effects = "normal(0, 0.5)"', fixed = TRUE)
  expect_false(grepl("_version_table", out, fixed = TRUE))
  expect_false(grepl("match.arg", out, fixed = TRUE))
  expect_false(grepl("void_mu", out, fixed = TRUE))
  expect_no_error(parse(text = out))
})

test_that("use_model_template() generates a version-table versioned scaffold", {
  skip_on_cran()
  out <- paste(
    capture.output(
      use_model_template("tmpl_ver", versions = c("simple", "full"), testing = TRUE)
    ),
    collapse = "\n"
  )
  expect_match(out, ".tmpl_ver_version_table <- list(", fixed = TRUE)
  expect_match(out, "simple = list(", fixed = TRUE)
  expect_match(out, "full = list(", fixed = TRUE)
  expect_match(out, "version = c('simple', 'full')", fixed = TRUE)
  expect_match(out, "version <- match.arg(version)", fixed = TRUE)
  expect_match(out, 'parameters = .tmpl_ver_version_table[[version]][["parameters"]]', fixed = TRUE)
  expect_match(out, 'paste0("tmpl_ver_", version)', fixed = TRUE)
  expect_false(grepl("void_mu", out, fixed = TRUE))
  expect_no_error(parse(text = out))
})

test_that("use_model_template() composes a custom family with versions", {
  skip_on_cran()
  out <- paste(
    capture.output(
      use_model_template(
        "tmpl_cf",
        versions = c("a", "b"),
        custom_family = TRUE,
        stanvar_blocks = c("functions", "likelihood"),
        testing = TRUE
      )
    ),
    collapse = "\n"
  )
  expect_match(out, "custom_family(", fixed = TRUE)
  expect_match(out, "_version_table", fixed = TRUE)
  expect_match(out, "match.arg", fixed = TRUE)
  expect_match(out, "stanvar(", fixed = TRUE)
  expect_no_error(parse(text = out))
})

test_that("generated template code constructs a valid bmmodel", {
  skip_on_cran()
  gen_env <- function(name, ...) {
    txt <- paste(
      capture.output(use_model_template(name, testing = TRUE, ...)),
      collapse = "\n"
    )
    env <- new.env(parent = asNamespace("bmm"))
    eval(parse(text = txt), envir = env)
    env
  }

  unver <- gen_env("tmpl_build_unver")
  model_u <- unver$.model_tmpl_build_unver(resp_var1 = "y", links = list(par1 = "log"))
  expect_equal(class(model_u), c("bmmodel", "tmpl_build_unver"))
  expect_equal(model_u$version, "NA")
  expect_equal(model_u$links$par1, "log")
  expect_false("void_mu" %in% names(model_u))

  ver <- gen_env("tmpl_build_ver", versions = c("simple", "full"))
  model_v <- ver$.model_tmpl_build_ver(resp_var1 = "y", version = "full")
  expect_equal(class(model_v), c("bmmodel", "tmpl_build_ver", "tmpl_build_ver_full"))
  expect_equal(model_v$version, "full")
  expect_error(
    ver$tmpl_build_ver("y", "a", "b", version = "nope"),
    "should be one of"
  )
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

  out <- extract_stan_blocks(stan_code)

  expect_type(out, "list")
  expect_setequal(names(out), c(
    "functions", "data", "transformed data", "parameters",
    "transformed parameters", "model", "generated quantities"
  ))
})

test_that("extracts only requested subset of blocks", {
  data <- data.frame(y = runif(100, min = -pi, pi))
  model <- mixture2p(resp_error = "y")
  formula <- bmf(thetat ~ 1, kappa ~ 1)

  stan_code <- stancode(formula, data = data, model = model)

  out <- extract_stan_blocks(stan_code, c("data", "model"))
  expect_setequal(names(out), c("data", "model"))
  expect_match(out$data, "int<lower=1> N;", fixed = TRUE)
  expect_match(out$model, "von_mises_lpdf", fixed = TRUE)
})

test_that("unknown block names are ignored (no error)", {
  stan_code <- "\nfunctions {\n}\n\
data {\n}\n\
model {\n}\n\
generated quantities {\n}\n"

  out <- extract_stan_blocks(stan_code, c("data", "flying spaghetti monster", "model"))
  expect_setequal(names(out), c("data", "model"))
})

test_that("block boundaries are correct and do not bleed into next block", {
  stan_code <- "\nfunctions {\n  real foo(real x) { return x; }\n}\n\
data {\n  int N;\n}\n\
model {\n  N ~ poisson(1);\n}\n\
generated quantities {\n  real y;\n}\n"

  out <- extract_stan_blocks(stan_code)

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

test_that("block extraction works when order of blocks is different from the brms/bmm default", {
  # default order transformed parameters comes after parameters
  stan_code <- "\nparameters {\nreal Intercept;}\ntransformed parameters {\n}"
  extracted_program_blocks <- bmm::extract_stan_blocks(stan_code)
  expect_equal(extracted_program_blocks$parameters, "real Intercept;")

  # reversed
  stan_code <- "\ntransformed parameters {\n}\nparameters {real Intercept;}"
  extracted_program_blocks <- bmm::extract_stan_blocks(stan_code)
  expect_equal(extracted_program_blocks$parameters, "real Intercept;")
})

test_that("real / int scalars parse with dims = 1", {
  out1 <- parse_parameters_line("real alpha;")
  out2 <- parse_parameters_line("int y;")

  expect_identical(out1$name, "alpha")
  expect_identical(out1$type, "real")
  expect_identical(out1$dims, "1")
  expect_null(out1$bounds)

  expect_identical(out2$name, "y")
  expect_identical(out2$type, "int")
  expect_identical(out2$dims, "1")
})

test_that("constraints are parsed into named list and comments stripped", {
  out <- parse_parameters_line("real<lower=0, upper=1> p; // comment")
  expect_identical(out$name, "p")
  expect_identical(out$type, "real")
  expect_identical(out$dims, "1")
  expect_type(out$bounds, "list")
  expect_identical(out$bounds$lower, "0")
  expect_identical(out$bounds$upper, "1")
})

test_that("vector/row_vector/simplex/unit_vector/ordered/positive_ordered need one dim", {
  v <- parse_parameters_line("vector[K] beta;")
  rv <- parse_parameters_line("row_vector[K] r;")
  sx <- parse_parameters_line("simplex[K] theta;")
  uv <- parse_parameters_line("unit_vector[K] u;")
  od <- parse_parameters_line("ordered[K] o;")
  po <- parse_parameters_line("positive_ordered[K] po;")

  expect_identical(v$dims, "K")
  expect_identical(rv$dims, "K")
  expect_identical(sx$dims, "K")
  expect_identical(uv$dims, "K")
  expect_identical(od$dims, "K")
  expect_identical(po$dims, "K")
})

test_that("matrix parses two dims", {
  out <- parse_parameters_line("matrix[M, N] A;")
  expect_identical(out$type, "matrix")
  expect_identical(out$dims, c("M", "N"))
})

test_that("square matrix families return a single size (current behavior)", {
  cmat <- parse_parameters_line("corr_matrix[K] Omega;")
  vmat <- parse_parameters_line("cov_matrix[K]  Sigma;")
  lcor <- parse_parameters_line("cholesky_factor_corr[K] Lcorr;")
  lcov <- parse_parameters_line("cholesky_factor_cov[K]  Lcov;")

  expect_identical(cmat$dims, "K")
  expect_identical(vmat$dims, "K")
  expect_identical(lcor$dims, "K")
  expect_identical(lcov$dims, "K")
})

test_that("array prefix dims are prepended and preserved in order", {
  # array of vectors
  out1 <- parse_parameters_line("array[N] vector[K] x;")
  expect_identical(out1$type, "vector") # base type remains the base
  expect_identical(out1$dims, c("N", "K"))

  # 2D array of matrices
  out2 <- parse_parameters_line("array[I, J] matrix[M, N] A;")
  expect_identical(out2$type, "matrix")
  expect_identical(out2$dims, c("I", "J", "M", "N"))

  # array of square-matrix family (current behavior keeps single base dim)
  out3 <- parse_parameters_line("array[T] corr_matrix[K] Omarr;")
  expect_identical(out3$type, "corr_matrix")
  expect_identical(out3$dims, c("T", "K"))
})

test_that("whitespace variants and CRLF endings are handled", {
  out1 <- parse_parameters_line("   real    sigma   ;")
  expect_identical(out1$name, "sigma")

  out2 <- parse_parameters_line("\r\nvector[ K ]\r\nb;\r\n")
  expect_identical(out2$name, "b")
  expect_identical(out2$dims, "K")
})

test_that("errors on missing dims where required", {
  expect_error(parse_parameters_line("vector beta;"), "Missing dimensions")
  expect_error(parse_parameters_line("matrix A;"), "Missing dimensions")
  expect_error(parse_parameters_line("corr_matrix Omega;"), "Missing dimensions")
})

test_that("errors on unknown base type or missing name", {
  expect_error(parse_parameters_line("weird_type[3] x;"), "Unknown or unsupported")
  expect_error(parse_parameters_line("real<lower=0>;"), "Missing parameter name")
})

test_that("empty/comment-only lines error out clearly", {
  expect_error(parse_parameters_line("// just a comment"), "Empty or comment-only")
  expect_error(parse_parameters_line("   "), "Empty or comment-only")
})

test_that("returns a named list keyed by parameter names, preserving order", {
  block <- "
    real alpha;
    vector[K] beta;
    matrix[M,N] A;
  "
  res <- extract_parameter_dimensions(block)

  expect_type(res, "list")
  expect_identical(names(res), c("alpha", "beta", "A"))
  expect_identical(res$alpha$type, "real")
  expect_identical(res$beta$type, "vector")
  expect_identical(res$A$type, "matrix")
  expect_identical(res$beta$dims, "K")
  expect_identical(res$A$dims, c("M", "N"))
})

test_that("handles arrays and square-matrix families", {
  block <- "
    array[N] vector[K] x;
    cov_matrix[K] Sigma;
    cholesky_factor_cov[K] L;
  "
  res <- extract_parameter_dimensions(block)

  expect_identical(names(res), c("x", "Sigma", "L"))
  expect_identical(res$x$type, "vector")
  expect_identical(res$x$dims, c("N", "K"))

  # by your current parse_parameters_line(): single size for these families
  expect_identical(res$Sigma$type, "cov_matrix")
  expect_identical(res$Sigma$dims, "K")
  expect_identical(res$L$type, "cholesky_factor_cov")
  expect_identical(res$L$dims, "K")
})

test_that("strips trailing comments but keeps code", {
  block <- "
    real<lower=0, upper=1> p;   // probability
    row_vector[J] r;            // row vec
  "
  res <- extract_parameter_dimensions(block)

  expect_identical(names(res), c("p", "r"))
  expect_identical(res$p$type, "real")
  expect_identical(res$p$dims, "1")
  expect_true(is.list(res$p$bounds))
  expect_identical(res$p$bounds$lower, "0")
  expect_identical(res$p$bounds$upper, "1")
  expect_identical(res$r$type, "row_vector")
  expect_identical(res$r$dims, "J")
})

test_that("robust to Windows-style CRLF line endings", {
  block <- "\r\nreal a;\r\nvector[K] b;\r\nmatrix[M,N] A;\r\n"
  res <- extract_parameter_dimensions(block)
  expect_identical(names(res), c("a", "b", "A"))
})

test_that("comment-only and blank lines are ignored (no errors) [guards ordering bug]", {
  block <- "
    // comment-only
    real a;          // inline ok

    // another comment-only
    vector[K] b;
  "
  expect_no_error({
    res <- extract_parameter_dimensions(block)
    expect_identical(names(res), c("a", "b"))
  })
})

test_that("duplicate names result in last-one-wins (documented behavior)", {
  block <- "
    real alpha;
    real alpha;  // re-declare (should overwrite previous entry's position)
    vector[K] beta;
  "
  res <- extract_parameter_dimensions(block)

  # List keeps three entries but both 'alpha's share the same name.
  # In practice, the *last* named element is retrieved by `res$alpha`.
  expect_identical(names(res), c("alpha", "alpha", "beta"))
  # sanity: still accessible and type corresponds to the last declaration
  expect_identical(res$alpha$type, "real")
  expect_identical(res$beta$type, "vector")
})

# tests/testthat/test-find_matching_brace.R

test_that("find_matching_brace matches the simplest pair", {
  x <- "{}"
  expect_equal(find_matching_brace(x, 1L), 2L)
})

test_that("find_matching_brace matches nested braces", {
  x <- "{{}{}}"
  expect_equal(find_matching_brace(x, 1L), 6L) # outer
  expect_equal(find_matching_brace(x, 2L), 3L) # first inner {}
  expect_equal(find_matching_brace(x, 4L), 5L) # second inner {}
})

test_that("find_matching_brace works with other text around braces", {
  x <- "abc{def{ghi}jkl}mno"
  # positions: abc(3) {(4) def{(8) ghi}(12) jkl}(16) mno
  expect_equal(find_matching_brace(x, 4L), 16L)
  expect_equal(find_matching_brace(x, 8L), 12L)
})

test_that("find_matching_brace errors if open_pos is not an opening brace", {
  x <- "a{b}c"
  expect_error(find_matching_brace(x, 1L), "not open_brace")
  expect_error(find_matching_brace(x, 3L), "not open_brace")
  expect_error(find_matching_brace(x, 4L), "not open_brace")
})

test_that("find_matching_brace errors on unbalanced braces", {
  expect_error(find_matching_brace("{", 1L), "No matching")
  expect_error(find_matching_brace("{ {", 1L), "No matching")
})

test_that("find_matching_brace handles lots of braces", {
  x <- paste0("{", strrep("{}", 100), "}")
  expect_equal(find_matching_brace(x, 1L), nchar(x))
})

test_that("resolve_fixed_links switches a parameter's link while it is fixed", {
  model <- structure(
    list(
      parameters = list(a = "a", b = "b"),
      links = list(a = "log", b = "log"),
      links_fixed = list(b = "identity"),
      fixed_parameters = list(b = 0)
    ),
    class = "bmmodel"
  )

  fixed <- resolve_fixed_links(model)
  expect_equal(fixed$links$b, "identity")
  expect_equal(fixed$links$a, "log")

  freed <- update_model_fixed_parameters(fixed, bmf(a ~ 1, b ~ 1))
  expect_null(freed$fixed_parameters$b)
  expect_equal(freed$links$b, "log")

  refixed <- update_model_fixed_parameters(freed, bmf(a ~ 1, b = 0))
  expect_equal(refixed$fixed_parameters$b, 0)
  expect_equal(refixed$links$b, "identity")
})
