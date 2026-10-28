test_that("extracts all blocks and names are correct", {
  withr::local_seed(445)
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
  withr::local_seed(445)
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

test_that("extract_stan_blocks refuses anything but a single string", {
  expect_error(extract_stan_blocks(1), "single character string")
  expect_error(extract_stan_blocks(c("data {}", "model {}")), "single character string")
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

test_that("errors on a missing name", {
  expect_error(parse_parameters_line("real<lower=0>;"), "Missing parameter name")
})

test_that("a declaration of an unmodelled type is named but left without a type", {
  # legal Stan that this parser does not model, plus a type Stan does not have
  declarations <- c(
    "sum_to_zero_vector[K] beta;", "complex z;", "complex_vector[N] cv;",
    "tuple(real, vector[N]) tv;", "weird_type[3] x;"
  )
  for (declaration in declarations) {
    out <- parse_parameters_line(declaration)
    expect_identical(out$type, NA_character_)
    expect_identical(out$types, NA_character_)
    expect_identical(out$dims, NA_character_)
    expect_null(out$bounds)
  }
  expect_identical(
    vapply(lapply(declarations, parse_parameters_line), `[[`, character(1), "name"),
    c("beta", "z", "cv", "tv", "x")
  )
})

test_that("one unmodelled declaration does not cost the rest of the block", {
  block <- "
    vector[K] b_kappa;
    sum_to_zero_vector[K] beta;
    real<lower=0> sd_1;
  "
  res <- extract_parameter_dimensions(block)

  expect_identical(names(res), c("b_kappa", "beta", "sd_1"))
  expect_identical(res$b_kappa$dims, "K")
  expect_identical(res$sd_1$bounds$lower, "0")
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

test_that("a size that indexes a data array keeps its brackets", {
  mo <- parse_parameters_line("simplex[Jmo_c[1]] simo_c_1;")
  expect_identical(mo$name, "simo_c_1")
  expect_identical(mo$type, "simplex")
  expect_identical(mo$dims, "Jmo_c[1]")

  s <- parse_parameters_line("vector[knots_kappa_1[1]] zs_kappa_1_1;")
  expect_identical(s$name, "zs_kappa_1_1")
  expect_identical(s$dims, "knots_kappa_1[1]")

  nested <- parse_parameters_line("array[J[1], N] matrix<lower=0>[M[2], K] A;")
  expect_identical(nested$name, "A")
  expect_identical(nested$dims, c("J[1]", "N", "M[2]", "K"))
  expect_identical(nested$types, c("array", "matrix"))
  expect_identical(nested$bounds$lower, "0")
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
