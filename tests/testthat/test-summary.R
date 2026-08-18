test_that("summary has reasonable outputs", {
  skip_on_cran()
  path <- test_path("assets/bmmfit_example1.rds")
  skip_if_not(file.exists(path), "SDM fixture not available (excluded by .Rbuildignore)")
  fit <- readRDS(path)
  summary1 <- suppressWarnings(summary(fit))
  expect_true(is.data.frame(summary1$fixed))
  expect_equal(
    rownames(summary1$fixed),
    c("mu_Intercept", "kappa_Intercept", "c_set_size1", "c_set_size2", "c_set_size3", "c_set_size4")
  )
  expect_equal(
    colnames(summary1$fixed),
    c("Estimate", "Est.Error", "l-95% CI", "u-95% CI", "Rhat", "Bulk_ESS", "Tail_ESS")
  )
  expect_output(print(summary1), "Constant Parameters:")
  expect_output(print(summary1), "Model: sdm")
  expect_output(print(summary1), "Links: mu = tan_half; c = log; kappa = log")
  expect_output(print(summary1), "Formula: mu = 0")
})

# minimal bmmsummary around a real model so the print helpers run
make_bmmsummary <- function(model, formula, fixed) {
  structure(
    list(
      fixed = fixed, random = list(), ngrps = list(),
      formula = formula, model = model,
      data = structure(data.frame(y = 0), data_name = "d"),
      iter = 100, warmup = 50, thin = 1, chains = 1,
      sampler = "NUTS", algorithm = "sampling"
    ),
    class = "bmmsummary"
  )
}

make_fixed <- function(rows) {
  n <- length(rows)
  data.frame(
    Estimate = seq_len(n), Est.Error = rep(0.1, n),
    "l-95% CI" = seq_len(n) - 1, "u-95% CI" = seq_len(n) + 1,
    Rhat = rep(1, n), Bulk_ESS = rep(500, n), Tail_ESS = rep(500, n),
    check.names = FALSE, row.names = rows
  )
}

test_that("print.bmmsummary selects rows by exact parameter prefix (#379)", {
  # imm has parameters a and kappa; "a_" is a substring of "kappa_", so an
  # unanchored grepl kept rows whose true parameter is not among those printed.
  model <- imm(
    resp_error = "y", nt_features = "nt", nt_distances = "d",
    set_size = "ss", version = "abc"
  )
  formula <- bmf(kappa ~ 1, a ~ 1, c ~ 1)
  fixed <- make_fixed(c("kappa_Intercept", "a_Intercept", "c_Intercept", "Xa_decoy"))
  out <- capture.output(print(make_bmmsummary(model, formula, fixed), color = FALSE))

  expect_true(any(grepl("kappa_Intercept", out)))
  expect_true(any(grepl("a_Intercept", out)))
  expect_false(any(grepl("Xa_decoy", out)))
})

test_that("print.bmmsummary handles a single regression coefficient row (#369)", {
  model <- sdm(resp_error = "y")
  formula <- bmf(c ~ 1, kappa ~ 1)
  fixed <- make_fixed("kappa_Intercept")
  out <- capture.output(print(make_bmmsummary(model, formula, fixed), color = FALSE))
  expect_true(any(grepl("kappa_Intercept", out)))
})

test_that(".summary_fixed_rows keeps a single fixed-effect row", {
  # gumbel-min sdt_ranking has one population coefficient (d_Intercept) but
  # two printed parameters (d, sdratio); the old sapply+apply errored here.
  one_row <- data.frame(Estimate = 0.6, Rhat = 1, row.names = "d_Intercept")
  out <- .summary_fixed_rows(one_row, c("d", "sdratio"))
  expect_s3_class(out, "data.frame")
  expect_identical(rownames(out), "d_Intercept")
})

test_that(".summary_fixed_rows selects all rows matching the printed parameters", {
  fixed <- data.frame(
    Estimate = 1:3, Rhat = c(1, NA, 1),
    row.names = c("d_Intercept", "sdratio_Intercept", "nuisance_Intercept")
  )
  out <- .summary_fixed_rows(fixed, c("d", "sdratio"))
  expect_identical(sort(rownames(out)), c("d_Intercept", "sdratio_Intercept"))
})
