test_that("bmm version is added to mock model", {
  fit <- bmm(
    formula = bmmformula(c ~ 1, kappa ~ 1),
    data = data.frame(y = rsdm(n = 10)),
    model = sdm(resp_error = "y"),
    backend = "mock",
    rename = F,
    mock = 1
  )
  expect_true("bmm" %in% names(fit$version))
})


test_that("get_mu_pars works", {
  a <- brms::brm(y ~ a, data.frame(y = c(1, 2, 3), a = c("A", "B", "C")),
    backend = "mock", mock_fit = 1, rename = F
  )
  mus <- get_mu_pars(a)
  expect_equal(mus, c("Intercept", "aB", "aC"))
})


# posterior_epred (#475) ------------------------------------------------------

test_that("posterior_epred_undefined() builds a refusal that names the model", {
  refuse <- posterior_epred_undefined("demo")
  expect_identical(names(formals(refuse))[1], "prep")
  expect_error(refuse(NULL), "not defined for the demo model; use native_parameters()")
})

test_that("posterior_epred() refuses the circular mixture models by name", {
  fit <- readRDS(test_path("assets/mock_bmmfit_mixture2p.rds"))
  expect_error(brms::posterior_epred(fit),
               "The expected response is not defined for the mixture2p model")

  skip_on_cran()
  dat <- oberauer_lin_2017
  fit <- bmm(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1), dat,
             mixture3p("dev_rad", nt_features = paste0("col_nt", 1:7), set_size = "set_size"),
             backend = "mock", mock = 1, rename = FALSE)
  expect_error(brms::posterior_epred(fit), "not defined for the mixture3p model")
  fit <- bmm(bmf(kappa ~ 1, c ~ 1, s ~ 1), dat,
             imm("dev_rad", nt_features = paste0("col_nt", 1:7),
                 nt_distances = paste0("dist_nt", 1:7), set_size = "set_size",
                 version = "bsc"),
             backend = "mock", mock = 1, rename = FALSE)
  expect_error(brms::posterior_epred(fit), "not defined for the imm model")
})

test_that("posterior_epred() of a circular model still returns its parameters", {
  fit <- load_fixture_fit("bmmfit_example1.rds")
  class(fit$bmm$model) <- c("bmmodel", "circular", "sdm", "sdm_simple")
  expect_error(brms::posterior_epred(fit, ndraws = 5), "not defined for the sdm model")
  expect_equal(dim(brms::posterior_epred(fit, dpar = "kappa", ndraws = 5)),
               c(5L, nrow(fit$data)))
})
