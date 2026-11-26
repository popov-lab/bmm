# Fits the throwaway extension in tests/local/bmmspike: a shifted log-normal
# response-time model that another package defines and registers with
# register_bmm_model() when it loads (#442). Nothing is installed; both
# packages load from source. Run by hand from the package root:
#   Rscript tests/local/check-extension.R
# The real fit compiles two Stan programs and samples each with four chains
# (about two minutes with CmdStan); it is skipped when CmdStan is not
# available. 2000 draws per chain, because with 1000 the participant SD's bulk
# ESS fell below 400 in one of the two fits. The script exits with an error when a test fails.
Sys.setenv(NOT_CRAN = "true")
suppressMessages(devtools::load_all(quiet = TRUE, export_all = FALSE)) # bmm's startup banner
suppressMessages(pkgload::load_all("tests/local/bmmspike", quiet = TRUE, export_all = FALSE)) # pkgload's dependency notes
library(testthat)
local_edition(3)

# 20 participants x 2 conditions x 40 trials, drawn with brms's own generator
# rather than the extension's rslnorm(), so a bug in the extension's
# distribution functions cannot cancel out between simulation and fit
spike_data <- function() {
  withr::local_seed(442)
  n_id <- 20
  u_id <- stats::rnorm(n_id, 0, 0.2)
  data <- expand.grid(trial = seq_len(40), cond = c("a", "b"), id = seq_len(n_id))
  mu <- -0.6 + 0.3 * (data$cond == "b") + u_id[data$id]
  data$rt <- brms::rshifted_lnorm(nrow(data), meanlog = mu, sdlog = 0.35, shift = 0.15)
  data$id <- factor(data$id)
  data
}

cmdstan_available <- function() {
  requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL))
}

data <- spike_data()
formula <- bmf(mu ~ 1 + cond + (1 | id), sigma ~ 1, ndt ~ 1)

test_that("loading the extension registers its model under its package", {
  expect_true("slnorm" %in% bmm_models())
  expect_identical(bmm:::model_package("slnorm"), "bmmspike")
  expect_identical(bmm:::registry_entry("slnorm")$version, utils::packageVersion("bmmspike"))
  expect_true("- slnorm() (bmmspike): Shifted log-normal response-time model (throwaway spike)" %in%
                capture.output(print(bmm_models())))
})

test_that("the extractors accept the extension's model", {
  code <- stancode(formula, data, slnorm("rt"))
  expect_match(code, paste0(" and bmmspike ", utils::packageVersion("bmmspike"), "\n"), fixed = TRUE)
  expect_match(code, "slnorm_lpdf", fixed = TRUE)
  prior <- default_prior(formula, data, slnorm("rt"))
  expect_true(any(prior$prior == "normal(-2, 1)" & prior$dpar == "ndt"))
  expect_equal(standata(formula, data, slnorm("rt"))$N, nrow(data))
})

test_that("a mock fit records the extension's package and version", {
  fit <- bmm(formula, data, slnorm("rt"), backend = "mock", mock_fit = 1, rename = FALSE)
  expect_s3_class(fit, "bmmfit")
  expect_identical(fit$bmm$package, "bmmspike")
  expect_identical(fit$version$bmmspike, utils::packageVersion("bmmspike"))
  expect_true(isTRUE(attr(fit, "slnorm_postprocess_ran")))
})

test_that("a cmdstanr fit of the extension's model matches brms's shifted_lognormal()", {
  skip_if_not(cmdstan_available(), "CmdStan is not available")
  fit <- bmm(formula, data, slnorm("rt"), backend = "cmdstanr",
             chains = 4, iter = 3000, warmup = 1000, cores = 4, seed = 1, refresh = 0)
  expect_identical(fit$bmm$package, "bmmspike")

  # the extension's default priors, written out for brms
  prior <- c(
    brms::prior(normal(-0.5, 1), class = "Intercept"),
    brms::prior(normal(0, 0.5), class = "b"),
    brms::prior(exponential(2), class = "sd"),
    brms::prior(normal(-1, 0.5), class = "Intercept", dpar = "sigma"),
    brms::prior(normal(-2, 1), class = "Intercept", dpar = "ndt")
  )
  reference <- brms::brm(
    brms::bf(rt ~ 1 + cond + (1 | id), sigma ~ 1, ndt ~ 1), data = data,
    family = brms::shifted_lognormal(), prior = prior, backend = "cmdstanr",
    chains = 4, iter = 3000, warmup = 1000, cores = 4, seed = 1, refresh = 0
  )

  variables <- c("b_Intercept", "b_condb", "sd_id__Intercept", "Intercept_sigma", "Intercept_ndt")
  summarise <- function(f) {
    posterior::summarise_draws(posterior::subset_draws(brms::as_draws_df(f), variable = variables))
  }
  ours <- summarise(fit)
  theirs <- summarise(reference)
  comparison <- data.frame(
    variable = ours$variable, bmm_mean = ours$mean, brms_mean = theirs$mean,
    bmm_sd = ours$sd, brms_sd = theirs$sd,
    difference_in_sd = (ours$mean - theirs$mean) / theirs$sd,
    bmm_rhat = ours$rhat, bmm_ess_bulk = ours$ess_bulk
  )
  print(comparison, digits = 3)

  # a comparison of means is only read once both fits have converged
  expect_true(all(c(ours$rhat, theirs$rhat) < 1.01))
  expect_true(all(c(ours$ess_bulk, theirs$ess_bulk) > 400))
  # the same likelihood and priors: means agree within Monte-Carlo error
  expect_true(all(abs(comparison$difference_in_sd) < 0.3))
})
