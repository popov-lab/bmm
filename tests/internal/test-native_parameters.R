# Model-fitting integration tests for native_parameters().
# These fit real models and are too slow for tests/testthat. Run manually with:
#   devtools::load_all(); testthat::test_file("tests/internal/test-native_parameters.R")

library(bmm)
library(testthat)

fit_mixture3p <- function() {
  bmm(
    bmf(thetat ~ 1, thetant ~ 1, kappa ~ 1),
    data = subset(oberauer_lin_2017, ID %in% 1:4),
    model = mixture3p(
      resp_error = "dev_rad",
      nt_features = paste0("col_nt", 1:7),
      set_size = "set_size"
    ),
    chains = 1, iter = 400, warmup = 200, refresh = 0, silent = 2
  )
}

fit_ddm <- function() {
  data <- subset(
    data_color_judgement_task,
    ID %in% unique(data_color_judgement_task$ID)[1:3] & rt > 0.2 & rt < 3
  )
  data$resp <- as.integer(data$response_correct)
  bmm(
    bmf(drift ~ 1, bound ~ 1, ndt ~ 1),
    data = data,
    model = ddm(rt = "rt", response = "resp"),
    chains = 1, iter = 300, warmup = 150, refresh = 0, silent = 2
  )
}

test_that("native_parameters returns mixture3p weights as probabilities", {
  fit <- fit_mixture3p()

  native <- native_parameters(fit, re_formula = NA)
  weights <- native[native$parameter %in% c("thetat", "thetant"), ]
  expect_true(all(weights$value > 0 & weights$value < 1))

  # the guessing component is the implicit reference, so the two weights
  # returned must leave room for it
  totals <- tapply(weights$value, weights$.draw, sum)
  expect_true(all(totals < 1))

  sampling <- native_parameters(fit, re_formula = NA, scale = "sampling")
  reference <- bmm:::.np_softmax(list(
    thetat = matrix(sampling$value[sampling$parameter == "thetat"]),
    thetant = matrix(sampling$value[sampling$parameter == "thetant"])
  ))
  expect_equal(
    weights$value[weights$parameter == "thetat"],
    as.vector(reference$thetat)
  )
  expect_equal(
    native$value[native$parameter == "kappa"],
    exp(sampling$value[sampling$parameter == "kappa"])
  )
  expect_true(all(native$value[native$parameter == "mu1"] == 0))
})

test_that("native_parameters agrees with conditional_effects on the native scale", {
  fit <- fit_mixture3p()

  native <- native_parameters(fit, re_formula = NA, pars = "thetat")
  ce <- conditional_effects(fit, "thetat", scale = "native", re_formula = NA)
  expect_equal(mean(native$value), ce[[1]]$estimate__[1], tolerance = 1e-8)
})

test_that("native_parameters reports link-scale constants on the native scale", {
  fit <- fit_ddm()

  # zr is fixed at 0 on the logit scale, which is an unbiased starting point
  expect_equal(parameters(fit)$value[parameters(fit)$parameter == "zr"], "0")

  out <- native_parameters(fit, summary = TRUE)
  expect_equal(out$Estimate[out$parameter == "zr"], 0.5)
  expect_equal(out$Est.Error[out$parameter == "zr"], 0)
  expect_true(all(out$Estimate[out$parameter %in% c("bound", "ndt")] > 0))

  sampling <- native_parameters(fit, scale = "sampling")
  native <- native_parameters(fit)
  expect_equal(
    native$value[native$parameter == "bound"],
    exp(sampling$value[sampling$parameter == "bound"])
  )
})

test_that("native_parameters excludes internal variables for models with non-targets", {
  fit <- bmm(
    bmf(c ~ 1, a ~ 1, kappa ~ 1),
    data = subset(oberauer_lin_2017, ID %in% 1:4),
    model = imm(
      resp_error = "dev_rad",
      nt_features = paste0("col_nt", 1:7),
      set_size = "set_size",
      version = "abc"
    ),
    chains = 1, iter = 400, warmup = 200, refresh = 0, silent = 2
  )

  out <- native_parameters(fit)
  expect_setequal(unique(out$parameter), c("mu1", "kappa", "a", "c"))
  expect_named(out, c(".chain", ".iteration", ".draw", "parameter", "value"))

  sampling <- native_parameters(fit, scale = "sampling")
  expect_equal(
    out$value[out$parameter == "a"],
    exp(sampling$value[sampling$parameter == "a"])
  )
})
