# =============================================================================
# Tests for rdm() default priors: ndt intercept, group-level sd rates, and
# that both reach report_priors()/fit$prior and the generated Stan code.
# =============================================================================

test_that("rdm() default ndt prior is the ddm's normal(-1.5, 0.5)", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

  expect_equal(model$default_priors$ndt$main, "normal(-1.5, 0.5)")
})

test_that("rdm() default sd priors follow the racing-model rate", {
  model <- rdm(rt = "rt", response = "response", n_choices = 2)

  expect_equal(model$default_priors$driftc$sd, "exponential(2)")
  expect_equal(model$default_priors$drifte$sd, "exponential(2)")
  expect_equal(model$default_priors$gap$sd, "exponential(2)")
  expect_equal(model$default_priors$sp$sd, "exponential(2)")
  expect_equal(model$default_priors$ndt$sd, "exponential(4)")
  expect_equal(model$default_priors$s$sd, "exponential(4)")
})

test_that("custom version inherits the same sd rates for its accumulator parameters", {
  model <- rdm(rt = "rt", response = "resp", version = "custom")

  expect_equal(model$default_priors$gap$sd, "exponential(2)")
  expect_equal(model$default_priors$ndt$sd, "exponential(4)")
  expect_equal(model$default_priors$s$sd, "exponential(4)")
  expect_equal(model$default_priors$sp$sd, "exponential(2)")
})

test_that("a hierarchical mock fit reports the sd default priors in report_priors() and fit$prior", {
  skip_on_cran()

  withr::local_seed(1)
  dat <- do.call(rbind, lapply(1:5, function(i) {
    d <- rrdm(n = 30, drift = c(3, 1.5), gap = 1, ndt = 0.2)
    d$id <- i
    d
  }))
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(
    driftc ~ 1 + (1 | id), drifte ~ 1 + (1 | id),
    gap ~ 1 + (1 | id), ndt ~ 1 + (1 | id)
  )
  fit <- bmm(f, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)

  sd_rows <- fit$prior[fit$prior$class == "sd" & fit$prior$group == "", ]
  expect_equal(sd_rows$prior[sd_rows$dpar == "driftc"], "exponential(2)")
  expect_equal(sd_rows$prior[sd_rows$dpar == "drifte"], "exponential(2)")
  expect_equal(sd_rows$prior[sd_rows$dpar == "gap"], "exponential(2)")
  expect_equal(sd_rows$prior[sd_rows$dpar == "ndt"], "exponential(4)")

  out <- report_priors(fit)
  sd_out <- out[out$class == "sd", ]
  expect_equal(sd_out$prior[sd_out$parameter == "driftc"], "exponential(2)")
  expect_equal(sd_out$prior[sd_out$parameter == "drifte"], "exponential(2)")
  expect_equal(sd_out$prior[sd_out$parameter == "gap"], "exponential(2)")
  expect_equal(sd_out$prior[sd_out$parameter == "ndt"], "exponential(4)")
  expect_true(all(sd_out$source == "bmm default"))
})

test_that("the sd default priors reach the generated Stan code", {
  dat <- data.frame(
    rt = c(0.5, 0.6, 0.7, 0.8), response = c(1, 2, 1, 2),
    id = c(1, 1, 2, 2)
  )
  model <- rdm(rt = "rt", response = "response", n_choices = 2)
  f <- bmf(
    driftc ~ 1 + (1 | id), drifte ~ 1 + (1 | id),
    gap ~ 1 + (1 | id), ndt ~ 1 + (1 | id)
  )
  code <- suppressWarnings(stancode(f, dat, model))

  expect_true(grepl("lprior \\+= exponential_lpdf\\(sd_1 \\| 2\\);", code))
  expect_true(grepl("lprior \\+= exponential_lpdf\\(sd_2 \\| 2\\);", code))
  expect_true(grepl("lprior \\+= exponential_lpdf\\(sd_3 \\| 2\\);", code))
  expect_true(grepl("lprior \\+= exponential_lpdf\\(sd_4 \\| 4\\);", code))
})
