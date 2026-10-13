latent_test_spec <- function(n = 12) {
  withr::local_seed(1)
  vwm <- data.frame(
    id = factor(rep(seq_len(n), each = 10)),
    error = runif(n * 10, -pi, pi),
    set_size = rep(c(2, 4), length.out = n * 10)
  )
  rt <- data.frame(id = factor(rep(seq_len(n), each = 5)), rt = rlnorm(n * 5))
  bmm_component(
    bmf(thetat ~ 1 + (1 | p | id), kappa ~ 1 + set_size + (1 + set_size | p | id)),
    model = mixture2p(resp_error = "error"), data = vwm
  ) +
    bmm_component(
      bmf(rt ~ 1 + (1 | p | id), sigma ~ 1 + (1 | p | id)),
      family = brms::lognormal(), data = rt
    )
}

latent_persons <- function(n = 12) {
  withr::local_seed(2)
  data.frame(id = seq_len(n), span = rnorm(n), age = rnorm(n, 40, 10))
}

test_that("latent_indicators() names random effects by response, parameter and coefficient", {
  ind <- latent_indicators(latent_test_spec())
  expect_setequal(
    ind$name,
    c("error_thetat", "error_kappa", "error_kappa_set_size", "rt_mu", "rt_sigma")
  )
  expect_equal(ind$link[ind$name == "error_kappa"], "log")
  expect_equal(ind$link[ind$name == "rt_sigma"], "log")
})

test_that("latent_indicators() lists the person-level variables of an attached structure", {
  spec <- latent_test_spec() +
    bmm_latent("WM =~ error_kappa + error_thetat + span", group = "id", data = latent_persons())
  ind <- latent_indicators(spec)
  expect_true(all(c("span", "age") %in% ind$name))
  expect_true(all(is.na(ind$component[ind$name %in% c("span", "age")])))
})

test_that("a latent structure attaches to the specification from either side of +", {
  latent <- bmm_latent("WM =~ error_kappa + error_thetat + rt_mu", group = "id")
  spec <- latent_test_spec() + latent
  expect_s3_class(attr(spec, "latent"), "bmm_latent")
  expect_length(spec, 2)

  comps <- latent_test_spec()
  spec2 <- latent + comps[[1]] + comps[[2]]
  expect_s3_class(attr(spec2, "latent"), "bmm_latent")
  expect_length(spec2, 2)

  expect_error(spec + latent, "only one latent structure")
})

test_that("names that are neither latent, random effects nor person data are rejected", {
  spec <- latent_test_spec() + bmm_latent("WM =~ error_kappa + error_thetat + rt_drift", group = "id")
  expect_error(configure_fit(spec, init = FALSE), "rt_drift.*latent_indicators")
})

test_that("only the terms of named parameters move into the uncorrelated latent block", {
  spec <- latent_test_spec() + bmm_latent("WM =~ error_kappa + error_thetat + rt_mu", group = "id")
  cfg <- suppressMessages(configure_fit(spec, init = FALSE))
  re <- brms:::frame_re(brms::brmsterms(cfg$config_args$formula), cfg$config_args$data)
  key <- paste(re$resp, re$nlpar, re$dpar, re$coef)

  latent_rows <- key %in% c(
    "error kappa  Intercept", "error kappa  set_size", "error thetat  Intercept", "rt   Intercept"
  )
  expect_true(all(!re$cor[latent_rows]))
  expect_length(unique(re$id[latent_rows]), 1)
  expect_false(re$id[key == "rt  sigma Intercept"] %in% re$id[latent_rows])
})

test_that("generated Stan code adds the latent structure to declared random-effect vectors", {
  spec <- latent_test_spec() + bmm_latent("
    WM =~ error_kappa + error_thetat + span
    speed =~ rt_mu + rt_sigma
    WM ~ age
    error_kappa_set_size ~ age
    error_kappa ~~ error_thetat
  ", group = "id", data = latent_persons())
  code <- suppressMessages(stancode(spec))
  updated <- regmatches(code, gregexpr("r_[0-9]+_[A-Za-z_]+_[0-9]+(?= \\+=)", code, perl = TRUE))[[1]]
  expect_length(unique(updated), 5)
  for (r in unique(updated)) {
    expect_match(code, paste0("vector\\[N_[0-9]+\\] ", r, ";"))
  }
  expect_match(code, "lv_eta\\[, 2\\] = lv_beta\\[1\\] \\* lv_X\\[, 1\\]")
  expect_match(code, "normal_lpdf\\(lv_Y\\[, 1\\]")
})

test_that("person-level data are aligned with the brms levels of the grouping variable", {
  persons <- latent_persons()
  shuffled <- persons[c(12:1), ]
  spec <- latent_test_spec() +
    bmm_latent("WM =~ error_kappa + error_thetat + span", group = "id", data = shuffled)
  cfg <- suppressMessages(configure_fit(spec, init = FALSE))
  latent <- attr(cfg$model, "latent")
  expect_equal(
    as.numeric(latent$person_values[, "span"]),
    as.numeric(scale(persons$span[as.integer(latent$levels)]))
  )
})

test_that("equality labels share one Stan parameter", {
  spec <- latent_test_spec() +
    bmm_latent("WM =~ error_kappa + a*error_thetat + a*rt_mu", group = "id")
  cfg <- suppressMessages(configure_fit(spec, init = FALSE))
  pt <- attr(cfg$model, "latent")$partable
  expect_equal(
    pt$stan[pt$rhs == "error_thetat" & pt$op == "=~"],
    pt$stan[pt$rhs == "rt_mu" & pt$op == "=~"]
  )
})

test_that("unsupported latent structures stop with an informative error", {
  higher_order <- latent_test_spec() + bmm_latent("
    A =~ error_kappa + error_thetat
    B =~ rt_mu + rt_sigma
    G =~ A + B
  ", group = "id")
  expect_error(configure_fit(higher_order, init = FALSE), "not supported yet")

  endo_cov <- latent_test_spec() + bmm_latent("
    A =~ error_kappa + error_thetat
    B =~ rt_mu + rt_sigma
    A ~ age
    A ~~ B
  ", group = "id", data = latent_persons())
  expect_error(configure_fit(endo_cov, init = FALSE), "endogenous")
})

test_that("a latent model runs through the full pipeline", {
  spec <- latent_test_spec() +
    bmm_latent("WM =~ error_kappa + error_thetat + rt_mu", group = "id")
  fit <- suppressMessages(bmm(spec, backend = "mock", mock_fit = 1, rename = FALSE))
  expect_equal(fit$bmm$latent$latents, "WM")
})
