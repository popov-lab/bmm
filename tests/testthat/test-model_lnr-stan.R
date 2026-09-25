# =============================================================================
# R vs Stan parity for the generated lnr likelihood
#
# cmdstanr's expose_functions()/compile_model_methods build fails on some
# toolchains (TBB symbol mismatch), so the generated function is evaluated
# through a generated quantities block under fixed_param instead.
# =============================================================================

skip_unless_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not(
    !is.null(tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)),
    "CmdStan is not installed"
  )
}

# the winner contributes its density and its own n - 1 survivors, every other
# category n_j survivors; this is the likelihood the model is defined by, not a
# transcription of either implementation
lnr_reference <- function(rt, response, m, s, ndt, counts) {
  t <- rt - ndt
  log(counts[response]) +
    stats::dlnorm(t, m[response], s, log = TRUE) +
    sum(vapply(seq_along(m), function(j) {
      reps <- if (j == response) counts[j] - 1 else counts[j]
      reps * stats::plnorm(t, m[j], s, lower.tail = FALSE, log.p = TRUE)
    }, numeric(1)))
}

lnr_stan_values <- function(grid_simple, grid_custom, cat_names) {
  code <- paste0(
    "functions {\n",
    bmm:::.lnr_stan_code("lnr_simple", c("correct", "error")), "\n",
    bmm:::.lnr_stan_code("lnr_custom", cat_names), "\n",
    "}\n",
    "data {\n",
    "  int S; vector[S] s_rt; vector[S] s_c; vector[S] s_e; vector[S] s_s;\n",
    "  vector[S] s_ndt; array[S] int s_resp; array[S] int s_n1; array[S] int s_n2;\n",
    "  int C; vector[C] c_rt; matrix[C, 4] c_m; vector[C] c_s; vector[C] c_ndt;\n",
    "  array[C] int c_resp; array[C, 4] int c_n;\n",
    "}\n",
    "generated quantities {\n",
    "  vector[S] lp_simple; vector[C] lp_custom;\n",
    "  for (i in 1:S) {\n",
    "    lp_simple[i] = lnr_simple_lpdf([s_rt[i]]' | [0]', [s_c[i]]', [s_e[i]]',\n",
    "      [s_ndt[i]]', [s_s[i]]', {s_resp[i]}, {s_n1[i]}, {s_n2[i]});\n",
    "  }\n",
    "  for (i in 1:C) {\n",
    "    lp_custom[i] = lnr_custom_lpdf([c_rt[i]]' | [0]', [c_m[i, 1]]',\n",
    "      [c_m[i, 2]]', [c_m[i, 3]]', [c_m[i, 4]]', [c_ndt[i]]', [c_s[i]]',\n",
    "      {c_resp[i]}, {c_n[i, 1]}, {c_n[i, 2]}, {c_n[i, 3]}, {c_n[i, 4]});\n",
    "  }\n",
    "}\n"
  )
  f <- file.path(tempdir(), "lnr_parity.stan")
  writeLines(code, f)
  mod <- cmdstanr::cmdstan_model(f)

  sdata <- list(
    S = nrow(grid_simple), s_rt = grid_simple$rt, s_c = grid_simple$correct,
    s_e = grid_simple$error, s_s = grid_simple$s, s_ndt = grid_simple$ndt,
    s_resp = as.integer(grid_simple$response),
    s_n1 = as.integer(grid_simple$n1), s_n2 = as.integer(grid_simple$n2),
    C = nrow(grid_custom),
    c_rt = grid_custom$rt,
    c_m = as.matrix(grid_custom[, paste0("m", 1:4)]),
    c_s = grid_custom$s, c_ndt = grid_custom$ndt,
    c_resp = as.integer(grid_custom$response),
    c_n = as.matrix(grid_custom[, paste0("n", 1:4)])
  )
  draws <- mod$sample(sdata, fixed_param = TRUE, chains = 1, iter_warmup = 0,
                      iter_sampling = 1, refresh = 0,
                      show_messages = FALSE)$draws(format = "df")
  list(
    simple = as.numeric(draws[1, paste0("lp_simple[", seq_len(nrow(grid_simple)), "]")]),
    custom = as.numeric(draws[1, paste0("lp_custom[", seq_len(nrow(grid_custom)), "]")])
  )
}

test_that("the generated Stan likelihood equals the R reference", {
  skip_unless_cmdstan()

  # decision times down to 1e-3 s and sdlogs down to 0.1 drive the losers'
  # survivors to log S well past -30, which is where a naive lccdf gives up
  grid_simple <- expand.grid(
    t = c(1e-3, 0.01, 0.1, 0.3, 1, 3, 10, 100),
    s = c(0.1, 0.3, 0.5, 1),
    response = 1:2,
    n2 = c(1L, 3L),
    ndt = 0.2
  )
  grid_simple$correct <- -1
  grid_simple$error <- 0
  grid_simple$n1 <- 1L
  grid_simple$rt <- grid_simple$t + grid_simple$ndt

  grid_custom <- expand.grid(
    t = c(1e-3, 0.01, 0.1, 0.3, 1, 3, 10, 100),
    s = c(0.1, 0.3, 0.5, 1),
    response = 1:4,
    ndt = 0.2
  )
  grid_custom$m1 <- -1.0
  grid_custom$m2 <- -0.4
  grid_custom$m3 <- 0.2
  grid_custom$m4 <- 0.8
  grid_custom$n1 <- 1L
  grid_custom$n2 <- 2L
  grid_custom$n3 <- 3L
  grid_custom$n4 <- 5L
  grid_custom$rt <- grid_custom$t + grid_custom$ndt

  ref_simple <- vapply(seq_len(nrow(grid_simple)), function(i) {
    lnr_reference(
      grid_simple$rt[i], grid_simple$response[i],
      c(grid_simple$correct[i], grid_simple$error[i]), grid_simple$s[i],
      grid_simple$ndt[i], c(grid_simple$n1[i], grid_simple$n2[i])
    )
  }, numeric(1))
  ref_custom <- vapply(seq_len(nrow(grid_custom)), function(i) {
    lnr_reference(
      grid_custom$rt[i], grid_custom$response[i],
      unlist(grid_custom[i, paste0("m", 1:4)]), grid_custom$s[i],
      grid_custom$ndt[i], unlist(grid_custom[i, paste0("n", 1:4)])
    )
  }, numeric(1))

  loser_log_s <- min(stats::plnorm(
    grid_simple$t, grid_simple$error, grid_simple$s,
    lower.tail = FALSE, log.p = TRUE
  ))
  expect_lt(loser_log_s, -30)

  stan <- suppressWarnings(
    lnr_stan_values(grid_simple, grid_custom, c("aa", "bb", "cc", "dd"))
  )

  # a log density below -700 underflows to zero probability, so no sampler
  # visits it; that is where Stan's lognormal_lccdf gives up and returns -Inf.
  # The comparison runs over the region a fit can reach, and the unreachable
  # rows are pinned separately so the boundary cannot drift inwards unnoticed.
  # The tolerance is set by the precision Stan writes its CSV with.
  for (part in c("simple", "custom")) {
    ref <- if (part == "simple") ref_simple else ref_custom
    got <- stan[[part]]
    reachable <- is.finite(ref) & ref > -700
    expect_gt(sum(reachable), 50)
    expect_true(all(is.finite(got[reachable])), label = part)
    expect_equal(got[reachable], ref[reachable], tolerance = 1e-7,
                 label = part)
    expect_true(all(ref[is.finite(ref) & !is.finite(got)] < -700),
                label = part)
  }
})

test_that("log_lik for the simple version carries the category constant", {
  skip_on_cran()
  path <- test_path("assets", "bmmfit_lnr_ppcheck.rds")
  skip_if_not(file.exists(path), "fixture not available")
  fit <- readRDS(path)

  ll <- brms::log_lik(fit, draw_ids = 1)[1, ]
  prep <- brms::prepare_predictions(fit, draw_ids = 1)
  counts <- c(prep$data$vint2[1], prep$data$vint3[1])

  reference <- vapply(seq_along(ll), function(i) {
    lnr_reference(
      prep$data$Y[i], prep$data$vint1[i],
      c(brms::get_dpar(prep, "correct", i = i),
        brms::get_dpar(prep, "error", i = i)),
      brms::get_dpar(prep, "s", i = i),
      brms::get_dpar(prep, "ndt", i = i),
      c(prep$data$vint2[i], prep$data$vint3[i])
    )
  }, numeric(1))
  expect_equal(ll, reference, tolerance = 1e-10)

  # an error trial is the probability that SOME error accumulator won, which
  # for K alternatives is log(K - 1) above the likelihood of one named error
  # accumulator winning while the others race. Constant in the parameters, so
  # posteriors are unaffected, but loo() and waic() are shifted by it; pinned
  # here so that a change of convention is deliberate (see ?lnr)
  expect_equal(counts[2], 2L)
  errors <- which(prep$data$vint1 == 2L)
  expect_gt(length(errors), 0)

  # built from first principles: one NAMED error accumulator wins while the
  # correct accumulator and the other error accumulator race and lose
  per_accumulator <- vapply(errors, function(i) {
    t <- prep$data$Y[i] - brms::get_dpar(prep, "ndt", i = i)
    sd_log <- brms::get_dpar(prep, "s", i = i)
    m_c <- brms::get_dpar(prep, "correct", i = i)
    m_e <- brms::get_dpar(prep, "error", i = i)
    stats::dlnorm(t, m_e, sd_log, log = TRUE) +
      stats::plnorm(t, m_c, sd_log, lower.tail = FALSE, log.p = TRUE) +
      stats::plnorm(t, m_e, sd_log, lower.tail = FALSE, log.p = TRUE)
  }, numeric(1))
  expect_equal(ll[errors] - per_accumulator, rep(log(2), length(errors)),
               tolerance = 1e-10)
})
