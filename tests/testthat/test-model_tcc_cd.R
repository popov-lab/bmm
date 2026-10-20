# Target confusability competition for single-probe change detection:
# constructor, check_data(), dtcc_cd() and rtcc_cd(), the generated Stan code,
# the post-processing functions, and the Stan likelihood against R

# exponential stand-in for a measured colour similarity, f(0) = 1, f(180 deg) = 0
tcc_cd_test_similarity <- function(n_steps = 360) {
  delta <- seq(0, 180, length.out = n_steps / 2 + 1)
  (exp(-delta / 20) - exp(-9)) / (1 - exp(-9))
}

tcc_cd_test_model <- function(...) {
  tcc_cd("change", "dist", similarity = tcc_cd_test_similarity(), ...)
}

# P("change") from the model's statement, without the package's binning or SDT
# helpers: the probe's familiarity is normal(d * f(probe), 1), and the observer
# says "same" when it exceeds the criterion. `deg` is the probe's distance from
# the target in whole degrees on a 1-degree similarity
ref_p_change <- function(deg, d, criterion, similarity = tcc_cd_test_similarity()) {
  pnorm(criterion - d * similarity[abs(deg) + 1])
}

tcc_cd_mock_fit <- function(model = tcc_cd_test_model(), formula = bmf(d ~ 1), ...) {
  dat <- data.frame(change = c(0, 1, 1, 0, 1, 0), dist = c(0, pi, 0.3, -0.3, 2, 0),
                    cond = rep(c("a", "b"), 3))
  suppressMessages(bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE, ...))
}

############################################################################# !
# CONSTRUCTOR                                                            ####
############################################################################# !

test_that("tcc_cd() builds an independent model with tcc()'s d and sdt_yn()'s criterion and bias", {
  model <- tcc_cd_test_model()
  expect_identical(class(model), c("bmmodel", "tcc_cd"))
  expect_named(model$parameters, c("d", "criterion", "bias"))
  expect_identical(model$alternative_parameters, list(c("criterion", "bias")))
  expect_identical(model$links, list(d = "log", criterion = "identity", bias = "identity"))
  # the priors of tcc()'s d and of sdt_yn()'s criterion and bias
  expect_identical(model$default_priors, list(
    d = list(main = "normal(0.5, 1)", effects = "normal(0, 1)", sd = "exponential(1)"),
    criterion = list(main = "normal(0.5, 1.5)", effects = "normal(0, 0.5)", sd = "exponential(1)"),
    bias = list(main = "normal(0, 1.5)", effects = "normal(0, 0.5)", sd = "exponential(2)")
  ))
  expect_identical(model$other_vars$units, "radians")
  expect_identical(tcc_cd_test_model(units = "bins")$other_vars$units, "bins")
  expect_warning(identity <- tcc_cd_test_model(links = list(d = "identity")), "allow values")
  expect_identical(identity$links$d, "identity")
})

test_that("tcc_cd() refuses a similarity vector the model cannot read", {
  f <- tcc_cd_test_similarity()
  expect_error(tcc_cd("change", "dist", similarity = c(0.9, f[-1])), "1 at distance 0")
  expect_error(tcc_cd("change", "dist", similarity = c(1, f[-c(1, 2)], 1.2)), "between 0 and 1")
  expect_error(tcc_cd("change", "dist", similarity = replace(f, 5, NA)), "between 0 and 1")
  expect_error(tcc_cd("change", "dist", similarity = as.character(f)), "numeric")
  expect_error(tcc_cd("change", "dist", similarity = c(1, 0.5)), "at least 3 values")
  expect_error(tcc_cd("change", "dist", similarity = f, units = "degrees"), "should be one of")
  expect_error(tcc_cd("change", "dist"), "missing.*similarity")
})

############################################################################# !
# CHECK_DATA                                                             ####
############################################################################# !

test_that("check_data() bins the probe distance to the similarity grid, either sign, +-pi together", {
  w <- 2 * pi / 360
  dat <- data.frame(change = 0, dist = c(0, pi, -pi, 10 * w, -10 * w, 0.4999 * w, 0.5001 * w,
                                         3 * pi / 2, NA))
  out <- check_data(tcc_cd_test_model(), dat, bmf(d ~ 1))
  expect_identical(out$tcc_cd_dist, c(0L, 180L, 180L, 10L, 10L, 0L, 1L, 90L, NA))

  # a 91-value similarity is on a grid of 180 steps, 2 degrees each
  model <- tcc_cd("change", "step", similarity = tcc_cd_test_similarity(180), units = "bins")
  out <- check_data(model, data.frame(change = 1, step = c(0, 3, -7, 90, 91, 181)), bmf(d ~ 1))
  expect_identical(out$tcc_cd_dist, c(0L, 3L, 7L, 90L, 89L, 1L))
})

test_that("check_data() codes the response 1 = change and accepts a logical one", {
  dat <- data.frame(change = c(TRUE, FALSE, NA), dist = c(0, pi, 1))
  out <- check_data(tcc_cd_test_model(), dat, bmf(d ~ 1))
  expect_identical(out$change, c(1L, 0L, NA))
})

test_that("check_data() refuses data the likelihood cannot read", {
  model <- tcc_cd_test_model()
  dat <- data.frame(change = c(0, 1), dist = c(0, pi))
  expect_error(check_data(model, dat[, "dist", drop = FALSE], bmf(d ~ 1)),
               "response variable 'change' is not present")
  expect_error(check_data(model, transform(dat, change = c(0, 2)), bmf(d ~ 1)),
               "coded 0 = 'same' and 1 = 'change'")
  expect_error(check_data(model, transform(dat, change = c("same", "change")), bmf(d ~ 1)),
               "coded 0 = 'same' and 1 = 'change'")
  expect_error(check_data(model, dat[, "change", drop = FALSE], bmf(d ~ 1)),
               "distance variable 'dist' is not present")
  expect_error(check_data(model, transform(dat, dist = c("0", "pi")), bmf(d ~ 1)),
               "'dist' must be numeric")
  expect_error(check_data(model, transform(dat, dist = c(0, 180)), bmf(d ~ 1)),
               "beyond 2\\*pi")
  bins <- tcc_cd_test_model(units = "bins")
  expect_error(check_data(bins, transform(dat, dist = c(0, 1.5)), bmf(d ~ 1)), "whole numbers")
})

test_that("check_data() warns when every probe is equally similar to its target", {
  model <- tcc_cd_test_model()
  # one similarity gives one rate for d and the criterion together
  expect_warning(check_data(model, data.frame(change = c(0, 1), dist = c(0, 0)), bmf(d ~ 1)),
                 "only one similarity")
  expect_silent(check_data(model, data.frame(change = c(0, 1), dist = c(0, pi)), bmf(d ~ 1)))
})

test_that("check_data() warns before overwriting its reserved column", {
  dat <- data.frame(change = c(0, 1), dist = c(0, pi), tcc_cd_dist = 3)
  expect_warning(out <- check_data(tcc_cd_test_model(), dat, bmf(d ~ 1)), "reserved")
  expect_identical(out$tcc_cd_dist, c(0L, 180L))
})

############################################################################# !
# DISTRIBUTION FUNCTIONS                                                 ####
############################################################################# !

test_that("dtcc_cd() is the probability that the probe's familiarity falls short of the criterion", {
  grid <- expand.grid(deg = c(0, 1, 7, 30, 90, 179, 180), d = c(0.2, 1.5, 6),
                      criterion = c(-1, 0.4, 2.5))
  rad <- grid$deg * pi / 180
  f <- tcc_cd_test_similarity()
  p <- ref_p_change(grid$deg, grid$d, grid$criterion)
  expect_equal(dtcc_cd(1, rad, grid$d, f, criterion = grid$criterion), p, tolerance = 1e-12)
  expect_equal(dtcc_cd(0, rad, grid$d, f, criterion = grid$criterion), 1 - p, tolerance = 1e-12)
  expect_equal(dtcc_cd(1, -rad, grid$d, f, criterion = grid$criterion, log = TRUE), log(p),
               tolerance = 1e-12)
  # bias is the criterion measured from the midpoint between 0 and d
  expect_equal(dtcc_cd(1, rad, grid$d, f, bias = grid$criterion),
               ref_p_change(grid$deg, grid$d, grid$criterion + grid$d / 2), tolerance = 1e-12)
})

test_that("dtcc_cd() stays finite on the log scale far into either tail", {
  f <- tcc_cd_test_similarity()
  expect_equal(dtcc_cd(1, 0, 40, f, criterion = 0, log = TRUE), pnorm(-40, log.p = TRUE))
  expect_equal(dtcc_cd(0, pi, 1, f, criterion = 40, log = TRUE), pnorm(-40, log.p = TRUE))
})

test_that("at similarities 0 and 1, dtcc_cd() is dsdt_yn() on single trials", {
  # the yes/no SDT kernel with 'same' as 'old': a probe at the target is a
  # signal trial, a probe opposite it (f = 0) a noise trial
  f <- tcc_cd_test_similarity()
  grid <- expand.grid(stimulus = c(0L, 1L), same = c(0L, 1L), d = c(0.3, 2), crit = c(-0.5, 1.2))
  dist <- ifelse(grid$stimulus == 1L, 0, pi)
  expect_equal(dtcc_cd(1 - grid$same, dist, grid$d, f, criterion = grid$crit, log = TRUE),
               dsdt_yn(grid$same, 1, grid$stimulus, grid$d, criterion = grid$crit, log = TRUE),
               tolerance = 1e-12)
  expect_equal(dtcc_cd(1 - grid$same, dist, grid$d, f, bias = grid$crit, log = TRUE),
               dsdt_yn(grid$same, 1, grid$stimulus, grid$d, bias = grid$crit, log = TRUE),
               tolerance = 1e-12)
})

test_that("rtcc_cd() matches a simulation of the probe's familiarity against the criterion", {
  withr::local_seed(566)
  f <- tcc_cd_test_similarity()
  n <- 40000
  for (case in list(c(deg = 0, d = 2, crit = 0.8), c(deg = 12, d = 2, crit = 0.8),
                    c(deg = 45, d = 3.5, crit = 0.2), c(deg = 180, d = 1, crit = -0.4))) {
    sampled <- rtcc_cd(n, case[["deg"]] * pi / 180, case[["d"]], f, criterion = case[["crit"]])
    expect_true(all(sampled %in% 0:1))
    familiarity <- stats::rnorm(n, case[["d"]] * f[case[["deg"]] + 1], 1)
    simulated <- mean(familiarity <= case[["crit"]])
    se <- sqrt(2 * simulated * (1 - simulated) / n)
    expect_lt(abs(mean(sampled) - simulated), 4 * se + 1e-9)
  }
})

test_that("dtcc_cd() and rtcc_cd() recycle their arguments and refuse what they cannot read", {
  f <- tcc_cd_test_similarity()
  expect_length(dtcc_cd(c(0, 1, 1), 0, 1, f, criterion = 0.5), 3)
  expect_length(rtcc_cd(5, c(0, pi), 1, f, criterion = 0.5), 5)
  expect_identical(dtcc_cd(numeric(0), 0, 1, f, criterion = 0.5), numeric(0))
  expect_identical(is.na(dtcc_cd(c(1, NA, 0), 0, 1, f, criterion = 0.5)), c(FALSE, TRUE, FALSE))
  expect_error(dtcc_cd(2, 0, 1, f, criterion = 0.5), "0 = 'same' and 1 = 'change'")
  expect_error(dtcc_cd(1, 0, 1, f), "exactly one of 'criterion'")
  expect_error(dtcc_cd(1, 0, 1, f, criterion = 0, bias = 0), "exactly one of 'criterion'")
  expect_error(rtcc_cd(5, 0, 1, f), "exactly one of 'criterion'")
  expect_error(dtcc_cd(1, 0, 1, f[-1], criterion = 0), "1 at distance 0")
  expect_error(rtcc_cd(0, 0, 1, f, criterion = 0), "single positive integer")
})

############################################################################# !
# CONFIGURE_MODEL AND STAN CODE                                          ####
############################################################################# !

test_that("the likelihood reads the probe's similarity through its binned distance", {
  code <- brms::stancode(tcc_cd_mock_fit())
  expect_match(code, "target += tcc_cd_lpmf(Y[n] | mu[n], d[n], criterion[n], vint1[n], tcc_cd_similarity);",
               fixed = TRUE)
  expect_match(code, "vector[181] tcc_cd_similarity;", fixed = TRUE)
  expect_match(code, "d = exp(d);", fixed = TRUE)
  # the yes/no SDT distribution functions are pasted in, not rewritten
  expect_match(code, "real sdt_log_cumprob(real eta, int dist_type)", fixed = TRUE)
  expect_match(code, "real tcc_cd_lpmf(int y, real mu, real d, real criterion, int dist,",
               fixed = TRUE)
  expect_identical(as.numeric(brms::standata(tcc_cd_mock_fit())$tcc_cd_similarity),
                   tcc_cd_test_similarity())
})

test_that("a formula with bias selects the bias family, and one with both is refused", {
  code <- brms::stancode(tcc_cd_mock_fit(formula = bmf(d ~ cond, bias ~ 1)))
  expect_match(code, "target += tcc_cd_bias_lpmf(Y[n] | mu[n], d[n], bias[n], vint1[n], tcc_cd_similarity);",
               fixed = TRUE)
  expect_error(tcc_cd_mock_fit(formula = bmf(d ~ 1, criterion ~ 1, bias ~ 1)),
               "gives both 'criterion', 'bias'")
})

test_that("threading passes the similarity to the sliced likelihood as data", {
  code <- brms::stancode(tcc_cd_mock_fit(threads = brms::threading(2)))
  expect_match(code, "reduce_sum(partial_log_lik_lpmf", fixed = TRUE)
  expect_match(code, "tcc_cd_lpmf(Y[nn] | mu[n], d[n], criterion[n], vint1[nn], tcc_cd_similarity)",
               fixed = TRUE)
  expect_match(code, "data vector tcc_cd_similarity", fixed = TRUE)
})

############################################################################# !
# LOG_LIK, POSTERIOR_PREDICT AND POSTERIOR_EPRED                         ####
############################################################################# !

tcc_cd_prep <- function(dpars, change, dist, similarity = tcc_cd_test_similarity()) {
  epred_prep(dpars, list(Y = change, vint1 = dist, tcc_cd_similarity = similarity))
}

test_that("log_lik_tcc_cd() is dtcc_cd() at each draw, with criterion or bias", {
  d <- matrix(c(0.5, 2, 4.5), 3, 4)
  crit <- matrix(c(-0.2, 0.7, 1.9), 3, 4)
  change <- c(1L, 0L, 1L, 0L)
  dist <- c(0L, 15L, 90L, 180L)
  f <- tcc_cd_test_similarity()
  for (anchor in c("criterion", "bias")) {
    prep <- tcc_cd_prep(stats::setNames(list(d, crit), c("d", anchor)), change, dist)
    for (i in 1:4) {
      args <- stats::setNames(list(crit[, i]), anchor)
      expected <- do.call(dtcc_cd, c(list(change[i], dist[i] * pi / 180, d[, i], f, log = TRUE), args))
      expect_equal(log_lik_tcc_cd(i, prep), expected, tolerance = 1e-12)
    }
  }
})

test_that("posterior_predict_tcc_cd() answers 'change' with the model's probability", {
  withr::local_seed(566)
  n_draws <- 40000
  prep <- tcc_cd_prep(list(d = matrix(2, n_draws, 2), criterion = matrix(0.9, n_draws, 2)),
                      c(0L, 1L), c(0L, 20L))
  for (i in 1:2) {
    draws <- posterior_predict_tcc_cd(i, prep)
    p <- ref_p_change(c(0, 20)[i], 2, 0.9)
    expect_true(all(draws %in% 0:1))
    expect_lt(abs(mean(draws) - p), 4 * sqrt(p * (1 - p) / n_draws))
  }
})

test_that("posterior_epred_tcc_cd() is the probability of 'change', observation by observation", {
  d <- matrix(c(0.5, 2, 4.5, 1, 3, 0.2), 2, 3)
  crit <- matrix(c(-0.2, 0.7, 1.9, 0.1, 2, 0), 2, 3)
  dist <- c(0L, 15L, 180L)
  out <- posterior_epred_tcc_cd(tcc_cd_prep(list(d = d, criterion = crit), c(0L, 1L, 0L), dist))
  # draws in rows, observations in columns, each column at its own distance
  deg <- matrix(dist[col(d)], 2, 3)
  expect_equal(out, ref_p_change(deg, d, crit), tolerance = 1e-12)
  out_bias <- posterior_epred_tcc_cd(tcc_cd_prep(list(d = d, bias = crit - d / 2), c(0L, 1L, 0L), dist))
  expect_equal(out_bias, out, tolerance = 1e-12)
})

test_that("a tcc_cd fit stores its post-processing functions", {
  family <- tcc_cd_mock_fit()$family
  expect_identical(family$log_lik, log_lik_tcc_cd)
  expect_identical(family$posterior_predict, posterior_predict_tcc_cd)
  expect_identical(family$posterior_epred, posterior_epred_tcc_cd)
})

test_that("the stored frame rebuilds the probe distance at the centre of its bin", {
  model <- tcc_cd_test_model()
  stored <- data.frame(change = c(0L, 1L), tcc_cd_dist = c(0L, 29L))
  out <- revert_check_data(model, stored)
  expect_equal(out$dist, c(0, 29 * 2 * pi / 360))
  expect_false("tcc_cd_dist" %in% colnames(out))
  expect_identical(attr(out, "rebuilt"), "dist")
  bins <- revert_check_data(tcc_cd_test_model(units = "bins"), stored)
  expect_identical(bins$dist, c(0L, 29L))
})

############################################################################# !
# STAN AGAINST R                                                         ####
############################################################################# !

# the two lpmfs in a fixed_param program over the stanvars bmm generates, as
# the other parity tests do (expose_functions() fails inside the suite)
tcc_cd_stan_lp <- function(model, rows) {
  stanvars <- bmm:::.tcc_cd_stanvars(model)
  funs <- paste(vapply(stanvars, function(s) if (s$block == "functions") s$scode else "", ""),
                collapse = "\n")
  data_decl <- paste(vapply(stanvars, function(s) if (s$block == "data") s$scode else "", ""),
                     collapse = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n", data_decl, "\n  int<lower=1> N;\n  array[N] int y;\n  array[N] int dist;\n",
    "  vector[N] d;\n  vector[N] crit;\n}\n",
    "generated quantities {\n",
    "  vector[N] lp_crit;\n  vector[N] lp_bias;\n",
    "  for (i in 1:N) {\n",
    "    lp_crit[i] = tcc_cd_lpmf(y[i] | 0, d[i], crit[i], dist[i], tcc_cd_similarity);\n",
    "    lp_bias[i] = tcc_cd_bias_lpmf(y[i] | 0, d[i], crit[i], dist[i], tcc_cd_similarity);\n",
    "  }\n",
    "}\n"
  )
  sdata <- lapply(stanvars, function(s) s$sdata)
  sdata <- sdata[!vapply(sdata, is.null, TRUE)]
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(sdata, list(N = nrow(rows), y = rows$y, dist = rows$dist, d = rows$d, crit = rows$crit)),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  idx <- seq_len(nrow(rows))
  list(crit = as.numeric(csv[1, paste0("lp_crit.", idx)]),
       bias = as.numeric(csv[1, paste0("lp_bias.", idx)]))
}

test_that("the Stan likelihood matches dtcc_cd() with criterion or bias, into the tails", {
  skip_on_cran()
  skip_unless_stan_tests()
  model <- tcc_cd_test_model()
  f <- tcc_cd_test_similarity()
  # every distance class (target, near, far, opposite), both responses, and
  # rows whose familiarity lies 30 to 40 SDs from the criterion
  rows <- data.frame(
    y = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    dist = c(0L, 0L, 3L, 17L, 60L, 120L, 180L, 180L, 0L, 180L),
    d = c(1.7, 0.4, 2.2, 3.1, 5, 0.9, 1.2, 2.5, 40, 0.5),
    crit = c(0.8, -0.3, 1.1, 0.2, 2.4, 1.5, -0.6, 0.9, 0, 30)
  )
  out <- tcc_cd_stan_lp(model, rows)
  rad <- rows$dist * pi / 180
  r_crit <- dtcc_cd(rows$y, rad, rows$d, f, criterion = rows$crit, log = TRUE)
  r_bias <- dtcc_cd(rows$y, rad, rows$d, f, bias = rows$crit, log = TRUE)
  expect_equal(out$crit, r_crit, tolerance = 1e-10)
  expect_equal(out$bias, r_bias, tolerance = 1e-10)
  prep <- tcc_cd_prep(list(d = matrix(rows$d, 1), criterion = matrix(rows$crit, 1)), rows$y, rows$dist)
  expect_equal(vapply(seq_len(nrow(rows)), log_lik_tcc_cd, numeric(1), prep = prep), out$crit,
               tolerance = 1e-10)
})
