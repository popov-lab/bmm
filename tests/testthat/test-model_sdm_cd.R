sdm_cd_data <- function(n = 100, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(target = stats::runif(n, -pi, pi), probe = stats::runif(n, -pi, pi))
    dat$change <- rsdm_cd(n, wrap(dat$probe - dat$target), c = 5, kappa = 4)
    dat
  })
}

sdm_cd_model <- function() sdm_cd("change", probe = "probe", target = "target")

# P("same") from stats::integrate() and uniroot(), sharing no code with
# .sdm_cd_psame(): the arc where the SDM density centred on the probe exceeds
# exp(-criterion) / (2 pi), integrated under the density centred on mu
sdm_cd_reference <- function(probe, mu, c, kappa, criterion) {
  eta <- function(u) c * exp(kappa * (cos(u) - 1)) * sqrt(kappa / (2 * pi))
  peak <- eta(0)
  f <- function(x) exp(eta(x - mu) - peak)
  z <- stats::integrate(f, mu - pi, mu + pi, rel.tol = 1e-13, subdivisions = 5000)$value
  g <- function(d) exp(eta(d) - peak) / z - exp(-criterion) / (2 * pi)
  if (g(0) <= 0) return(0)
  if (g(pi) >= 0) return(1)
  hw <- stats::uniroot(g, c(0, pi), tol = 1e-14)$root
  cuts <- sort(unique(c(probe - hw, probe + hw, mu + c(-2, 0, 2) * pi)))
  cuts <- cuts[cuts >= probe - hw & cuts <= probe + hw]
  sum(vapply(seq_len(length(cuts) - 1), function(i) {
    stats::integrate(f, cuts[i], cuts[i + 1], rel.tol = 1e-13, subdivisions = 5000)$value
  }, numeric(1))) / z
}

test_that("sdm_cd() is a change-detection model independent of sdm()", {
  m <- sdm_cd_model()
  expect_s3_class(m, "change_detection")
  expect_equal(class(m), c("bmmodel", "change_detection", "sdm_cd", "sdm_cd_simple"))
  expect_false(inherits(m, "sdm"))
  expect_false(inherits(m, "circular"))
  expect_equal(m$fixed_parameters, list(mu = 0, criterion = 0))
  expect_equal(m$resp_vars, list(response = "change", probe = "probe", target = "target"))
  de <- sdm("y")
  expect_equal(class(de), c("bmmodel", "circular", "sdm", "sdm_simple"))
  expect_false(inherits(de, "change_detection"))
})

test_that("sdm_cd() names the parameters of sdm() plus criterion", {
  # joint fits of both tasks (#394) rely on the shared names
  expect_equal(names(sdm_cd_model()$parameters), c(names(sdm("y")$parameters), "criterion"))
  expect_equal(sdm_cd_model()$links, c(sdm("y")$links, list(criterion = "identity")))
})

test_that("sdm_cd() requires its three variables", {
  expect_error(sdm_cd("change", probe = "probe"), "target")
})

test_that("dsdm_cd() returns a probability mass function", {
  p <- dsdm_cd(c(0, 1), probe = 0.7, c = 5, kappa = 4)
  expect_equal(sum(p), 1)
  expect_equal(dsdm_cd(c(0, 1), probe = 0.7, c = 5, kappa = 4, log = TRUE), log(p))
})

test_that("dsdm_cd() and rsdm_cd() refuse invalid arguments", {
  expect_error(dsdm_cd(1, 0, kappa = -1), "kappa must be non-negative")
  expect_error(dsdm_cd(1, 0, c = -1), "c must be non-negative")
  expect_error(dsdm_cd(2, 0), "coded 0")
  expect_error(rsdm_cd(c(1, 2), 0), "single integer")
})

test_that("P('same') matches integration with stats::integrate()", {
  grid <- expand.grid(
    probe = c(0, 0.3, 1.5, 3), c = c(0.5, 5, 40), kappa = c(1, 8, 100),
    criterion = c(0, 0.5, -0.4)
  )
  ours <- .sdm_cd_psame(grid$probe, 0.2, grid$c, grid$kappa, grid$criterion)
  ref <- mapply(sdm_cd_reference, grid$probe, 0.2, grid$c, grid$kappa, grid$criterion)
  expect_lt(max(abs(ours - ref)), 1e-10)
})

test_that("the arc integral holds where one Gauss-Legendre rule fails", {
  reference <- function(lo, hi, c, kappa) {
    eta <- function(u) c * exp(kappa * (cos(u) - 1)) * sqrt(kappa / (2 * pi))
    pieces <- unique(sort(c(lo, 0, hi)))
    eta(0) + log(sum(vapply(seq_len(length(pieces) - 1L), function(i) {
      stats::integrate(function(x) exp(eta(x) - eta(0)), pieces[i], pieces[i + 1L],
        rel.tol = 1e-12, subdivisions = 5000
      )$value
    }, numeric(1))))
  }
  grid <- expand.grid(c = c(1, 10, 40), kappa = c(1, 20, 100), hw = c(0.2, 1.5, pi))
  ours <- .sdm_cd_log_int(-grid$hw, grid$hw, 0, grid$c, grid$kappa)
  ref <- mapply(function(c, k, h) reference(-h, h, c, k), grid$c, grid$kappa, grid$hw)
  expect_lt(max(abs(ours - ref)), 1e-8)
  # the normalising constant does not depend on where the circle starts
  expect_equal(.sdm_cd_log_int(-pi, pi, 0, 6, 12), .sdm_cd_log_int(-pi + 1.3, pi + 1.3, 1.3, 6, 12))
})

test_that("an identical probe is called 'same' far above chance", {
  expect_lt(dsdm_cd(1, probe = 0, c = 5, kappa = 4), 0.3)
  expect_gt(dsdm_cd(1, probe = pi, c = 5, kappa = 4), 0.9)
})

test_that("P('change') rises with the distance of the probe", {
  p <- dsdm_cd(1, seq(0, pi, length.out = 25), c = 5, kappa = 4)
  expect_true(all(diff(p) > 0))
})

test_that("stronger memory improves change detection", {
  expect_gt(dsdm_cd(1, pi / 2, c = 15, kappa = 8), dsdm_cd(1, pi / 2, c = 1, kappa = 8))
  expect_lt(dsdm_cd(1, 0, c = 15, kappa = 8), dsdm_cd(1, 0, c = 1, kappa = 8))
})

test_that("a larger criterion makes 'change' responses less likely", {
  p <- dsdm_cd(1, probe = 0.5, c = 4, kappa = 6, criterion = c(-0.5, 0, 0.5))
  expect_true(all(diff(p) < 0))
})

test_that("rsdm_cd() draws 'change' at the rate dsdm_cd() gives", {
  p <- withr::with_seed(4, mean(rsdm_cd(1e5, probe = 1, c = 5, kappa = 4)))
  expect_equal(p, dsdm_cd(1, 1, c = 5, kappa = 4), tolerance = 0.02)
})

test_that("the sdm_cd pipeline wires its Stan code and data", {
  dat <- sdm_cd_data()
  m <- sdm_cd_model()
  f <- bmf(c ~ 1, kappa ~ 1)
  expect_equal(bmm(f, dat, m, backend = "mock", mock_fit = 1, rename = FALSE)$fit, 1)

  sc <- stancode(f, dat, model = m)
  expect_match(
    sc,
    "sdm_cd_lpmf(Y | mu, c, kappa, criterion, vreal1, G_sdm_runs, sdm_run_start, sdm_run_count, cd_gl_x, cd_gl_w)",
    fixed = TRUE
  )
  expect_match(sc, "Intercept_criterion = 0", fixed = TRUE)
  sd <- standata(f, dat, model = m)
  expect_length(sd$cd_gl_x, 32)
  expect_equal(as.numeric(sd$vreal1), wrap(dat$probe - dat$target))
  expect_equal(c(sd$G_sdm_runs, sd$sdm_run_start, sd$sdm_run_count), c(1, 1, nrow(dat)))
})

test_that("sdm_cd() computes the normalising constant once per run of predictor values", {
  dat <- sdm_cd_data()
  dat$cond <- rep(c("a", "b"), each = 50)
  sd <- standata(bmf(c ~ 0 + cond, kappa ~ 1), dat, model = sdm_cd_model())
  expect_equal(sd$G_sdm_runs, 2)
  expect_equal(as.vector(sd$sdm_run_count), c(50, 50))
})

test_that("sdm_cd() refuses within-chain threading", {
  dat <- sdm_cd_data()
  expect_error(
    bmm(bmf(c ~ 1, kappa ~ 1), dat, sdm_cd_model(),
      threads = brms::threading(2), backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "does not support within-chain threading"
  )
  withr::with_options(list(brms.threads = 2), expect_error(
    stancode(bmf(c ~ 1, kappa ~ 1), dat, model = sdm_cd_model()),
    "does not support within-chain threading"
  ))
})

test_that("a fitted sdm_cd reports c on its natural scale", {
  fit <- bmm(bmf(c ~ 1, kappa ~ 1), sdm_cd_data(), sdm_cd_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_equal(fit$family$link_c, "log")
  expect_equal(fit$formula$family$link_c, "log")
  reverted <- revert_postprocess_brm(fit$bmm$model, fit)
  expect_equal(reverted$family$link_c, "identity")
  expect_equal(reverted$formula$family$link_c, "identity")
})

test_that("posterior_epred_sdm_cd() is P('change') and matches posterior_predict", {
  fit <- bmm(bmf(c ~ 1, kappa ~ 1), sdm_cd_data(), sdm_cd_model(),
    backend = "mock", mock_fit = 1, rename = FALSE
  )
  expect_identical(fit$formula$family$posterior_epred, posterior_epred_sdm_cd)

  sets <- data.frame(mu = 0, c = c(5, 2, 20), kappa = c(4, 10, 2), criterion = c(0, 0.5, -0.3))
  probe <- c(0.3, 1.5, 3)
  res <- epred_vs_predict(sets, posterior_epred_sdm_cd,
    function(i, prep) posterior_predict_sdm_cd(i, prep),
    data = list(vreal1 = probe), n_draws = 20000
  )
  expect_equal(res[, "epred"], dsdm_cd(1, probe, 0, sets$c, sets$kappa, sets$criterion))
  expect_true(all(abs(res[, "mc"] - res[, "epred"]) < 4 * sqrt(0.25 / 20000)))

  dpars <- list(
    mu = matrix(0, 3, 2), c = matrix(c(5, 6, 7, 2, 3, 4), 3),
    kappa = matrix(c(4, 5, 6, 9, 8, 7), 3), criterion = matrix(0.2, 3, 2)
  )
  expect_epred_by_cell(posterior_epred_sdm_cd, dpars, list(vreal1 = c(0.4, -2)))
})

test_that("log_lik_sdm_cd() scores the response with dsdm_cd()", {
  prep <- epred_prep(
    list(mu = matrix(0, 2, 1), c = matrix(c(5, 3), 2), kappa = matrix(c(4, 9), 2),
         criterion = matrix(c(0, 0.4), 2)),
    data = list(Y = 1, vreal1 = 0.8)
  )
  expect_equal(
    log_lik_sdm_cd(1, prep),
    dsdm_cd(1, 0.8, 0, c(5, 3), c(4, 9), c(0, 0.4), log = TRUE)
  )
})

# Evaluates sdm_cd_funs.stan through a fixed_param run. c reaches Stan on the
# log scale, as brms passes it.
sdm_cd_stan_values <- function(data) {
  sc_path <- system.file("stan_chunks", package = "bmm")
  code <- paste0(
    "functions {\n",
    read_lines2(file.path(sc_path, "cd_funs.stan")), "\n",
    read_lines2(file.path(sc_path, "sdm_cd_funs.stan")), "\n}\n",
    "data {
      int N; array[N] int y; vector[N] mu; vector[N] c; vector[N] kappa;
      vector[N] criterion; array[N] real probe; int G; array[G] int run_start;
      array[G] int run_count; vector[32] gl_x; vector[32] gl_w;
    }
    generated quantities {
      vector[N] log_z; vector[N] hw; vector[N] log_arc; real lpmf;
      for (n in 1:N) {
        log_z[n] = sdm_cd_log_int(-pi(), pi(), 0, c[n], kappa[n], gl_x, gl_w);
        hw[n] = sdm_cd_crit_angle(c[n], kappa[n], criterion[n], log_z[n]);
        log_arc[n] = sdm_cd_log_int(probe[n] - 0.7, probe[n] + 0.7, mu[n], c[n],
                                    kappa[n], gl_x, gl_w);
      }
      lpmf = sdm_cd_lpmf(y | mu, c, kappa, criterion, probe, G, run_start,
                         run_count, gl_x, gl_w);
    }"
  )
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = data, fixed_param = TRUE, iter_sampling = 1, chains = 1,
    refresh = 0, show_messages = FALSE, sig_figs = 18
  )
  function(variable) as.numeric(fit$draws(variable))
}

test_that("the Stan and R sdm_cd likelihoods agree", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the sdm_cd functions"
  )
  # three runs, each with its own c, kappa and criterion, covering a narrow
  # spike, a flat density, an arc that spans the whole circle and hw = 0
  par <- data.frame(
    c = c(5, 5, 5, 0.5, 0.5, 40, 40, 40, 3),
    kappa = c(4, 4, 4, 1, 1, 100, 100, 100, 0.2),
    criterion = c(0, 0, 0, 0.4, 0.4, -0.3, -0.3, -0.3, -3)
  )
  run_count <- c(3L, 2L, 3L, 1L)
  probe <- c(0, 1.2, -3, 0.4, 2.8, 0.05, -0.1, 2, 1)
  mu <- c(0, 0, 0, 0.3, 0.3, 0, 0, 0, 0)
  y <- c(0L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 0L)
  gl <- .cd_gl_rule()
  stan <- sdm_cd_stan_values(list(
    N = length(y), y = y, mu = mu, c = log(par$c), kappa = par$kappa,
    criterion = par$criterion, probe = probe, G = length(run_count),
    run_start = as.integer(cumsum(c(1, head(run_count, -1)))), run_count = run_count,
    gl_x = gl$x, gl_w = gl$w
  ))

  log_z <- .sdm_cd_log_int(-pi, pi, 0, par$c, par$kappa)
  expect_equal(stan("log_z"), log_z, tolerance = 1e-12)
  expect_equal(stan("hw"), .sdm_cd_crit_angle(par$c, par$kappa, par$criterion, log_z),
    tolerance = 1e-12
  )
  expect_equal(stan("log_arc"), .sdm_cd_log_int(probe - 0.7, probe + 0.7, mu, par$c, par$kappa),
    tolerance = 1e-12
  )
  expect_equal(
    stan("lpmf"),
    # the last row has hw = 0, where dsdm_cd() gives -Inf and the clamp a finite value
    sum(.cd_bernoulli_ld(y, .sdm_cd_psame(probe, mu, par$c, par$kappa, par$criterion))),
    tolerance = 1e-12
  )
})

test_that("a sampled sdm_cd fit returns c on its natural scale", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to fit the model"
  )
  dat <- sdm_cd_data(200)
  fit <- suppressWarnings(suppressMessages(bmm(
    bmf(c ~ 1, kappa ~ 1), dat, sdm_cd_model(),
    backend = "cmdstanr", chains = 1, iter = 60, warmup = 30, refresh = 0,
    seed = 3, rename = FALSE
  )))
  draws <- brms::as_draws_df(fit)
  c_nat <- brms::posterior_linpred(fit, dpar = "c", transform = TRUE)
  expect_equal(c_nat[, 1], exp(draws$b_c_Intercept), ignore_attr = TRUE)

  ll <- brms::log_lik(fit)
  kappa <- exp(draws$b_kappa_Intercept)
  expect_equal(
    ll[, 5],
    dsdm_cd(dat$change[5], wrap(dat$probe[5] - dat$target[5]), 0,
      exp(draws$b_c_Intercept), kappa, 0, log = TRUE
    ),
    ignore_attr = TRUE
  )
  # the expected response is P("change"), which log_lik gives for y = 1
  changed <- matrix(fit$data$change == 1, nrow(ll), ncol(ll), byrow = TRUE)
  expect_equal(
    brms::posterior_epred(fit),
    ifelse(changed, exp(ll), 1 - exp(ll)),
    ignore_attr = TRUE
  )
})
