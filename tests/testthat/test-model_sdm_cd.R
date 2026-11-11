# P("same") from stats::integrate() and uniroot(), sharing no code with
# .sdm_cd_psame(): the arc where the log-likelihood ratio of "change" stays
# at or below the criterion, bounded by tanh to the trial's range of evidence,
# integrated under the SDM density centred on mu
sdm_cd_reference <- function(probe, mu, c, kappa, criterion) {
  eta <- function(u) c * sqrt(kappa / (2 * pi)) * exp(kappa * (cos(u) - 1))
  peak <- eta(0)
  f <- function(x) exp(eta(x - mu) - peak)
  width <- 1 / sqrt(kappa * max(peak, 1))
  splits <- function(lo, hi) {
    at <- outer(c(-2, 0, 2) * pi + mu, c(-16, -4, -1, 0, 1, 4, 16) * width, "+")
    sort(unique(c(lo, at[at > lo & at < hi], hi)))
  }
  area <- function(lo, hi) {
    cuts <- splits(lo, hi)
    sum(vapply(seq_len(length(cuts) - 1), function(i) {
      stats::integrate(f, cuts[i], cuts[i + 1], rel.tol = 1e-13, subdivisions = 5000)$value
    }, numeric(1)))
  }
  log_z <- peak + log(area(mu - pi, mu + pi))
  llr <- function(d) -(log(2 * pi) + eta(d) - log_z)
  bounded <- if (criterion > 0) {
    llr(pi) * tanh(criterion / llr(pi))
  } else if (criterion < 0) {
    -llr(0) * tanh(criterion / -llr(0))
  } else {
    0
  }
  g <- function(d) llr(d) - bounded
  if (g(0) > 0) return(0)
  if (g(pi) <= 0) return(1)
  hw <- stats::uniroot(g, c(0, pi), tol = 1e-14)$root
  area(probe - hw, probe + hw) * exp(peak - log_z)
}

test_that("sdm_cd() validates its arguments", {
  expect_error(sdm_cd("r", "p"), "target")
  expect_error(sdm_cd("r", "p", "t", version = "full"), "should be")
  expect_error(sdm_cd("r", "p", "t", links = list(c = "identity")), "cannot be changed")
  # arguments copied from the circular mixtures are reported, not silently dropped
  expect_warning(sdm_cd("r", "p", "t", knowledge = "rich"), "knowledge")
  expect_warning(sdm_cd("r", "p", "t", variable_precision = TRUE), "variable_precision")
})

test_that("sdm_cd() takes sdm's parameters, priors, links and inits, plus criterion", {
  cd <- sdm_cd("r", "p", "t")
  de <- sdm("y")
  # a joint fit of both tasks relies on the shared names
  expect_identical(names(cd$parameters), c(names(de$parameters), "criterion"))
  expect_identical(cd$links, c(de$links, list(criterion = "identity")))
  expect_identical(cd$default_priors[names(de$default_priors)], de$default_priors)
  expect_identical(cd$init_ranges[names(de$init_ranges)], de$init_ranges)
  expect_identical(cd$fixed_parameters, list(mu = 0, criterion = 0))
  expect_s3_class(cd, c("bmmodel", "change_detection", "sdm_cd", "sdm_cd_simple"), exact = TRUE)
  expect_false(inherits(cd, "circular"))
  expect_null(cd$knowledge)
  expect_identical(cd$task, "Detection, recognition and confidence judgments")
  expect_identical(cd$domain, "Working memory")
  expect_true("sdm_cd" %in% bmm_models())
})

test_that("criterion's init range starts on the linear part of its bound", {
  # at the corners of the inits of c and kappa, the ends of criterion's init
  # range keep more than half of their value
  init <- sdm_cd("r", "p", "t")$init_ranges
  corners <- expand.grid(c = init$c, kappa = init$kappa, criterion = init$criterion)
  log_z <- .sdm_cd_log_int(-pi, pi, 0, corners$c, corners$kappa)
  bounded <- .sdm_cd_bounded_criterion(corners$criterion, corners$c, corners$kappa, log_z)
  expect_gt(min(bounded / corners$criterion), 0.5)
})

test_that("dsdm_cd() and rsdm_cd() refuse what sdm's functions refuse", {
  expect_error(dsdm_cd(1, 0, kappa = -1), "kappa must be non-negative")
  expect_error(dsdm_cd(1, 0, c = -1), "c must be non-negative")
  expect_error(rsdm_cd(5, 0, c = -1), "c must be non-negative")
  expect_error(dsdm_cd(2, 0), "response must be 0")
  expect_error(dsdm_cd(NA, 0), "response must be 0")
  expect_error(rsdm_cd(c(1, 2), 0), "single positive whole number")
  expect_error(rsdm_cd(2.5, 0), "single positive whole number")
  expect_error(rsdm_cd(0, 0), "single positive whole number")
})

test_that("dsdm_cd() integrates the SDM density over the arc that integrate() finds", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  # a flat memory, the default inits and a sharp one, against both sides of
  # the criterion; probes inside, across and opposite the peak
  grid <- expand.grid(
    probe = c(0, 0.3, 1.5, 3), c = c(0.5, 5, 40), kappa = c(1, 8, 100),
    criterion = c(0, 0.5, -0.4)
  )
  ours <- dsdm_cd(0, grid$probe, 0.2, grid$c, grid$kappa, grid$criterion)
  ref <- mapply(sdm_cd_reference, grid$probe, 0.2, grid$c, grid$kappa, grid$criterion)
  expect_lt(max(abs(ours - ref)), 1e-9)
  expect_equal(
    dsdm_cd(c(0, 1), 0.4, mu = 0.3, c = 5, kappa = 4, log = TRUE),
    log(dsdm_cd(c(0, 1), 0.4, mu = 0.3, c = 5, kappa = 4))
  )
})

test_that("the normalising constant holds at the corner cells of the sampled range (E5)", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  # integrate() on [0, pi] split at multiples of the spike width; the rule's
  # error peaks at log c 4.5-5 and kappa 100-300 (1.6e-6 at log c 5, kappa 200)
  reference <- function(lc, kappa) {
    log_peak <- lc + 0.5 * log(kappa / (2 * pi))
    peak <- exp(log_peak)
    w <- min(pi, sqrt(2 * max(log_peak, 1) / kappa))
    breaks <- sort(unique(c(0, pmin(pi, w * c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16)), pi)))
    f <- function(x) exp(peak * expm1(kappa * (cos(x) - 1)))
    peak + log(2 * sum(vapply(seq_len(length(breaks) - 1), function(i) {
      stats::integrate(f, breaks[i], breaks[i + 1], rel.tol = 1e-13, abs.tol = 0,
                       subdivisions = 5000L)$value
    }, numeric(1))))
  }
  corners <- expand.grid(lc = c(4.5, 5), kappa = c(100, 200, 300))
  ours <- .sdm_cd_log_int(-pi, pi, 0, exp(corners$lc), corners$kappa)
  expect_lt(max(abs(ours - mapply(reference, corners$lc, corners$kappa))), 5e-6)
  inner <- expand.grid(lc = c(-1, 1, 3), kappa = c(0.5, 5, 30))
  ours <- .sdm_cd_log_int(-pi, pi, 0, exp(inner$lc), inner$kappa)
  expect_lt(max(abs(ours - mapply(reference, inner$lc, inner$kappa))), 1e-9)
  # the constant does not depend on where the circle starts
  expect_equal(.sdm_cd_log_int(-pi, pi, 0, 6, 12), .sdm_cd_log_int(-pi + 1.3, pi + 1.3, 1.3, 6, 12))
})

test_that("P('change') rises with the probe's distance and falls with the criterion", {
  p <- dsdm_cd(1, seq(0, pi, length.out = 25), c = 5, kappa = 4)
  expect_true(all(diff(p) > 0))
  expect_lt(p[1], 0.3)
  expect_gt(p[25], 0.9)
  expect_true(all(diff(dsdm_cd(1, 0.5, c = 4, kappa = 6, criterion = c(-0.5, 0, 0.5))) < 0))
  # a stronger memory detects a change better and false-alarms less
  expect_gt(dsdm_cd(1, pi / 2, c = 15, kappa = 8), dsdm_cd(1, pi / 2, c = 1, kappa = 8))
  expect_lt(dsdm_cd(1, 0, c = 15, kappa = 8), dsdm_cd(1, 0, c = 1, kappa = 8))
})

test_that("rsdm_cd() draws the observer that dsdm_cd() integrates", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(136)
  n <- 1e5
  cells <- expand.grid(c = c(1.5, 5), kappa = c(2, 8), criterion = c(0, 0.6, -0.4))
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    expected <- dsdm_cd(1, 0.9, mu = 0.1, c = cell$c, kappa = cell$kappa, criterion = cell$criterion)
    drawn <- mean(rsdm_cd(n, 0.9, mu = 0.1, c = cell$c, kappa = cell$kappa, criterion = cell$criterion))
    expect_lt(abs(drawn - expected), 5 * sqrt(expected * (1 - expected) / n), label = i)
  }
})

test_that("check_data() refuses missing values in every column the model reads", {
  data <- sdm_cd_data()
  model <- sdm_cd("change", "probe", "target")
  formula <- bmf(c ~ 1, kappa ~ 1)
  for (column in c("change", "probe", "target")) {
    broken <- data
    broken[[column]][3] <- NA
    expect_error(check_data(model, broken, formula), "missing", label = column)
  }
  checked <- check_data(model, data, formula)
  expect_equal(checked$probe_centered, wrap(data$probe - data$target))
  expect_identical(built_data_columns(model), "probe_centered")
})

test_that("an estimated criterion is announced with the bound", {
  data <- sdm_cd_data()
  formula <- bmf(c ~ 1, kappa ~ 1, criterion ~ 1)
  model <- check_model(sdm_cd("change", "probe", "target"), data, formula)
  expect_message(check_formula(model, data, formula), "Decision criterion")
  fixed <- bmf(c ~ 1, kappa ~ 1)
  model <- check_model(sdm_cd("change", "probe", "target"), data, fixed)
  expect_no_message(check_formula(model, data, fixed), message = "criterion is estimated")
})

test_that("the likelihood is one call over the runs, with c through its log link", {
  code <- stancode(bmf(c ~ 1, kappa ~ 1), sdm_cd_data(), sdm_cd("change", "probe", "target"))
  expect_match(
    code,
    paste0("sdm_cd_simple_lpmf(Y | mu, c, kappa, criterion, vreal1, G_sdm_runs, ",
           "sdm_run_start, sdm_run_count, cd_gl_x, cd_gl_w)"),
    fixed = TRUE
  )
  expect_match(code, "c = exp(c);", fixed = TRUE)
  expect_match(code, "Intercept_criterion = 0", fixed = TRUE)
  expect_match(code, "sdm_run_total != N", fixed = TRUE)
})

test_that("the runs are the cells of every predictor, criterion's included", {
  data <- data.frame(
    change = rep(c(0L, 1L, 1L), 10), target = 0, probe = rep(c(0, 1.5, -2), 10),
    cond = rep(c("a", "b", "c"), c(12, 10, 8))
  )
  model <- sdm_cd("change", "probe", "target")
  runs <- function(formula, data) {
    sdata <- suppressMessages(standata(formula, data, model))
    c(G = sdata$G_sdm_runs, start = sdata$sdm_run_start, count = sdata$sdm_run_count)
  }
  expect_equal(runs(bmf(c ~ 1, kappa ~ 1), data), c(G = 1, start = 1, count = 30))
  expect_equal(
    runs(bmf(c ~ cond, kappa ~ 1), data),
    c(G = 3, start1 = 1, start2 = 13, start3 = 23, count1 = 12, count2 = 10, count3 = 8)
  )
  expect_equal(
    runs(bmf(c ~ 1, kappa ~ 1, criterion ~ cond), data),
    c(G = 3, start1 = 1, start2 = 13, start3 = 23, count1 = 12, count2 = 10, count3 = 8)
  )
  # unsorted data are fitted as they are, one run per block of equal values,
  # and the user is told that sorting would be faster
  interleaved <- data[c(1:6, 13:18, 7:12, 19:30), ]
  messages <- capture_messages(check_data(model, interleaved, bmf(c ~ cond, kappa ~ 1)))
  expect_match(paste(messages, collapse = ""), "not ordered by predictors")
  expect_equal(
    runs(bmf(c ~ cond, kappa ~ 1), interleaved),
    c(G = 5, start1 = 1, start2 = 7, start3 = 13, start4 = 19, start5 = 23,
      count1 = 6, count2 = 6, count3 = 6, count4 = 4, count5 = 8)
  )
})

test_that("within-chain threading is refused", {
  # brms would slice Y and the dpars per thread but pass the runs whole
  model <- sdm_cd("change", "probe", "target")
  formula <- bmf(c ~ 1, kappa ~ 1)
  expect_error(
    bmm(formula, sdm_cd_data(), model, backend = "mock", mock_fit = 1,
        threads = brms::threading(2)),
    "does not support within-chain threading"
  )
  withr::local_options(brms.threads = 2)
  expect_error(
    bmm(formula, sdm_cd_data(), model, backend = "mock", mock_fit = 1),
    "does not support within-chain threading"
  )
})

test_that("log_lik(), posterior_predict() and posterior_epred() evaluate the same P('same')", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(552)
  n <- 4e4
  draws <- list(mu = 0.1, c = c(4, 1.5), kappa = c(3, 9), criterion = c(0, 0.3))
  dpars <- lapply(draws, function(x) matrix(rep_len(x, n), n, 1))
  prep <- structure(
    list(dpars = dpars, data = list(Y = 1L, vreal1 = 0.8), ndraws = n, nobs = 1),
    class = "brmsprep"
  )
  p_change <- dsdm_cd(1, 0.8, mu = 0.1, c = c(4, 1.5), kappa = c(3, 9), criterion = c(0, 0.3))
  expect_equal(log_lik_sdm_cd_simple(1, prep)[1:2], log(p_change), tolerance = 1e-12)
  expect_equal(as.vector(posterior_epred_sdm_cd_simple(prep))[1:2], p_change, tolerance = 1e-12)
  drawn <- posterior_predict_sdm_cd_simple(1, prep)
  expect_lt(abs(mean(drawn) - mean(p_change)), 0.015)
})

test_that("simulate_data() draws the observer's decision on the design's probes", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(649)
  # c e^1.5, kappa e, the probe 1 from the target: P("change") from the
  # integrate() reference
  design <- data.frame(target = 2, probe = 3)[rep(1, 1e5), ]
  model <- sdm_cd("change", "probe", "target")
  data <- simulate_data(model, list(c = 1.5, kappa = 1), design)
  expected <- 1 - sdm_cd_reference(1, 0, exp(1.5), exp(1), 0)
  expect_lt(abs(mean(data$change) - expected), 5 * sqrt(0.25 / nrow(data)))
  expect_error(simulate_data(model, list(c = 1.5, kappa = 1), n_trials = 5), "probe")
})

test_that("the likelihood computes log Z once per run and matches R row by row", {
  skip_unless_stan_tests()
  skip_on_cran()
  # three runs of lengths 3, 1 and 2, c and kappa constant within a run and
  # distinct between runs, as check_data() lays them out. Every parameter
  # takes distinct values, so swapping c and kappa, or mu and criterion, in
  # the signature changes every row. Within run 1 the rows differ in mu,
  # criterion and probe, which share no log Z with c and kappa. The length-1
  # run 2 sits between runs whose c and kappa differ from its own, so log Z
  # read from a run start off by one, or an inner loop one row too short or
  # too long, changes run 2 (and run 1 or 3). Run 3 holds the corner of E5
  # (log c 5, kappa 200).
  cases <- data.frame(
    y = c(0, 1, 1, 0, 1, 0),
    mu = c(0.1, -0.2, 0.3, 0.05, 0, -0.1),
    c = c(3.5, 3.5, 3.5, 0.5, exp(5), exp(5)),
    kappa = c(6, 6, 6, 1, 200, 200),
    criterion = c(0.35, -0.25, 0, 0.4, 0, -0.3),
    probe = c(0.9, 2.8, -0.4, 1.7, 0.15, 0.02)
  )
  runs <- list(G = 3L, start = c(1L, 4L, 5L), count = c(3L, 1L, 2L))
  chunks <- c("circmix_functions.stan", "cd_functions.stan", "circmix_cd_functions.stan",
              "sdm_cd_functions.stan")
  program <- paste0(
    "functions {\n",
    paste(vapply(chunks, function(f) read_lines2(system.file("stan_chunks", f, package = "bmm")), ""),
          collapse = "\n"),
    "}\n",
    "data { int N; array[N] int y; vector[N] mu; vector[N] c; vector[N] kappa;",
    " vector[N] criterion; array[N] real probe; int G; array[G] int run_start;",
    " array[G] int run_count; vector[32] gl_x; vector[32] gl_w; }\n",
    "generated quantities {\n",
    "  real total = sdm_cd_simple_lpmf(y | mu, c, kappa, criterion, probe, G, run_start,",
    " run_count, gl_x, gl_w);\n",
    "  vector[G] per_run;\n",
    "  vector[G] log_z;\n",
    "  for (g in 1:G) {\n",
    "    int s = run_start[g];\n",
    "    int e = s + run_count[g] - 1;\n",
    "    per_run[g] = sdm_cd_simple_lpmf(y[s:e] | mu[s:e], c[s:e], kappa[s:e], criterion[s:e],",
    " probe[s:e], 1, {1}, {run_count[g]}, gl_x, gl_w);\n",
    "    log_z[g] = sdm_cd_log_int(-pi(), pi(), 0, c[s], kappa[s], gl_x, gl_w);\n",
    "  }\n",
    "}\n"
  )
  gl <- .cd_gl_rule()
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(
      list(N = nrow(cases), G = runs$G, run_start = runs$start, run_count = runs$count,
           gl_x = gl$x, gl_w = gl$w),
      as.list(cases)
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  r_rows <- with(cases, dsdm_cd(y, probe, mu, c, kappa, criterion, log = TRUE))
  r_runs <- vapply(seq_len(runs$G), function(g) {
    sum(r_rows[runs$start[g] + seq_len(runs$count[g]) - 1])
  }, numeric(1))
  expect_equal(as.numeric(fit$draws("total")), sum(r_rows), tolerance = 1e-10)
  expect_equal(as.numeric(fit$draws("per_run")), r_runs, tolerance = 1e-10)
  expect_equal(
    as.numeric(fit$draws("log_z")),
    .sdm_cd_log_int(-pi, pi, 0, cases$c[runs$start], cases$kappa[runs$start]),
    tolerance = 1e-10
  )
})
