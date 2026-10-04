imm_cd_data <- function(n = 60) {
  withr::with_seed(136, {
    dat <- data.frame(
      target = runif(n, -pi, pi), probe = runif(n, -pi, pi),
      nt1 = runif(n, -pi, pi), nt2 = runif(n, -pi, pi),
      d1 = runif(n, 0.2, 2), d2 = runif(n, 0.2, 2),
      ss = rep(1:3, length.out = n), resp = rbinom(n, 1, 0.5)
    )
  })
  dat[dat$ss < 3, c("nt2", "d2")] <- NA
  dat[dat$ss < 2, c("nt1", "d1")] <- NA
  dat
}

imm_cd_model <- function(version = "full", ...) {
  imm_cd("resp", "probe", "target",
    nt_features = c("nt1", "nt2"),
    nt_distances = if (version != "abc") c("d1", "d2"),
    set_size = "ss", version = version, ...
  )
}

imm_cd_formula <- function(version = "full") {
  switch(version,
    full = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
    bsc = bmf(kappa ~ 1, c ~ 1, s ~ 1),
    abc = bmf(kappa ~ 1, c ~ 1, a ~ 1)
  )
}

test_that("imm_cd() has imm()'s parameters plus the criterion", {
  # joint fits of both tasks (#394) pair the parameters by name
  for (version in c("full", "bsc", "abc")) {
    for (vp in c(FALSE, TRUE)) {
      de <- imm("y", "nt", if (version != "abc") "d", "ss",
        version = version,
        variable_precision = vp
      )
      cd <- imm_cd_model(version, variable_precision = vp)
      expect_identical(names(cd$parameters), c(names(de$parameters), "criterion"))
      expect_identical(cd$links[names(de$links)], de$links)
      expect_identical(cd$default_priors[names(de$default_priors)], de$default_priors)
    }
  }
})

test_that("imm_cd() is its own model, not a version of imm()", {
  model <- imm_cd_model()
  expect_s3_class(model, c("bmmodel", "change_detection", "non_targets", "imm_cd", "imm_cd_full"),
    exact = TRUE
  )
  expect_equal(model$fixed_parameters, list(mu = 0, b = 0, criterion = 0))
  expect_equal(model$knowledge, "limited")
  expect_equal(imm_cd_model(knowledge = "rich")$knowledge, "rich")
  expect_error(imm_cd_model(version = "xyz"), "should be one of")
  expect_error(imm_cd_model(knowledge = "full"), "should be one of")
  expect_error(imm_cd_model(vp_nodes = 21), "at least 41")
  expect_null(imm_cd_model("abc")$other_vars$nt_distances)
  expect_equal(model_versions("imm_cd"), c("full", "bsc", "abc"))
})

test_that("imm_cd() cites Lin and Oberauer (2022), and van den Berg et al. with variable precision", {
  expect_length(model_citation(imm_cd_model()), 2)
  expect_match(model_citation(imm_cd_model())[2], "10.1016/j.cogpsych.2022.101463", fixed = TRUE)
  expect_length(model_citation(imm_cd_model(variable_precision = TRUE)), 3)
})

# ---- the decision rule on the imm retrieval mixture -------------------------

test_that("dimm_cd() is a probability mass function", {
  args <- list(probe = 0.7, mu = c(0, 1.5, -2), dist = c(0, 0.5, 1.2), c = 5, a = 0.6, s = 1.35, kappa = 6)
  for (tau in c(0, 1)) {
    p <- do.call(dimm_cd, c(list(x = c(0, 1), tau = tau), args))
    expect_equal(sum(p), 1)
    expect_true(all(p > 0 & p < 1))
  }
  expect_equal(
    do.call(dimm_cd, c(list(x = 1, log = TRUE), args)),
    log(do.call(dimm_cd, c(list(x = 1), args)))
  )
})

test_that("an identical probe is rejected far above chance and a distant one is detected", {
  args <- list(mu = c(0, 2, -2.5), dist = c(0, 1.5, 2), kappa = 8, c = 20, a = 0.3, s = 1.5)
  expect_lt(do.call(dimm_cd, c(list(x = 1, probe = 0), args)), 0.35)
  expect_gt(do.call(dimm_cd, c(list(x = 1, probe = pi), args)), 0.85)
})

test_that("nearer non-targets interfere more", {
  args <- list(x = 1, probe = 2, mu = c(0, 2, -2.5), kappa = 8, c = 6, a = 0.4, s = 1.2)
  near <- do.call(dimm_cd, c(args, list(dist = c(0, 0.1, 2))))
  far <- do.call(dimm_cd, c(args, list(dist = c(0, 3, 2))))
  # a probe matching a nearby non-target is harder to call a change
  expect_lt(near, far)
})

test_that("stronger context activation improves change detection", {
  args <- list(x = 1, probe = pi / 2, mu = c(0, 2, -2.5), dist = c(0, 1, 1.5), kappa = 8, a = 0.4, s = 1.2)
  expect_gt(do.call(dimm_cd, c(args, c = 20)), do.call(dimm_cd, c(args, c = 0.5)))
})

test_that("at criterion 0, P('same') is the retrieval mass inside the closed-form arc", {
  # the observer's prior of storage cancels at criterion 0 (Lin & Oberauer,
  # 2022, Appendix B), so the boundary depends on kappa alone
  probe <- 0.9
  activation <- c(5 + 0.6, 5 * exp(-1.35 * 0.5) + 0.6, 1)
  w <- activation / sum(activation)
  hw <- acos(.circmix_log_besselI0(6) / 6)
  vm <- function(x, mu) exp(6 * cos(x - mu) - .circmix_log_besselI0(6)) / (2 * pi)
  mass <- vapply(c(0, 1.5), function(mu) {
    stats::integrate(vm, probe - hw, probe + hw, mu = mu, rel.tol = 1e-12)$value
  }, numeric(1))
  expect_equal(
    dimm_cd(0, probe, mu = c(0, 1.5), dist = c(0, 0.5), c = 5, a = 0.6, s = 1.35, kappa = 6),
    sum(w[1:2] * mass) + w[3] * hw / pi,
    tolerance = 1e-8
  )
})

test_that("the version weights reproduce what dimm_cd() computes", {
  nt <- c(2, -1.5, 1)
  dist <- c(0.5, 1.5, 2)
  probe <- seq(-pi, pi, length.out = 9)
  for (tau in c(0, 0.8)) {
    for (rich in c(FALSE, TRUE)) {
      expect_equal(
        .imm_cd_psame(probe, 0, 6,
          c = 5, a = 2, s = 2, b = 1, criterion = 0.3, set_size = 4,
          nt = nt, dist = dist, tau = tau, nodes = 41L, rich = rich, version = "full"
        ),
        dimm_cd(0, probe,
          mu = c(0, nt), dist = c(0, dist), c = 5, a = 2, b = 1, s = 2,
          kappa = 6, criterion = 0.3, tau = tau, knowledge = if (rich) "rich" else "limited"
        ),
        tolerance = 1e-12
      )
    }
  }
})

test_that("the reduced versions drop the term they are named for", {
  nt <- c(2, -1.5, 1)
  dist <- c(0.5, 1.5, 2)
  probe <- seq(-pi, pi, length.out = 9)
  psame <- function(version, a, s, ss = 4) {
    .imm_cd_psame(probe, 0, 6,
      c = 5, a = a, s = s, b = 1, criterion = 0, set_size = ss,
      nt = nt, dist = dist, tau = 0, nodes = 41L, rich = FALSE, version = version
    )
  }
  expect_equal(psame("bsc", NULL, 2), psame("full", 1e-12, 2), tolerance = 1e-9)
  # with s -> infinity, every non-target keeps only its cue-independent a
  expect_equal(psame("abc", 2, NULL), psame("full", 2, 1e6), tolerance = 1e-9)
  # a set size of one leaves the target, weighted c + a, and the background
  expect_equal(
    psame("full", 2, 2, ss = 1),
    .circmix_cd_vp_psame(
      matrix(probe), matrix(log(7 / 8), length(probe)), rep(log(1 / 8), length(probe)),
      rep(6, length(probe)), rep(0, length(probe)), 0, 7 / 8, FALSE
    ),
    tolerance = 1e-12
  )
})

test_that("rimm_cd() draws from dimm_cd()", {
  withr::local_seed(2022)
  args <- list(probe = 1.6, mu = c(0, 2, -1.5), dist = c(0, 0.5, 2), c = 4, a = 0.5, s = 2, kappa = 8)
  for (tau in c(0, 1)) {
    draws <- do.call(rimm_cd, c(list(n = 4e4, tau = tau), args))
    expect_setequal(unique(draws), c(0, 1))
    # binomial SE at 4e4 draws is below 0.0025
    expect_equal(mean(draws), do.call(dimm_cd, c(list(x = 1, tau = tau), args)), tolerance = 0.01)
  }
})

# ---- the bmm pipeline ---------------------------------------------------------

test_that("every version runs through the bmm() pipeline", {
  dat <- imm_cd_data()
  for (version in c("full", "bsc", "abc")) {
    expect_silent(bmm(imm_cd_formula(version), dat, imm_cd_model(version),
      backend = "mock", mock_fit = 1, rename = FALSE
    ))
  }
  expect_silent(bmm(bmf(kappa ~ 1, tau ~ 1, c ~ 1, a ~ 1, s ~ 1, criterion ~ 1), dat,
    imm_cd_model(variable_precision = TRUE, knowledge = "rich"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
})

test_that("the likelihood receives the centred probe, the set size and the non-targets", {
  dat <- imm_cd_data()
  code <- stancode(imm_cd_formula(), dat, imm_cd_model())
  expect_match(code, paste0(
    "imm_cd_full_lpmf(Y[n] | mu[n], kappa[n], c[n], a[n], s[n], b[n], criterion[n], ",
    "vint1[n], vreal1[n], vreal2[n], vreal3[n], vreal4[n], vreal5[n]"
  ), fixed = TRUE)
  # criterion is fixed by a constant(0) intercept, which reaches Stan as exactly 0
  expect_match(code, "Intercept_criterion = 0", fixed = TRUE)
  # without variable precision tau is the literal 0, and knowledge "limited" is 0
  expect_match(code, "return imm_cd_full_core(y, mu, kappa, 0.0, c, a, s, b, criterion, ss, probe, to_vector({nt1, nt2}), to_vector({dist1, dist2}), 41, 0,", fixed = TRUE)
  expect_match(
    stancode(imm_cd_formula("abc"), dat, imm_cd_model("abc", knowledge = "rich")),
    "to_vector({nt1, nt2}), 41, 1,", fixed = TRUE
  )

  sd <- standata(imm_cd_formula(), dat, imm_cd_model())
  expect_length(sd$cd_gl_x, 32)
  expect_equal(as.numeric(sd$vreal1), wrap(dat$probe - dat$target))
  expect_equal(as.numeric(sd$vint1), dat$ss)
  expect_equal(as.numeric(sd$vreal2), ifelse(is.na(dat$nt1), 0, dat$nt1))
  expect_equal(as.numeric(sd$vreal4), ifelse(is.na(dat$d1), 999, dat$d1))
})

test_that("the generated Stan code of every version parses", {
  skip_on_cran()
  skip_if_not(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr is required")
  dat <- imm_cd_data()
  for (version in c("full", "bsc", "abc")) {
    code <- suppressMessages(stancode(
      imm_cd_formula(version), dat,
      imm_cd_model(version, variable_precision = version == "bsc")
    ))
    model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
    expect_true(model$check_syntax(quiet = TRUE))
  }
})

test_that("a set-size-1 factor level pins the parameters that trial cannot inform", {
  dat <- imm_cd_data()
  dat$ss_f <- factor(dat$ss)
  model <- imm_cd("resp", "probe", "target", c("nt1", "nt2"), c("d1", "d2"), "ss_f")
  pr <- default_prior(bmf(kappa ~ 1, c ~ 1, a ~ 0 + ss_f, s ~ 0 + ss_f), dat, model)
  pinned <- pr[pr$class == "b" & pr$coef == "ss_f1", ]
  expect_setequal(pinned$dpar, c("a", "s"))
  expect_true(all(pinned$prior == "constant(0)"))
  expect_error(
    bmm(bmf(kappa ~ 1, c ~ 1, a ~ 1 + ss_f, s ~ 0 + ss_f), dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "contains \\s*an intercept"
  )
})

test_that("check_data() validates the non-target distances", {
  dat <- imm_cd_data()
  dat$d1[dat$ss == 3] <- -1
  expect_error(
    check_data(imm_cd_model(), dat, imm_cd_formula()),
    "distances to the target need to be postive"
  )
})

# ---- post-processing ------------------------------------------------------------

imm_cd_prep_data <- list(
  Y = c(0, 1, 1), vint1 = c(1L, 2L, 3L), vreal1 = c(0.2, -1.1, 2.4),
  vreal2 = c(0, 1.2, -2), vreal3 = c(0, 0, 0.7),
  vreal4 = c(999, 0.4, 1.5), vreal5 = c(999, 999, 0.6)
)

imm_cd_prep_sets <- data.frame(
  mu = 0, kappa = c(4, 8, 15), c = c(3, 5, 2), a = c(0.3, 0.6, 1),
  s = c(2, 1, 0.5), b = 1, criterion = c(0, 0.3, -0.2)
)

test_that("log_lik_imm_cd_full() reads the probe, set size and non-targets of its row", {
  # one response against a vector of draws; cd-core recycles y in .cd_bernoulli_ld()
  skip_if(
    length(.cd_bernoulli_ld(1, c(0.5, 0.5))) == 1,
    "needs the cd-core fix that recycles y in .cd_bernoulli_ld()"
  )
  prep <- epred_prep(lapply(imm_cd_prep_sets, matrix, nrow = 2, ncol = 3, byrow = TRUE),
    data = imm_cd_prep_data
  )
  d <- imm_cd_prep_data
  for (i in 1:3) {
    ss <- d$vint1[i]
    nt <- c(d$vreal2[i], d$vreal3[i])[seq_len(ss - 1)]
    dist <- c(d$vreal4[i], d$vreal5[i])[seq_len(ss - 1)]
    par <- imm_cd_prep_sets[i, ]
    expected <- dimm_cd(d$Y[i], d$vreal1[i],
      mu = c(0, nt), dist = c(0, dist), c = par$c, a = par$a, b = par$b,
      s = par$s, kappa = par$kappa, criterion = par$criterion, log = TRUE
    )
    expect_equal(log_lik_imm_cd_full(i, prep), rep(expected, 2), tolerance = 1e-12)
    p_change <- posterior_epred_imm_cd_full(prep)[, i]
    expect_equal(p_change, rep(dimm_cd(1, d$vreal1[i],
      mu = c(0, nt), dist = c(0, dist), c = par$c, a = par$a, b = par$b,
      s = par$s, kappa = par$kappa, criterion = par$criterion
    ), 2), tolerance = 1e-12)
  }
})

test_that("posterior_epred_imm_cd_full() is the rate posterior_predict_imm_cd_full() simulates", {
  skip_on_cran()
  withr::local_seed(475)
  res <- epred_vs_predict(imm_cd_prep_sets, posterior_epred_imm_cd_full,
    posterior_predict_imm_cd_full,
    data = imm_cd_prep_data
  )
  expect_lt(max(abs(res[, "mc"] - res[, "epred"])), 0.015)
})

test_that("posterior_epred_imm_cd_full() returns one column per observation", {
  dpars <- lapply(imm_cd_prep_sets, matrix, nrow = 2, ncol = 3)
  expect_epred_by_cell(posterior_epred_imm_cd_full, dpars, imm_cd_prep_data)
})

test_that("each imm_cd family stores its functions under its own name", {
  dat <- imm_cd_data()
  for (version in c("full", "bsc", "abc")) {
    fit <- bmm(imm_cd_formula(version), dat, imm_cd_model(version),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
    family <- fit$formula$family
    for (fun in c("log_lik", "posterior_predict", "posterior_epred")) {
      # add_posterior_epred() finds the function for old fits by the family name
      expect_identical(
        family[[fun]],
        get(paste0(fun, "_", family$name), envir = asNamespace("bmm"))
      )
    }
    expect_equal(family$knowledge, "limited")
  }
})

# ---- Stan and R agree -----------------------------------------------------------

test_that("the Stan and R likelihoods agree for every version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the imm_cd likelihood"
  )
  chunks <- c(
    "circmix_funs.stan", "cd_funs.stan", "circmix_cd_funs.stan", "imm_funs.stan",
    "imm_cd_funs.stan"
  )
  code <- paste0(
    "functions {\n",
    paste(vapply(chunks, function(f) {
      read_lines2(system.file("stan_chunks", f, package = "bmm"))
    }, character(1)), collapse = "\n"),
    "\n}\n",
    "data {
      int NK; vector[NK] logk; vector[NK] dlogk; real logJ_min; real dlogJ;
      vector[32] gl_x; vector[32] gl_w;
      int N; int NT; array[N] int y; vector[N] probe; vector[N] kappa;
      vector[N] c; vector[N] a; vector[N] s; vector[N] b;
      array[N] int ss; matrix[N, NT] nt; matrix[N, NT] dist;
      int NTAU; vector[NTAU] tau; int NCRIT; vector[NCRIT] crit;
    }
    generated quantities {
      array[2, NTAU, NCRIT] vector[N] out_full;
      array[2, NTAU, NCRIT] vector[N] out_bsc;
      array[2, NTAU, NCRIT] vector[N] out_abc;
      for (r in 1:2) for (t in 1:NTAU) for (k in 1:NCRIT) for (i in 1:N) {
        out_full[r, t, k, i] = imm_cd_full_core(y[i], 0.0, kappa[i], tau[t], c[i],
          a[i], s[i], b[i], crit[k], ss[i], probe[i], to_vector(nt[i]),
          to_vector(dist[i]), 41, r - 1, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
        out_bsc[r, t, k, i] = imm_cd_bsc_core(y[i], 0.0, kappa[i], tau[t], c[i],
          s[i], b[i], crit[k], ss[i], probe[i], to_vector(nt[i]),
          to_vector(dist[i]), 41, r - 1, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
        out_abc[r, t, k, i] = imm_cd_abc_core(y[i], 0.0, kappa[i], tau[t], c[i],
          a[i], b[i], crit[k], ss[i], probe[i], to_vector(nt[i]), 41, r - 1,
          logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
      }
    }"
  )

  y <- c(0L, 1L, 1L, 0L, 1L, 0L)
  probe <- c(-2.5, -0.8, 0, 0.4, 1.9, 3.0)
  kappa <- c(6, 7, 9, 12, 20, 80)
  cc <- c(0.5, 1, 3, 5, 8, 12)
  a <- c(0.1, 0.3, 1, 2, 0.5, 4)
  s_grad <- c(0.2, 1, 2, 3, 0.5, 5)
  b <- c(1, 1, 1, 2, 0.5, 1)
  ss <- c(1L, 2L, 3L, 4L, 6L, 8L)
  nt <- matrix(seq(-3, 3, length.out = 6 * 7), nrow = 6)
  dist <- matrix(seq(0.1, 3, length.out = 6 * 7), nrow = 6)
  tau <- c(0, 0.5, 2)
  crit <- c(0, 0.4, -0.3)
  tab <- .circmix_kappa_table()
  gl <- .cd_gl_rule()

  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, gl_x = gl$x, gl_w = gl$w,
      N = length(y), NT = ncol(nt), y = y, probe = probe, kappa = kappa,
      c = cc, a = a, s = s_grad, b = b, ss = ss, nt = nt, dist = dist,
      NTAU = length(tau), tau = tau, NCRIT = length(crit), crit = crit
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  stan <- function(variable) as.numeric(fit$draws(variable))

  r_side <- function(version) {
    out <- array(NA_real_, c(2, length(tau), length(crit), length(y)))
    for (r in 1:2) for (t in seq_along(tau)) for (k in seq_along(crit)) {
      out[r, t, k, ] <- vapply(seq_along(y), function(i) {
        .cd_bernoulli_ld(y[i], .imm_cd_psame(probe[i], 0, kappa[i],
          c = cc[i], a = if (version != "bsc") a[i],
          s = if (version != "abc") s_grad[i], b = b[i], criterion = crit[k],
          set_size = ss[i], nt = nt[i, ], dist = dist[i, ], tau = tau[t],
          nodes = 41L, rich = r == 2, version = version
        ))
      }, numeric(1))
    }
    as.vector(out)
  }

  # cmdstanr flattens arrays with the first index fastest, as R does
  expect_equal(stan("out_full"), r_side("full"), tolerance = 1e-12)
  expect_equal(stan("out_bsc"), r_side("bsc"), tolerance = 1e-12)
  expect_equal(stan("out_abc"), r_side("abc"), tolerance = 1e-12)
})
