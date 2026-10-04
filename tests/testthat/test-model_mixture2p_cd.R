cd_data <- function(n = 200, seed = 1) {
  withr::with_seed(seed, {
    dat <- data.frame(target = runif(n, -pi, pi), set_size = rep(c(2, 4, 6, 8), length.out = n))
    changed <- rep(c(FALSE, TRUE), length.out = n)
    dat$probe <- wrap(dat$target + ifelse(changed, runif(n, -pi, pi), 0))
    dat$change <- rmixture2p_cd(n, wrap(dat$probe - dat$target), kappa = 6, p_mem = 0.7)
    dat
  })
}

cd_model <- function(...) {
  mixture2p_cd("change", "probe", "target", set_size = "set_size", ...)
}

cd_formula <- function(version) {
  if (version == "simple") bmf(kappa ~ 1, thetat ~ 1) else bmf(kappa ~ 1, K ~ 1)
}

versions <- c("simple", "slot", "slot_averaging")

test_that("mixture2p_cd() validates its arguments", {
  expect_error(mixture2p_cd("change", "probe"), "target")
  expect_error(mixture2p_cd("change", "probe", "target", version = "xyz"), "should be one of")
  expect_error(mixture2p_cd("change", "probe", "target", knowledge = "xyz"), "should be one of")
  expect_error(
    mixture2p_cd("change", "probe", "target", version = "slot"),
    "set_size argument is required"
  )
  expect_error(cd_model(vp_nodes = 40), "odd number")
})

test_that("mixture2p_cd is its own model, not a version of mixture2p", {
  model <- cd_model(version = "slot")
  expect_equal(
    class(model),
    c("bmmodel", "change_detection", "mixture2p_cd", "mixture2p_cd_slot")
  )
  expect_false(inherits(model, "circular"))
  expect_equal(model$resp_vars, list(response = "change", probe = "probe", target = "target"))
  expect_equal(model$knowledge, "limited")
  expect_equal(cd_model(knowledge = "rich")$knowledge, "rich")
  expect_output(print(model), "1 = 'change'")
})

test_that("mixture2p_cd has mixture2p's parameters plus criterion", {
  # joint fits of both tasks rely on the shared names
  for (version in versions) {
    for (vp in c(FALSE, TRUE)) {
      cd <- cd_model(version = version, variable_precision = vp)
      de <- mixture2p("y", set_size = "set_size", version = version, variable_precision = vp)
      expect_equal(names(cd$parameters), c(names(de$parameters), "criterion"))
      expect_equal(cd$links[names(de$links)], de$links)
    }
  }
})

test_that("criterion is fixed to 0 unless a formula frees it", {
  model <- cd_model()
  expect_equal(model$fixed_parameters, list(mu = 0, criterion = 0))
  freed <- check_model(model, cd_data(), bmf(kappa ~ 1, thetat ~ 1, criterion ~ 1))
  expect_false("criterion" %in% names(freed$fixed_parameters))
})

test_that("mixture2p_cd cites Lin & Oberauer, and van den Berg et al. only with variable precision", {
  refs <- model_citation(cd_model())
  expect_length(refs, 2)
  expect_match(refs[2], "^Lin, H.-Y., & Oberauer, K. \\(2022\\).*10.1016/j.cogpsych.2022.101463$")
  expect_length(model_citation(cd_model(variable_precision = TRUE)), 3)
})

test_that("dmixture2p_cd() is a probability mass function for every version", {
  for (version in versions) {
    for (tau in c(0, 1)) {
      p <- dmixture2p_cd(c(0, 1), 0.7, kappa = 5, p_mem = 0.8, K = 2.5, set_size = 4,
        tau = tau, version = version
      )
      expect_true(all(p > 0 & p < 1))
      expect_equal(sum(p), 1)
    }
  }
  expect_equal(
    dmixture2p_cd(1, 0.7, kappa = 5, log = TRUE),
    log(dmixture2p_cd(1, 0.7, kappa = 5))
  )
})

test_that("dmixture2p_cd() rejects invalid arguments", {
  expect_error(dmixture2p_cd(1, 0, kappa = -1), "kappa")
  expect_error(dmixture2p_cd(1, 0, p_mem = 1.2), "p_mem")
  expect_error(dmixture2p_cd(1, 0, K = 0, version = "slot"), "K must be positive")
  expect_error(dmixture2p_cd(2, 0), "response must be 0")
  expect_error(dmixture2p_cd(1, 0, knowledge = "some"), "should be one of")
})

test_that("an identical probe is called 'same' well above chance", {
  # under the decision rule of Lin & Oberauer (2022), not exactly 0.5
  expect_lt(dmixture2p_cd(1, 0, kappa = 8, p_mem = 0.8), 0.35)
  p <- dmixture2p_cd(1, seq(0, pi, length.out = 25), kappa = 6, p_mem = 0.75)
  expect_true(all(diff(p) > 0))
})

test_that("better memory and a higher criterion lower P('change') as expected", {
  expect_gt(
    dmixture2p_cd(1, pi / 2, kappa = 20, p_mem = 0.8),
    dmixture2p_cd(1, pi / 2, kappa = 3, p_mem = 0.8)
  )
  p <- vapply(c(-1, 0, 1), function(b) {
    dmixture2p_cd(1, 1, kappa = 6, p_mem = 0.7, criterion = b)
  }, numeric(1))
  expect_true(all(diff(p) < 0))
  # with nothing in memory, "change" is the share of the circle outside the arc
  expect_equal(
    dmixture2p_cd(1, 0, kappa = 6, p_mem = 0),
    1 - acos(.circmix_log_besselI0(6) / 6) / pi
  )
})

test_that("dmixture2p_cd() matches a fine trapezoid rule over the decision arc", {
  vm <- function(x, kappa) exp(kappa * cos(x) - .circmix_log_besselI0(kappa)) / (2 * pi)
  grid <- expand.grid(probe = c(0, 0.8, 2.2), kappa = c(1, 5, 25), p_mem = c(0.3, 0.9))
  reference <- mapply(function(probe, kappa, p_mem) {
    hw <- acos(.circmix_log_besselI0(kappa) / kappa)
    x <- seq(probe - hw, probe + hw, length.out = 2e5 + 1)
    f <- p_mem * vm(x, kappa) + (1 - p_mem) / (2 * pi)
    sum((f[-1] + f[-length(f)]) / 2) * (x[2] - x[1])
  }, grid$probe, grid$kappa, grid$p_mem)
  expect_equal(
    dmixture2p_cd(0, grid$probe, kappa = grid$kappa, p_mem = grid$p_mem),
    reference,
    tolerance = 1e-8
  )
})

test_that("the capacity versions reduce to the simple version where they should", {
  probe <- c(0, 0.5, 1.7)
  # slot: an item is held with probability min(1, K / set_size)
  expect_equal(
    dmixture2p_cd(1, probe, kappa = 6, K = 3, set_size = 5, version = "slot"),
    dmixture2p_cd(1, probe, kappa = 6, p_mem = 0.6)
  )
  expect_equal(
    dmixture2p_cd(1, probe, kappa = 6, K = 7, set_size = 5, version = "slot"),
    dmixture2p_cd(1, probe, kappa = 6, p_mem = 1)
  )
  # slot averaging with K a multiple of the set size: two slots per item
  two_slots <- .circmix_kappa(2 * .circmix_J(6))
  for (knowledge in c("limited", "rich")) {
    expect_equal(
      dmixture2p_cd(1, probe, kappa = 6, K = 4, set_size = 2, version = "slot_averaging",
        knowledge = knowledge
      ),
      dmixture2p_cd(1, probe, kappa = two_slots, p_mem = 1),
      tolerance = 1e-10
    )
  }
})

test_that("the observer's knowledge matters only when precision varies", {
  probe <- c(0.2, 1.1)
  for (version in versions) {
    # K / set_size = 1.25: an item holds one slot or two
    args <- list(1, probe, kappa = 6, p_mem = 0.8, K = 2.5, set_size = 2, criterion = 0.3,
      version = version
    )
    constant <- lapply(c("limited", "rich"), function(k) {
      do.call(dmixture2p_cd, c(args, knowledge = k))
    })
    if (version == "slot_averaging") {
      # two slot counts are two precision states even at constant precision
      expect_false(isTRUE(all.equal(constant[[1]], constant[[2]])))
    } else {
      expect_equal(constant[[1]], constant[[2]])
    }
    varying <- lapply(c("limited", "rich"), function(k) {
      do.call(dmixture2p_cd, c(args, tau = 1, knowledge = k))
    })
    expect_false(isTRUE(all.equal(varying[[1]], varying[[2]])))
  }
})

test_that("dmixture2p_cd() approaches constant precision as tau goes to 0", {
  for (version in versions) {
    args <- list(1, c(0.2, 1.1), kappa = 6, p_mem = 0.8, K = 2.5, set_size = 4,
      version = version
    )
    expect_equal(
      do.call(dmixture2p_cd, c(args, tau = 1e-3)),
      do.call(dmixture2p_cd, args),
      tolerance = 1e-3
    )
  }
})

test_that("rmixture2p_cd() draws from dmixture2p_cd()", {
  cases <- list(
    list(probe = 1.2, kappa = 6, p_mem = 0.7),
    list(probe = 0.4, kappa = 6, p_mem = 0.7, tau = 1),
    list(probe = 0.4, kappa = 6, K = 3, set_size = 5, version = "slot_averaging")
  )
  for (case in cases) {
    p <- withr::with_seed(3, mean(do.call(rmixture2p_cd, c(n = 1e5, case))))
    expect_equal(p, do.call(dmixture2p_cd, c(response = 1, case)), tolerance = 0.02)
  }
})

test_that("every version runs through the bmm() pipeline and passes the probe to Stan", {
  dat <- cd_data()
  for (version in versions) {
    for (vp in c(FALSE, TRUE)) {
      model <- cd_model(version = version, variable_precision = vp)
      formula <- cd_formula(version)
      if (vp) formula <- formula + bmf(tau ~ 1)
      expect_silent(bmm(formula, dat, model, backend = "mock", mock_fit = 1, rename = FALSE))
      code <- stancode(formula, dat, model)
      expect_match(
        code,
        glue::glue(
          "mixture2p_cd_{version}_lpmf(Y[n] | mu[n], kappa[n], {if (vp) 'tau[n], ' else ''}",
          "{if (version == 'simple') 'thetat' else 'K'}[n], criterion[n], ",
          "{if (version == 'simple') '' else 'vint1[n], '}vreal1[n]"
        ),
        fixed = TRUE
      )
      expect_match(code, "Intercept_criterion = 0", fixed = TRUE)
    }
  }
  sdata <- standata(bmf(kappa ~ 1, K ~ 1), dat, cd_model(version = "slot"))
  expect_length(sdata$cd_gl_x, 32)
  expect_equal(as.numeric(sdata$vreal1), wrap(dat$probe - dat$target))
  expect_equal(as.numeric(sdata$vint1), dat$set_size)
})

test_that("a freed criterion is sampled rather than held at 0", {
  code <- stancode(bmf(kappa ~ 1, thetat ~ 1, criterion ~ 1), cd_data(), cd_model())
  expect_no_match(code, "Intercept_criterion = 0", fixed = TRUE)
})

test_that("the knowledge setting reaches the Stan likelihood and the family", {
  dat <- cd_data()
  for (knowledge in c("limited", "rich")) {
    model <- cd_model(variable_precision = TRUE, knowledge = knowledge)
    code <- stancode(bmf(kappa ~ 1, tau ~ 1, thetat ~ 1), dat, model)
    expect_match(
      code,
      glue::glue("criterion, probe, 41, {as.integer(knowledge == 'rich')}, circmix_logk"),
      fixed = TRUE
    )
    fit <- bmm(bmf(kappa ~ 1, tau ~ 1, thetat ~ 1), dat, model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )
    expect_equal(fit$formula$family$knowledge, knowledge)
  }
})

test_that("check_data() reads the set size only for the capacity versions", {
  dat <- cd_data()
  checked <- check_data(cd_model(version = "slot"), dat, cd_formula("slot"))
  expect_equal(checked$ss_numeric, dat$set_size)
  expect_null(check_data(cd_model(), dat, cd_formula("simple"))$ss_numeric)
  expect_warning(
    check_data(cd_model(version = "slot"), transform(dat, set_size = 4), cd_formula("slot")),
    "not identified by the data"
  )
})

test_that("posterior_predict and posterior_epred agree with dmixture2p_cd()", {
  sets <- data.frame(
    mu = c(0, 0.1), kappa = c(6, 12), thetat = c(0.7, 0.9), criterion = c(0, 0.3)
  )
  probe <- c(0.4, 1.5)
  res <- epred_vs_predict(
    sets, posterior_epred_mixture2p_cd_simple, posterior_predict_mixture2p_cd_simple,
    data = list(vreal1 = probe)
  )
  expected <- dmixture2p_cd(1, probe - sets$mu, kappa = sets$kappa, p_mem = sets$thetat,
    criterion = sets$criterion
  )
  expect_equal(unname(res[, "epred"]), expected)
  expect_lt(max(abs(res[, "mc"] - res[, "epred"])), 0.015)

  withr::with_seed(5, {
    dpars <- list(
      mu = matrix(0, 3, 4), kappa = matrix(runif(12, 2, 20), 3, 4),
      K = matrix(runif(12, 1, 6), 3, 4), criterion = matrix(runif(12, -0.3, 0.3), 3, 4),
      tau = matrix(runif(12, 0.2, 1), 3, 4)
    )
  })
  data <- list(vreal1 = c(0, 0.5, 1.5, -2), vint1 = c(2L, 3L, 5L, 8L))
  expect_epred_by_cell(posterior_epred_mixture2p_cd_slot_averaging, dpars, data)
  expect_epred_by_cell(posterior_epred_mixture2p_cd_slot, dpars, data)
})

test_that("log_lik returns one value per draw, matching dmixture2p_cd()", {
  # needs .cd_bernoulli_ld() to recycle y over the draws (cd-core fix)
  skip_if(
    length(.cd_bernoulli_ld(0, c(0.3, 0.4))) == 1,
    ".cd_bernoulli_ld() does not yet recycle the response over the draws"
  )
  sets <- data.frame(
    mu = c(0, 0.1), kappa = c(6, 12), thetat = c(0.7, 0.9), criterion = c(0, 0.3)
  )
  prep <- epred_prep(lapply(sets, function(x) matrix(x, nrow = 2, ncol = 1)),
    data = list(Y = 0L, vreal1 = 0.4)
  )
  expect_equal(
    log_lik_mixture2p_cd_simple(1, prep),
    dmixture2p_cd(0, 0.4 - sets$mu, kappa = sets$kappa, p_mem = sets$thetat,
      criterion = sets$criterion, log = TRUE
    )
  )
})

test_that("a mixture2p_cd fit stores its posterior_epred and returns P('change')", {
  dat <- cd_data()
  for (version in versions) {
    fit <- bmm(cd_formula(version), dat, cd_model(version = version),
      backend = "mock", mock_fit = 1, rename = FALSE
    )
    expect_identical(
      fit$formula$family$posterior_epred,
      get(paste0("posterior_epred_mixture2p_cd_", version), envir = asNamespace("bmm"))
    )
    expect_identical(
      fit$formula$family$log_lik,
      get(paste0("log_lik_mixture2p_cd_", version), envir = asNamespace("bmm"))
    )
  }
})

test_that("the Stan and R likelihoods agree for every version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the mixture2p_cd likelihood"
  )

  chunks <- c(
    "circmix_funs.stan", "cd_funs.stan", "circmix_cd_funs.stan",
    "mixture2p_funs.stan", "mixture2p_cd_funs.stan"
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
      int N; array[N] int y; vector[N] mu; vector[N] probe; vector[N] kappa;
      vector[N] tau; vector[N] thetat; vector[N] K; array[N] int ss;
      vector[N] criterion;
    }
    generated quantities {
      array[2] vector[N] out_simple; array[2] vector[N] out_slot;
      array[2] vector[N] out_sa;
      for (r in 1:2) {
        for (i in 1:N) {
          out_simple[r, i] = mixture2p_cd_simple_core(y[i], mu[i], kappa[i], tau[i],
            thetat[i], criterion[i], probe[i], 41, r - 1, logk, dlogk, logJ_min, dlogJ,
            gl_x, gl_w);
          out_slot[r, i] = mixture2p_cd_slot_core(y[i], mu[i], kappa[i], tau[i], K[i],
            criterion[i], ss[i], probe[i], 41, r - 1, logk, dlogk, logJ_min, dlogJ,
            gl_x, gl_w);
          out_sa[r, i] = mixture2p_cd_slot_averaging_core(y[i], mu[i], kappa[i], tau[i],
            K[i], criterion[i], ss[i], probe[i], 41, r - 1, logk, dlogk, logJ_min,
            dlogJ, gl_x, gl_w);
        }
      }
    }"
  )

  base <- data.frame(
    mu = c(0, 0.2, -0.1), probe = c(0.3, 2.4, -1.1), kappa = c(4, 8, 30),
    thetat = c(0.5, 0.8, 0.95), K = c(1.5, 3, 4.5), ss = c(4L, 2L, 6L)
  )
  grid <- expand.grid(row = seq_len(nrow(base)), tau = c(0, 0.5, 2), criterion = c(0, 0.4, -0.3))
  cases <- cbind(base[grid$row, ], grid[c("tau", "criterion")])
  cases$y <- rep(c(0L, 1L), length.out = nrow(cases))

  tab <- .circmix_kappa_table()
  gl <- .cd_gl_rule()
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = c(
      list(
        NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
        logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, gl_x = gl$x, gl_w = gl$w,
        N = nrow(cases)
      ),
      as.list(cases)
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )

  r_ld <- function(version, rich) {
    with(cases, .cd_bernoulli_ld(y, .mixture2p_cd_psame(
      probe, mu, kappa, thetat, criterion, tau, K, ss, version, rich
    )))
  }
  # cmdstanr flattens arrays with the first index fastest
  outputs <- c(simple = "out_simple", slot = "out_slot", slot_averaging = "out_sa")
  for (version in names(outputs)) {
    stan <- matrix(as.numeric(fit$draws(outputs[[version]])), nrow = 2)
    expect_equal(stan[1, ], r_ld(version, FALSE), tolerance = 1e-12)
    expect_equal(stan[2, ], r_ld(version, TRUE), tolerance = 1e-12)
  }
})
