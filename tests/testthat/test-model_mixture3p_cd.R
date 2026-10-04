cd3p_model <- function(...) {
  mixture3p_cd("resp", "probe", "target", nt_features = c("nt1", "nt2"),
    set_size = "ss", ...
  )
}

cd3p_data <- function(n = 60) {
  withr::with_seed(136, {
    dat <- data.frame(
      target = stats::runif(n, -pi, pi), nt1 = stats::runif(n, -pi, pi),
      nt2 = stats::runif(n, -pi, pi), ss = rep(c(2, 3), length.out = n)
    )
    dat$probe <- wrap(dat$target + c(0, 1.5, -2)[(seq_len(n) %% 3) + 1])
    dat$resp <- stats::rbinom(n, 1, 0.5)
    dat
  })
}

cd3p_formulas <- list(
  simple = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
  slot = bmf(kappa ~ 1, K ~ 1, pnt ~ 1),
  slot_averaging = bmf(kappa ~ 1, K ~ 1, pnt ~ 1)
)

vm_arc <- function(lo, hi, mu, kappa) {
  stats::integrate(
    function(x) exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, 0)),
    lo, hi, rel.tol = 1e-12
  )$value
}

test_that("mixture3p_cd() is its own model with mixture3p's parameters plus criterion", {
  for (version in c("simple", "slot", "slot_averaging")) {
    for (vp in c(FALSE, TRUE)) {
      cd <- cd3p_model(version = version, variable_precision = vp)
      de <- mixture3p("y", c("nt1", "nt2"), "ss", version = version, variable_precision = vp)
      expect_equal(names(cd$parameters), c(names(de$parameters), "criterion"))
      expect_equal(cd$links[names(de$links)], de$links)
      expect_equal(
        class(cd),
        c("bmmodel", "change_detection", "non_targets", "mixture3p_cd",
          paste0("mixture3p_cd_", version))
      )
    }
  }
  model <- cd3p_model()
  expect_false(inherits(model, "mixture3p"))
  expect_false(inherits(model, "circular"))
  expect_equal(model$fixed_parameters, list(mu = 0, criterion = 0))
  expect_equal(model$knowledge, "limited")
  expect_equal(cd3p_model(knowledge = "rich")$knowledge, "rich")
  expect_error(cd3p_model(version = "slots"), "should be one of")
  expect_error(cd3p_model(knowledge = "full"), "should be one of")
  expect_error(cd3p_model(vp_nodes = 21), "at least 41")
  expect_error(mixture3p_cd("resp", "probe", nt_features = "nt1", set_size = 2), "target")
})

test_that("mixture3p_cd() takes links except on its softmax weights", {
  expect_equal(settable_links(cd3p_model()), c("mu", "kappa", "criterion"))
  expect_error(cd3p_model(links = list(thetat = "logit")), "cannot be changed")
  dat <- cd3p_data()
  model <- cd3p_model(version = "slot", links = list(kappa = "softplus", pnt = "probit"))
  ff <- cd3p_formulas$slot
  family <- configure_model(model, check_data(model, dat, ff), ff)$formula$family
  expect_equal(c(family$link_kappa, family$link_pnt), c("softplus", "probit"))
  model <- cd3p_model()
  ff <- cd3p_formulas$simple
  family <- configure_model(model, check_data(model, dat, ff), ff)$formula$family
  expect_equal(c(family$link_thetat, family$link_thetant), c("identity", "identity"))
})

test_that("mixture3p_cd() cites the mixture model, the decision rule and its options", {
  cites <- model_citation(cd3p_model())
  expect_length(cites, 2)
  expect_match(cites[2], "Lin, H.-Y., & Oberauer", fixed = TRUE)
  vp <- model_citation(cd3p_model(version = "slot", variable_precision = TRUE))
  expect_length(vp, 4)
  expect_match(vp[2], "Zhang, W., & Luck", fixed = TRUE)
  expect_match(vp[4], "van den Berg", fixed = TRUE)
})

test_that("P('same') is the retrieval mass inside the unbiased arc", {
  # Lin & Oberauer (2022), Appendix B: at criterion 0 the arc half-width solves
  # vM(hw | 0, kappa) = 1 / (2 pi), whatever p_s
  kappa <- 6
  hw <- acos(log(besselI(kappa, 0)) / kappa)
  mu <- c(0, 2, -1.5)
  for (probe in c(0, 0.7, 2, -2.8)) {
    expected <- 0.6 * vm_arc(probe - hw, probe + hw, 0, kappa) +
      0.1 * vm_arc(probe - hw, probe + hw, 2, kappa) +
      0.1 * vm_arc(probe - hw, probe + hw, -1.5, kappa) +
      0.2 * hw / pi
    expect_equal(
      dmixture3p_cd(0, probe, mu = mu, kappa = kappa, p_mem = 0.6, p_nt = 0.2),
      expected,
      tolerance = 1e-9
    )
  }
})

test_that("dmixture3p_cd() is a probability mass function and rmixture3p_cd() samples it", {
  p <- dmixture3p_cd(c(0, 1), probe = 0.7, mu = c(0, 1.5, -2), kappa = 5)
  expect_equal(sum(p), 1)
  expect_equal(
    dmixture3p_cd(1, 0.7, mu = c(0, 1.5, -2), kappa = 5, log = TRUE), log(p[2])
  )
  withr::local_seed(136)
  draws <- rmixture3p_cd(4e4, probe = 0.7, mu = c(0, 1.5, -2), kappa = 5, tau = 1)
  expect_equal(
    mean(draws),
    dmixture3p_cd(1, 0.7, mu = c(0, 1.5, -2), kappa = 5, tau = 1),
    tolerance = 0.01
  )
  expect_error(dmixture3p_cd(2, 0), "coded 0")
})

test_that("an identical probe is rejected far above chance", {
  p <- dmixture3p_cd(1, probe = 0, mu = c(0, 2, -2.5), kappa = 8, p_mem = 0.8, p_nt = 0.1)
  expect_lt(p, 0.35)
})

test_that("a probe at a non-target is harder to reject than a new one (intrusion cost)", {
  mu <- c(0, 2, -0.9)
  intrusion <- dmixture3p_cd(1, probe = 2, mu = mu, kappa = 8, p_mem = 0.6, p_nt = 0.3)
  new_probe <- dmixture3p_cd(1, probe = -2, mu = mu, kappa = 8, p_mem = 0.6, p_nt = 0.3)
  expect_lt(intrusion, new_probe)
  # without swaps the non-targets play no role
  expect_equal(
    dmixture3p_cd(1, probe = 2, mu = mu, kappa = 8, p_mem = 0.6, p_nt = 0),
    dmixture3p_cd(1, probe = -2, mu = mu, kappa = 8, p_mem = 0.6, p_nt = 0),
    tolerance = 1e-12
  )
})

test_that("P('change') grows with the distance of the probe from the target", {
  for (knowledge in c("limited", "rich")) {
    # one item, so no non-target near the far probes pulls P("change") back
    p <- dmixture3p_cd(1, probe = seq(0, pi, length.out = 12), mu = 0,
      kappa = 6, p_mem = 0.7, p_nt = 0, tau = 1, knowledge = knowledge
    )
    expect_true(all(diff(p) > 0))
  }
})

test_that("non-targets beyond a trial's set size drop out", {
  args <- list(
    probe = 1.3, mu = 0, kappa = 6, thetat = 1, thetant = 0.4, criterion = 0.2,
    tau = 0.5, rich = FALSE, nodes = 41L
  )
  short <- do.call(.mixture3p_cd_psame_simple, c(args, list(set_size = 2, nt = 1.3)))
  padded <- do.call(.mixture3p_cd_psame_simple, c(args, list(set_size = 2, nt = c(1.3, 2))))
  expect_equal(padded, short)
  expect_false(isTRUE(all.equal(
    do.call(.mixture3p_cd_psame_simple, c(args, list(set_size = 3, nt = c(1.3, 2)))),
    short
  )))
})

test_that("the capacity versions reduce to the simple one where they should", {
  common <- list(probe = 0.9, mu = 0, kappa = 7, criterion = 0.3, set_size = 3,
                 nt = c(0.5, -2), tau = 0, rich = FALSE, nodes = 41L)
  # below capacity nothing is guessed: p_mem = 1, the swap probability splits
  # the rest. thetat/thetant on the softmax scale with guessing at -Inf
  # is approximated by a large reference gap
  slot <- do.call(.mixture3p_cd_psame_slot, c(common, list(K = 5, p_nt = 0.2)))
  simple <- do.call(.mixture3p_cd_psame_simple, c(common, list(
    thetat = log(0.8) + 40, thetant = log(0.2) + 40
  )))
  expect_equal(slot, simple, tolerance = 1e-12)
  # slot averaging with K a multiple of the set size holds every item in the
  # same number of slots, so it is the slot model at that precision. Its
  # observer knows the item is stored (p_s = 1) where the slot observer takes
  # the target weight, which only matters away from criterion 0
  # the precision goes through the tabulated J -> kappa inverse there
  # (relative error 1.7e-11), hence the looser tolerance
  common$criterion <- 0
  expect_equal(
    do.call(.mixture3p_cd_psame_slot_averaging, c(common, list(K = 3, p_nt = 0.2))),
    do.call(.mixture3p_cd_psame_slot, c(common, list(K = 3, p_nt = 0.2))),
    tolerance = 1e-9
  )
})

test_that("the swap parameter is flagged as prior-only at set size 1", {
  dat <- cd3p_data()
  dat$ss[1:10] <- 1
  expect_warning(
    bmm(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1 + ss), dat,
      cd3p_model(), backend = "mock", mock_fit = 1, rename = FALSE
    ),
    "set-size-1 level"
  )
})

test_that("every version runs through the bmm() pipeline and wires its Stan data", {
  dat <- cd3p_data()
  for (version in names(cd3p_formulas)) {
    for (vp in c(FALSE, TRUE)) {
      model <- cd3p_model(version = version, variable_precision = vp)
      ff <- cd3p_formulas[[version]]
      if (vp) ff <- ff + bmf(tau ~ 1)
      expect_s3_class(
        suppressMessages(bmm(ff, dat, model, backend = "mock", mock_fit = 1, rename = FALSE)),
        "bmmfit"
      )
    }
  }
  code <- suppressMessages(stancode(cd3p_formulas$simple, dat, cd3p_model()))
  expect_match(code, paste0(
    "mixture3p_cd_simple_lpmf(Y[n] | mu[n], kappa[n], thetat[n], thetant[n], ",
    "criterion[n], vint1[n], vreal1[n], vreal2[n], vreal3[n]"
  ), fixed = TRUE)
  expect_match(code, "Intercept_criterion = 0", fixed = TRUE)
  standata <- suppressMessages(standata(cd3p_formulas$simple, dat, cd3p_model()))
  expect_length(standata$cd_gl_x, 32)
  expect_equal(as.numeric(standata$vreal1), wrap(dat$probe - dat$target))
  expect_equal(as.numeric(standata$vint1), dat$ss)
  # a freed criterion is sampled rather than pinned
  code <- suppressMessages(stancode(
    cd3p_formulas$simple + bmf(criterion ~ 1), dat, cd3p_model(knowledge = "rich")
  ))
  expect_no_match(code, "Intercept_criterion = 0", fixed = TRUE)
  expect_match(code, "to_vector({nt1, nt2}), 41, 1, circmix_logk", fixed = TRUE)
})

test_that("every version's family carries its own epred, log_lik and predict", {
  dat <- cd3p_data()
  for (version in names(cd3p_formulas)) {
    fit <- suppressMessages(bmm(cd3p_formulas[[version]], dat, cd3p_model(version = version),
      backend = "mock", mock_fit = 1, rename = FALSE
    ))
    family <- fit$formula$family
    for (method in c("posterior_epred", "log_lik", "posterior_predict")) {
      expect_identical(
        family[[method]],
        get(paste0(method, "_mixture3p_cd_", version), envir = asNamespace("bmm"))
      )
    }
    expect_equal(family$knowledge, "limited")
    fit$formula$family$posterior_epred <- NULL
    expect_identical(add_posterior_epred(fit)$formula$family$posterior_epred, family$posterior_epred)
  }
})

cd3p_prep_data <- list(
  Y = c(1, 0, 1), vint1 = c(3, 2, 3), vreal1 = c(0.4, -1, 2.5),
  vreal2 = c(2.5, 1, -0.5), vreal3 = c(-1.5, 0, 1)
)

test_that("posterior_epred is the probability of a 'change' that posterior_predict simulates", {
  skip_on_cran()
  withr::local_seed(136)
  sets <- data.frame(mu = 0, kappa = c(4, 8, 15), thetat = c(1, 0.5, 2),
                     thetant = c(-1, 0, 0.5), criterion = c(0, 0.4, -0.3))
  res <- epred_vs_predict(sets, posterior_epred_mixture3p_cd_simple,
    posterior_predict_mixture3p_cd_simple, data = cd3p_prep_data
  )
  # Monte-Carlo SE of a proportion near 0.5 at 20000 draws is 0.0035
  expect_lt(max(abs(res[, "mc"] - res[, "epred"])), 0.015)
  expect_epred_by_cell(
    posterior_epred_mixture3p_cd_slot,
    list(mu = matrix(0, 2, 3), kappa = matrix(c(5, 9), 2, 3), K = matrix(c(1.5, 2.5), 2, 3),
         pnt = matrix(0.2, 2, 3), criterion = matrix(c(0, 0.3), 2, 3)),
    cd3p_prep_data
  )
})

test_that("log_lik scores the observed response with the R twin", {
  # depends on the cd-core fix that recycles y in .cd_bernoulli_ld() to the
  # draws; until it is merged, log_lik returns one value per observation
  skip_if(
    length(.cd_bernoulli_ld(1, c(0.2, 0.3))) == 1,
    "needs the cd-core fix of .cd_bernoulli_ld()"
  )
  prep <- epred_prep(
    list(mu = matrix(0, 2, 3), kappa = matrix(c(5, 9), 2, 3), K = matrix(2.5, 2, 3),
         pnt = matrix(c(0.1, 0.3), 2, 3), criterion = matrix(0, 2, 3)),
    cd3p_prep_data
  )
  p_same <- .mixture3p_cd_psame_slot_averaging(
    probe = 2.5, mu = 0, kappa = c(5, 9), K = 2.5, p_nt = c(0.1, 0.3),
    criterion = 0, set_size = 3, nt = c(-0.5, 1), tau = 0, rich = FALSE, nodes = 41L
  )
  expect_equal(log_lik_mixture3p_cd_slot_averaging(3, prep), log1p(-p_same))
})

test_that("the Stan and R change-detection likelihoods agree for every version", {
  skip_on_cran()
  skip_if_not(
    requireNamespace("cmdstanr", quietly = TRUE),
    "cmdstanr is required to compile the mixture3p_cd likelihood"
  )
  chunks <- c(
    "circmix_funs.stan", "cd_funs.stan", "circmix_cd_funs.stan",
    "mixture3p_funs.stan", "mixture3p_cd_funs.stan"
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
      int N; int NT; array[N] int y; vector[N] probe; vector[N] kappa; vector[N] tau;
      vector[N] thetat; vector[N] thetant; vector[N] K; vector[N] pnt;
      vector[N] crit; array[N] int ss; matrix[N, NT] nt;
    }
    generated quantities {
      array[2] vector[N] out_simple; array[2] vector[N] out_slot; array[2] vector[N] out_sa;
      for (r in 1:2) {
        for (i in 1:N) {
          out_simple[r, i] = mixture3p_cd_simple_core(y[i], 0.0, kappa[i], tau[i],
            thetat[i], thetant[i], crit[i], ss[i], probe[i], to_vector(nt[i]), 41,
            r - 1, logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
          out_slot[r, i] = mixture3p_cd_slot_core(y[i], 0.0, kappa[i], tau[i], K[i],
            pnt[i], crit[i], ss[i], probe[i], to_vector(nt[i]), 41, r - 1, logk,
            dlogk, logJ_min, dlogJ, gl_x, gl_w);
          out_sa[r, i] = mixture3p_cd_slot_averaging_core(y[i], 0.0, kappa[i], tau[i],
            K[i], pnt[i], crit[i], ss[i], probe[i], to_vector(nt[i]), 41, r - 1,
            logk, dlogk, logJ_min, dlogJ, gl_x, gl_w);
        }
      }
    }"
  )

  # every tau crossed with every criterion, over set sizes 1 to 4 and both
  # sides of the slot-averaging zero-slot branch
  grid <- expand.grid(tau = c(0, 0.5, 2), crit = c(0, 0.4, -0.3))
  n <- nrow(grid)
  y <- rep(c(0L, 1L), length.out = n)
  probe <- seq(-2.8, 2.9, length.out = n)
  kappa <- rep(c(2, 6, 15), length.out = n)
  thetat <- rep(c(-0.5, 0.8, 1.5), length.out = n)
  thetant <- rep(c(-1, 0, 0.5), length.out = n)
  K <- rep(c(1.5, 2.5, 3.5, 0.8), length.out = n)
  pnt <- rep(c(0.1, 0.3), length.out = n)
  ss <- rep(c(1L, 2L, 3L, 4L), length.out = n)
  nt <- matrix(seq(-3, 3, length.out = n * 3), nrow = n)

  tab <- .circmix_kappa_table()
  gl <- .cd_gl_rule()
  model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code))
  fit <- model$sample(
    data = list(
      NK = length(tab$logkappa), logk = tab$logkappa, dlogk = tab$dlogkappa,
      logJ_min = tab$logJ_min, dlogJ = tab$dlogJ, gl_x = gl$x, gl_w = gl$w,
      N = n, NT = ncol(nt), y = y, probe = probe, kappa = kappa, tau = grid$tau,
      thetat = thetat, thetant = thetant, K = K, pnt = pnt, crit = grid$crit,
      ss = ss, nt = nt
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  stan <- function(variable) as.numeric(fit$draws(variable))

  r_side <- function(psame) {
    # cmdstanr flattens arrays with the first index fastest
    as.vector(t(vapply(c(FALSE, TRUE), function(rich) {
      vapply(seq_len(n), function(i) .cd_bernoulli_ld(y[i], psame(i, rich)), numeric(1))
    }, numeric(n))))
  }
  common <- function(i, rich) {
    list(probe = probe[i], mu = 0, kappa = kappa[i], criterion = grid$crit[i],
         set_size = ss[i], nt = nt[i, ], tau = grid$tau[i], rich = rich, nodes = 41L)
  }
  expect_equal(stan("out_simple"), r_side(function(i, rich) {
    do.call(.mixture3p_cd_psame_simple, c(common(i, rich), list(
      thetat = thetat[i], thetant = thetant[i]
    )))
  }), tolerance = 1e-12)
  expect_equal(stan("out_slot"), r_side(function(i, rich) {
    do.call(.mixture3p_cd_psame_slot, c(common(i, rich), list(K = K[i], p_nt = pnt[i])))
  }), tolerance = 1e-12)
  expect_equal(stan("out_sa"), r_side(function(i, rich) {
    do.call(.mixture3p_cd_psame_slot_averaging, c(common(i, rich), list(K = K[i], p_nt = pnt[i])))
  }), tolerance = 1e-12)
})
