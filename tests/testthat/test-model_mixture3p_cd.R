cd3p_model <- function(version = "simple", ...) {
  mixture3p_cd("change", "probe", nt_features = paste0("col_nt", 1:7),
               set_size = "set_size", version = version, ...)
}

cd3p_formulas <- list(
  simple = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1),
  slot = bmf(kappa ~ 1, K ~ 1, ptarget ~ 1),
  slot_averaging = bmf(kappa ~ 1, K ~ 1, ptarget ~ 1)
)

# von Mises mass on [lo, hi] by adaptive integration, independent of the
# Gauss-Legendre arc of the package
vm_mass <- function(lo, hi, mu, kappa) {
  stats::integrate(
    function(x) exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, 0)),
    lo, hi, rel.tol = 1e-12
  )$value
}

test_that("mixture3p_cd() validates its arguments", {
  expect_error(mixture3p_cd("r", "p"), "nt_features")
  expect_error(cd3p_model(version = "full"), "should be one of")
  expect_error(cd3p_model(knowledge = "partial"), "should be one of")
  expect_error(cd3p_model(variable_precision = TRUE), "planned for bmm 1.6.0")
  expect_error(cd3p_model(variable_precision = NA), "TRUE or FALSE")
  # a vp_nodes copied from mixture3p() is reported, not silently dropped
  expect_warning(cd3p_model(vp_nodes = 81L), "vp_nodes")
})

test_that("mixture3p_cd() takes mixture3p's parameters, priors, links and inits, plus criterion", {
  for (version in names(cd3p_formulas)) {
    cd <- cd3p_model(version)
    de <- mixture3p("y", nt_features = paste0("col_nt", 1:7), set_size = "set_size", version = version)
    expect_identical(names(cd$parameters), c(names(de$parameters), "criterion"), label = version)
    for (field in c("links", "default_priors", "init_ranges")) {
      expect_identical(cd[[field]][names(de[[field]])], de[[field]], label = paste(version, field))
    }
    expect_identical(cd$fixed_parameters, list(mu = 0, criterion = 0), label = version)
    expect_identical(cd$links$criterion, "identity")
    expect_s3_class(
      cd, c("bmmodel", "change_detection", "non_targets", "mixture3p_cd", paste0("mixture3p_cd_", version)),
      exact = TRUE
    )
  }
  expect_identical(cd3p_model("slot_averaging", knowledge = "rich")$knowledge, "rich")
  expect_identical(cd3p_model()$knowledge, "limited")
  expect_identical(cd3p_model()$task, "Detection, recognition and confidence judgments")
  expect_identical(cd3p_model()$domain, "Working memory")
  expect_true("mixture3p_cd" %in% bmm_models())
  # the constructor's arguments shared with the reproduction twin keep its
  # order; knowledge sits after variable_precision, before links (CD-D1 F1)
  expect_identical(
    names(formals(mixture3p_cd)),
    c("response", "probe", "nt_features", "set_size", "regex", "version",
      "variable_precision", "knowledge", "links", "...")
  )
})

test_that("knowledge = 'rich' warns in the versions with one precision state", {
  expect_warning(cd3p_model(knowledge = "rich"), "no effect in the 'simple' version")
  expect_warning(cd3p_model("slot", knowledge = "rich"), "no effect in the 'slot' version")
  expect_no_warning(cd3p_model("slot_averaging", knowledge = "rich"))
})

test_that("the link of criterion cannot be set", {
  expect_error(cd3p_model(links = list(criterion = "log")), "'criterion' cannot be changed")
})

test_that("criterion's init range starts on the linear part of its bound", {
  # the bound is per trial; at the corners of the inits of the memory
  # parameters and set sizes 1 to 8, the ends of criterion's init range keep
  # more than half of their value under the capacity versions. p_s is the
  # target's retrieval weight; at the simple corner thetat -1.1, thetant 1.1 it
  # is 0.077 at every set size above 1, and criterion 0.2 keeps 0.39 of its
  # value there
  simple <- cd3p_model()$init_ranges
  p_s_simple <- unlist(lapply(1:8, function(ss) {
    corners <- expand.grid(thetat = simple$thetat, thetant = simple$thetant)
    exp(.mixture3p_softmax_weights(corners$thetat, corners$thetant, ss)$logw[, 1])
  }))
  capacity <- cd3p_model("slot")$init_ranges
  p_s_slot <- as.vector(outer(pmin(1, outer(capacity$K, 1:8, "/")), capacity$ptarget))
  kept <- function(p_s) {
    corners <- expand.grid(p_s = p_s, kappa = simple$kappa, criterion = simple$criterion)
    min(.circmix_cd_bounded_criterion(corners$criterion, corners$p_s, corners$kappa, 0) /
      corners$criterion)
  }
  expect_gt(kept(c(p_s_slot, 1)), 0.5)
  expect_equal(kept(p_s_simple), 0.393, tolerance = 1e-3)
})

test_that("dmixture3p_cd() and rmixture3p_cd() refuse what they cannot evaluate", {
  # `log` passed by position, after the memory parameters, lands on criterion
  expect_error(dmixture3p_cd(1, 0.5, c(0, 2, -1.5), 5, 0.6, 0.2, TRUE), "Did you mean log = TRUE")
  expect_error(rmixture3p_cd(5, 0.5, criterion = TRUE), "criterion must be numeric")
  expect_true(is.na(dmixture3p_cd(1, 0.5, criterion = NA)))
  for (n in list(2.5, c(2, 3), 0, "3")) {
    expect_error(rmixture3p_cd(n, 0.5), "single positive whole number")
  }
  expect_error(dmixture3p_cd(2, 0.5), "response must be 0")
  expect_error(dmixture3p_cd(NA, 0.5), "response must be 0")
  expect_error(dmixture3p_cd(1, 0.5, kappa = -1), "non-negative")
  expect_error(dmixture3p_cd(1, 0.5, p_mem = 0.9, p_nt = 0.2), "smaller than one")
  # dmixture3p() would let p_nt go missing from the mixture with one item
  expect_error(dmixture3p_cd(1, 0.5, mu = 0, p_nt = 0.1), "p_nt must be 0")
  expect_error(rmixture3p_cd(5, 0.5, mu = 0, p_nt = 0.1), "p_nt must be 0")
  expect_no_error(dmixture3p_cd(1, 0.5, mu = 0, p_nt = 0))
})

test_that("dmixture3p_cd() integrates the retrieval mixture over the unbiased arc", {
  # at criterion 0 the arc half-width solves vM(hw | 0, kappa) = 1 / 2pi; each
  # component contributes its von Mises mass in the arc, guessing hw / pi
  kappa <- 6
  mu <- c(0.3, 2, -2.4)
  p_mem <- 0.6
  p_nt <- 0.3
  hw <- acos(log(besselI(kappa, 0)) / kappa)
  for (probe in c(0.3, 1.8, -2.9)) {
    p_same <- p_mem * vm_mass(probe - hw, probe + hw, mu[1], kappa) +
      p_nt / 2 * (vm_mass(probe - hw, probe + hw, mu[2], kappa) +
        vm_mass(probe - hw, probe + hw, mu[3], kappa)) +
      (1 - p_mem - p_nt) * hw / pi
    expect_equal(dmixture3p_cd(0, probe, mu = mu, kappa = kappa, p_mem = p_mem, p_nt = p_nt),
                 p_same, tolerance = 1e-10, label = paste("probe", probe))
    expect_equal(dmixture3p_cd(1, probe, mu = mu, kappa = kappa, p_mem = p_mem, p_nt = p_nt,
                               log = TRUE),
                 log1p(-p_same), tolerance = 1e-10)
  }
})

test_that("a probe at a non-target is harder to reject than a new one (intrusion cost)", {
  # the same distance from the target, once onto a non-target and once onto
  # empty space
  args <- list(mu = c(0, 2, 2.9), kappa = 8, p_mem = 0.6, p_nt = 0.3)
  intrusion <- do.call(dmixture3p_cd, c(list(1, probe = 2), args))
  new <- do.call(dmixture3p_cd, c(list(1, probe = -2), args))
  expect_lt(intrusion, new - 0.05)
  # without non-target retrievals the cost vanishes
  args$p_nt <- 0
  expect_equal(do.call(dmixture3p_cd, c(list(1, probe = 2), args)),
               do.call(dmixture3p_cd, c(list(1, probe = -2), args)), tolerance = 1e-12)
})

test_that("the versions reduce to one another where they should", {
  nt <- c(1.7, -2.2)
  args <- list(probe = c(-0.4, 0.9, 2.5), mu = 0.1, kappa = 6, criterion = 0.35, nt = nt)
  # slot at K >= set size and ptarget is the simple version with the same
  # target and swap probabilities and no guessing
  slot <- do.call(.mixture3p_cd_psame, c(args, list(
    w1 = 5, w2 = 0.7, set_size = 3, version = "slot", rich = FALSE
  )))
  marginal <- dmixture3p_cd(0, args$probe, mu = c(0.1, nt), kappa = 6, p_mem = 0.7, p_nt = 0.3,
                            criterion = 0.35)
  expect_equal(slot, marginal, tolerance = 1e-12)
  # slot averaging below one slot per item retrieves as the slot model; the
  # observers' p_s differ (K / set_size against K / set_size * ptarget), which
  # drops out at criterion 0
  expect_equal(
    do.call(.mixture3p_cd_psame, c(replace(args, "criterion", 0), list(
      w1 = 2.4, w2 = 0.7, set_size = 3, version = "slot_averaging", rich = FALSE
    ))),
    do.call(.mixture3p_cd_psame, c(replace(args, "criterion", 0), list(
      w1 = 2.4, w2 = 0.7, set_size = 3, version = "slot", rich = FALSE
    ))),
    tolerance = 1e-10
  )
})

test_that("knowledge changes the predictions of slot averaging only", {
  args <- list(probe = c(0.2, 1.1), mu = 0, kappa = 5, criterion = 0.4, set_size = 2, nt = 2)
  psame <- function(version, rich, w1, w2) {
    do.call(.mixture3p_cd_psame, c(args, list(w1 = w1, w2 = w2, version = version, rich = rich)))
  }
  expect_identical(psame("simple", TRUE, 0.5, -0.3), psame("simple", FALSE, 0.5, -0.3))
  expect_identical(psame("slot", TRUE, 1.3, 0.8), psame("slot", FALSE, 1.3, 0.8))
  # K / set_size = 1.4: one slot held, a second with probability 0.4
  diff <- psame("slot_averaging", TRUE, 2.8, 0.8) - psame("slot_averaging", FALSE, 2.8, 0.8)
  expect_true(all(abs(diff) > 1e-4))
})

test_that("the observer that simulates decides as the likelihood integrates", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(136)
  n <- 4e4
  nt <- c(1.7, -2.2)
  cases <- list(
    simple = list(w1 = 0.6, w2 = -0.2, set_size = 3, rich = FALSE),
    simple_ss1 = list(w1 = 0.6, w2 = -0.2, set_size = 1, rich = FALSE),
    slot = list(w1 = 2.2, w2 = 0.75, set_size = 3, rich = FALSE),
    sa_empty = list(w1 = 2.2, w2 = 0.75, set_size = 3, rich = FALSE),
    sa_limited = list(w1 = 2.7, w2 = 0.75, set_size = 2, rich = FALSE),
    sa_rich = list(w1 = 2.7, w2 = 0.75, set_size = 2, rich = TRUE)
  )
  versions <- c(simple = "simple", simple_ss1 = "simple", slot = "slot", sa_empty = "slot_averaging",
                sa_limited = "slot_averaging", sa_rich = "slot_averaging")
  for (name in names(cases)) {
    case <- cases[[name]]
    for (criterion in c(0, 0.6, -0.5)) {
      for (probe in c(0.15, 1.7, -1)) {
        p_change <- 1 - .mixture3p_cd_psame(
          probe, 0.1, 6, case$w1, case$w2, criterion, case$set_size, nt, versions[[name]], case$rich
        )
        drawn <- .mixture3p_cd_observe(
          rep(probe, n), 0.1, matrix(nt, n, 2, byrow = TRUE), case$set_size, 6,
          case$w1, case$w2, criterion, versions[[name]], case$rich
        )
        expect_lt(abs(mean(drawn) - p_change), 5 * sqrt(p_change * (1 - p_change) / n) + 1e-3,
                  label = paste(name, criterion, probe))
      }
    }
  }
  # rmixture3p_cd() draws through the same observer
  draws <- rmixture3p_cd(n, probe = 2, mu = c(0, 2, -2.5), kappa = 8, p_mem = 0.7, p_nt = 0.1,
                         criterion = 0.3)
  p <- dmixture3p_cd(1, probe = 2, mu = c(0, 2, -2.5), kappa = 8, p_mem = 0.7, p_nt = 0.1,
                     criterion = 0.3)
  expect_lt(abs(mean(draws) - p), 5 * sqrt(p * (1 - p) / n))
})

test_that("check_data() refuses a missing set size and warns where K cannot be tested", {
  data <- mixture3p_cd_data()
  data$set_size <- as.numeric(as.character(data$set_size))
  formula <- cd3p_formulas$slot
  missing_ss <- data
  missing_ss$set_size[3] <- NA
  expect_error(check_data(cd3p_model("slot"), missing_ss, formula), "contains missing values")
  one_ss <- data[data$set_size == 4, ]
  four <- function(version) {
    mixture3p_cd("change", "probe", nt_features = paste0("col_nt", 1:3),
                 set_size = "set_size", version = version)
  }
  expect_warning(check_data(four("slot"), one_ss, formula), "does not vary")
  expect_no_warning(check_data(four("simple"), one_ss, cd3p_formulas$simple))
  attr(one_ss, "prediction_newdata") <- TRUE
  expect_no_warning(check_data(four("slot"), one_ss, formula))
})

test_that("the swap parameter is flagged as undetermined at set size 1", {
  data <- mixture3p_cd_data()
  expect_warning(
    bmm(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size), data, cd3p_model(),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "thetant"
  )
  expect_warning(
    bmm(bmf(kappa ~ 1, K ~ 1, ptarget ~ 0 + set_size), data, cd3p_model("slot"),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "ptarget"
  )
})

test_that("an estimated criterion is announced with the bound", {
  expect_message(
    bmm(bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1, criterion ~ 1), mixture3p_cd_data(), cd3p_model(),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "criterion is estimated"
  )
})

test_that("the generated program passes the set size as an integer, the probe and the non-targets as reals", {
  data <- mixture3p_cd_data()
  code <- stancode(cd3p_formulas$simple, data, cd3p_model())
  expect_match(code, "array[N] int vint1", fixed = TRUE)
  expect_match(
    code,
    "mixture3p_cd_simple_lpmf(Y[n] | mu[n], kappa[n], thetat[n], thetant[n], criterion[n], vint1[n], vreal1[n], vreal2[n]",
    fixed = TRUE
  )
  code <- stancode(cd3p_formulas$slot_averaging, data, cd3p_model("slot_averaging", knowledge = "rich"))
  # rich reaches the core as the literal 1, after the packed non-targets
  expect_match(code, "ss, probe, to_vector({nt1, nt2, nt3, nt4, nt5, nt6, nt7}), 1, circmix_logk",
               fixed = TRUE)
})

test_that("log_lik(), posterior_predict() and posterior_epred() evaluate the same P('same')", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(552)
  n <- 4e4
  nt <- c(1.7, -2.2)
  draws <- list(mu = 0.1, kappa = c(4, 7), thetat = c(0.6, -0.2), thetant = c(-0.5, 0.4),
                K = c(2.2, 3.4), ptarget = c(0.6, 0.85), criterion = c(0, 0.3))
  for (version in names(cd3p_formulas)) {
    weights <- if (version == "simple") c("thetat", "thetant") else c("K", "ptarget")
    dpars <- lapply(draws[c("mu", "kappa", weights, "criterion")], function(x) matrix(rep_len(x, n), n, 1))
    prep <- structure(
      list(dpars = dpars, data = list(Y = 1L, vint1 = 3L, vreal1 = 0.8, vreal2 = nt[1], vreal3 = nt[2]),
           ndraws = n, nobs = 1, family = list(knowledge = "rich")),
      class = "brmsprep"
    )
    p_change <- 1 - .mixture3p_cd_psame(
      0.8, 0.1, c(4, 7), draws[[weights[1]]], draws[[weights[2]]], c(0, 0.3), 3, nt, version, TRUE
    )
    expect_equal(.mixture3p_cd_log_lik(version)(1, prep)[1:2], log(p_change), tolerance = 1e-12)
    expect_equal(as.vector(.mixture3p_cd_posterior_epred(version)(prep))[1:2], p_change,
                 tolerance = 1e-12)
    drawn <- .mixture3p_cd_posterior_predict(version)(1, prep)
    expect_lt(abs(mean(drawn) - mean(p_change)), 0.015)
  }
})

test_that("simulate_data() draws the observer's decision on the design's probes", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(649)
  # two items, the probe on the non-target: P("change") by integrate() over
  # the unbiased arc, with thetat = thetant = 0 (a third each)
  n <- 6e4
  design <- data.frame(probe = 2, nt1 = 2, set_size = 2)[rep(1, n), ]
  model <- mixture3p_cd("change", "probe", nt_features = "nt1", set_size = "set_size")
  data <- simulate_data(model, list(kappa = log(6), thetat = 0, thetant = 0), design)
  hw <- acos(log(besselI(6, 0)) / 6)
  p_same <- (vm_mass(2 - hw, 2 + hw, 0, 6) + vm_mass(2 - hw, 2 + hw, 2, 6) + hw / pi) / 3
  expect_lt(abs(mean(data$change) - (1 - p_same)), 5 * sqrt(0.25 / n))
  expect_error(simulate_data(model, list(kappa = 1, thetat = 0, thetant = 0), n_trials = 5), "probe")
})

test_that("the generated likelihood passes each parameter to its place in the core, and Stan matches R", {
  skip_unless_stan_tests()
  skip_on_cran()
  data <- mixture3p_cd_data()
  data <- data[as.numeric(as.character(data$set_size)) <= 3, ]
  model_for <- function(version, knowledge) {
    mixture3p_cd("change", "probe", nt_features = c("col_nt1", "col_nt2"),
                 set_size = "set_size", version = version, knowledge = knowledge)
  }
  # one item and three; nothing held (K / ss < 1) and a fractional slot held
  # under each knowledge. Every value differs, so a swapped slot would show
  cases <- list(
    simple = list(version = "simple", ss = 3),
    simple_ss1 = list(version = "simple", ss = 1),
    slot = list(version = "slot", ss = 3),
    sa_empty = list(version = "slot_averaging", ss = 3),
    sa_limited = list(version = "slot_averaging", ss = 2),
    sa_rich = list(version = "slot_averaging", ss = 2, knowledge = "rich")
  )
  values <- c(mu = 0.1, kappa = 6, thetat = 0.4, thetant = -0.3, K = 2.6, ptarget = 0.7,
              criterion = 0.35, probe = 0.9, nt1 = 1.7, nt2 = -2.2)
  wrappers <- character(0)
  calls <- character(0)
  r_side <- numeric(0)
  for (name in names(cases)) {
    case <- cases[[name]]
    knowledge <- case$knowledge %||% "limited"
    model <- model_for(case$version, knowledge)
    formula <- cd3p_formulas[[case$version]]
    family <- paste0("mixture3p_cd_", case$version)
    configured <- suppressWarnings(configure_model(model, check_data(model, data, formula), formula))
    code <- configured$stanvars[[paste0(family, "_lpmf")]]$scode
    code <- gsub(paste0(family, "_lpmf"), paste0("w_", name, "_lpmf"), code, fixed = TRUE)
    wrappers <- c(wrappers, code)
    signature <- sub(".*_lpmf\\((.*?)\\).*", "\\1", strsplit(code, "\n")[[1]][1])
    arg_names <- sub(".* ", "", strsplit(signature, ", ")[[1]])
    arg_values <- ifelse(arg_names %in% names(values), format(values[arg_names], digits = 17), arg_names)
    arg_values[arg_names == "y"] <- "0"
    arg_values[arg_names == "ss"] <- as.character(case$ss)
    calls <- c(calls, glue::glue(
      "out[{length(calls) + 1}] = w_{name}_lpmf({arg_values[1]} | {paste(arg_values[-1], collapse = ', ')});"
    ))
    weights <- if (case$version == "simple") values[c("thetat", "thetant")] else values[c("K", "ptarget")]
    r_side[name] <- .cd_bernoulli_ld(0, .mixture3p_cd_psame(
      values[["probe"]], values[["mu"]], values[["kappa"]], weights[[1]], weights[[2]],
      values[["criterion"]], case$ss, values[c("nt1", "nt2")], case$version, knowledge == "rich"
    ))
  }
  chunks <- c("circmix_functions.stan", "cd_functions.stan", "circmix_cd_functions.stan",
              "mixture3p_functions.stan", "mixture3p_cd_functions.stan")
  program <- paste0(
    "functions {\n",
    paste(vapply(chunks, function(f) read_lines2(system.file("stan_chunks", f, package = "bmm")), ""),
          collapse = "\n"),
    "\n", paste(wrappers, collapse = "\n"), "}\n",
    "data { int NK; vector[NK] circmix_logk; vector[NK] circmix_dlogk;",
    " real circmix_logJ_min; real circmix_dlogJ; vector[32] cd_gl_x; vector[32] cd_gl_w; }\n",
    "generated quantities { vector[", length(calls), "] out;\n",
    paste(calls, collapse = "\n"), "\n}\n"
  )
  tab <- .circmix_kappa_table()
  gl <- .cd_gl_rule()
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = list(
      NK = length(tab$logkappa), circmix_logk = tab$logkappa, circmix_dlogk = tab$dlogkappa,
      circmix_logJ_min = tab$logJ_min, circmix_dlogJ = tab$dlogJ, cd_gl_x = gl$x, cd_gl_w = gl$w
    ),
    fixed_param = TRUE, iter_sampling = 1, chains = 1, refresh = 0,
    show_messages = FALSE, sig_figs = 18
  )
  stan_side <- as.numeric(fit$draws("out"))
  expect_equal(stan_side, unname(r_side), tolerance = 1e-10)
  expect_identical(length(unique(round(stan_side, 6))), length(cases))
})

test_that("the threaded program parses", {
  skip_unless_stan_tests()
  skip_on_cran()
  data <- mixture3p_cd_data()
  # the capacity versions name a parameter ptarget, which brms declares itself
  # in a threaded program, so they fail to parse there, as mixture3p()'s do;
  # they join this loop once the circmix stack renames the parameter
  for (version in "simple") {
    code <- suppressWarnings(stancode(cd3p_formulas[[version]], data, cd3p_model(version),
                                      threads = brms::threading(2)))
    expect_match(code, "partial_log_lik", fixed = TRUE)
    stan_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
    expect_true(stan_model$check_syntax(quiet = TRUE), label = version)
  }
})
