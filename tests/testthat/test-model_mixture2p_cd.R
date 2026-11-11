test_that("mixture2p_cd() validates its arguments", {
  expect_error(mixture2p_cd("r"), "probe")
  expect_error(mixture2p_cd("r", "p", version = "full"), "should be one of")
  expect_error(mixture2p_cd("r", "p", knowledge = "partial"), "should be one of")
  expect_error(mixture2p_cd("r", "p", version = "slot"), "set_size argument is required")
  expect_error(mixture2p_cd("r", "p", version = "slot_averaging"), "set_size argument is required")
  expect_error(mixture2p_cd("r", "p", variable_precision = TRUE), "planned for bmm 1.6.0")
  expect_error(mixture2p_cd("r", "p", variable_precision = NA), "TRUE or FALSE")
  # a vp_nodes copied from mixture2p() is reported, not silently dropped
  expect_warning(mixture2p_cd("r", "p", vp_nodes = 81L), "vp_nodes")
  # the argument of the first drafts of this model
  expect_warning(mixture2p_cd("r", "p", target = "t"), "`target` is not an argument")
  expect_error(mixture2p_cd("r", "p", links = list(criterion = "log")), "'criterion' cannot be changed")
})

test_that("knowledge = 'rich' warns in the versions with one precision state", {
  expect_warning(mixture2p_cd("r", "p", knowledge = "rich"), "no effect in the 'simple' version")
  expect_warning(
    mixture2p_cd("r", "p", set_size = "s", version = "slot", knowledge = "rich"),
    "no effect in the 'slot' version"
  )
  expect_no_warning(
    mixture2p_cd("r", "p", set_size = "s", version = "slot_averaging", knowledge = "rich")
  )
  expect_warning(
    dmixture2p_cd(1, 0.5, K = 2, set_size = 4, version = "slot", knowledge = "rich"),
    "no effect in the 'slot' version"
  )
  expect_warning(rmixture2p_cd(5, 0.5, knowledge = "rich"), "no effect in the 'simple' version")
})

test_that("mixture2p_cd() takes mixture2p's parameters, priors, links and inits, plus criterion", {
  for (version in c("simple", "slot", "slot_averaging")) {
    set_size <- if (version != "simple") "ss"
    cd <- mixture2p_cd("r", "p", set_size = set_size, version = version)
    de <- mixture2p("y", set_size = set_size, version = version)
    expect_identical(names(cd$parameters), c(names(de$parameters), "criterion"), label = version)
    for (field in c("links", "default_priors", "init_ranges")) {
      expect_identical(cd[[field]][names(de[[field]])], de[[field]], label = paste(version, field))
    }
    expect_identical(cd$fixed_parameters, list(mu = 0, criterion = 0), label = version)
    expect_identical(cd$links$criterion, "identity")
    expect_s3_class(cd, c("bmmodel", "change_detection", "mixture2p_cd", paste0("mixture2p_cd_", version)),
                    exact = TRUE)
  }
  model <- mixture2p_cd("r", "p", set_size = "s", version = "slot_averaging", knowledge = "rich")
  expect_identical(model$knowledge, "rich")
  expect_identical(model$task, "Detection, recognition and confidence judgments")
  expect_identical(model$domain, "Working memory")
  expect_true("mixture2p_cd" %in% bmm_models())
})

test_that("criterion's init range starts on the linear part of its bound", {
  # the bound is per trial; at the corners of the inits of the memory
  # parameters and set sizes 1 to 8, the ends of criterion's init range keep
  # more than half of their value
  init <- mixture2p_cd("r", "p", "t")$init_ranges
  p_s <- c(init$thetat, pmin(1, outer(c(1.9, 3.8), 1:8, "/")))
  corners <- expand.grid(p_s = p_s, kappa = init$kappa, criterion = init$criterion)
  bounded <- .circmix_cd_bounded_criterion(corners$criterion, corners$p_s, corners$kappa, 0)
  expect_gt(min(bounded / corners$criterion), 0.5)
})

test_that("dmixture2p_cd() refuses what the constructor refuses", {
  expect_error(dmixture2p_cd(1, 0.5, version = "slot"), "set_size argument is required")
  expect_error(rmixture2p_cd(5, 0.5, version = "slot_averaging"), "set_size argument is required")
  # `log` passed by position, after the memory parameters, lands on criterion
  expect_error(dmixture2p_cd(1, 0.5, 0, 5, 0.6, TRUE), "Did you mean log = TRUE")
  expect_error(rmixture2p_cd(5, 0.5, criterion = TRUE), "criterion must be numeric")
  expect_true(is.na(dmixture2p_cd(1, 0.5, criterion = NA)))
  for (n in list(2.5, c(2, 3), 0, "3")) {
    expect_error(rmixture2p_cd(n, 0.5), "single positive whole number")
  }
  expect_error(dmixture2p_cd(1, 0.5, knowledge = "partial"), "should be one of")
  expect_error(dmixture2p_cd(1, 0.5, version = "full"), "should be one of")
  expect_error(dmixture2p_cd(2, 0.5), "response must be 0")
  expect_error(dmixture2p_cd(NA, 0.5), "response must be 0")
  expect_error(dmixture2p_cd(1, 0.5, K = 2, set_size = 2.5, version = "slot"), "whole numbers")
  expect_error(rmixture2p_cd(5, 0.5, K = 2, set_size = c(2, 0), version = "slot"), "whole numbers")
  expect_error(dmixture2p_cd(1, 0.5, kappa = -1), "non-negative")
  expect_error(dmixture2p_cd(1, 0.5, p_mem = 1.2), "smaller than one")
  expect_warning(dmixture2p_cd(1, 0.5, K = 2), "ignores 'K'")
  expect_warning(dmixture2p_cd(1, 0.5, p_mem = 0.5, K = 2, set_size = 3, version = "slot"), "ignores 'p_mem'")
})

test_that("dmixture2p_cd() integrates the retrieval mixture over the unbiased arc", {
  # P("same") = p_mem times the von Mises mass within hw of the probe plus
  # (1 - p_mem) hw / pi, with cos(hw) = log I0(kappa) / kappa, by integrate()
  expect_equal(
    dmixture2p_cd(0, c(0, 1.2, -2.5), kappa = 6, p_mem = 0.7),
    c(0.730333741913013, 0.192678499287274, 0.075905318617164),
    tolerance = 1e-10
  )
  expect_equal(
    dmixture2p_cd(1, 0.4, kappa = 6, K = 2.6, set_size = 4, version = "slot"),
    1 - 0.622355441908816,
    tolerance = 1e-10
  )
  expect_equal(
    dmixture2p_cd(c(0, 1), 0.4, mu = 0.3, kappa = 6, p_mem = 0.7, log = TRUE),
    log(dmixture2p_cd(c(0, 1), 0.4, mu = 0.3, kappa = 6, p_mem = 0.7))
  )
  # an item held in every version is the slot version at K >= set_size
  expect_equal(
    dmixture2p_cd(0, 0.4, kappa = 6, K = 4, set_size = 4, version = "slot"),
    dmixture2p_cd(0, 0.4, kappa = 6, p_mem = 1)
  )
})

test_that("knowledge changes the predictions of slot averaging only", {
  probe <- c(0.2, 1, 2)
  for (version in c("simple", "slot")) {
    args <- list(response = 0, probe = probe, kappa = 4, criterion = 0.4, K = 2.6, set_size = 2,
                 version = version)
    if (version == "simple") args[c("K", "set_size")] <- NULL
    expect_warning(rich <- do.call(dmixture2p_cd, c(args, knowledge = "rich")), "no effect")
    expect_identical(do.call(dmixture2p_cd, c(args, knowledge = "limited")), rich, label = version)
  }
  # 2.6 slots over 2 items: one or two slots per item, two precision states
  limited <- dmixture2p_cd(0, probe, kappa = 4, K = 2.6, set_size = 2, version = "slot_averaging")
  rich <- dmixture2p_cd(0, probe, kappa = 4, K = 2.6, set_size = 2, version = "slot_averaging",
                        knowledge = "rich")
  expect_gt(min(abs(limited - rich)), 1e-3)
})

test_that("rmixture2p_cd() draws the observer that dmixture2p_cd() integrates", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(136)
  n <- 1e5
  # nothing held, a fractional slot held and a whole K / set size, under
  # either knowledge, and both sides of the criterion
  cells <- expand.grid(
    version = c("simple", "slot", "slot_averaging"), knowledge = c("limited", "rich"),
    criterion = c(0, 0.6, -0.4), set_size = c(5, 2, 1), stringsAsFactors = FALSE
  )
  # the rich observer differs, and is accepted without a warning, only where
  # precision varies
  cells <- subset(cells, knowledge == "limited" | version == "slot_averaging")
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    args <- list(probe = 0.9, mu = 0.1, kappa = 5, criterion = cell$criterion,
                 version = cell$version, knowledge = cell$knowledge)
    args <- c(args, if (cell$version == "simple") list(p_mem = 0.7) else list(K = 2.6, set_size = cell$set_size))
    expected <- do.call(dmixture2p_cd, c(list(response = 1), args))
    drawn <- mean(do.call(rmixture2p_cd, c(list(n = n), args)))
    expect_lt(abs(drawn - expected), 5 * sqrt(expected * (1 - expected) / n))
  }
})

test_that("check_data() refuses missing values and reads the set size as mixture2p does", {
  data <- mixture2p_cd_data()
  model <- mixture2p_cd("change", "probe", set_size = "set_size", version = "slot")
  formula <- bmf(kappa ~ 1, K ~ 1)
  for (column in c("change", "probe", "set_size")) {
    broken <- data
    broken[[column]][3] <- NA
    expect_error(check_data(model, broken, formula), "missing|NA", label = column)
  }
  fractional <- transform(data, set_size = as.numeric(as.character(set_size)) + 0.5)
  expect_error(check_data(model, fractional, formula), "whole")
  checked <- check_data(model, data, formula)
  expect_identical(checked$ss_numeric, as.numeric(as.character(data$set_size)))
  expect_identical(built_data_columns(model), "ss_numeric")
  expect_identical(built_data_columns(mixture2p_cd("change", "probe")), character(0))

  single <- subset(data, set_size == 4)
  expect_warning(check_data(model, single, formula), "cannot test the capacity rule")
  attr(single, "prediction_newdata") <- TRUE
  expect_no_warning(check_data(model, single, formula))
})

test_that("an estimated criterion is announced with the bound", {
  data <- mixture2p_cd_data()
  formula <- bmf(kappa ~ 1, thetat ~ 1, criterion ~ 1)
  model <- check_model(mixture2p_cd("change", "probe"), data, formula)
  expect_message(check_formula(model, data, formula), "Decision criterion")
  fixed <- bmf(kappa ~ 1, thetat ~ 1)
  model <- check_model(mixture2p_cd("change", "probe"), data, fixed)
  expect_no_message(check_formula(model, data, fixed), message = "criterion is estimated")
})

test_that("the capacity versions pass the set size as an integer and the probe as a real", {
  data <- mixture2p_cd_data()
  code <- stancode(bmf(kappa ~ 1, thetat ~ 1), data, mixture2p_cd("change", "probe"))
  expect_match(code, "mixture2p_cd_simple_lpmf(Y[n] | mu[n], kappa[n], thetat[n], criterion[n], vreal1[n]",
               fixed = TRUE)
  model <- mixture2p_cd("change", "probe", set_size = "set_size", version = "slot_averaging",
                        knowledge = "rich")
  code <- stancode(bmf(kappa ~ 1, K ~ 1), data, model)
  expect_match(code, "array[N] int vint1", fixed = TRUE)
  expect_match(code, "mixture2p_cd_slot_averaging_lpmf(Y[n] | mu[n], kappa[n], K[n], criterion[n], vint1[n], vreal1[n]",
               fixed = TRUE)
  # rich reaches the core as the literal 1
  expect_match(code, "K, criterion, ss, probe, 1, circmix_logk", fixed = TRUE)
})

test_that("log_lik(), posterior_predict() and posterior_epred() evaluate the same P('same')", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(552)
  n <- 4e4
  draws <- list(mu = 0.1, kappa = c(4, 7), thetat = c(0.6, 0.8), K = c(2.2, 3.4), criterion = c(0, 0.3))
  for (version in c("simple", "slot", "slot_averaging")) {
    weight <- if (version == "simple") "thetat" else "K"
    knowledge <- if (version == "slot_averaging") "rich" else "limited"
    dpars <- lapply(draws[c("mu", "kappa", weight, "criterion")], function(x) matrix(rep_len(x, n), n, 1))
    prep <- structure(
      list(dpars = dpars, data = list(Y = 1L, vint1 = 3L, vreal1 = 0.8), ndraws = n, nobs = 1,
           family = list(knowledge = knowledge)),
      class = "brmsprep"
    )
    weight_args <- if (version == "simple") list(p_mem = c(0.6, 0.8)) else list(K = c(2.2, 3.4), set_size = 3)
    p_change <- do.call(dmixture2p_cd, c(
      list(1, 0.8, mu = 0.1, kappa = c(4, 7), criterion = c(0, 0.3), version = version, knowledge = knowledge),
      weight_args
    ))
    expect_equal(.mixture2p_cd_log_lik(version)(1, prep)[1:2], log(p_change), tolerance = 1e-12)
    expect_equal(as.vector(.mixture2p_cd_posterior_epred(version)(prep))[1:2], p_change, tolerance = 1e-12)
    drawn <- .mixture2p_cd_posterior_predict(version)(1, prep)
    expect_lt(abs(mean(drawn) - mean(p_change)), 0.015)
  }
})

test_that("simulate_data() draws the observer's decision on the design's probes", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(649)
  # kappa e, P(memory) 0.5, the probe 1 from the target: P("change") by
  # integrate() over the unbiased arc
  design <- data.frame(probe = rep(1, 1e5))
  data <- simulate_data(mixture2p_cd("change", "probe"), list(kappa = 1, thetat = 0), design)
  expect_lt(abs(mean(data$change) - 0.571702760207019), 5 * sqrt(0.25 / nrow(data)))
  expect_error(simulate_data(mixture2p_cd("change", "probe"), list(kappa = 1, thetat = 0),
                             n_trials = 5), "probe")
})

test_that("the generated likelihood passes each parameter to its place in the core, and Stan matches R", {
  skip_unless_stan_tests()
  skip_on_cran()
  data <- mixture2p_cd_data()
  # nothing held (K / ss < 1) and a fractional slot held, under each knowledge;
  # below K / ss = 1 slot averaging is the slot model, so the two get
  # different set sizes
  cases <- list(
    simple = list(version = "simple", ss = 4, formula = bmf(kappa ~ 1, thetat ~ 1)),
    slot = list(version = "slot", ss = 3, formula = bmf(kappa ~ 1, K ~ 1)),
    sa_empty = list(version = "slot_averaging", ss = 4, formula = bmf(kappa ~ 1, K ~ 1)),
    sa_limited = list(version = "slot_averaging", ss = 2, formula = bmf(kappa ~ 1, K ~ 1)),
    sa_rich = list(version = "slot_averaging", ss = 2, knowledge = "rich", formula = bmf(kappa ~ 1, K ~ 1))
  )
  values <- c(y = 0, mu = 0.1, kappa = 6, thetat = 0.75, K = 2.6, criterion = 0.35, probe = 0.9)
  wrappers <- character(0)
  calls <- character(0)
  r_side <- numeric(0)
  for (name in names(cases)) {
    case <- cases[[name]]
    knowledge <- case$knowledge %||% "limited"
    model <- mixture2p_cd("change", "probe",
      set_size = if (case$version != "simple") "set_size", version = case$version, knowledge = knowledge
    )
    family <- paste0("mixture2p_cd_", case$version)
    configured <- configure_model(model, check_data(model, data, case$formula), case$formula)
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
    weight_args <- if (case$version == "simple") {
      list(p_mem = values[["thetat"]])
    } else {
      list(K = values[["K"]], set_size = case$ss)
    }
    r_side[name] <- do.call(dmixture2p_cd, c(list(
      0, values[["probe"]], mu = values[["mu"]], kappa = values[["kappa"]],
      criterion = values[["criterion"]], version = case$version, knowledge = knowledge, log = TRUE
    ), weight_args))
  }
  chunks <- c("circmix_functions.stan", "cd_functions.stan", "circmix_cd_functions.stan",
              "mixture2p_functions.stan", "mixture2p_cd_functions.stan")
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
  # every case evaluates a different P("same"), so a swapped slot would show
  expect_identical(length(unique(round(stan_side, 6))), length(cases))
})

test_that("the threaded program of every version parses", {
  skip_unless_stan_tests()
  skip_on_cran()
  data <- mixture2p_cd_data()
  for (version in c("simple", "slot", "slot_averaging")) {
    model <- mixture2p_cd("change", "probe",
      set_size = if (version != "simple") "set_size", version = version
    )
    formula <- if (version == "simple") bmf(kappa ~ 1, thetat ~ 1) else bmf(kappa ~ 1, K ~ 1)
    code <- stancode(formula, data, model, threads = brms::threading(2))
    expect_match(code, "partial_log_lik", fixed = TRUE)
    stan_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
    expect_true(stan_model$check_syntax(quiet = TRUE), label = version)
  }
})
