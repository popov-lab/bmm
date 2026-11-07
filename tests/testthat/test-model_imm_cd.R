cdimm_model <- function(version = "full", ...) {
  imm_cd("change", "probe", "target", nt_features = paste0("col_nt", 1:7),
         nt_distances = paste0("dist_nt", 1:7), set_size = "set_size", version = version, ...)
}

cdimm_formulas <- list(
  full = bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1),
  bsc = bmf(kappa ~ 1, c ~ 1, s ~ 1),
  abc = bmf(kappa ~ 1, c ~ 1, a ~ 1)
)

# von Mises mass on [lo, hi] by adaptive integration, independent of the
# Gauss-Legendre arc of the package
vm_mass <- function(lo, hi, mu, kappa) {
  stats::integrate(
    function(x) exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, 0)),
    lo, hi, rel.tol = 1e-12
  )$value
}

test_that("imm_cd() validates its arguments", {
  expect_error(imm_cd("r", "p", "t", nt_features = "nt", set_size = 2), "nt_distances")
  expect_no_error(imm_cd("r", "p", "t", nt_features = "nt", set_size = 2, version = "abc"))
  expect_error(cdimm_model(version = "slot"), "should be one of")
  expect_error(cdimm_model(knowledge = "partial"), "should be one of")
  expect_error(cdimm_model(variable_precision = TRUE), "planned for bmm 1.6.0")
  expect_error(cdimm_model(variable_precision = NA), "TRUE or FALSE")
  # a vp_nodes copied from imm() is reported, not silently dropped
  expect_warning(cdimm_model(vp_nodes = 81L), "vp_nodes")
  # b is fixed on its log link, as in imm()
  expect_error(cdimm_model(links = list(b = "identity")), "Links can be set for")
})

test_that("imm_cd() takes imm's parameters, priors, links and inits, plus criterion", {
  for (version in names(cdimm_formulas)) {
    cd <- cdimm_model(version)
    de <- imm("y", nt_features = paste0("col_nt", 1:7), nt_distances = paste0("dist_nt", 1:7),
              set_size = "set_size", version = version)
    expect_identical(names(cd$parameters), c(names(de$parameters), "criterion"), label = version)
    for (field in c("links", "default_priors", "init_ranges")) {
      expect_identical(cd[[field]][names(de[[field]])], de[[field]], label = paste(version, field))
    }
    expect_identical(cd$fixed_parameters, list(mu = 0, b = 0, criterion = 0), label = version)
    expect_identical(cd$links$criterion, "identity")
    expect_identical(cd$other_vars$nt_distances, if (version != "abc") paste0("dist_nt", 1:7))
    expect_s3_class(
      cd, c("bmmodel", "change_detection", "non_targets", "imm_cd", paste0("imm_cd_", version)),
      exact = TRUE
    )
  }
  expect_identical(cdimm_model(knowledge = "rich")$knowledge, "rich")
  expect_identical(cdimm_model()$knowledge, "limited")
  expect_identical(cdimm_model()$task, "Detection, recognition and confidence judgments")
  expect_identical(cdimm_model()$domain, "Working memory")
  expect_true("imm_cd" %in% bmm_models())
  # the constructor's arguments shared with the reproduction twin keep its
  # order; knowledge sits after variable_precision, before links (CD-D1 F1)
  expect_identical(
    names(formals(imm_cd)),
    c("response", "probe", "target", "nt_features", "nt_distances", "set_size", "regex",
      "version", "variable_precision", "knowledge", "links", "...")
  )
})

test_that("criterion's init range starts on the linear part of its bound", {
  # the bound is per trial and narrows as p_s, the target's retrieval weight,
  # falls; over the corners of the inits of c, a and s, set sizes 1 to 8 and
  # every non-target at the target's distance 0 (the smallest p_s), the ends
  # of criterion's init range keep more than half of their value. The
  # narrowest corner is "bsc" at set size 8: p_s 0.100 keeps 0.504
  inits <- cdimm_model()$init_ranges
  corners <- expand.grid(c = inits$c, a = inits$a, s = inits$s, ss = 1:8)
  kept <- vapply(names(cdimm_formulas), function(version) {
    p_s <- unlist(lapply(seq_len(nrow(corners)), function(j) {
      exp(.imm_log_weights(corners$c[j], corners$a[j], corners$s[j], 1, corners$ss[j], rep(0, 7),
                           version)$logw[, 1])
    }))
    grid <- expand.grid(p_s = p_s, kappa = inits$kappa, criterion = inits$criterion)
    min(.circmix_cd_bounded_criterion(grid$criterion, grid$p_s, grid$kappa, 0) / grid$criterion)
  }, numeric(1))
  expect_true(all(kept > 0.5))
  expect_equal(kept[["bsc"]], 0.504, tolerance = 1e-3)
})

test_that("dimm_cd() and rimm_cd() refuse what they cannot evaluate", {
  expect_error(dimm_cd(1, 0.5, variable_precision = TRUE), "planned for bmm 1.6.0")
  expect_error(rimm_cd(5, 0.5, variable_precision = TRUE), "planned for bmm 1.6.0")
  expect_error(dimm_cd(1, 0.5, variable_precision = NA), "TRUE or FALSE")
  expect_error(dimm_cd(2, 0.5), "response must be 0")
  expect_error(dimm_cd(NA, 0.5), "response must be 0")
  expect_error(dimm_cd(1, 0.5, kappa = -1), "non-negative")
  expect_error(dimm_cd(1, 0.5, s = -1), "non-negative")
  expect_error(rimm_cd(5, 0.5, dist = c(0, -0.5, 2)), "positive")
  expect_error(dimm_cd(1, 0.5, mu = c(0, 2)), "does not match the distances")
})

test_that("dimm_cd() integrates the retrieval mixture over the unbiased arc", {
  # at criterion 0 the arc half-width solves vM(hw | 0, kappa) = 1 / 2pi; each
  # item contributes its von Mises mass in the arc in proportion to its
  # activation c exp(-s dist) + a, the background b contributes hw / pi
  kappa <- 6
  mu <- c(0.3, 2, -2.4)
  dist <- c(0, 0.5, 2)
  activation <- 4 * exp(-1.5 * dist) + 0.5
  weights <- activation / (sum(activation) + 1.2)
  hw <- acos(log(besselI(kappa, 0)) / kappa)
  for (probe in c(0.3, 1.8, -2.9)) {
    p_same <- sum(weights * vapply(mu, function(m) vm_mass(probe - hw, probe + hw, m, kappa), 0)) +
      (1 - sum(weights)) * hw / pi
    expect_equal(dimm_cd(0, probe, mu = mu, dist = dist, c = 4, a = 0.5, s = 1.5, b = 1.2,
                         kappa = kappa),
                 p_same, tolerance = 1e-10, label = paste("probe", probe))
    expect_equal(dimm_cd(1, probe, mu = mu, dist = dist, c = 4, a = 0.5, s = 1.5, b = 1.2,
                         kappa = kappa, log = TRUE),
                 log1p(-p_same), tolerance = 1e-10)
  }
})

test_that("a probe at a non-target is harder to reject the nearer that non-target was", {
  # the same probe and the same non-target feature, once near the target and
  # once far from it, against a probe of a new feature
  p_change <- function(probe, dist2) {
    dimm_cd(1, probe, mu = c(0, 2, -2.5), dist = c(0, dist2, 3), c = 4, a = 0.3, s = 2, kappa = 8)
  }
  expect_lt(p_change(2, 0.2), p_change(2, 1.5) - 0.05)
  expect_lt(p_change(2, 1.5), p_change(-1.2, 1.5))
  # without any non-target activation the cost vanishes
  args <- list(mu = c(0, 2, -2.5), dist = c(0, 0.2, 3), c = 4, a = 0, s = 1e3, kappa = 8)
  expect_equal(do.call(dimm_cd, c(list(1, probe = 2), args)),
               do.call(dimm_cd, c(list(1, probe = -2), args)), tolerance = 1e-12)
})

test_that("the versions reduce to the full model where they should", {
  nt <- c(1.7, -2.2)
  dist <- c(0.6, 1.8)
  args <- list(probe = c(-0.4, 0.9, 2.5), mu = 0.1, kappa = 6, b = 1, criterion = 0.35,
               set_size = 3, nt = nt, dist = dist)
  psame <- function(version, c, a, s) {
    do.call(.imm_cd_psame, c(args, list(c = c, a = a, s = s, version = version)))
  }
  # dimm_cd() takes the full model with the target first
  expect_equal(
    psame("full", 3, 0.7, 1.2),
    dimm_cd(0, args$probe, mu = c(0.1, nt), dist = c(0, dist), c = 3, a = 0.7, s = 1.2, b = 1,
            kappa = 6, criterion = 0.35),
    tolerance = 1e-12
  )
  # bsc is the full model without cue-independent activation, abc the full
  # model whose non-targets receive none of the cue's activation
  expect_equal(psame("bsc", 3, NULL, 1.2), psame("full", 3, 1e-300, 1.2), tolerance = 1e-12)
  expect_equal(psame("abc", 3, 0.7, NULL), psame("full", 3, 0.7, 1e4), tolerance = 1e-12)
})

test_that("the observer that simulates decides as the likelihood integrates", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(136)
  n <- 4e4
  nt <- c(1.7, -2.2)
  dist <- c(0.6, 1.8)
  for (version in names(cdimm_formulas)) {
    a <- if (version != "bsc") 0.7
    s <- if (version != "abc") 1.2
    for (set_size in c(1, 3)) {
      for (criterion in c(0, 0.6, -0.5)) {
        for (probe in c(0.15, 1.7, -1)) {
          p_change <- 1 - .imm_cd_psame(probe, 0.1, 6, 3, a, s, 1, criterion, set_size, nt, dist,
                                        version)
          drawn <- .imm_cd_observe(
            rep(probe, n), 0.1, matrix(nt, n, 2, byrow = TRUE), matrix(dist, n, 2, byrow = TRUE),
            set_size, 6, 3, a, s, 1, criterion
          )
          expect_lt(abs(mean(drawn) - p_change), 5 * sqrt(p_change * (1 - p_change) / n) + 1e-3,
                    label = paste(version, set_size, criterion, probe))
        }
      }
    }
  }
  # rimm_cd() draws through the same observer, with mu[1] off 0
  draws <- rimm_cd(n, probe = 2, mu = c(0.4, 2, -1.5), dist = c(0, 0.5, 2), c = 4, a = 0.5, s = 2,
                   kappa = 8, criterion = 0.3)
  p <- dimm_cd(1, probe = 2, mu = c(0.4, 2, -1.5), dist = c(0, 0.5, 2), c = 4, a = 0.5, s = 2,
               kappa = 8, criterion = 0.3)
  expect_lt(abs(mean(draws) - p), 5 * sqrt(p * (1 - p) / n))
})

test_that("check_data() refuses a missing set size and checks the distances as imm() does", {
  data <- imm_cd_data()
  data$set_size <- as.numeric(as.character(data$set_size))
  missing_ss <- data
  missing_ss$set_size[3] <- NA
  expect_error(check_data(cdimm_model(), missing_ss, cdimm_formulas$full), "contains missing values")
  negative <- data
  negative$dist_nt1[negative$set_size > 1][1] <- -1
  expect_error(check_data(cdimm_model(), negative, cdimm_formulas$full), "postive")
  few <- imm_cd("change", "probe", "target", nt_features = paste0("col_nt", 1:7),
                nt_distances = paste0("dist_nt", 1:6), set_size = "set_size")
  expect_error(check_data(few, data, cdimm_formulas$full), "max\\(set_size\\)-1")
  # abc reads no distances
  data$dist_nt1 <- NULL
  expect_s3_class(check_data(cdimm_model("abc"), data, cdimm_formulas$abc), "data.frame")
})

test_that("the set-size-1 levels of s, a and c are flagged as undetermined", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  data <- imm_cd_data()
  expect_warning(
    bmm(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 0 + set_size), data, cdimm_model(),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "'s'"
  )
  expect_warning(
    bmm(bmf(kappa ~ 1, c ~ 0 + set_size, a ~ 0 + set_size), data, cdimm_model("abc"),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "'a' and 'c'"
  )
})

test_that("an estimated criterion is announced with the bound", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  expect_message(
    bmm(bmf(kappa ~ 1, c ~ 1, a ~ 1, s ~ 1, criterion ~ 1), imm_cd_data(), cdimm_model(),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "criterion is estimated"
  )
})

test_that("the generated program passes the set size as an integer, the probe, non-targets and distances as reals", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  data <- imm_cd_data()
  code <- stancode(cdimm_formulas$full, data, cdimm_model())
  expect_match(code, "array[N] int vint1", fixed = TRUE)
  expect_match(
    code,
    "imm_cd_full_lpmf(Y[n] | mu[n], kappa[n], c[n], a[n], s[n], b[n], criterion[n], vint1[n], vreal1[n], vreal2[n]",
    fixed = TRUE
  )
  # features first, then distances, as for imm(); rich reaches the core as
  # the literal 1
  expect_match(
    code,
    paste0("ss, probe, to_vector({nt1, nt2, nt3, nt4, nt5, nt6, nt7}), ",
           "to_vector({dist1, dist2, dist3, dist4, dist5, dist6, dist7}), 0, circmix_logk"),
    fixed = TRUE
  )
  code <- stancode(cdimm_formulas$abc, data, cdimm_model("abc", knowledge = "rich"))
  expect_match(code, "ss, probe, to_vector({nt1, nt2, nt3, nt4, nt5, nt6, nt7}), 1, circmix_logk",
               fixed = TRUE)
})

test_that("log_lik(), posterior_predict() and posterior_epred() evaluate the same P('same')", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(552)
  n <- 4e4
  nt <- c(1.7, -2.2)
  dist <- c(0.6, 1.8)
  draws <- list(mu = 0.1, kappa = c(4, 7), c = c(3, 1.5), a = c(0.6, 0.3), s = c(1.2, 0.4), b = 1,
                criterion = c(0, 0.3))
  for (version in names(cdimm_formulas)) {
    pars <- c("mu", "kappa", .imm_version_spec(version)$weight_parameters, "criterion")
    dpars <- lapply(draws[pars], function(x) matrix(rep_len(x, n), n, 1))
    covariates <- c(list(vreal1 = 0.8, vreal2 = nt[1], vreal3 = nt[2]),
                    if (version != "abc") list(vreal4 = dist[1], vreal5 = dist[2]))
    prep <- structure(
      list(dpars = dpars, data = c(list(Y = 1L, vint1 = 3L), covariates),
           ndraws = n, nobs = 1, family = list(knowledge = "limited")),
      class = "brmsprep"
    )
    p_change <- 1 - .imm_cd_psame(
      0.8, 0.1, c(4, 7), c(3, 1.5), if (version != "bsc") c(0.6, 0.3),
      if (version != "abc") c(1.2, 0.4), 1, c(0, 0.3), 3, nt, if (version != "abc") dist, version
    )
    expect_equal(.imm_cd_log_lik(version)(1, prep)[1:2], log(p_change), tolerance = 1e-12,
                 label = version)
    expect_equal(as.vector(.imm_cd_posterior_epred(version)(prep))[1:2], p_change,
                 tolerance = 1e-12)
    drawn <- .imm_cd_posterior_predict(version)(1, prep)
    expect_lt(abs(mean(drawn) - mean(p_change)), 0.015)
  }
})

test_that("simulate_data() draws the observer's decision on the design's probes", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(649)
  # two items, the probe on the non-target: P("change") by integrate() over
  # the unbiased arc, with c = 2, a = 0.5, s = 1 at distance 0.7 and b = 1
  n <- 6e4
  design <- data.frame(target = 0.5, probe = 2.5, nt1 = 2, dist1 = 0.7, set_size = 2)[rep(1, n), ]
  model <- imm_cd("change", "probe", "target", nt_features = "nt1", nt_distances = "dist1",
                  set_size = "set_size")
  data <- simulate_data(model, list(kappa = log(6), c = log(2), a = log(0.5), s = 0), design)
  activation <- c(2.5, 2 * exp(-0.7) + 0.5)
  weights <- activation / (sum(activation) + 1)
  hw <- acos(log(besselI(6, 0)) / 6)
  p_same <- weights[1] * vm_mass(2 - hw, 2 + hw, 0, 6) + weights[2] * vm_mass(2 - hw, 2 + hw, 2, 6) +
    (1 - sum(weights)) * hw / pi
  expect_lt(abs(mean(data$change) - (1 - p_same)), 5 * sqrt(0.25 / n))
  expect_error(simulate_data(model, list(kappa = 1, c = 0, a = 0, s = 0), n_trials = 5), "probe")
})

test_that("the generated likelihood passes each parameter to its place in the core, and Stan matches R", {
  skip_unless_stan_tests()
  skip_on_cran()
  data <- imm_cd_data()
  data <- data[as.numeric(as.character(data$set_size)) <= 3, ]
  model_for <- function(version, knowledge) {
    imm_cd("change", "probe", "target", nt_features = c("col_nt1", "col_nt2"),
           nt_distances = c("dist_nt1", "dist_nt2"), set_size = "set_size", version = version,
           knowledge = knowledge)
  }
  # every version at one item and three, the full model under both knowledge
  # settings (one precision state, so they agree). Every value differs, so a
  # swapped slot would show
  cases <- list(
    full = list(version = "full", ss = 3),
    full_rich = list(version = "full", ss = 3, knowledge = "rich"),
    full_ss1 = list(version = "full", ss = 1),
    bsc = list(version = "bsc", ss = 3),
    abc = list(version = "abc", ss = 3)
  )
  values <- c(mu = 0.1, kappa = 6, c = 2.7, a = 0.45, s = 1.3, b = 1.15, criterion = 0.35,
              probe = 0.9, nt1 = 1.7, nt2 = -2.2, dist1 = 0.6, dist2 = 1.8)
  wrappers <- character(0)
  calls <- character(0)
  r_side <- numeric(0)
  for (name in names(cases)) {
    case <- cases[[name]]
    knowledge <- case$knowledge %||% "limited"
    model <- model_for(case$version, knowledge)
    formula <- cdimm_formulas[[case$version]]
    family <- paste0("imm_cd_", case$version)
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
    r_side[name] <- .cd_bernoulli_ld(0, .imm_cd_psame(
      values[["probe"]], values[["mu"]], values[["kappa"]], values[["c"]],
      if (case$version != "bsc") values[["a"]], if (case$version != "abc") values[["s"]],
      values[["b"]], values[["criterion"]], case$ss, values[c("nt1", "nt2")],
      values[c("dist1", "dist2")], case$version
    ))
  }
  chunks <- c("circmix_functions.stan", "cd_functions.stan", "circmix_cd_functions.stan",
              "imm_functions.stan", "imm_cd_functions.stan")
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
  stan_side <- stats::setNames(as.numeric(fit$draws("out")), names(cases))
  expect_equal(stan_side, r_side, tolerance = 1e-10)
  expect_identical(stan_side[["full_rich"]], stan_side[["full"]])
  # the four other cases differ from one another
  expect_identical(length(unique(round(stan_side[names(stan_side) != "full_rich"], 6))), 4L)
})

test_that("the threaded program parses", {
  skip_unless_stan_tests()
  skip_on_cran()
  data <- imm_cd_data()
  for (version in names(cdimm_formulas)) {
    code <- suppressWarnings(stancode(cdimm_formulas[[version]], data, cdimm_model(version),
                                      threads = brms::threading(2)))
    expect_match(code, "partial_log_lik", fixed = TRUE)
    stan_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), compile = FALSE)
    expect_true(stan_model$check_syntax(quiet = TRUE), label = version)
  }
})
