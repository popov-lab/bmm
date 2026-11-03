# Test m-AFC SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_mafc reads a near-miss link target as the parameter, with a warning", {
  expect_warning(
    sdt_mafc("n_correct", "n_trials", m = 4, links = list(dd = "log")),
    "read as 'd'"
  )
})

test_that("sdt_mafc model stores m", {
  model <- sdt_mafc("n_correct", "n_trials", m = 6, dist = "gumbel_min")
  expect_equal(model$other_vars$m, 6L)
})

test_that("sdt_mafc accepts m as a constant or a column name", {
  m_const <- sdt_mafc("n_correct", "n_trials", m = 4)
  expect_equal(m_const$other_vars$m, 4L)

  m_col <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  expect_identical(m_col$other_vars$m, "set_size")
})

test_that("sdt_mafc requires response and m, and reads trial rows without n_trials", {
  expect_error(sdt_mafc("n_correct", "n_trials"), "missing in sdt_mafc\\(\\): m")
  trial <- sdt_mafc("correct", m = 4)
  expect_true(trial$other_vars$trial_level)
  expect_null(trial$other_vars$n_trials)
  expect_match(trial$requirements, "one row per trial")
  expect_false(sdt_mafc("n_correct", "n_trials", m = 4)$other_vars$trial_level)
})

test_that("sdt_mafc rejects m < 2 and invalid distributions", {
  expect_error(sdt_mafc("n_correct", "n_trials", m = 1), "m must be")
  expect_error(sdt_mafc("n_correct", "n_trials", m = 4, dist = "foo"))
})

test_that("sdt_mafc warns and truncates a non-integer m", {
  expect_warning(model <- sdt_mafc("n_correct", "n_trials", m = 2.7), "integer")
  expect_equal(model$other_vars$m, 2L)
})

test_that("sdt_mafc rejects NA and non-finite m with its own message", {
  expect_error(sdt_mafc("n_correct", "n_trials", m = NA), "m must be")
  expect_error(sdt_mafc("n_correct", "n_trials", m = NA_real_), "m must be")
  expect_error(sdt_mafc("n_correct", "n_trials", m = Inf), "m must be")
  expect_error(sdt_mafc("n_correct", "n_trials", m = NaN), "m must be")
})


############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

test_that("sdt_mafc check_data adds m_afc and dist_type columns", {
  model <- sdt_mafc("n_correct", "n_trials", m = 4, dist = "logistic")
  formula <- bmf(d ~ 1)
  dat <- data.frame(n_correct = c(70, 65), n_trials = c(100, 100))

  result <- check_data(model, dat, formula)
  expect_true(all(c("m_afc", "dist_type") %in% colnames(result)))
  expect_equal(unique(result$m_afc), 4L)
  expect_equal(unique(result$dist_type), 4L)
})

test_that("sdt_mafc check_data resolves per-row set size from a data column", {
  model <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  dat <- data.frame(
    set_size = c(2L, 4L, 6L),
    n_correct = c(70, 55, 40),
    n_trials = c(100, 100, 100)
  )

  result <- check_data(model, dat, bmf(d ~ 1))
  expect_equal(result$m_afc, c(2L, 4L, 6L))
})

test_that("sdt_mafc requirements list the set-size column only when m names one", {
  line <- "Number of alternatives (m): one integer >= 2 per row"
  expect_match(sdt_mafc("c", m = "k")$requirements, line, fixed = TRUE)
  expect_match(sdt_mafc("c", "n", m = "k")$requirements, line, fixed = TRUE)
  expect_no_match(sdt_mafc("c", m = 4)$requirements, line, fixed = TRUE)
})

test_that("sdt_mafc check_data errors when the set-size column is missing", {
  model <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  dat <- data.frame(n_correct = c(70, 65), n_trials = c(100, 100))
  expect_error(check_data(model, dat, bmf(d ~ 1)), "Set-size column")
})

test_that("trial-level sdt_mafc check_data() names the responses that are not 0 or 1", {
  model <- sdt_mafc("correct", m = 4)
  expect_error(
    check_data(model, data.frame(correct = c(1, 0, 4, 0.5)), bmf(d ~ 1)),
    "'correct' must hold one of '0', '1' in every row. Found: '4', '0.5'"
  )
  checked <- check_data(model, data.frame(correct = c(1, 0, 1)), bmf(d ~ 1))
  expect_identical(checked$correct, c(1L, 0L, 1L))
  expect_identical(checked$nTrials, c(1L, 1L, 1L))
})

test_that("trial rows and the cells aggregated from them give the same sdt_mafc Stan program", {
  skip_on_cran()
  # 2 participants x 20 trials (R3)
  trials <- data.frame(
    id = rep(1:2, each = 20),
    correct = c(1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1,
                1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0)
  )
  cells <- data.frame(id = 1:2, correct = c(15, 12), n_trials = 20)
  mock <- function(data, model) suppressMessages(
    bmm(bmf(d ~ 1 + (1 | id)), data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
  trial_fit <- mock(trials, sdt_mafc("correct", m = 4))
  cell_fit <- mock(cells, sdt_mafc("correct", "n_trials", m = 4))
  expect_identical(brms::stancode(trial_fit), brms::stancode(cell_fit))
  trial_data <- brms::standata(trial_fit)
  expect_equal(as.vector(trial_data$trials), rep(1, 40))
  expect_equal(as.vector(tapply(trial_data$Y, trials$id, sum)), c(15, 12))
})

# the family's log_lik on the Stan data of trial rows sums to its log_lik on
# the Stan data of the cell less log choose(n, y): 7 correct responses in 10
# trials, choose(10, 7) = 120
test_that("log_lik_sdt_mafc() on trial rows sums to the cell value less log choose(n, y)", {
  skip_on_cran()
  trials <- data.frame(correct = c(1, 1, 0, 1, 1, 0, 1, 1, 0, 1))
  cells <- data.frame(correct = 7, n = 10)
  mock <- function(data, model) suppressMessages(
    bmm(bmf(d ~ 1), data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
  ll <- function(fit) {
    data <- brms::standata(fit)
    prep <- epred_prep(list(d = matrix(1.1, 1, data$N)), data)
    vapply(seq_len(data$N), function(i) log_lik_sdt_mafc(i, prep), numeric(1))
  }
  by_trial <- ll(mock(trials, sdt_mafc("correct", m = 4)))
  by_cell <- ll(mock(cells, sdt_mafc("correct", "n", m = 4)))
  expect_equal(sum(by_trial), by_cell - log(120), tolerance = 1e-10)
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_mafc recycles vectorized parameters per observation", {
  withr::local_seed(445)
  n_correct <- rsdt_mafc(12, 100, m = 4, d = rnorm(12, 1.5, 0.3))
  expect_length(n_correct, 12)

  n2 <- rsdt_mafc(3, 100, m = c(2, 4, 8), d = 1.5)
  expect_length(n2, 3)
})

test_that("rsdt_mafc stays inside the binomial's domain in the far tail", {
  # sampled off the probability scale these two cells are the two ways the
  # gumbel_min Gamma ratio fails: at d' = -36.4 it has cancelled to 7.9e13,
  # which rbinom() answers with NA, and from d' = -40 it returns exactly 1,
  # which is silent and simulates all-correct where P(correct) is near 4e-18.
  withr::local_seed(445)
  expect_equal(rsdt_mafc(5, 100L, m = 2L, d = -36.4, dist = "gumbel_min"),
               rep(0L, 5))
  expect_equal(rsdt_mafc(5, 100L, m = 2L, d = -40, dist = "gumbel_min"),
               rep(0L, 5))
})

test_that("m-AFC probability correct is monotone increasing in d'", {
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    pcs <- vapply(c(0, 0.5, 1, 2, 3), .mafc_pc_r, numeric(1), m = 4L, dist = di)
    expect_true(all(diff(pcs) > 0), info = di)
    expect_true(all(pcs > 0 & pcs < 1), info = di)
  }
})

test_that(".mafc_pc_r vectorized path matches elementwise evaluation", {
  d <- c(0.3, 1.1, 2.2)
  mm <- c(2L, 4L, 6L)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    vec <- .mafc_pc_r(d, mm, di)
    ref <- vapply(seq_along(d),
                  function(i) .mafc_pc_r(d[i], mm[i], di), numeric(1))
    expect_equal(vec, ref, tolerance = 1e-12, info = di)
  }
})

test_that("m-AFC probability correct matches independent numerical oracles", {
  # Independent oracles: latent-space integration / binomial-sum closed form,
  # none of which reuse the implementation in .mafc_pc_r.
  o_norm <- function(d, m) {
    stats::integrate(function(z) dnorm(z) * pnorm(z + d)^(m - 1),
                     -Inf, Inf, rel.tol = 1e-12)$value
  }
  o_logis <- function(d, m) {
    stats::integrate(function(z) dlogis(z) * plogis(z + d)^(m - 1),
                     -Inf, Inf, rel.tol = 1e-12)$value
  }
  # largest extreme value: pdf exp(-z - exp(-z)), cdf exp(-exp(-z))
  o_gmax <- function(d, m) {
    stats::integrate(function(z) exp(-z - exp(-z)) * exp(-exp(-(z + d)))^(m - 1),
                     -Inf, Inf, rel.tol = 1e-12)$value
  }
  # smallest extreme value: binomial-sum closed form of the same integral
  o_gmin <- function(d, m) {
    a <- exp(d)
    k <- 0:(m - 1)
    sum(choose(m - 1, k) * (-1)^k / (a * k + 1))
  }
  for (m in c(3L, 4L, 6L)) {
    for (d in c(0.5, 1, 2)) {
      expect_equal(.mafc_pc_r(d, m, "normal"), o_norm(d, m), tolerance = 1e-5)
      expect_equal(.mafc_pc_r(d, m, "logistic"), o_logis(d, m), tolerance = 1e-7)
      expect_equal(.mafc_pc_r(d, m, "gumbel_max"), o_gmax(d, m), tolerance = 1e-6)
      expect_equal(.mafc_pc_r(d, m, "gumbel_min"), o_gmin(d, m), tolerance = 1e-9)
    }
  }
})


############################################################################# !
# FORMULA AND CONFIGURE_MODEL TESTS                                      ####
############################################################################# !

test_that("sdt_mafc produces valid stancode including its custom functions", {
  dat <- data.frame(n_correct = c(70, 65, 80, 60), n_trials = rep(100, 4),
                    id = c(1, 1, 2, 2))
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  code <- stancode(bmf(d ~ 1), data = dat, model = model)
  expect_true(grepl("sdt_mafc", code))
  expect_true(grepl("mafc_logit_pc", code))
  expect_true(grepl("sdt_quantile", code))
})

test_that("sdt_mafc integrates with the bmm pipeline via mock backend", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(445)
  dat <- data.frame(n_trials = rep(100L, 8))
  dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4, d = 1.2)
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_mafc("n_correct", "n_trials", m = 4, dist = d)
    expect_silent(
      bmm(bmf(d ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_mafc fits mixed set sizes through the pipeline (mock)", {
  withr::local_seed(445)
  dat <- data.frame(set_size = rep(c(2L, 4L, 6L), each = 8), n_trials = 100L)
  dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = dat$set_size,
                             d = 1.2)
  model <- sdt_mafc("n_correct", "n_trials", m = "set_size")
  expect_silent(
    bmm(bmf(d ~ 1), dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

############################################################################# !
# R <-> STAN QUADRATURE AND BRANCH CONTRACT                              ####
############################################################################# !

test_that("the Stan quadrature tables are the R tables, digit for digit", {
  # .mafc_pc_r and mafc_pc() integrate with the same nodes and weights, but the
  # 208 numbers are typed out twice. The Stan parity test below compares with
  # dsdt_mafc() only to 1e-8, so a corrupted digit of a small weight could pass
  # it while shifting every normal-noise and logistic m-AFC likelihood.
  # Text-level, as in test-model_cswald.R, so it costs no compilation.
  sc_path <- system.file("stan_chunks", package = "bmm")
  tdata <- read_lines2(file.path(sc_path, "sdt_mafc_tdata.stan"))
  decls <- regmatches(
    tdata,
    gregexpr("vector\\[[0-9]+\\] [a-z_]+ = to_vector\\(\\{[^}]*\\}\\)", tdata)
  )[[1]]
  stan_dim <- as.integer(sub("vector\\[([0-9]+)\\].*", "\\1", decls))
  stan_tab <- lapply(decls, function(decl) {
    as.numeric(regmatches(decl, gregexpr("-?[0-9][.][0-9]+e[+-][0-9]+", decl))[[1]])
  })
  names(stan_tab) <- names(stan_dim) <- sub(".*\\] ([a-z_]+) = .*", "\\1", decls)

  r_tab <- list(gh_nodes = .mafc_gh_nodes, gh_weights = .mafc_gh_weights,
                gl_nodes = .mafc_gl_nodes, gl_weights = .mafc_gl_weights)
  expected_n <- c(gh_nodes = 40L, gh_weights = 40L, gl_nodes = 64L, gl_weights = 64L)

  expect_named(stan_tab, names(r_tab))
  for (nm in names(r_tab)) {
    expect_length(r_tab[[nm]], expected_n[[nm]])
    expect_length(stan_tab[[nm]], expected_n[[nm]])
    expect_equal(stan_dim[[nm]], expected_n[[nm]], info = nm)
    expect_equal(stan_tab[[nm]], r_tab[[nm]], tolerance = 0, info = nm)
  }

  # mafc_logit_pc() runs once per row, so a table declared in it is rebuilt per
  # row per gradient: ~22% more time per leapfrog step.
  expect_false(grepl("to_vector({", read_lines2(file.path(sc_path, "sdt_mafc_functions.stan")),
                     fixed = TRUE))
})

test_that("the Stan mafc_logit_pc branches match the registry order", {
  # The R side dispatches by position in .sdt_dists; mafc_logit_pc() hardcodes
  # the integers. Swapping two branches in the chunk would leave every test
  # green while turning every dist = "gumbel_max" fit into a gumbel_min model.
  pc <- read_lines2(file.path(system.file("stan_chunks", package = "bmm"),
                              "sdt_mafc_functions.stan"))
  pc <- sub("(?s)\n\\}.*", "",
            sub("(?s).*real mafc_logit_pc\\(", "", pc, perl = TRUE), perl = TRUE)
  pc <- gsub("[[:space:]]+", " ", gsub("//[^\n]*", "", pc))

  expect_equal(names(.sdt_dists),
               c("normal", "gumbel_min", "gumbel_max", "logistic"))

  closed_form <- c(
    normal     = paste("{ if (m == 2) return sdt_log_cumprob(d / sqrt(2.0), dist_type)",
                       "- sdt_log_one_minus_cumprob(d / sqrt(2.0), dist_type);"),
    gumbel_min = paste("{ real e = exp(-d); real log_pc = 0;",
                       "for (k in 1:(m - 1)) log_pc -= log1p(e / k);",
                       "return log_pc - log(-expm1(log_pc)); }"),
    gumbel_max = "return d - log(m - 1);"
  )
  for (nm in names(closed_form)) {
    expect_match(pc,
                 paste0("dist_type == ", which(names(.sdt_dists) == nm), ") ",
                        closed_form[[nm]]),
                 fixed = TRUE, info = nm)
  }

  # logistic is the fallthrough, so it must not appear as a guard at all
  guards <- as.integer(gsub("\\D", "", regmatches(pc, gregexpr("dist_type == [0-9]+", pc))[[1]]))
  expect_equal(sort(guards), c(1L, 2L, 3L))
  expect_false(which(names(.sdt_dists) == "logistic") %in% guards)

  # both quadrature branches read log P(correct) off the log-sum-exp and its
  # complement off the sum of -expm1 terms, from one sweep of the same nodes
  quadrature <- c(
    gh = paste("real log_cdf = (m - 1) * sdt_log_cumprob(gh_nodes[i] + d, dist_type);",
               "log_terms[i] = log(gh_weights[i]) + log_cdf;",
               "q += gh_weights[i] * (-expm1(log_cdf));"),
    gl = paste("real log_cdf = (m - 1) * sdt_log_cumprob(sdt_quantile(gl_nodes[i],",
               "dist_type) + d, dist_type); log_terms[i] = log(gl_weights[i]) + log_cdf;",
               "q += gl_weights[i] * (-expm1(log_cdf));")
  )
  for (nm in names(quadrature)) {
    expect_match(pc, quadrature[[nm]], fixed = TRUE, info = nm)
  }
  expect_equal(length(gregexpr("log_sum_exp(log_terms) - log(fmin(q, 1))", pc,
                               fixed = TRUE)[[1]]), 2L)
})


test_that("the quadrature tables reach sdt_mafc_lpmf in its declared order", {
  skip_on_cran() # four stancode() calls, over 0.4 s in CRAN mode (ADR 0004, R10)
  # The four tables are positional and each pair shares a Stan type, so a
  # transposed pair compiles and runs. Transposing gl_nodes and gl_weights is
  # silent: every logistic cell stays finite and moves, at m = 16 and d' = 0
  # from P(correct) = 0.0625 to 1.0e-23. (Transposing the Gauss-Hermite pair is
  # loud instead -- log() of a negative node makes every normal cell with
  # m > 2 non-finite.) Nothing else in the suite reads the call site, so a
  # swapped `vars` in configure_model.sdt_mafc() would pass everything.
  tables <- function(x) regmatches(x, gregexpr("g[hl]_(nodes|weights)", x))[[1]]
  chunk <- read_lines2(file.path(system.file("stan_chunks", package = "bmm"),
                                 "sdt_mafc_functions.stan"))
  declared <- tables(sub("(?s)\\).*", "",
                         sub("(?s).*real sdt_mafc_lpmf\\(", "", chunk, perl = TRUE),
                         perl = TRUE))
  expect_equal(declared, c("gh_nodes", "gh_weights", "gl_nodes", "gl_weights"))

  dat <- data.frame(n_correct = c(70, 65, 80, 60), n_trials = rep(100, 4),
                    id = c(1, 1, 2, 2))
  model <- sdt_mafc("n_correct", "n_trials", m = 4)
  # the likelihood takes all rows at once, so m and the trials are sliced
  # with the response exactly when brms slices it
  rows <- list(serial = "(Y | mu, d, vint1, vint2[1], trials, ",
               threaded = "(Y[start:end] | mu, d, vint1[start:end], vint2[1], trials[start:end], ",
               force = "(Y | mu, d, vint1, vint2[1], trials, ")
  threads <- list(serial = NULL, threaded = brms::threading(2),
                  force = brms::threading(2, force = TRUE))
  for (nm in names(threads)) {
    code <- stancode(bmf(d ~ 1), data = dat, model = model, threads = threads[[nm]])
    call <- regmatches(code, regexpr("sdt_mafc_lpmf\\(Y[^;]*?\\);", code))
    expect_equal(tables(call), declared, info = nm)
    expect_match(call, rows[[nm]], fixed = TRUE, info = nm)
  }
  # reduce_sum gives the partial log-likelihood its own scope, so the tables
  # have to arrive there in the same order too
  code <- stancode(bmf(d ~ 1), data = dat, model = model,
                   threads = brms::threading(2))
  for (site in c("real partial_log_lik_lpmf\\([^)]*\\)",
                 "reduce_sum\\(partial_log_lik_lpmf[^;]*\\)")) {
    expect_equal(tables(regmatches(code, regexpr(site, code))), declared,
                 info = site)
  }
})

test_that("the Stan sdt_mafc likelihood matches dsdt_mafc row by row and over blocks", {
  skip_unless_stan_tests()
  skip_on_cran()
  # The likelihood takes all rows at once, with m varying by row, and drops the
  # rows whose logit is infinite before one binomial_logit_lpmf() call. One-row
  # slices are compared with dsdt_mafc(); whole blocks with the sum of their
  # rows, which a row read at another row's index would break. A row with an
  # infinite logit adds 0 when its count is the one that logit forces, and makes
  # the whole block -Inf otherwise
  sc <- system.file("stan_chunks", package = "bmm")
  tables <- "gh_nodes, gh_weights, gl_nodes, gl_weights"
  program <- paste0(
    "functions {\n", read_lines2(file.path(sc, "sdt_dist_functions.stan")), "\n",
    read_lines2(file.path(sc, "sdt_mafc_functions.stan")), "\n}\n",
    "data {\n  int N; int K; array[N] int y; array[N] int trials; array[N] int m;\n",
    "  array[N] int dist; vector[N] d; array[K] int first; array[K] int last;\n",
    "  int B; array[B] int bad;\n}\n",
    "transformed data {\n", read_lines2(file.path(sc, "sdt_mafc_tdata.stan")), "\n}\n",
    "generated quantities {\n  vector[N] row; vector[K] blk; vector[B] blk_bad;\n",
    "  for (i in 1:N) row[i] = sdt_mafc_lpmf(y[i:i] | rep_vector(0, 1), d[i:i], m[i:i], ",
    "dist[i], trials[i:i], ", tables, ");\n",
    "  for (k in 1:K) {\n    int s = first[k];\n    int e = last[k];\n",
    "    blk[k] = sdt_mafc_lpmf(y[s:e] | rep_vector(0, e - s + 1), d[s:e], m[s:e], dist[s], ",
    "trials[s:e], ", tables, ");\n  }\n",
    # each block again, with one row that contradicts its infinite logit added
    "  for (b in 1:B) {\n    int k = dist[bad[b]];\n    int s = first[k];\n    int e = last[k];\n",
    "    blk_bad[b] = sdt_mafc_lpmf(append_array(y[s:e], {y[bad[b]]}) | rep_vector(0, e - s + 2), ",
    "append_row(d[s:e], d[bad[b]]), append_array(m[s:e], {m[bad[b]]}), dist[s], ",
    "append_array(trials[s:e], {trials[bad[b]]}), ", tables, ");\n  }\n}\n"
  )
  grid <- expand.grid(m = c(2L, 3L, 4L, 8L), d = c(-1, 0.3, 1.5, 3), cnt = 1:5,
                      dist = seq_along(bmm:::.sdt_dists))
  counts <- rbind(c(20L, 13L), c(20L, 0L), c(20L, 20L), c(1L, 0L), c(1L, 1L))
  grid$trials <- counts[grid$cnt, 1]
  grid$y <- counts[grid$cnt, 2]
  # infinite logits: +Inf where every quadrature term rounds to certainty
  # (normal and logistic at d' = 800), -Inf where exp(-d') overflows
  # (gumbel_min at -800). `forced` rows have the count that logit forces,
  # `contradicting` rows do not
  inf_d <- c(800, 800, -800)
  inf_dist <- c(1L, 4L, 2L)
  forced <- data.frame(m = 4L, d = inf_d, dist = inf_dist, trials = 20L, y = c(20L, 20L, 0L))
  contradicting <- data.frame(m = 4L, d = inf_d, dist = inf_dist, trials = 20L, y = c(19L, 19L, 1L))
  good <- rbind(grid[names(forced)], forced)
  good <- good[order(good$dist), ]
  grid <- rbind(good, contradicting)
  ok <- seq_len(nrow(good))
  N <- nrow(grid)
  K <- length(bmm:::.sdt_dists)
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(list(N = N, K = K, first = match(seq_len(K), good$dist),
                  last = nrow(good) + 1L - match(seq_len(K), rev(good$dist)),
                  B = nrow(contradicting), bad = nrow(good) + seq_len(nrow(contradicting))),
             as.list(grid[c("y", "trials", "m", "dist", "d")])),
    fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  stan <- function(name, n) as.numeric(csv[1, paste0(name, ".", seq_len(n))])
  row <- stan("row", N)
  dist_names <- names(bmm:::.sdt_dists)[grid$dist]
  r_row <- vapply(seq_len(N), function(i) with(grid[i, ], dsdt_mafc(
    y, trials, m, d, dist = dist_names[i], log = TRUE)), numeric(1))
  finite <- abs(grid$d) < 100
  expect_equal(row[finite], r_row[finite], tolerance = 1e-8)
  expect_identical(row[!finite], rep(c(0, -Inf), each = nrow(forced)))
  expect_equal(stan("blk", K), as.vector(tapply(row[ok], good$dist, sum)), tolerance = 1e-12)
  expect_identical(stan("blk_bad", nrow(contradicting)), rep(-Inf, nrow(contradicting)))
})


test_that("the m-AFC logit is chance at d' = 0 and analytic for gumbel_max", {
  # 1e-7 is set by one cell, normal at m = 8, where the 40-point Gauss-Hermite
  # rule sits 5.2e-8 from -log(m - 1) *absolute* (2.7e-8 relative, which is the
  # measure the sibling probability test quotes); the other fifteen are all
  # below 6e-14. The four m = 2 cells are weaker than they look: -log(m - 1) is
  # exactly 0 there, so expect_equal() drops to an absolute comparison and they
  # assert only |x| < 1e-7, where three of the four are 0 by construction. The
  # agreement test below is immune because expect_lt() on an explicit
  # max(abs(...)) names the measure instead of inferring it from the target.
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    for (m in c(2L, 3L, 4L, 8L)) {
      expect_equal(.mafc_logit_pc_r(0, m, di), -log(m - 1), tolerance = 1e-7,
                   info = paste(di, "m =", m))
    }
  }
  # taking the max of gumbel_max variates gives a softmax, whose logit is
  # d' - log(m - 1) with nothing left to cancel: it is exact at every d',
  # including where P(correct) itself has rounded to 1
  for (m in c(2L, 3L, 4L, 8L)) {
    for (d in c(0, 1, 5, 40, 200, 700)) {
      expect_equal(.mafc_logit_pc_r(d, m, "gumbel_max"), d - log(m - 1),
                   tolerance = 0, info = paste("m =", m, "d =", d))
    }
  }
})

test_that(".mafc_logit_pc_r agrees with the probability scale while that scale holds", {
  # the logit and .mafc_pc_r() are separate code paths, so this is what ties
  # them together over the range where the probability scale is still exact.
  # The grid reaches d' = -10 on purpose: reading log P(correct) off the
  # complement alone costs the P(correct) -> 0 tail, and that route misses the
  # normal cell at d' = -10 by 9.7 rather than by the 4.4e-11 of the worst cell
  # here (gumbel_min, m = 3, d' = -10). Outside [-10, 4] it is .mafc_pc_r()
  # that loses precision first, so the comparison stops being informative.
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    for (m in c(2L, 3L, 4L, 8L, 16L)) {
      d <- c(-10, -5, -3, -1, 0, 0.5, 1, 2, 3, 4)
      expect_lt(max(abs(.mafc_logit_pc_r(d, m, di) -
                        stats::qlogis(.mafc_pc_r(d, m, di)))),
                1e-9, label = paste(di, "m =", m))
    }
  }
})

test_that("the m-AFC density keeps responding to d' where P(correct) rounds to 1", {
  # binomial_lpmf on a natural-scale P(correct) stops seeing d' as soon as the
  # complement drops below the double epsilon. Measured on the probability
  # scale at y = 90 of 100 trials and m = 4, stepping d' by 0.25: the normal
  # branch froze at -310.4951 from d' = 12, and gumbel_min, gumbel_max and
  # logistic stopped decreasing at d' = 34.5, 36.25 and 37.75. A log link on
  # d puts that region within reach -- 6.9% of the default normal(1, 1) prior
  # mass sits above d' = 12 under a log link, against 1.9e-28 under identity.
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    ld <- vapply(c(12, 16, 20, 25, 30, 37, 40), function(dd)
      dsdt_mafc(90, 100, m = 4, d = dd, dist = di, log = TRUE), numeric(1))
    expect_true(all(is.finite(ld)), info = di)
    expect_true(all(diff(ld) < 0), info = di)
  }
  # the three branches with no quadrature limit left reach arbitrarily far
  for (di in c("gumbel_max", "gumbel_min", "logistic")) {
    ld <- vapply(c(60, 100, 200, 700), function(dd)
      dsdt_mafc(90, 100, m = 4, d = dd, dist = di, log = TRUE), numeric(1))
    expect_true(all(is.finite(ld)), info = di)
    expect_true(all(diff(ld) < 0), info = di)
  }
  # a cell with no errors, or nothing but errors, is a probability and not NaN
  expect_equal(dsdt_mafc(100, 100, m = 4, d = 200), 1, tolerance = 1e-12)
  expect_equal(dsdt_mafc(0, 100, m = 4, d = -200), 1, tolerance = 1e-12)
})

test_that("the gumbel_min logit outlives its Gamma ratio", {
  # Gamma(1 + e) Gamma(m) / Gamma(m + e) is a difference of lgammas that has
  # cancelled to exactly zero by e = exp(-36), so .mafc_pc_r() returns 1 there
  # and the logit taken off it is +Inf. Telescoping the ratio to
  # prod(k / (k + e)) moves that wall out past d' = 700. Without this the logit
  # scale buys gumbel_min nothing: both forms die at d' = 34.5.
  expect_equal(.mafc_pc_r(36, 4L, "gumbel_min"), 1, tolerance = 0)
  for (d in c(36, 40, 100, 700)) {
    expect_true(is.finite(.mafc_logit_pc_r(d, 4L, "gumbel_min")), info = d)
  }
  # the telescoped product is the same number where lgamma can still form it
  for (d in c(-5, 0, 1, 5, 10)) {
    expect_equal(.mafc_logit_pc_r(d, 4L, "gumbel_min"),
                 stats::qlogis(exp(lgamma(1 + exp(-d)) + lgamma(4) -
                                   lgamma(4 + exp(-d)))),
                 tolerance = 1e-10, info = paste("d =", d))
  }
})

test_that("the gumbel_min logit is -Inf, not NaN, once exp(-d') overflows", {
  # e = exp(-d') is Inf below d' = -709.78, so a cell zeroed by multiplication
  # would be 0 * Inf, and the NaN would poison its column's sum. Only a column
  # whose m is below max(m) carries such cells, which is why a single m never
  # showed it. The unmasked Inf has to survive: -Inf is the right answer.
  expect_equal(.mafc_logit_pc_r(-800, 4L, "gumbel_min"), -Inf)
  expect_equal(.mafc_logit_pc_r(c(-800, -800), c(2L, 8L), "gumbel_min"),
               c(-Inf, -Inf))
  # and zeroing by index must not shift which k a column keeps
  d <- c(0.3, 1.1, 2.2)
  mm <- c(2L, 4L, 6L)
  expect_equal(.mafc_logit_pc_r(d, mm, "gumbel_min"),
               vapply(seq_along(d), function(i) {
                 .mafc_logit_pc_r(d[i], mm[i], "gumbel_min")
               }, numeric(1)),
               tolerance = 1e-12)
})

# posterior_epred (#475) ------------------------------------------------------

test_that("posterior_epred_sdt_mafc() is the mean count posterior_predict_sdt_mafc() simulates", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  withr::local_seed(475)
  sets <- data.frame(d = c(1.5, 0.5, 2, 1))
  for (dist in .sdt_dist_names) {
    data <- list(vint1 = c(2, 4, 3, 6), vint2 = rep(.sdt_dist_id(dist), 4),
                 trials = c(50, 80, 20, 100))
    res <- epred_vs_predict(sets, posterior_epred_sdt_mafc, posterior_predict_sdt_mafc,
                            data = data)
    expect_lt(max(abs(res[, "rel_error"])), 0.02)
  }
})

test_that("posterior_epred_sdt_mafc() returns one column per observation", {
  dpars <- list(d = matrix(c(1.5, 0.5, 2, 1, 0.7, 1.2), 2))
  data <- list(vint1 = c(2, 4, 6), vint2 = rep(.sdt_dist_id("normal"), 3),
               trials = c(50, 80, 20))
  expect_epred_by_cell(posterior_epred_sdt_mafc, dpars, data)
})

test_that("an sdt_mafc fit stores posterior_epred_sdt_mafc() in its family", {
  withr::local_seed(445)
  dat <- data.frame(n_trials = rep(100L, 8))
  dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4, d = 1.2)
  fit <- bmm(bmf(d ~ 1), dat, sdt_mafc("n_correct", "n_trials", m = 4),
             backend = "mock", mock_fit = 1, rename = FALSE)
  expect_identical(fit$formula$family$posterior_epred, posterior_epred_sdt_mafc)
})

test_that("log_lik() of an sdt_mafc fit reads each observation's counts, m and distribution", {
  # two draws of d for two observations; the second observation has its own
  # number of alternatives and the gumbel_max distribution (id 3)
  prep <- epred_prep(list(d = matrix(c(0.8, 1.1, 1.6, 2.2), 2, 2)),
                     list(Y = c(30L, 55L), trials = c(40L, 60L), vint1 = c(2L, 4L),
                          vint2 = c(1L, 3L)))
  expect_equal(log_lik_sdt_mafc(1, prep),
               dsdt_mafc(30, 40, 2, c(0.8, 1.1), dist = "normal", log = TRUE))
  expect_equal(log_lik_sdt_mafc(2, prep),
               dsdt_mafc(55, 60, 4, c(1.6, 2.2), dist = "gumbel_max", log = TRUE))
})
