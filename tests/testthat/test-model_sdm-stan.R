# R <-> Stan checks of the sdm() likelihood (#615) ------------------------------

# The normalising constant Z = integral of exp(eta(u)) over the circle, with
# eta(u) = c sqrt(kappa / (2 pi)) exp(kappa (cos(u) - 1)). integrate() gets
# breakpoints at multiples of the two widths of the spike at u = 0: where eta
# falls to about 1, and the SD of exp(eta) near the peak, which is far smaller
# once the peak is large. The integrand is exp(eta - peak), with eta - peak
# written without cancellation, so that exp() neither overflows at large c and
# kappa nor loses digits near the peak. lc is log c.
sdm_log_z_reference <- function(lc, kappa) {
  log_peak <- lc + 0.5 * log(kappa / (2 * pi))
  width <- min(pi, sqrt(2 * max(log_peak, 1) / kappa))
  sd_peak <- min(pi, 1 / sqrt(exp(log_peak) * kappa))
  breaks <- sort(unique(pmin(pi, c(
    0, width * c(0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16),
    sd_peak * c(0.5, 1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64), pi
  ))))
  integrand <- function(u) exp(exp(log_peak) * expm1(-2 * kappa * sin(u / 2)^2))
  half <- vapply(seq_len(length(breaks) - 1), function(i) {
    stats::integrate(integrand, breaks[i], breaks[i + 1],
      rel.tol = 1e-13, abs.tol = 0, subdivisions = 5000L
    )$value
  }, numeric(1))
  exp(log_peak) + log(2 * sum(half))
}

# the grid of #615; the default priors put warm-up near its upper corner
sdm_issue_grid <- expand.grid(
  lc = seq(-1, 5, by = 0.5),
  kappa = c(0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300)
)

# where develop's error was 0.24 to 7.5
sdm_corner_cells <- data.frame(
  lc = c(4.5, 4.5, 4.5, 5, 5, 5),
  kappa = c(100, 200, 300, 100, 200, 300)
)

# beyond the grid: where a refit's warm-up went along the ridge of constant
# c kappa^1.5, and cutting only at multiples of the eta width was off by 8e-5
# to 1.4e-2; small c with large kappa, where cutting only at multiples of
# the SD near the peak was off by 7e-6; and large c with large kappa, where
# eta - peak formed as peak (exp(kappa (cos(u) - 1)) - 1) was off by 3e-7
sdm_extended_cells <- data.frame(
  lc = c(6, 7, 8, 10, 10, 12, -1, 12),
  kappa = c(300, 30, 5, 1, 10, 0.5, 3000, 3000)
)

# a pair on either side of the switch between the two rules at kappa (peak +
# 1) = 4000, at a large and at a small peak; and log peak 1.3 with kappa 2000,
# where the last cut at 3 instead of 5 eta widths was off by 3e-8
sdm_tight_cells <- data.frame(
  lc = c(2.25, 2.3, -1.5, -1.4, -1.6),
  kappa = c(100, 100, 1000, 1000, 2000)
)
sdm_stan_cells <- rbind(sdm_issue_grid, sdm_extended_cells, sdm_tight_cells)

# dsdm() integrates without scaling by the peak and stops with a non-finite
# value once log Z passes about 709 (log c 5 with kappa 200), so the
# likelihood cells stop at kappa 100; (5, 100) is a corner cell
sdm_parity_cells <- data.frame(lc = c(1, 3, 5), kappa = c(2, 30, 100))
sdm_parity_y <- c(-pi, -2, -1, -0.3, -0.05, -0.01, 0, 0.004, 0.02, 0.1, 0.6, 2.5)

# Runs the installed tdata and likelihood chunks inside generated quantities,
# once per cell, with the names they read from a brms program (N, the run
# metadata, c on the log scale, kappa) declared locally and `target +=` turned
# into `lp +=`. The model's other data stanvars come from configure_model(),
# so the program passes the chunks what a fit would. One program serves both
# tests: compiling takes most of the time of the first one.
sdm_stan_program <- function(data_scode) {
  chunk <- function(file) read_lines2(system.file("stan_chunks", file, package = "bmm"))
  funs <- gsub("#include[^\n]*", "", chunk("sdm_simple_funs.stan"))
  tdata <- chunk("sdm_simple_tdata.stan")
  likelihood <- gsub("target +=", "lp +=", chunk("sdm_simple_likelihood.stan"), fixed = TRUE)
  paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n",
    "  int<lower=1> K; vector[K] grid_lc; vector[K] grid_kappa;\n",
    "  int<lower=1> P; int<lower=1> n_y; vector[n_y] y_obs;\n",
    "  vector[P] par_lc; vector[P] par_kappa;\n",
    data_scode, "\n}\n",
    "generated quantities {\n",
    "  vector[K] log_z;\n",
    "  vector[P + 1] log_lik;\n",
    "  for (k in 1:K) {\n",
    "    int N = 1;\n",
    "    vector[N] c = rep_vector(grid_lc[k], N);\n",
    "    vector[N] kappa = rep_vector(grid_kappa[k], N);\n",
    "    int G_sdm_runs = 1;\n",
    "    array[1] int sdm_run_start = {1};\n",
    "    array[1] int sdm_run_count = {1};\n",
    "    real lp = 0;\n",
    tdata, "\n", likelihood, "\n",
    "    log_z[k] = -lp;\n",
    "  }\n",
    # p <= P: one run of n_y rows at parity cell p; p = P + 1: all P cells as
    # P runs of one design, which checks that each run gets its own constant
    "  for (p in 1:(P + 1)) {\n",
    "    int G_sdm_runs = p <= P ? 1 : P;\n",
    "    int N = G_sdm_runs * n_y;\n",
    "    vector[N] Y;\n",
    "    vector[N] mu = rep_vector(0, N);\n",
    "    vector[N] c;\n",
    "    vector[N] kappa;\n",
    "    array[G_sdm_runs] int sdm_run_start;\n",
    "    array[G_sdm_runs] int sdm_run_count;\n",
    "    real lp = 0;\n",
    "    for (g in 1:G_sdm_runs) {\n",
    "      int cell = p <= P ? p : g;\n",
    "      int first = (g - 1) * n_y + 1;\n",
    "      sdm_run_start[g] = first;\n",
    "      sdm_run_count[g] = n_y;\n",
    "      Y[first:(first + n_y - 1)] = y_obs;\n",
    "      c[first:(first + n_y - 1)] = rep_vector(par_lc[cell], n_y);\n",
    "      kappa[first:(first + n_y - 1)] = rep_vector(par_kappa[cell], n_y);\n",
    "    }\n",
    tdata, "\n", likelihood, "\n",
    "    lp += sdm_simple_lpdf(Y | mu, c, kappa);\n",
    "    log_lik[p] = lp;\n",
    "  }\n",
    "}\n"
  )
}

sdm_stan_values <- function() {
  formula <- bmf(c ~ 1, kappa ~ 1)
  model <- sdm(resp_error = "y")
  dat <- data.frame(y = sdm_parity_y)
  stanvars <- configure_model(model, check_data(model, dat, formula), formula)$stanvars
  run_metadata <- c("G_sdm_runs", "sdm_run_start", "sdm_run_count")
  model_data <- Filter(function(sv) sv$block == "data" && !sv$name %in% run_metadata, stanvars)
  data <- c(
    list(
      K = nrow(sdm_stan_cells), grid_lc = sdm_stan_cells$lc, grid_kappa = sdm_stan_cells$kappa,
      P = nrow(sdm_parity_cells), n_y = length(sdm_parity_y), y_obs = sdm_parity_y,
      par_lc = sdm_parity_cells$lc, par_kappa = sdm_parity_cells$kappa
    ),
    stats::setNames(lapply(model_data, `[[`, "sdata"), vapply(model_data, `[[`, "", "name"))
  )
  program <- sdm_stan_program(paste(vapply(model_data, `[[`, "", "scode"), collapse = "\n"))

  dir <- withr::local_tempdir()
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program, dir = dir), quiet = TRUE)$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0,
    refresh = 0, show_messages = FALSE, sig_figs = 17, output_dir = dir
  )
  # base-R read: cmdstanr's reader misparses long plain-decimal numbers
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  list(
    log_z = as.numeric(csv[1, paste0("log_z.", seq_len(data$K))]),
    log_lik = as.numeric(csv[1, paste0("log_lik.", seq_len(data$P + 1))])
  )
}

sdm_stan_cached <- local({
  values <- NULL
  function() {
    if (is.null(values)) values <<- sdm_stan_values()
    values
  }
})

test_that("the Stan normalising constant of sdm() matches integrate() across the parameter range", {
  skip_unless_stan_tests()
  skip_on_cran()

  # about 10 s: the first call compiles the program both tests read
  stan <- sdm_stan_cached()
  reference <- mapply(sdm_log_z_reference, sdm_stan_cells$lc, sdm_stan_cells$kappa)
  error <- abs(stan$log_z - reference)
  on_grid <- seq_len(nrow(sdm_issue_grid))

  # the rule's measured worst: 1.3e-13 on the grid, one ulp of log Z (4.7e-10
  # at log c 12, kappa 3000) at the cells beyond it and 7.1e-15 at the tight
  # cells, so a rule with fewer nodes in either regime fails here
  expect_lt(max(error[on_grid]), 1e-7)
  expect_lt(max(error[on_grid][sdm_issue_grid$kappa <= 30]), 1e-10)
  for (i in seq_len(nrow(sdm_corner_cells))) {
    cell <- sdm_corner_cells[i, ]
    at <- which(sdm_issue_grid$lc == cell$lc & sdm_issue_grid$kappa == cell$kappa)
    expect_lt(error[at], 1e-7,
      label = glue("|error in log Z| at log c {cell$lc}, kappa {cell$kappa}")
    )
  }
  for (i in seq_len(nrow(sdm_extended_cells))) {
    cell <- sdm_extended_cells[i, ]
    expect_lt(error[nrow(sdm_issue_grid) + i], 1e-7,
      label = glue("|error in log Z| at log c {cell$lc}, kappa {cell$kappa}")
    )
  }
  # within 1e-10 on both sides, so the value jumps by less than 2e-10 per run
  # where a sampler crosses the switch
  for (i in seq_len(nrow(sdm_tight_cells))) {
    cell <- sdm_tight_cells[i, ]
    expect_lt(error[nrow(sdm_issue_grid) + nrow(sdm_extended_cells) + i], 1e-10,
      label = glue("|error in log Z| at log c {cell$lc}, kappa {cell$kappa}")
    )
  }
})

test_that("the Stan log-likelihood of sdm() runs matches dsdm()", {
  skip_unless_stan_tests()
  skip_on_cran()

  stan <- sdm_stan_cached()
  r_log_lik <- vapply(seq_len(nrow(sdm_parity_cells)), function(i) {
    cell <- sdm_parity_cells[i, ]
    sum(dsdm(sdm_parity_y, 0, exp(cell$lc), cell$kappa, log = TRUE))
  }, numeric(1))

  for (i in seq_len(nrow(sdm_parity_cells))) {
    cell <- sdm_parity_cells[i, ]
    expect_lt(abs(stan$log_lik[i] - r_log_lik[i]), 1e-6,
      label = glue("|Stan - dsdm()| log-likelihood at log c {cell$lc}, kappa {cell$kappa}")
    )
  }
  expect_lt(abs(stan$log_lik[nrow(sdm_parity_cells) + 1] - sum(r_log_lik)), 1e-6,
    label = "|Stan - dsdm()| log-likelihood of the three cells as runs of one design"
  )
})

# The cells other than the default (#558) -------------------------------------

# the #615 cells taken as log c of the density scaling and translated to the
# log peak, and very small kappa with log c taken as the log peak. Stan gets
# lc and kappa and translates: cmdstanr writes data with 15 significant
# digits, which moves log Z by 1e-7 at log peak 15
sdm_laplace_cells <- rbind(
  data.frame(
    lc = sdm_stan_cells$lc, kappa = sdm_stan_cells$kappa, translate = 1,
    on_grid = seq_len(nrow(sdm_stan_cells)) <= nrow(sdm_issue_grid)
  ),
  data.frame(expand.grid(lc = c(-1, 0, 2, 5), kappa = c(0.05, 0.1)), translate = 0, on_grid = TRUE)
)
sdm_laplace_cells$log_peak <- with(sdm_laplace_cells, lc + translate * 0.5 * log(kappa / (2 * pi)))

# Like sdm_stan_program(), for the sdm_cont chunks: the Laplace normaliser on
# its cells, sdm_vm_ldenom_peak() against sdm_simple_ldenom() on the #615
# cells, and the likelihood of every cell at the parity cells, one run each
# and all as runs of one design
sdm_cont_stan_program <- function(data_scode) {
  chunk <- function(file) read_lines2(system.file("stan_chunks", file, package = "bmm"))
  funs <- gsub("#include[^\n]*", "", paste(
    chunk("sdm_simple_funs.stan"), chunk("sdm_kernel_funs.stan"), chunk("sdm_cont_funs.stan"),
    sep = "\n"
  ))
  tdata <- chunk("sdm_simple_tdata.stan")
  likelihood <- gsub("target +=", "lp +=", chunk("sdm_cont_likelihood.stan"), fixed = TRUE)
  paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n",
    "  int<lower=1> K; vector[K] grid_lc; vector[K] grid_kappa;\n",
    "  int<lower=1> L; vector[L] lap_lc; vector[L] lap_kappa; vector[L] lap_translate;\n",
    "  int<lower=1> P; int<lower=1> n_y; vector[n_y] y_obs;\n",
    "  vector[P] par_lc; vector[P] par_kappa;\n",
    data_scode, "\n}\n",
    "generated quantities {\n",
    "  vector[L] laplace_log_z;\n",
    "  vector[K] vm_peak_minus_simple;\n",
    "  array[2, 3] vector[P + 1] log_lik;\n",
    "  for (l in 1:L) {\n",
    "    real log_peak = lap_lc[l] + lap_translate[l] * 0.5 * log(lap_kappa[l] / (2 * pi()));\n",
    "    laplace_log_z[l] = log(2 * pi()) - sdm_laplace_ldenom(log_peak, lap_kappa[l], sdm_gl_x, sdm_gl_w);\n",
    "  }\n",
    "  for (k in 1:K) {\n",
    "    vm_peak_minus_simple[k] =\n",
    "      sdm_vm_ldenom_peak(grid_lc[k] + 0.5 * log(grid_kappa[k] / (2 * pi())), grid_kappa[k],\n",
    "                         sdm_gl_x, sdm_gl_w, sdm_wide_cosm1, sdm_wide_w)\n",
    "      - sdm_simple_ldenom(grid_lc[k], grid_kappa[k], sdm_gl_x, sdm_gl_w, sdm_wide_cosm1, sdm_wide_w);\n",
    "  }\n",
    "  for (sdm_activation in 1:2) for (sdm_scaling in 1:3) for (p in 1:(P + 1)) {\n",
    "    int G_sdm_runs = p <= P ? 1 : P;\n",
    "    int N = G_sdm_runs * n_y;\n",
    "    vector[N] Y;\n",
    "    vector[N] mu = rep_vector(0, N);\n",
    "    vector[N] c;\n",
    "    vector[N] kappa;\n",
    "    array[G_sdm_runs] int sdm_run_start;\n",
    "    array[G_sdm_runs] int sdm_run_count;\n",
    "    real lp = 0;\n",
    "    for (g in 1:G_sdm_runs) {\n",
    "      int cell = p <= P ? p : g;\n",
    "      int first = (g - 1) * n_y + 1;\n",
    "      sdm_run_start[g] = first;\n",
    "      sdm_run_count[g] = n_y;\n",
    "      Y[first:(first + n_y - 1)] = y_obs;\n",
    "      c[first:(first + n_y - 1)] = rep_vector(par_lc[cell], n_y);\n",
    "      kappa[first:(first + n_y - 1)] = rep_vector(par_kappa[cell], n_y);\n",
    "    }\n",
    tdata, "\n", likelihood, "\n",
    "    lp += sdm_cont_lpdf(Y | mu, c, kappa, sdm_activation, sdm_scaling);\n",
    "    log_lik[sdm_activation, sdm_scaling, p] = lp;\n",
    "  }\n",
    "}\n"
  )
}

sdm_cont_stan_values <- function() {
  formula <- bmf(c ~ 1, kappa ~ 1)
  model <- suppressMessages(sdm(resp_error = "y", activation = "laplace"))
  dat <- data.frame(y = sdm_parity_y)
  stanvars <- configure_model(model, check_data(model, dat, formula), formula)$stanvars
  not_needed <- c("G_sdm_runs", "sdm_run_start", "sdm_run_count", "sdm_activation", "sdm_scaling")
  model_data <- Filter(function(sv) sv$block == "data" && !sv$name %in% not_needed, stanvars)
  data <- c(
    list(
      K = nrow(sdm_stan_cells), grid_lc = sdm_stan_cells$lc, grid_kappa = sdm_stan_cells$kappa,
      L = nrow(sdm_laplace_cells), lap_lc = sdm_laplace_cells$lc,
      lap_kappa = sdm_laplace_cells$kappa, lap_translate = sdm_laplace_cells$translate,
      P = nrow(sdm_parity_cells), n_y = length(sdm_parity_y), y_obs = sdm_parity_y,
      par_lc = sdm_parity_cells$lc, par_kappa = sdm_parity_cells$kappa
    ),
    stats::setNames(lapply(model_data, `[[`, "sdata"), vapply(model_data, `[[`, "", "name"))
  )
  program <- sdm_cont_stan_program(paste(vapply(model_data, `[[`, "", "scode"), collapse = "\n"))

  dir <- withr::local_tempdir()
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program, dir = dir), quiet = TRUE)$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0,
    refresh = 0, show_messages = FALSE, sig_figs = 17, output_dir = dir
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  read <- function(names) as.numeric(csv[1, names])
  list(
    laplace_log_z = read(paste0("laplace_log_z.", seq_len(data$L))),
    vm_peak_minus_simple = read(paste0("vm_peak_minus_simple.", seq_len(data$K))),
    log_lik = function(a, s) read(paste0("log_lik.", a, ".", s, ".", seq_len(data$P + 1)))
  )
}

sdm_cont_stan_cached <- local({
  values <- NULL
  function() {
    if (is.null(values)) values <<- sdm_cont_stan_values()
    values
  }
})

test_that("the Stan normaliser of the Laplace activation matches integrate()", {
  skip_unless_stan_tests()
  skip_on_cran()

  # about 15 s: the first call compiles the program the three tests read
  stan <- sdm_cont_stan_cached()
  reference <- mapply(sdm_log_z_integrate, sdm_laplace_cells$log_peak, sdm_laplace_cells$kappa,
    MoreArgs = list(activation = "laplace")
  )
  error <- abs(stan$laplace_log_z - reference)
  # the rule's measured worst is 1.7e-12 against the exponential-integral
  # series; beyond the grid an ulp of the log peak is up to 1e-8 of log Z
  expect_lt(max(error), 1e-7)
  expect_lt(max(error[sdm_laplace_cells$on_grid & sdm_laplace_cells$kappa <= 30]), 1e-10)
})

test_that("sdm_vm_ldenom_peak() is sdm_simple_ldenom() taking the log peak", {
  skip_unless_stan_tests()
  skip_on_cran()

  stan <- sdm_cont_stan_cached()
  expect_lt(max(abs(stan$vm_peak_minus_simple)), 1e-14)
})

test_that("the Stan log-likelihood of every sdm() cell matches dsdm()", {
  skip_unless_stan_tests()
  skip_on_cran()

  stan <- sdm_cont_stan_cached()
  for (a in 1:2) for (s in 1:3) {
    activation <- c("vonmises", "laplace")[a]
    scaling <- c("density", "peak", "range")[s]
    r_log_lik <- vapply(seq_len(nrow(sdm_parity_cells)), function(i) {
      cell <- sdm_parity_cells[i, ]
      sum(dsdm(sdm_parity_y, 0, exp(cell$lc), cell$kappa,
        log = TRUE, activation = activation, scaling = scaling
      ))
    }, numeric(1))
    expect_lt(max(abs(stan$log_lik(a, s) - c(r_log_lik, sum(r_log_lik)))), 1e-6,
      label = glue("|Stan - dsdm()| log-likelihood of {activation}, {scaling}")
    )
  }
})

# Discrete response options (#558, K2) -----------------------------------------

# Runs the installed sdm_grid chunks inside generated quantities for every
# activation x scaling cell: the serial likelihood chunk, and the threaded one
# inside a function whose signature is built from the stanvars' pll_args, as
# brms builds partial_log_lik(), called once per slice. The data stanvars and
# the cell tables come from configure_model(); the parameters are data, one
# value per row, constant within each run of `run`
sdm_grid_stan_program <- function(stanvars) {
  chunk <- function(file) read_lines2(system.file("stan_chunks", file, package = "bmm"))
  funs <- gsub("#include[^\n]*", "", paste(chunk("sdm_kernel_funs.stan"), chunk("sdm_grid_funs.stan"),
                                            sep = "\n"))
  data_vars <- Filter(function(sv) sv$block == "data", stanvars)
  pll <- vapply(data_vars, `[[`, "", "pll_args")
  threaded <- gsub("target +=", "lp +=", chunk("sdm_grid_likelihood_threaded.stan"), fixed = TRUE)
  serial <- gsub("target +=", "lp +=", chunk("sdm_grid_likelihood.stan"), fixed = TRUE)
  paste0(
    "functions {\n", funs, "\n",
    "  real sdm_slice_loglik(int start, int end, data array[] int Y, vector mu_all, vector c_all,\n",
    "                    vector kappa_all, ", paste(pll, collapse = ", "), ") {\n",
    "    int N = end - start + 1;\n",
    "    vector[N] mu = mu_all[start:end];\n",
    "    vector[N] c = c_all[start:end];\n",
    "    vector[N] kappa = kappa_all[start:end];\n",
    "    real lp = sdm_grid_lpmf(Y[start:end] | mu, c, kappa, sdm_K, sdm_activation, sdm_scaling);\n",
    threaded, "\n",
    "    return lp;\n",
    "  }\n",
    "}\n",
    "data {\n",
    "  int<lower=1> N; array[N] int Y; vector[N] par_lc; vector[N] par_kappa; vector[N] par_mu;\n",
    "  int<lower=1> n_slices; array[n_slices] int slice_start; array[n_slices] int slice_end;\n",
    paste(vapply(Filter(function(sv) !sv$name %in% c("sdm_activation", "sdm_scaling"), data_vars),
                 `[[`, "", "scode"), collapse = "\n"), "\n",
    "}\n",
    "transformed data {\n", chunk("sdm_grid_tdata.stan"), "\n}\n",
    "generated quantities {\n",
    "  array[2, 3] real serial;\n",
    "  array[2, 3] real sliced;\n",
    "  for (sdm_activation in 1:2) for (sdm_scaling in 1:3) {\n",
    "    vector[N] mu = par_mu;\n",
    "    vector[N] c = par_lc;\n",
    "    vector[N] kappa = par_kappa;\n",
    "    real lp = sdm_grid_lpmf(Y | mu, c, kappa, sdm_K, sdm_activation, sdm_scaling);\n",
    serial, "\n",
    "    serial[sdm_activation, sdm_scaling] = lp;\n",
    "    sliced[sdm_activation, sdm_scaling] = 0;\n",
    "    for (k in 1:n_slices) {\n",
    "      sliced[sdm_activation, sdm_scaling] += sdm_slice_loglik(slice_start[k], slice_end[k], Y, mu, c,\n",
    "        kappa, ", paste(vapply(data_vars, `[[`, "", "name"), collapse = ", "), ");\n",
    "    }\n",
    "  }\n",
    "}\n"
  )
}

# both designs give the same program, which is compiled once per session
sdm_grid_stan_compiled <- local({
  programs <- list()
  function(program) {
    key <- rlang::hash(program)
    if (is.null(programs[[key]])) {
      dir <- tempfile("sdm_grid_stan")
      dir.create(dir)
      programs[[key]] <<- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program, dir = dir),
                                                  quiet = TRUE)
    }
    programs[[key]]
  }
})

# One fixed-param run over the design `dat` (columns y, lc, kappa, mu, run and
# the option columns of `model`); returns Stan's serial and sliced
# log-likelihood per cell
sdm_grid_stan_run <- function(model, dat, slice = 7L) {
  formula <- bmf(c ~ 0 + run, kappa ~ 0 + run)
  stanvars <- withr::with_options(list(bmm.sort_data = FALSE), {
    suppressMessages(configure_model(model, check_data(model, dat, formula), formula)$stanvars)
  })
  checked <- suppressMessages(check_data(model, dat, formula))
  data_vars <- Filter(function(sv) sv$block == "data", stanvars)
  starts <- seq(1L, nrow(dat), by = slice)
  data <- c(
    list(
      N = nrow(dat), Y = checked$sdm_bin, par_lc = dat$lc, par_kappa = dat$kappa, par_mu = dat$mu,
      n_slices = length(starts), slice_start = starts, slice_end = pmin(starts + slice - 1L, nrow(dat))
    ),
    stats::setNames(lapply(data_vars, `[[`, "sdata"), vapply(data_vars, `[[`, "", "name"))
  )
  data[c("sdm_activation", "sdm_scaling")] <- NULL
  dir <- withr::local_tempdir()
  fit <- sdm_grid_stan_compiled(sdm_grid_stan_program(stanvars))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1, iter_warmup = 0,
    refresh = 0, show_messages = FALSE, sig_figs = 17, output_dir = dir
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#", check.names = FALSE)
  read <- function(name) matrix(as.numeric(csv[1, paste0(name, ".", rep(1:2, 3), ".", rep(1:3, each = 2))]), 2, 3)
  list(serial = read("serial"), sliced = read("sliced"))
}

# dsdm() row by row: the discrete likelihood written without the cell tables
sdm_grid_r_log_lik <- function(dat, options_of) {
  out <- matrix(NA_real_, 2, 3)
  for (a in 1:2) for (s in 1:3) {
    out[a, s] <- sum(vapply(seq_len(nrow(dat)), function(i) {
      dsdm(dat$y[i], dat$mu[i], exp(dat$lc[i]), dat$kappa[i], n_options = 360, options = options_of(i),
           log = TRUE, activation = c("vonmises", "laplace")[a], scaling = c("density", "peak", "range")[s])
    }, numeric(1)))
  }
  out
}

# Runs of five rows at mu = 0, then rows with parameters of their own (the
# per-trial regime), some with mu off the grid, and with kappa on both sides
# of the switch to the continuous normaliser at w = 2.25 steps: log c* about
# 0 with kappa 600 (w 2.34 steps) and 700 (2.17 steps), and log c* 3 with
# kappa 30 and 400; and two rows with w about one step (log c* 0 with kappa
# 3300, log c* 3 with kappa 160), where the grid sum and the continuous
# normaliser differ by about 1e-5, so a switch below 2.25 steps fails here
sdm_grid_full_design <- function() {
  withr::local_seed(558)
  runs <- data.frame(lc = c(1, 3, 5, 0), kappa = c(2, 30, 100, 600), mu = 0)
  own <- data.frame(
    lc = c(0, 0, 3, 3, -1, 2, 4, 1.5, 0.5, 2.5, 0, 3),
    kappa = c(600, 700, 30, 400, 0.5, 5, 50, 10, 650, 3, 3300, 160),
    mu = c(0, 0, 0.1, -0.37, 0.25, 0, 1.2, -2.9, 0, 0.004, 0, 0.3)
  )
  pars <- rbind(runs[rep(1:4, each = 5), ], own)
  pars$run <- sprintf("r%02d", seq_len(nrow(pars)))
  pars$run[1:20] <- sprintf("r%02d", rep(1:4, each = 5))
  pars$y <- deg2rad(sample(-179:180, nrow(pars), replace = TRUE))
  pars$y[1:6] <- deg2rad(c(0, 1, -1, 180, 0, 2))
  pars
}

# Rows offering sets of 2 to 9 options at irregular offsets, different on
# almost every row, each with the target and the response; three runs, two of
# them with mu off the grid
sdm_grid_sets_design <- function() {
  withr::local_seed(436)
  n <- 30
  options <- lapply(seq_len(n), function(i) {
    m <- sample(2:9, 1)
    c(0, sample(setdiff(-179:180, 0), m - 1))
  })
  options[[7]] <- options[[3]][c(2, 1, seq_along(options[[3]])[-(1:2)])]
  width <- max(lengths(options))
  dat <- data.frame(
    lc = rep(c(1, 3, 0.5), each = 10), kappa = rep(c(3, 40, 0.8), each = 10),
    mu = rep(c(0, -0.2, 1.1), each = 10), run = rep(c("a", "b", "c"), each = 10)
  )
  dat$y <- deg2rad(vapply(options, function(o) o[sample.int(length(o), 1)], numeric(1)))
  for (j in seq_len(width)) {
    dat[[paste0("o", j)]] <- deg2rad(vapply(options, function(o) o[j], numeric(1)))
  }
  dat
}

test_that("the Stan likelihood of the full wheel matches dsdm(), serial and sliced", {
  skip_unless_stan_tests()
  skip_on_cran()

  # about 15 s: compiles the program both grid tests use
  dat <- sdm_grid_full_design()
  stan <- sdm_grid_stan_run(sdm("y", n_options = 360), dat)
  r <- sdm_grid_r_log_lik(dat, function(i) NULL)
  # the switch to the continuous normaliser errs by at most 5e-12 per row
  expect_lt(max(abs(stan$serial - r)), 1e-8)
  expect_lt(max(abs(stan$sliced - stan$serial)), 1e-10)
})

test_that("rows with a different number and spacing of options each match dsdm(options = )", {
  skip_unless_stan_tests()
  skip_on_cran()

  dat <- sdm_grid_sets_design()
  columns <- grep("^o[0-9]+$", names(dat), value = TRUE)
  model <- sdm("y", n_options = 360, options = columns)
  stan <- sdm_grid_stan_run(model, dat, slice = 4L)
  r <- sdm_grid_r_log_lik(dat, function(i) {
    o <- unlist(dat[i, columns])
    o[!is.na(o)]
  })
  expect_lt(max(abs(stan$serial - r)), 1e-8)
  expect_lt(max(abs(stan$sliced - stan$serial)), 1e-10)
})

# On a coarse wheel the grid sum and the continuous normaliser differ by 2.6e-6
# per row at w = 2.25 steps with K = 4 (measured from K = 90 up, the bound is
# 5e-12), so the switch is limited to K >= 90; rows near w = 2.25 steps, kappa
# about 0.08 with c* about 1 for K = 4
test_that("a coarse wheel keeps the exact grid sum, serial and sliced", {
  skip_unless_stan_tests()
  skip_on_cran()

  dat <- data.frame(
    lc = c(0, 0, -0.5, 0.2, 1, 0, 0, -1),
    kappa = c(0.08, 0.06, 0.1, 0.07, 0.03, 0.08, 0.5, 0.09),
    mu = c(0, 0, 0, 0, 0, 0.3, 0, 0)
  )
  dat$run <- sprintf("r%02d", seq_len(nrow(dat)))
  dat$y <- c(0, 1, 2, -1, 0, 2, 1, 0) * pi / 2
  stan <- sdm_grid_stan_run(sdm("y", n_options = 4), dat, slice = 3L)
  r <- matrix(NA_real_, 2, 3)
  for (a in 1:2) for (s in 1:3) {
    r[a, s] <- sum(vapply(seq_len(nrow(dat)), function(i) {
      dsdm(dat$y[i], dat$mu[i], exp(dat$lc[i]), dat$kappa[i], n_options = 4, log = TRUE,
           activation = c("vonmises", "laplace")[a], scaling = c("density", "peak", "range")[s])
    }, numeric(1)))
  }
  expect_lt(max(abs(stan$serial - r)), 1e-10)
  expect_lt(max(abs(stan$sliced - stan$serial)), 1e-10)
})
