# Generates every fixture in tests/testthat/assets/ (ADR 0004, R5: one chain,
# at most 300 KB and 200 draws). Run from the repo root, for all of them or
# for the ones named:
#   Rscript tests/local/make-fixtures.R
#   Rscript tests/local/make-fixtures.R bmmfit_m3_ppcheck.rds
# Fixtures are git-tracked in tests/testthat/assets/ but excluded from the
# built package via .Rbuildignore; tests skip when they are absent. Tests load
# them with load_fixture_fit(), and a model registers its stored fit as a
# `fixture` case in its provider (tests/testthat/helper-model-<model>.R).

devtools::load_all(quiet = TRUE)

wanted <- commandArgs(trailingOnly = TRUE)
make <- function(name) length(wanted) == 0 || name %in% wanted

fit_args <- list(
  backend = "cmdstanr", chains = 1, iter = 450, warmup = 300,
  refresh = 0, silent = 2, seed = 123
)

# the seed and the versions that made a fixture travel with it
save_fixture <- function(fit, name, seed) {
  attr(fit, "fixture") <- list(
    seed = seed, bmm = as.character(utils::packageVersion("bmm")),
    brms = as.character(utils::packageVersion("brms")), date = as.character(Sys.Date())
  )
  path <- file.path("tests/testthat/assets", name)
  saveRDS(fit, path)
  cat(name, ":", round(file.size(path) / 1024), "KB\n")
}

# ddm: two conditions so the group tests have a predictor column
if (make("bmmfit_ddm_ppcheck.rds")) {
  set.seed(1)
  ddm_data <- rbind(
    cbind(rddm(120, drift = 1.5, bound = 1.4, ndt = 0.25), cond = "easy"),
    cbind(rddm(120, drift = 0.5, bound = 1.4, ndt = 0.25), cond = "hard")
  )
  ddm_fit <- do.call(bmm, c(list(
    formula = bmf(drift ~ 0 + cond),
    data = ddm_data,
    model = ddm(rt = "rt", response = "response")
  ), fit_args))
  save_fixture(ddm_fit, "bmmfit_ddm_ppcheck.rds", seed = 1)
}

# cswald: two conditions so the group tests have a predictor column, at an
# error rate the 'simple' version is meant for (~10%) while still leaving
# enough lower-boundary responses for the response check
if (make("bmmfit_cswald_ppcheck.rds")) {
  set.seed(2)
  cswald_data <- rbind(
    cbind(rcswald(75, drift = 3.0, bound = 0.8, ndt = 0.2), cond = "easy"),
    cbind(rcswald(75, drift = 2.0, bound = 0.8, ndt = 0.2), cond = "hard")
  )
  cswald_fit <- do.call(bmm, c(list(
    formula = bmf(drift ~ 0 + cond),
    data = cswald_data,
    model = cswald(rt = "rt", response = "response", version = "simple")
  ), fit_args))
  save_fixture(cswald_fit, "bmmfit_cswald_ppcheck.rds", seed = 2)
}

if (make("bmmfit_ezdm3_ppcheck.rds")) {
  set.seed(3)
  ezdm3_data <- rezdm(10, n_trials = 40, drift = 1.5, bound = 1.2, ndt = 0.3,
                      version = "3par")
  ezdm3_fit <- do.call(bmm, c(list(
    formula = bmf(drift ~ 1),
    data = ezdm3_data,
    model = ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
                 n_trials = "n_trials", version = "3par")
  ), fit_args))
  save_fixture(ezdm3_fit, "bmmfit_ezdm3_ppcheck.rds", seed = 3)
}

# ezdm 4par: moderate drift and small cells so both boundaries carry
# statistics in the data, while simulated replicates occasionally produce
# undefined statistics (the NA path)
if (make("bmmfit_ezdm4_ppcheck.rds")) {
  set.seed(4)
  repeat {
    ezdm4_data <- rezdm(10, n_trials = 20, drift = 0.8, bound = 1.2, ndt = 0.3,
                        zr = 0.5, version = "4par")
    if (!anyNA(ezdm4_data)) break
  }
  ezdm4_fit <- do.call(bmm, c(list(
    formula = bmf(drift ~ 1),
    data = ezdm4_data,
    model = ezdm(mean_rt = c("mean_rt_upper", "mean_rt_lower"),
                 var_rt = c("var_rt_upper", "var_rt_lower"),
                 n_upper = "n_upper", n_trials = "n_trials", version = "4par")
  ), fit_args))
  save_fixture(ezdm4_fit, "bmmfit_ezdm4_ppcheck.rds", seed = 4)
}

# m3 custom with random effects, for the m3 pp_check() and the
# post-processing tests; ten participants keep the random-effect draws small
if (make("bmmfit_m3_ppcheck.rds")) {
  m3_fit <- bmm(
    bmf(
      corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b,
      c ~ 1 + cond + (1 + cond || ID), a ~ 1 + cond + (1 + cond || ID),
      d ~ 1 + (1 || ID)
    ),
    data = oberauer_lewandowsky_2019_e1[oberauer_lewandowsky_2019_e1$ID %in% 1:10, ],
    model = m3(
      resp_cats = c("corr", "other", "dist", "npl"),
      num_options = c(n_opt_corr = "n_corr", n_opt_other = "n_other",
                      n_opt_dist = "n_dist", n_opt_npl = "n_npl"),
      choice_rule = "simple", version = "custom",
      links = list(c = "log", a = "log", d = "log")
    ),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200,
    refresh = 0, silent = 2, seed = 5
  )
  save_fixture(m3_fit, "bmmfit_m3_ppcheck.rds", seed = 5)
}

# mpt: the two-high-threshold model with person effects on D and g, for the
# ppp() T2 test, which needs brms's own draw arrays, the tree column rebuilt
# from the tree indicators and participants read from (1 | id)
if (make("bmmfit_mpt_ppcheck.rds")) {
  set.seed(8)
  n_id <- 10
  D <- stats::plogis(stats::rnorm(n_id, 0.5, 0.7))
  g <- stats::plogis(stats::rnorm(n_id, 0, 0.5))
  mpt_data <- data.frame(
    id = factor(rep(seq_len(n_id), 2)), item_type = rep(c("old", "new"), each = n_id)
  )
  p_old <- ifelse(mpt_data$item_type == "old", D + (1 - D) * g, (1 - D) * g)
  mpt_data$old <- stats::rbinom(2 * n_id, 40, p_old)
  mpt_data$new <- 40L - mpt_data$old
  mpt_fit <- do.call(bmm, c(list(
    formula = bmf(D ~ 1 + (1 | id), g ~ 1 + (1 | id)),
    data = mpt_data,
    model = mpt(list(
      mpt_tree("old", old ~ D + (1 - D) * g, new ~ (1 - D) * (1 - g)),
      mpt_tree("new", old ~ (1 - D) * g, new ~ D + (1 - D) * (1 - g))
    ), tree_id = "item_type")
  ), fit_args))
  save_fixture(mpt_fit, "bmmfit_mpt_ppcheck.rds", seed = 8)
}

# sdm on rstan, so the tests can read what an rstan fit stores: no warmup
# saved, draws from the prior; two participants at set sizes 1 to 4
if (make("bmmfit_example1.rds")) {
  sdm_data <- oberauer_lin_2017[oberauer_lin_2017$ID %in% 1:2 &
                                  oberauer_lin_2017$set_size %in% 1:4, ]
  sdm_data$ID <- factor(sdm_data$ID)
  sdm_data$set_size <- factor(sdm_data$set_size)
  sdm_fit <- bmm(
    bmf(c ~ 0 + set_size, kappa ~ 1), data = sdm_data, model = sdm("dev_rad"),
    backend = "rstan", chains = 1, iter = 400, warmup = 200, init = 1,
    save_warmup = FALSE, sample_prior = "yes", refresh = 0, silent = 2, seed = 6
  )
  # the compiled model and the sampler instance are 90% of the file; the
  # C++ code stays, because stan_version() reads stanc's version from it
  compiled <- sdm_fit$fit@stanmodel
  sdm_fit$fit@stanmodel <- methods::new("stanmodel",
    model_name = compiled@model_name, model_code = compiled@model_code,
    model_cpp = compiled@model_cpp
  )
  sdm_fit$fit@.MISC <- new.env()
  save_fixture(sdm_fit, "bmmfit_example1.rds", seed = 6)
}

# mixture3p: the weights are softmaxed jointly, and set size 1 has no
# non-target, so its thetant is switched off; two participants at set sizes
# 1, 2 and 4
if (make("bmmfit_mixture3p.rds")) {
  m3p_data <- oberauer_lin_2017[oberauer_lin_2017$ID %in% 1:2 &
                                  oberauer_lin_2017$set_size %in% c(1, 2, 4), ]
  m3p_data$set_size <- factor(m3p_data$set_size)
  m3p_fit <- do.call(bmm, c(list(
    formula = bmf(kappa ~ 1, thetat ~ 0 + set_size, thetant ~ 0 + set_size),
    data = m3p_data,
    model = mixture3p("dev_rad", nt_features = paste0("col_nt", 1:3), set_size = "set_size")
  ), fit_args))
  save_fixture(m3p_fit, "bmmfit_mixture3p.rds", seed = 123)
}

# The two mock fixtures keep what the parameter routing reads (the formula,
# the data and the bmm slot of a mock fit), so they carry no draws
mock_fixture <- function(...) {
  fit <- bmm(..., backend = "mock", mock_fit = 1, rename = FALSE)
  structure(unclass(fit)[c("formula", "data", "bmm")], class = class(fit))
}

if (make("mock_bmmfit_mixture2p.rds")) {
  set.seed(7)
  mix_data <- data.frame(y = rmixture2p(20, kappa = 5, p_mem = 0.7),
                         x = rep(c("a", "b"), 10))
  save_fixture(
    mock_fixture(bmf(kappa ~ 1 + x, thetat ~ 1), mix_data, mixture2p("y")),
    "mock_bmmfit_mixture2p.rds", seed = 7
  )
}

if (make("mock_bmmfit_m3.rds")) {
  save_fixture(
    mock_fixture(
      bmf(c ~ 1 + cond, a ~ 1),
      oberauer_lewandowsky_2019_e1[oberauer_lewandowsky_2019_e1$ID %in% 1:2, ],
      m3(resp_cats = c("corr", "other", "npl"), num_options = c("n_corr", "n_other", "n_npl"),
         choice_rule = "softmax", version = "ss")
    ),
    "mock_bmmfit_m3.rds", seed = NA
  )
}

# A fit saved by bmm 0.2.1 (February 2024), for restructure() of old fits.
# It is cut from a reference fit that lived in tests/internal until #422,
# read from the git history so that this block still runs after their
# deletion: one chain of 200 draws, without the compiled model and the
# sampler instance. Its custom family's environment held the whole frame of
# the 0.2.1 fitting call; restructure() reads only the model and the
# formula from it, so only those stay, with the formulas detached from the
# frames they were written in
if (make("bmmfit_legacy_sdm_v0.2.1.rds")) {
  source_commit <- "6870f5d69ba94430192e6947cda6fb1b5d0fd744"
  source_path <- "tests/internal/ref_fits/20240215_v0.2.1_sdmSimple_seed-365_d406b23125.rds"
  blob <- tempfile(fileext = ".rds")
  stopifnot(system2("git", c("show", paste0(source_commit, ":", source_path)), stdout = blob) == 0)
  legacy <- readRDS(blob)

  sim <- legacy$fit@sim
  keep <- seq_len(200)
  # rstan 2.21 stores a chain as a data frame whose attributes carry the
  # sampler parameters; subsetting the rows drops them, so they are put back
  chain <- sim$samples[[1]]
  extra <- attributes(chain)[setdiff(names(attributes(chain)), c("names", "row.names", "class"))]
  extra$sampler_params <- lapply(extra$sampler_params, `[`, keep)
  cut <- chain[keep, , drop = FALSE]
  attributes(cut)[names(extra)] <- extra
  sim$samples <- list(cut)
  sim$chains <- 1L
  sim$n_save <- length(keep)
  sim$warmup2 <- 0L
  sim$iter <- sim$warmup + length(keep)
  sim$permutation <- list(keep)
  legacy$fit@sim <- sim
  legacy$fit@stan_args <- legacy$fit@stan_args[1]
  legacy$fit@inits <- legacy$fit@inits[1]
  legacy$fit@stanmodel <- methods::new("stanmodel")
  legacy$fit@.MISC <- new.env()

  detach_formulas <- function(x) {
    rapply(x, function(f) {
      environment(f) <- globalenv()
      f
    }, classes = "formula", how = "replace")
  }
  old_env <- legacy$family$env
  env <- new.env(parent = globalenv())
  env$model <- old_env$model
  env$formula <- detach_formulas(unclass(old_env$formula))
  class(env$formula) <- class(old_env$formula)
  env$formula$family$env <- env
  legacy$family$env <- env
  legacy$formula <- detach_formulas(unclass(legacy$formula))
  class(legacy$formula) <- class(old_env$formula)
  legacy$formula$family$env <- env
  attr(attr(legacy$data, "terms"), ".Environment") <- globalenv()
  save_fixture(legacy, "bmmfit_legacy_sdm_v0.2.1.rds", seed = 365)
}
