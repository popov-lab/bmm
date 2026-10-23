# Generates the fixtures of fits saved by bmm 1.4.1, the last release that
# fitted mixture2p, mixture3p and imm as brms mixture families with non-linear
# parameters (#434). test-circmix-migration.R reads them to check that such
# fits keep working after the move to custom families. The fits must come from
# 1.4.1 itself, so the script exports the v1.4.1 tag to a temporary directory
# and loads that tree instead of this one. Run from the repo root, for all of
# them or for the ones named (the real fit needs CmdStan):
#   Rscript tests/local/make-old-circmix-fixtures.R
#   Rscript tests/local/make-old-circmix-fixtures.R mock_v141_imm.rds
# Same rules as tests/local/make-fixtures.R (ADR 0004, R5).

wanted <- commandArgs(trailingOnly = TRUE)
make <- function(name) length(wanted) == 0 || name %in% wanted

assets <- normalizePath("tests/testthat/assets", mustWork = TRUE)
tree <- tempfile("bmm-v1.4.1-")
dir.create(tree)
stopifnot(system(paste("git archive v1.4.1 | tar -x -C", shQuote(tree))) == 0)
pkgload::load_all(tree, quiet = TRUE)
stopifnot(identical(as.character(utils::packageVersion("bmm")), "1.4.1"))
# bmm.silent stays at its default: update() reads it from the stored
# configure_opts, and a silent fit would hide the message the tests expect
options(bmm.sort_data = FALSE)

save_fixture <- function(fit, name, seed) {
  attr(fit, "fixture") <- list(
    seed = seed, bmm = as.character(utils::packageVersion("bmm")),
    brms = as.character(utils::packageVersion("brms")), date = as.character(Sys.Date())
  )
  path <- file.path(assets, name)
  saveRDS(fit, path)
  cat(name, ":", round(file.size(path) / 1024), "KB\n")
}

# mock fits store 1 as their stanfit: an empty stanfit object would add 90 KB
# to each file (the environments of its prototype), so the tests that refit
# put one in themselves
mock_fit <- function(...) {
  bmm(..., backend = "mock", mock_fit = 1, rename = FALSE)
}

# two participants at set sizes 1, 2 and 4: set size 1 is in the data, so 1.4
# added its constant priors for a parameter predicted by the set-size factor
old_data <- oberauer_lin_2017[oberauer_lin_2017$ID %in% 1:2 &
                                oberauer_lin_2017$set_size %in% c(1, 2, 4), ]
old_data$set_size <- droplevels(old_data$set_size)
nt_features <- paste0("col_nt", 1:3)
nt_distances <- paste0("dist_nt", 1:3)

# moved from make-fixtures.R: the parameter-routing tests keep only the
# formula, the data and the bmm slot of a mock fit, so it carries no draws.
# The tracked file predates this script and is kept as it is, so it is
# rebuilt only when named
if ("mock_bmmfit_mixture2p.rds" %in% wanted) {
  set.seed(7)
  mix_data <- data.frame(y = rmixture2p(20, kappa = 5, p_mem = 0.7),
                         x = rep(c("a", "b"), 10))
  fit <- mock_fit(bmf(kappa ~ 1 + x, thetat ~ 1), mix_data, mixture2p("y"))
  save_fixture(structure(unclass(fit)[c("formula", "data", "bmm")], class = class(fit)),
               "mock_bmmfit_mixture2p.rds", seed = 7)
}

# mock fits carry everything update() and prior_info() read, but no draws.
# mixture3p: a free mu1 with its 1.4 prior (mu1 was a dpar), a class-b prior on
# a parameter with an intercept and a slope (kappa; 1.4 refused set_size next
# to an intercept, so the slope is session) and on one without an intercept
# (thetat), and the set-size-1 constants on the fixed and the random effect of
# thetant
if (make("mock_v141_mixture3p.rds")) {
  fit <- mock_fit(
    bmf(mu1 ~ 1, kappa ~ 1 + session, thetat ~ 0 + set_size,
        thetant ~ 0 + set_size + (0 + set_size || ID)),
    old_data,
    mixture3p("dev_rad", nt_features = nt_features, set_size = "set_size"),
    prior = brms::set_prior("student_t(1, 0, 0.3)", class = "Intercept", dpar = "mu1") +
      brms::set_prior("normal(2, 0.5)", class = "b", nlpar = "kappa") +
      brms::set_prior("normal(0, 1)", class = "b", nlpar = "thetat")
  )
  save_fixture(fit, "mock_v141_mixture3p.rds", seed = NA)
}

# imm: the set-size-1 constants on the fixed and the random effect of a, a
# class-b prior on an intercept-only parameter (c) and an sd prior (a)
if (make("mock_v141_imm.rds")) {
  fit <- mock_fit(
    bmf(kappa ~ 1, c ~ 1, a ~ 0 + set_size + (0 + set_size || ID), s ~ 1),
    old_data,
    imm("dev_rad", nt_features = nt_features, nt_distances = nt_distances,
        set_size = "set_size", version = "full"),
    prior = brms::set_prior("normal(0, 0.7)", class = "b", nlpar = "c") +
      brms::set_prior("exponential(2)", class = "sd", nlpar = "a")
  )
  save_fixture(fit, "mock_v141_imm.rds", seed = NA)
}

# a sampled fit, for emmeans() and native_parameters(), with the three ways
# 1.4 users wrote priors on a non-linear parameter. Built after the mock fits:
# the link functions of a brms family keep the frame they were made in, and
# after a sampled fit in the same session that frame doubles the size of the
# mixture3p mock (282 KB instead of 132 KB)
if (make("bmmfit_v141_mixture2p.rds")) {
  fit <- bmm(
    bmf(kappa ~ 1 + set_size + (1 | ID), thetat ~ 1), old_data, mixture2p("dev_rad"),
    prior = brms::set_prior("normal(1, 0.5)", class = "b", coef = "Intercept", nlpar = "kappa") +
      brms::set_prior("normal(0, 0.3)", class = "b", nlpar = "kappa") +
      brms::set_prior("exponential(3)", class = "sd", nlpar = "kappa"),
    backend = "cmdstanr", chains = 1, iter = 400, warmup = 200, refresh = 0, seed = 21
  )
  save_fixture(fit, "bmmfit_v141_mixture2p.rds", seed = 21)
}
