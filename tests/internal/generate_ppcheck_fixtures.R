# Generates the small fitted-model fixtures used by tests/testthat/test-pp_check.R
# for the multi-observable pp_check() tests (#401). Run from the repo root:
#   Rscript tests/internal/generate_ppcheck_fixtures.R
# Fixtures are git-tracked in tests/testthat/assets/ but excluded from the
# built package via .Rbuildignore; tests skip when they are absent.

devtools::load_all(quiet = TRUE)

fit_args <- list(
  backend = "cmdstanr", chains = 1, iter = 450, warmup = 300,
  refresh = 0, silent = 2, seed = 123
)

save_fixture <- function(fit, name) {
  path <- file.path("tests/testthat/assets", name)
  saveRDS(fit, path)
  cat(name, ":", round(file.size(path) / 1024), "KB\n")
}

# ddm: two conditions so the group tests have a predictor column
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
save_fixture(ddm_fit, "bmmfit_ddm_ppcheck.rds")

set.seed(3)
ezdm3_data <- rezdm(10, n_trials = 40, drift = 1.5, bound = 1.2, ndt = 0.3,
                    version = "3par")
ezdm3_fit <- do.call(bmm, c(list(
  formula = bmf(drift ~ 1),
  data = ezdm3_data,
  model = ezdm(mean_rt = "mean_rt", var_rt = "var_rt", n_upper = "n_upper",
               n_trials = "n_trials", version = "3par")
), fit_args))
save_fixture(ezdm3_fit, "bmmfit_ezdm3_ppcheck.rds")

# ezdm 4par: moderate drift and small cells so both boundaries carry
# statistics in the data, while simulated replicates occasionally produce
# undefined statistics (the NA path)
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
save_fixture(ezdm4_fit, "bmmfit_ezdm4_ppcheck.rds")
