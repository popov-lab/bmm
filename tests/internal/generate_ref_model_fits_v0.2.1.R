###############################################################################!
# GENERATE REFERENCE MODEL FITS WITH FIXED VERSION AND SEEDS FOR TESTING  ####
###############################################################################!

# this function will be built into the next version of the package
install_and_load_bmm_version <- function(version) {
  if ("package:bmm" %in% search()) {
    detach("package:bmm", unload=TRUE)
  }
  path <- paste0(.libPaths()[1], "/bmm-", version)
  if (!dir.exists(path) || length(list.files(path)) == 0 || length(list.files(paste0(path, "/bmm"))) == 0) {
    dir.create(path)
    remotes::install_github(paste0("venpopov/bmm@",version), lib=path)
  }
  library(bmm, lib.loc=path)
}

# Load data for circular models (Participant 1-10, SetSize 1-4, from oberauer_lin_2017)
# TODO: generalize in the future to allow for different datasets for different models
ref_data <- function() {
  withr::local_package('dplyr')
  dat <- oberauer_lin_2017
  dat <- dat |>
    mutate(ID = as.factor(ID),
           SetSize = as.factor(SetSize)) |>
    filter(ID %in% c(1,2,3,4,5,6,7,8,9,10),
           SetSize %in% c(1,2,3,4)) |>
    arrange(ID, SetSize)
  dat
}

run_sdmSimple <- function(...) {
  cat(paste0("\n##########################################################################\n",
             "# Running sdmSimple\n",
             "##########################################################################\n\n"))
  date <- format(Sys.time(), "%Y%m%d")
  hash <- paste0(sample(c(letters[1:6],0:9), 10, replace = TRUE), collapse="")
  file <- paste0(date,"_", bmm_version, "_sdmSimple_seed-", seed, "_",hash,".rds")
  dir  <- here::here("tests/internal/ref_fits")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- paste0(dir, "/", file)
  dat <- ref_data()
  formula <- brms::bf(dev_rad ~ 1,
                      c ~ 0 + SetSize,
                      kappa ~ 1)
  model <- sdmSimple()
  fit <- fit_model(formula, dat, model,
                   cores = parallel::detectCores(),
                   refresh = 500,
                   init = 1,
                   silent = 1,
                   backend = 'cmdstanr',
                   seed = seed)
  saveRDS(fit, path)
}

run_mixture3p <- function(...) {
  cat(paste0("\n##########################################################################\n",
             "# Running mixture3p \n",
             "##########################################################################\n\n"))
  date <- format(Sys.time(), "%Y%m%d")
  hash <- paste0(sample(c(letters[1:6],0:9), 10, replace = TRUE), collapse="")
  file <- paste0(date,"_", bmm_version, "_mixture3p_seed-", seed, "_",hash,".rds")
  dir  <- here::here("tests/internal/ref_fits")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- paste0(dir, "/", file)
  dat <- ref_data()
  formula <- brms::bf(dev_rad ~ 1,
                      thetat ~ 0 + SetSize,
                      thetant ~ 0 + SetSize,
                      kappa ~ 0 + SetSize)
  model <- mixture3p(non_targets = paste0("Item", 2:4,"_Col_rad"),
                     set_size = "SetSize")
  fit <- fit_model(formula, dat, model,
                   cores = parallel::detectCores(),
                   iter = 1000,
                   refresh = 500,
                   init = 1,
                   silent = 1,
                   backend = 'cmdstanr',
                   seed = seed)
  saveRDS(fit, path)
}

run_mixture2p <- function(...) {
  cat(paste0("\n##########################################################################\n",
             "# Running mixture2p \n",
             "##########################################################################\n\n"))
  date <- format(Sys.time(), "%Y%m%d")
  hash <- paste0(sample(c(letters[1:6],0:9), 10, replace = TRUE), collapse="")
  file <- paste0(date,"_", bmm_version, "_mixture2p_seed-", seed, "_",hash,".rds")
  dir  <- here::here("tests/internal/ref_fits")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- paste0(dir, "/", file)
  dat <- ref_data()
  formula <- brms::bf(dev_rad ~ 1,
                      thetat ~ 0 + SetSize,
                      kappa ~ 0 + SetSize)
  model <- mixture2p()
  fit <- fit_model(formula, dat, model,
                   cores = parallel::detectCores(),
                   iter = 1000,
                   refresh = 500,
                   init = 1,
                   silent = 1,
                   backend = 'cmdstanr',
                   seed = seed)
  saveRDS(fit, path)
}

run_IMMabc <- function(...) {
  cat(paste0("\n##########################################################################\n",
             "# Running IMMabc \n",
             "##########################################################################\n\n"))
  date <- format(Sys.time(), "%Y%m%d")
  hash <- paste0(sample(c(letters[1:6],0:9), 10, replace = TRUE), collapse="")
  file <- paste0(date,"_", bmm_version, "_IMMabc_seed-", seed, "_",hash,".rds")
  dir  <- here::here("tests/internal/ref_fits")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- paste0(dir, "/", file)
  dat <- ref_data()
  formula <- brms::bf(dev_rad ~ 1,
                      a ~ 0 + SetSize,
                      c ~ 0 + SetSize,
                      kappa ~ 0 + SetSize)
  model <- IMMabc(non_targets = paste0("Item", 2:4,"_Col_rad"),
                  set_size = "SetSize")
  fit <- fit_model(formula, dat, model,
                   cores = parallel::detectCores(),
                   iter = 1000,
                   refresh = 500,
                   init = 1,
                   silent = 1,
                   backend = 'cmdstanr',
                   seed = seed)
  saveRDS(fit, path)
}


run_IMMbsc <- function(...) {
  cat(paste0("\n##########################################################################\n",
             "# Running IMMbsc \n",
             "##########################################################################\n\n"))
  date <- format(Sys.time(), "%Y%m%d")
  hash <- paste0(sample(c(letters[1:6],0:9), 10, replace = TRUE), collapse="")
  file <- paste0(date,"_", bmm_version, "_IMMbsc_seed-", seed, "_",hash,".rds")
  dir  <- here::here("tests/internal/ref_fits")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- paste0(dir, "/", file)
  dat <- ref_data()
  formula <- brms::bf(dev_rad ~ 1,
                      c ~ 0 + SetSize,
                      s ~ 0 + SetSize,
                      kappa ~ 0 + SetSize)
  model <- IMMbsc(non_targets = paste0("Item", 2:4,"_Col_rad"),
                  spaPos = paste0("Item", 2:4,"_Pos_rad"),
                  set_size = "SetSize")
  fit <- fit_model(formula, dat, model,
                   cores = parallel::detectCores(),
                   iter = 1000,
                   refresh = 500,
                   init = 1,
                   silent = 1,
                   backend = 'cmdstanr',
                   seed = seed)
  saveRDS(fit, path)
}

run_IMMfull <- function(...) {
  cat(paste0("\n##########################################################################\n",
             "# Running IMMfull \n",
             "##########################################################################\n\n"))
  date <- format(Sys.time(), "%Y%m%d")
  hash <- paste0(sample(c(letters[1:6],0:9), 10, replace = TRUE), collapse="")
  file <- paste0(date,"_", bmm_version, "_IMMfull_seed-", seed, "_",hash,".rds")
  dir  <- here::here("tests/internal/ref_fits")
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  path <- paste0(dir, "/", file)
  dat <- ref_data()
  formula <- brms::bf(dev_rad ~ 1,
                      a ~ 0 + SetSize,
                      c ~ 0 + SetSize,
                      s ~ 0 + SetSize,
                      kappa ~ 0 + SetSize)
  model <- IMMfull(non_targets = paste0("Item", 2:4,"_Col_rad"),
                   spaPos = paste0("Item", 2:4,"_Pos_rad"),
                   set_size = "SetSize")
  fit <- fit_model(formula, dat, model,
                   cores = parallel::detectCores(),
                   chains = 4,
                   iter = 1000,
                   refresh = 500,
                   init = 1,
                   silent = 1,
                   backend = 'cmdstanr',
                   seed = seed)
  saveRDS(fit, path)
}

###############################################################################!
# SETUP                                                                    ####
###############################################################################!

library(brms)
bmm_version <- "v0.2.1"
install_and_load_bmm_version(bmm_version)
seed <- 365


###############################################################################!
# RUN MODELS                                                               ####
###############################################################################!


run_sdmSimple(seed=seed, bmm_version=bmm_version)
run_mixture2p(seed=seed, bmm_version=bmm_version)
run_mixture3p(seed=seed, bmm_version=bmm_version)
run_IMMabc(seed=seed, bmm_version=bmm_version)
run_IMMbsc(seed=seed, bmm_version=bmm_version)
run_IMMfull(seed=seed, bmm_version=bmm_version)
