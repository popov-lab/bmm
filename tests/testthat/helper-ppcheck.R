# a saved fit for the pp_check() tests; a helper so that a model's own test
# file can load its fixture too
load_ppcheck_fit <- function(name) {
  path <- test_path("assets", name)
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}
