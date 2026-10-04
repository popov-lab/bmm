# a saved fit for the pp_check() tests; a helper so that a model's own test
# file can load its fixture too
load_ppcheck_fit <- function(name) {
  path <- test_path("assets", name)
  skip_if_not(file.exists(path), "fixture not available (excluded by .Rbuildignore)")
  readRDS(path)
}

# a brmsprep with the given dpars (scalars, or vectors filled draw-major into
# ndraws x nobs matrices) and data, for pp_simulate() without a fit
fake_prep <- function(ndraws, nobs, dpars, data = list()) {
  structure(
    list(ndraws = ndraws, nobs = nobs,
         dpars = lapply(dpars, function(v) {
           if (length(v) == 1L) v else matrix(v, ndraws, nobs)
         }),
         data = data),
    class = "brmsprep"
  )
}
