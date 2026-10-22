# Every use of an unexported bmm object is marked "bmm:::".
.onLoad <- function(libname, pkgname) {
  # check_formula() is not exported by bmm, so its method is registered into
  # bmm's own S3 table. `S3method(bmm::check_formula, slnorm)` in NAMESPACE
  # also works; the port to bmm's exported author surface replaces this
  registerS3method("check_formula", "slnorm", check_formula.slnorm,
                   envir = asNamespace("bmm"))   # bmm::: (unexported generic)
  bmm::register_bmm_model("slnorm", .model_slnorm, slnorm, package = pkgname)
}
