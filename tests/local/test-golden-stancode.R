# The golden Stan code snapshots of tests/testthat/test-model_mpt-stancode.R
# skip under any brms other than golden_brms, so after a brms release no PR
# checks them any more. This weekly test fails until they are re-recorded. It
# reads golden_brms from the parsed file rather than sourcing it, which would
# run that file's tests here.
test_that("the golden mpt Stan code was recorded under the installed brms", {
  is_golden_brms <- function(e) {
    is.call(e) && identical(e[[1]], as.name("<-")) && identical(e[[2]], as.name("golden_brms"))
  }
  assignment <- Filter(is_golden_brms, as.list(parse("../testthat/test-model_mpt-stancode.R")))
  expect_length(assignment, 1)
  golden_brms <- eval(assignment[[1]][[3]], baseenv())
  installed <- as.character(utils::packageVersion("brms"))
  expect(
    identical(golden_brms, installed),
    paste0(
      "The golden Stan code in tests/testthat/_snaps/model_mpt-stancode.md was ",
      "recorded under brms ", golden_brms, ", but brms ", installed, " is ",
      "installed, so its snapshot tests skip. Re-record it: set golden_brms <- \"",
      installed, "\" in tests/testthat/test-model_mpt-stancode.R, run that file, ",
      "read the diff and accept it with testthat::snapshot_accept(\"model_mpt-stancode\")."
    )
  )
})
