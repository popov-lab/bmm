# plyr registers an as.data.frame() method for functions that returns a
# function instead of failing. A pp_check() loads plyr (through bayesplot), so
# the regression showed only when an earlier test on the same worker had
# loaded it. A method of the same kind is registered here instead of loading
# plyr, which DESCRIPTION does not declare.
test_that("a function is not coercible data when as.data.frame() returns it", {
  local_mocked_s3_method("as.data.frame", "function", function(x, ...) x)
  expect_error(
    bmm(bmf(kappa ~ 1, thetat ~ 1), sum, mixture2p("y"), backend = "mock", mock_fit = 1),
    "must be coercible to a data.frame"
  )
})
