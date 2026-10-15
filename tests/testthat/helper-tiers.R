# Tier gates (ADR 0004). A new test that compiles a Stan program (a parity or
# gradient check, the bmm_setup() smoke test) opens with
# skip_unless_stan_tests(). The Windows and macOS pull-request jobs
# set BMM_TEST_STAN=false and ship no CmdStan; the Ubuntu jobs and every push
# run them. A new test that needs a stored fit loads it through
# load_fixture_fit(); one that needs only the assets folder opens with
# skip_unless_fixtures(). Tests that sample live in tests/local and need no
# gate here.
skip_unless_stan_tests <- function() {
  skip_if(
    identical(Sys.getenv("BMM_TEST_STAN"), "false"),
    "Stan compile tests are off on this job (BMM_TEST_STAN=false)"
  )
  skip_if_not_installed("cmdstanr")
  skip_if(
    is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)),
    "CmdStan is not installed"
  )
}

skip_unless_fixtures <- function() {
  skip_if_not(
    dir.exists(test_path("assets")),
    "fixtures are not in the built package (tests/testthat/assets is build-ignored)"
  )
}
