local_machine <- function(build_tools = TRUE,
                          cmdstanr = list(version = "0.9.0", error = NULL),
                          cmdstan = "2.40.0",
                          cmdstan_toolchain = NULL,
                          rstan = list(version = "2.32.7", error = NULL),
                          smoke = function(backend) list(error = NULL, seconds = 42),
                          os = "macos",
                          env = parent.frame()) {
  local_mocked_bindings(
    probe_build_tools = function() build_tools,
    probe_package = function(pkg) switch(pkg, cmdstanr = cmdstanr, rstan = rstan),
    probe_cmdstan_version = function() cmdstan,
    probe_cmdstan_toolchain = function() cmdstan_toolchain,
    probe_smoke_test = smoke,
    probe_os = function() os,
    .env = env
  )
}

not_installed <- list(version = NULL, error = "not installed")
smoke_must_not_run <- function(backend) stop("the smoke test ran")

row_of <- function(report, check) report[report$check == check, ]

test_that("a working machine passes every check and has no fixes", {
  local_machine()
  report <- suppressMessages(bmm_setup(backend = "cmdstanr"))

  expect_s3_class(report, "bmm_setup")
  expect_named(report, c("check", "status", "detail", "fix"))
  expect_equal(
    report$check,
    c("C++ toolchain", "cmdstanr", "CmdStan", "CmdStan toolchain", "rstan", "Backend", "Smoke test")
  )
  expect_equal(unique(report$status), "pass")
  expect_true(all(is.na(report$fix)))
})

test_that("the smoke test runs on the backend bmm() will use and reports its duration", {
  used <- NULL
  local_machine(smoke = function(backend) {
    used <<- backend
    list(error = NULL, seconds = 42.3)
  })

  report <- suppressMessages(bmm_setup(backend = "rstan"))

  expect_equal(used, "rstan")
  expect_match(row_of(report, "Smoke test")$detail, "42 s")
  expect_match(row_of(report, "Backend")$detail, "rstan")
})

test_that("the smoke test announces how long it takes before it starts", {
  local_machine()
  expect_message(bmm_setup(backend = "cmdstanr"), "small test model")
})

test_that("a missing C++ toolchain fails both backends and offers no switch", {
  local_machine(build_tools = FALSE, smoke = smoke_must_not_run)

  for (backend in c("cmdstanr", "rstan")) {
    report <- bmm_setup(backend = backend)
    expect_equal(row_of(report, "C++ toolchain")$status, "fail")
    expect_match(row_of(report, "C++ toolchain")$fix, "xcode-select --install", fixed = TRUE)
    expect_equal(row_of(report, "Backend")$status, "fail")
    expect_no_match(row_of(report, "Backend")$fix, "brms.backend", fixed = TRUE)
    expect_equal(row_of(report, "Smoke test")$status, "skip")
  }
})

test_that("without pkgbuild the toolchain check is skipped, not failed", {
  local_machine(build_tools = NA)
  report <- suppressMessages(bmm_setup(backend = "rstan"))

  expect_equal(row_of(report, "C++ toolchain")$status, "skip")
  expect_match(row_of(report, "C++ toolchain")$detail, "pkgbuild")
  expect_equal(row_of(report, "Backend")$status, "pass")
})

test_that("a missing cmdstanr is optional with rstan and fatal with cmdstanr", {
  local_machine(cmdstanr = not_installed, cmdstan = stop("not reached"))

  with_rstan <- suppressMessages(bmm_setup(backend = "rstan"))
  expect_equal(
    with_rstan$status,
    c("pass", "skip", "skip", "skip", "pass", "pass", "pass")
  )
  expect_true(all(is.na(with_rstan$fix)))

  local_mocked_bindings(probe_smoke_test = smoke_must_not_run)
  with_cmdstanr <- bmm_setup(backend = "cmdstanr")
  expect_equal(row_of(with_cmdstanr, "cmdstanr")$status, "fail")
  expect_match(row_of(with_cmdstanr, "cmdstanr")$fix, "https://stan-dev.r-universe.dev", fixed = TRUE)
  expect_equal(row_of(with_cmdstanr, "CmdStan")$status, "skip")
  expect_equal(row_of(with_cmdstanr, "Backend")$status, "fail")
  expect_equal(row_of(with_cmdstanr, "Smoke test")$status, "skip")
})

test_that("a missing CmdStan fails the cmdstanr backend and offers rstan instead", {
  local_machine(cmdstan = NULL, smoke = smoke_must_not_run)
  report <- bmm_setup(backend = "cmdstanr")

  expect_equal(row_of(report, "CmdStan")$status, "fail")
  expect_match(row_of(report, "CmdStan")$fix, "cmdstanr::install_cmdstan()", fixed = TRUE)
  expect_equal(row_of(report, "Backend")$status, "fail")
  expect_match(row_of(report, "Backend")$fix, 'options(brms.backend = "rstan")', fixed = TRUE)
  expect_equal(row_of(report, "Smoke test")$status, "skip")
})

test_that("a failing CmdStan toolchain check reports cmdstanr's reason", {
  reason <- "A C++ compiler was not found. Please install the command line tools."
  local_machine(cmdstan_toolchain = reason, smoke = smoke_must_not_run)
  report <- bmm_setup(backend = "cmdstanr")

  expect_equal(row_of(report, "CmdStan toolchain")$status, "fail")
  expect_match(row_of(report, "CmdStan toolchain")$detail, "C++ compiler was not found", fixed = TRUE)
  expect_match(row_of(report, "CmdStan toolchain")$fix, "xcode-select --install", fixed = TRUE)
  expect_equal(row_of(report, "Backend")$status, "fail")
})

test_that("a broken rstan fails both backends, because brms needs it after every fit", {
  broken <- list(version = NULL, error = "package or namespace load failed for 'rstan'")
  local_machine(rstan = broken, smoke = smoke_must_not_run)

  for (backend in c("cmdstanr", "rstan")) {
    report <- bmm_setup(backend = backend)
    expect_equal(row_of(report, "rstan")$status, "fail")
    expect_match(row_of(report, "rstan")$fix, 'remove.packages(c("rstan", "StanHeaders"))', fixed = TRUE)
    expect_equal(row_of(report, "Backend")$status, "fail")
    expect_no_match(row_of(report, "Backend")$fix, "brms.backend", fixed = TRUE)
    expect_equal(row_of(report, "Smoke test")$status, "skip")
  }
})

test_that("a failing smoke test reports the first line of its error and where to report it", {
  local_machine(smoke = function(backend) {
    list(error = "Stan program failed to compile\nline 2 of the compiler output", seconds = 3)
  })
  report <- suppressMessages(bmm_setup(backend = "cmdstanr"))
  smoke <- row_of(report, "Smoke test")

  expect_equal(smoke$status, "fail")
  expect_equal(smoke$detail, "Stan program failed to compile")
  expect_match(smoke$fix, "~/.R/Makevars", fixed = TRUE)
  expect_match(smoke$fix, "https://github.com/popov-lab/bmm/issues", fixed = TRUE)
  expect_no_match(smoke$fix, "\\\\")
  expect_equal(
    smoke$fix,
    paste(
      "look for a compiler set in ~/.R/Makevars or in the CXX or CXX17 environment",
      "variables; if there is none, please report this at",
      "https://github.com/popov-lab/bmm/issues, with the output of bmm_setup()"
    )
  )
})

test_that("smoke_test = FALSE skips the fit", {
  local_machine(smoke = smoke_must_not_run)
  report <- expect_silent(bmm_setup(smoke_test = FALSE, backend = "cmdstanr"))

  expect_equal(row_of(report, "Smoke test")$status, "skip")
  expect_true(is.na(row_of(report, "Smoke test")$fix))
})

test_that("backend = NULL checks the backend bmm() would choose", {
  withr::local_options(brms.backend = NULL)
  local_machine(smoke = smoke_must_not_run)
  report <- bmm_setup(smoke_test = FALSE, backend = NULL)

  chosen <- if (requireNamespace("cmdstanr", quietly = TRUE)) "cmdstanr" else "rstan"
  expect_equal(attr(report, "backend"), chosen)
  expect_match(row_of(report, "Backend")$detail, chosen)
})

test_that("the default backend is the brms.backend option", {
  withr::local_options(brms.backend = "rstan")
  local_machine(smoke = smoke_must_not_run)

  expect_equal(attr(bmm_setup(smoke_test = FALSE), "backend"), "rstan")
})

test_that("the Backend detail names why each backend was chosen", {
  local_machine(smoke = smoke_must_not_run)

  withr::with_options(list(brms.backend = "rstan"), {
    option_only <- bmm_setup(smoke_test = FALSE)
    expect_equal(
      row_of(option_only, "Backend")$detail,
      "bmm() will use rstan (set by options(brms.backend))"
    )

    option_and_explicit <- bmm_setup(smoke_test = FALSE, backend = "cmdstanr")
    expect_equal(row_of(option_and_explicit, "Backend")$detail, "cmdstanr (requested)")

    option_and_null <- bmm_setup(smoke_test = FALSE, backend = NULL)
    expect_equal(
      row_of(option_and_null, "Backend")$detail,
      "bmm() will use cmdstanr (the cmdstanr package is installed)"
    )
  })

  withr::local_options(brms.backend = NULL)
  neither <- bmm_setup(smoke_test = FALSE)
  expect_equal(
    row_of(neither, "Backend")$detail,
    "bmm() will use cmdstanr (the cmdstanr package is installed)"
  )
})

test_that("probe_package() reports a package that is not installed", {
  result <- probe_package("nonexistentpkgxyz")

  expect_equal(result$error, "not installed")
  expect_null(result$version)
})

test_that("bmm_setup() validates its arguments", {
  expect_error(bmm_setup(smoke_test = "yes"), "smoke_test")
  expect_error(bmm_setup(smoke_test = NA), "smoke_test")
  expect_error(bmm_setup(backend = "stan"), "should be one of")
})

test_that("the toolchain fix names the right tools for each operating system", {
  windows <- setup_fixes("windows", numeric_version("4.6.1"))
  macos <- setup_fixes("macos", numeric_version("4.6.1"))
  linux <- setup_fixes("linux", numeric_version("4.6.1"))

  expect_match(windows$toolchain, "Rtools", fixed = TRUE)
  expect_match(windows$toolchain, "R 4.6", fixed = TRUE)
  expect_match(windows$toolchain, "https://cran.r-project.org/bin/windows/Rtools/", fixed = TRUE)
  expect_match(macos$toolchain, "xcode-select --install", fixed = TRUE)
  expect_match(linux$toolchain, "build-essential", fixed = TRUE)
  expect_equal(length(unique(c(windows$toolchain, macos$toolchain, linux$toolchain))), 3)

  for (fix in c("cmdstanr", "cmdstan", "rstan")) {
    expect_equal(windows[[fix]], macos[[fix]])
    expect_equal(macos[[fix]], linux[[fix]])
  }
})

test_that("the report uses the fixes of the machine's operating system", {
  local_machine(build_tools = FALSE, os = "linux", smoke = smoke_must_not_run)
  report <- bmm_setup(backend = "rstan")

  expect_match(row_of(report, "C++ toolchain")$fix, "build-essential", fixed = TRUE)
})

test_that("print() shows one line per check and a fix under each failure only", {
  local_machine(cmdstan = NULL, smoke = smoke_must_not_run)
  report <- bmm_setup(backend = "cmdstanr")
  out <- capture.output(print(report, color = FALSE))

  expect_length(grep("^(PASS|FAIL|SKIP) ", out), 7)
  expect_length(grep("^FAIL ", out), 2)
  fixes <- grep("^ +Fix: ", out)
  expect_length(fixes, 2)
  expect_match(out[fixes - 1], "^FAIL ")
  expect_match(out[grep("^FAIL +CmdStan ", out) + 1], "install_cmdstan", fixed = TRUE)
  expect_match(out[length(out)], "2 checks failed")
})

test_that("print() says bmm is ready when everything passes", {
  local_machine()
  report <- suppressMessages(bmm_setup(backend = "cmdstanr"))
  out <- capture.output(print(report, color = FALSE))

  expect_match(out[1], "macOS")
  expect_match(out[length(out)], "can fit models with cmdstanr")
})

test_that("print() says which failures do not matter once the smoke test passed", {
  local_machine(cmdstan = NULL)
  report <- suppressMessages(bmm_setup(backend = "rstan"))
  out <- capture.output(print(report, color = FALSE))

  expect_equal(row_of(report, "CmdStan")$status, "fail")
  expect_match(out[length(out) - 1], "can fit models with rstan")
  expect_match(out[length(out)], "only for the cmdstanr backend")
})

test_that("print() without a smoke test does not claim bmm can fit models", {
  local_machine(smoke = smoke_must_not_run)
  report <- bmm_setup(smoke_test = FALSE, backend = "cmdstanr")
  out <- capture.output(print(report, color = FALSE))

  expect_no_match(out[length(out)], "can fit models")
  expect_match(out[length(out)], "run bmm_setup() to compile", fixed = TRUE)
})

test_that("a cmdstanr that fails to load under a CMDSTAN variable points to the variable", {
  withr::local_envvar(CMDSTAN = "/path/without/cmdstan")
  broken <- list(version = NULL, error = ".onLoad failed in loadNamespace() for 'cmdstanr'")
  local_machine(cmdstanr = broken, smoke = smoke_must_not_run)
  report <- bmm_setup(backend = "cmdstanr")

  expect_match(row_of(report, "cmdstanr")$fix, "/path/without/cmdstan", fixed = TRUE)
  expect_no_match(row_of(report, "cmdstanr")$fix, "r-universe", fixed = TRUE)
})

test_that("a missing cmdstanr is installed, whatever CMDSTAN says", {
  withr::local_envvar(CMDSTAN = "/path/without/cmdstan")
  local_machine(cmdstanr = not_installed, smoke = smoke_must_not_run)
  report <- bmm_setup(backend = "cmdstanr")

  expect_match(row_of(report, "cmdstanr")$fix, "r-universe", fixed = TRUE)
})

test_that("a subset of the report still prints", {
  local_machine(smoke = smoke_must_not_run)
  report <- bmm_setup(smoke_test = FALSE, backend = "cmdstanr")

  expect_no_error(capture.output(print(head(report, 3), color = FALSE)))
})

test_that("a column subset of the report falls back to data.frame printing", {
  local_machine(smoke = smoke_must_not_run)
  report <- bmm_setup(smoke_test = FALSE, backend = "cmdstanr")

  expect_no_error(capture.output(print(report[, c("check", "status")])))
})

test_that("print() shows the last messages of a failed smoke test", {
  local_machine(smoke = function(backend) {
    list(
      error = "An error occured during compilation! See the message above for more information.",
      seconds = 2,
      output = c("Compiling Stan program...", "make: /nonexistent/clang++: No such file or directory")
    )
  })
  report <- suppressMessages(bmm_setup(backend = "cmdstanr"))
  out <- capture.output(print(report, color = FALSE))

  expect_true(any(grepl("make: /nonexistent/clang++: No such file or directory", out, fixed = TRUE)))
})

# verbatim texts, provoked on macOS in throwaway processes (local/431_bmm_setup/error_texts.md)
setup_error_texts <- list(
  cmdstan_missing = simpleError("CmdStan path has not been set yet. See ?set_cmdstan_path."),
  cmdstanr_missing = simpleError("Please install the 'cmdstanr' package.", quote(require_package("cmdstanr"))),
  cmdstanr_compiler = simpleError("An error occured during compilation! See the message above for more information."),
  rstan_compiler = simpleError("invalid connection", quote(sink(type = "output"))),
  rstan_after_cmdstanr_fit = simpleError("unable to load required package ‘rstan’", quote(.requirePackage(package))),
  rstan_on_load = simpleError(paste(
    ".onLoad failed in loadNamespace() for 'rstan', details:",
    "  call: fun(libname, pkgname)",
    "  error: unable to load shared object '/fake/rstan.so'",
    sep = "\n"
  )),
  cmdstan_toolchain = simpleError(paste(
    "A suitable C++ compiler was not found. Please install the command line tools for Mac with",
    "'xcode-select --install' or install Xcode from the app store. Then restart R and run",
    "cmdstanr::check_cmdstan_toolchain()."
  )),
  rstan_compile_code = simpleError("Compilation ERROR, function(s)/method(s) not created!"),
  rstan_after_cmdstanr_fit_german = simpleError(
    "kann benötigtes Paket ‘rstan’ nicht laden", quote(.requirePackage(package))
  ),
  rstan_on_load_german = simpleError(".onLoad in loadNamespace() für 'rstan' fehlgeschlagen"),
  rstan_compiler_german = simpleError("ungültige Verbindung", quote(sink(type = "output"))),
  rstan_removed = simpleError("there is no package called ‘rstan’", quote(loadNamespace(x))),
  rstan_removed_german = simpleError("es gibt kein Paket namens ‘rstan’", quote(loadNamespace(x))),
  rstan_shared_object_broken = simpleError(
    "unable to load shared object '/Library/Frameworks/R.framework/Versions/4.6/Resources/library/rstan/libs/rstan.so': dlopen(...)",
    quote(dyn.load(file, DLLpath = DLLpath, ...))
  ),
  rstan_stanheaders_version = simpleError(
    "namespace ‘StanHeaders’ 2.26.28 is already loaded, but >= 2.32.0 is required",
    quote(loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]))
  )
)

other_error_texts <- list(
  prior = simpleError(paste(
    "The following priors do not correspond to any model parameter: ",
    "b_doesnotexist ~ normal(0, 1)",
    "Function 'default_prior' might be helpful to you.",
    sep = "\n"
  ), quote(.validate_prior(prior, bframe = bframe, sample_prior = sample_prior))),
  data = simpleError("The response variable 'y' is not present in the data."),
  init = simpleError("Assertion on 'init' failed: File does not exist: 'abc'.", quote(validate_init(self$init, num_inits))),
  connection_elsewhere = simpleError("invalid connection", quote(readLines(con))),
  sampling = simpleError("Fitting failed. Unable to retrieve the draws.")
)

test_that("is_setup_error() recognises toolchain and backend failures", {
  for (name in names(setup_error_texts)) {
    expect_true(is_setup_error(setup_error_texts[[name]]), label = name)
  }
})

test_that("is_setup_error() leaves prior, data, init and sampling errors alone", {
  for (name in names(other_error_texts)) {
    expect_false(is_setup_error(other_error_texts[[name]]), label = name)
  }
})

mock_bmm_error <- function(condition) {
  bmm(
    bmf(kappa ~ 1, thetat ~ 1), data.frame(y = c(-0.5, 0, 0.5)), mixture2p(resp_error = "y"),
    backend = "mock", mock_fit = function() stop(condition), rename = FALSE
  )
}

test_that("bmm() points a toolchain failure to bmm_setup() and keeps the error's class", {
  classed <- setup_error_texts$cmdstanr_missing
  class(classed) <- c("brms_error", class(classed))
  err <- expect_error(mock_bmm_error(classed), class = "brms_error")

  expect_match(conditionMessage(err), "^Please install the 'cmdstanr' package\\.")
  expect_match(conditionMessage(err), "Run bmm_setup()", fixed = TRUE)
  expect_match(conditionMessage(err), "https://github.com/popov-lab/bmm/issues", fixed = TRUE)
})

test_that("bmm() passes other errors from brm() through unchanged", {
  original <- other_error_texts$prior
  err <- expect_error(mock_bmm_error(original))

  expect_identical(conditionMessage(err), conditionMessage(original))
})

test_that("the toolchain probe ignores the compile flags rstan leaves behind", {
  skip_on_cran()
  skip_if_not_installed("pkgbuild")
  withr::local_envvar(PKG_CPPFLAGS = "-include nonexistent_header.hpp")
  expect_true(probe_build_tools())
})

test_that("probe_smoke_test() restores PKG_CPPFLAGS, PKG_LIBS and USE_CXX17 a user had set", {
  withr::local_envvar(PKG_CPPFLAGS = "-Ikeepme", PKG_LIBS = "-lkeepme", USE_CXX17 = "keepme")
  local_mocked_bindings(bmm = function(...) {
    Sys.setenv(PKG_CPPFLAGS = "-Irstan", PKG_LIBS = "-lrstan", USE_CXX17 = "1")
    stop("smoke ran")
  })
  probe_smoke_test("cmdstanr")

  expect_equal(Sys.getenv("PKG_CPPFLAGS"), "-Ikeepme")
  expect_equal(Sys.getenv("PKG_LIBS"), "-lkeepme")
  expect_equal(Sys.getenv("USE_CXX17"), "keepme")
})

test_that("a real machine check without a smoke test finds no failure", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  report <- bmm_setup(smoke_test = FALSE, backend = "cmdstanr")
  expect_false(
    any(report$status == "fail"),
    info = paste(capture.output(print(report, color = FALSE)), collapse = "\n")
  )
})

test_that("a real smoke test compiles and samples", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if(is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE)))

  report <- suppressMessages(bmm_setup(backend = "cmdstanr"))
  expect_equal(
    row_of(report, "Smoke test")$status, "pass",
    info = paste(capture.output(print(report, color = FALSE)), collapse = "\n")
  )
})
