#' @title Check whether this machine can fit bmm models
#' @description Fitting a model with [bmm()] needs a C++ toolchain (a compiler
#'   and `make`) and a Stan backend, the software that turns the model into a
#'   program and samples from it. There are two backends: the R package
#'   `cmdstanr`, which runs CmdStan, a separate program installed with
#'   `cmdstanr::install_cmdstan()`; and the R package `rstan`. `rstan` must
#'   load with either backend, because `brms` uses it to store every fit.
#'   `bmm_setup()` checks each of these, prints whether it passed, and gives
#'   one fix for every check that failed. It does not install or change
#'   anything.
#'
#'   Two checks concern the toolchain. "C++ toolchain" lets R compile a small
#'   test file; "CmdStan toolchain" is `cmdstanr`'s own check, which only looks
#'   for `make` and a compiler on the `PATH`. The second can pass while the
#'   first fails, e.g. on a Mac whose command line tools are missing, and the
#'   first is the one that decides.
#'
#'   With `smoke_test = TRUE`, it finally compiles and samples a small
#'   two-parameter mixture model ([mixture2p()]) through [bmm()], on the
#'   backend [bmm()] will use. This is the only check that shows the whole
#'   chain works. On an Apple Silicon Mac it takes about 10 seconds with
#'   `cmdstanr` and 45 seconds with `rstan`, and it is skipped when a check it
#'   needs has failed. Once every check passes, rerun the [bmm()] call that
#'   failed.
#'
#' @param smoke_test Logical. Compile and sample a small test model? Defaults
#'   to `TRUE`.
#' @param backend The backend to check, `"cmdstanr"` or `"rstan"`. The default
#'   checks the backend [bmm()] would choose: the `brms.backend` option if it
#'   is set, otherwise `cmdstanr` whenever the `cmdstanr` package is
#'   installed, even when CmdStan is not, and `rstan` otherwise.
#'
#' @return A data frame of class `bmm_setup` with one row per check and the
#'   columns `check`, `status` (`"pass"`, `"fail"` or `"skip"`), `detail` and
#'   `fix` (`NA` unless the check failed). It prints as a report.
#'
#' @seealso [bmm()], [bmm_data_check()]
#' @export
#'
#' @examplesIf interactive()
#' # checks without compiling anything
#' bmm_setup(smoke_test = FALSE)
#'
#' # compiles and samples a small model as well
#' bmm_setup()
bmm_setup <- function(smoke_test = TRUE, backend = getOption("brms.backend", NULL)) {
  from_option <- missing(backend) && !is.null(backend)
  stopif(!isTRUE(smoke_test) && !isFALSE(smoke_test), "'smoke_test' must be TRUE or FALSE.")
  if (!is.null(backend)) {
    backend <- match.arg(backend, c("cmdstanr", "rstan"))
  }
  chosen <- configure_options(list(backend = backend))$backend %||% "rstan"
  os <- probe_os()
  fixes <- setup_fixes(os, getRversion())

  report <- rbind(
    toolchain_check(fixes),
    cmdstanr_checks(chosen, fixes),
    rstan_check(fixes)
  )
  report <- rbind(report, backend_check(report, chosen, backend_reason(backend, chosen, from_option)))
  smoke <- smoke_check(report, chosen, smoke_test)

  structure(
    rbind(report, smoke),
    class = c("bmm_setup", "data.frame"),
    backend = chosen,
    os = os,
    r_version = format(getRversion()),
    smoke_output = attr(smoke, "output")
  )
}

setup_row <- function(check, status, detail, fix = NA_character_) {
  data.frame(check = check, status = status, detail = as.character(detail), fix = as.character(fix))
}

toolchain_check <- function(fixes) {
  found <- probe_build_tools()
  if (is.na(found)) {
    return(setup_row("C++ toolchain", "skip", "not checked: needs the pkgbuild package"))
  }
  if (found) {
    return(setup_row("C++ toolchain", "pass", "R compiled a small test file"))
  }
  setup_row("C++ toolchain", "fail", "R could not compile a small test file", fixes$toolchain)
}

cmdstanr_checks <- function(chosen, fixes) {
  cmdstanr <- probe_package("cmdstanr")
  if (!is.null(cmdstanr$error)) {
    needs_cmdstanr <- setup_row(c("CmdStan", "CmdStan toolchain"), "skip", "not checked: needs cmdstanr")
    if (chosen == "rstan" && cmdstanr$error == "not installed") {
      return(rbind(setup_row("cmdstanr", "skip", "not installed; bmm() will use rstan"), needs_cmdstanr))
    }
    # cmdstanr 0.9.0 cannot load when CMDSTAN names a directory without
    # CmdStan, and brms then asks to install cmdstanr
    fix <- if (cmdstanr$error != "not installed" && nzchar(Sys.getenv("CMDSTAN"))) {
      glue("point the CMDSTAN environment variable to a CmdStan installation, or unset it \\
      (it is '{Sys.getenv('CMDSTAN')}')")
    } else {
      fixes$cmdstanr
    }
    return(rbind(setup_row("cmdstanr", "fail", cmdstanr$error, fix), needs_cmdstanr))
  }
  version <- probe_cmdstan_version()
  toolchain_error <- probe_cmdstan_toolchain()
  rbind(
    setup_row("cmdstanr", "pass", cmdstanr$version),
    if (is.null(version)) {
      setup_row("CmdStan", "fail", "not installed", fixes$cmdstan)
    } else {
      setup_row("CmdStan", "pass", version)
    },
    if (is.null(toolchain_error)) {
      setup_row("CmdStan toolchain", "pass", "make and a C++ compiler are on the PATH")
    } else {
      setup_row("CmdStan toolchain", "fail", toolchain_error, fixes$toolchain)
    }
  )
}

rstan_check <- function(fixes) {
  rstan <- probe_package("rstan")
  if (is.null(rstan$error)) {
    return(setup_row("rstan", "pass", rstan$version))
  }
  setup_row("rstan", "fail", rstan$error, fixes$rstan)
}

# brms turns the draws of every fit into a stanfit with rstan, also under
# cmdstanr, so a broken rstan fails both backends after sampling. cmdstanr's
# own toolchain check only looks for make and clang++ on the PATH, which the
# macOS stubs satisfy without the command line tools; pkgbuild's test compile
# is the check that notices
backend_requirements <- list(
  cmdstanr = c("C++ toolchain", "cmdstanr", "CmdStan", "CmdStan toolchain", "rstan"),
  rstan = c("C++ toolchain", "rstan")
)

backend_works <- function(report, backend) {
  !any(report$status[report$check %in% backend_requirements[[backend]]] == "fail")
}

backend_reason <- function(backend, chosen, from_option) {
  if (from_option) {
    return("set by options(brms.backend)")
  }
  if (!is.null(backend)) {
    return("requested")
  }
  if (chosen == "cmdstanr") "the cmdstanr package is installed" else "the brms default"
}

backend_check <- function(report, chosen, reason) {
  detail <- if (reason == "requested") glue("{chosen} (requested)") else glue("bmm() will use {chosen} ({reason})")
  if (backend_works(report, chosen)) {
    return(setup_row("Backend", "pass", detail))
  }
  fix <- if (chosen == "cmdstanr" && backend_works(report, "rstan")) {
    'fix the failures above, or switch to rstan with options(brms.backend = "rstan")'
  } else {
    "nothing separate; this passes once the failures above are fixed"
  }
  setup_row("Backend", "fail", detail, fix)
}

# measured on an M-series Mac with CmdStan 2.40 and rstan 2.32
smoke_seconds <- c(cmdstanr = 10, rstan = 45)

smoke_check <- function(report, chosen, smoke_test) {
  if (!smoke_test) {
    return(setup_row("Smoke test", "skip", "not run (smoke_test = FALSE)"))
  }
  if (!backend_works(report, chosen)) {
    return(setup_row("Smoke test", "skip", glue("not run: {chosen} does not work yet")))
  }
  message2(
    "Compiling and sampling a small test model with {chosen}. This takes about \\
    {smoke_seconds[[chosen]]} s on an Apple Silicon Mac, longer on slower machines."
  )
  smoke <- probe_smoke_test(chosen)
  if (is.null(smoke$error)) {
    return(setup_row("Smoke test", "pass", glue("compiled and sampled a mixture2p model in {round(smoke$seconds)} s")))
  }
  structure(
    setup_row(
      "Smoke test", "fail", strsplit(smoke$error, "\n", fixed = TRUE)[[1]][1],
      glue("look for a compiler set in ~/.R/Makevars or in the CXX or CXX17 environment \\
      variables; if there is none, please report this at \\
      https://github.com/popov-lab/bmm/issues, with the output of bmm_setup()")
    ),
    output = smoke$output
  )
}

setup_fixes <- function(os, r_version) {
  list(
    toolchain = switch(os,
      windows = glue(
        "install the Rtools version for R {format(r_version[, 1:2])} from \\
        https://cran.r-project.org/bin/windows/Rtools/"
      ),
      macos = "run xcode-select --install in the Terminal",
      linux = "install g++ and make with your package manager, e.g. sudo apt install build-essential"
    ),
    cmdstanr = 'install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))',
    cmdstan = "cmdstanr::install_cmdstan()",
    rstan = 'remove.packages(c("rstan", "StanHeaders")), then in a fresh R session install.packages("rstan")'
  )
}

probe_os <- function() {
  if (.Platform$OS.type == "windows") {
    return("windows")
  }
  if (Sys.info()[["sysname"]] == "Darwin") "macos" else "linux"
}

# debug = TRUE skips pkgbuild's per-session cache, so a rerun after a fix sees
# the fix; buildtools.check = NULL keeps RStudio from offering an install
probe_build_tools <- function() {
  if (!requireNamespace("pkgbuild", quietly = TRUE)) {
    return(NA)
  }
  withr::local_options(buildtools.check = NULL, pkgbuild.has_compiler = NULL)
  # rstan leaves its compile flags in the session environment; they force a
  # C++ header into pkgbuild's C test file and fail an intact toolchain
  withr::local_envvar(PKG_CPPFLAGS = NA, PKG_LIBS = NA, USE_CXX17 = NA)
  utils::capture.output(found <- suppressMessages(pkgbuild::has_build_tools(debug = TRUE)))
  found
}

probe_package <- function(pkg) {
  if (!nzchar(system.file(package = pkg))) {
    return(list(version = NULL, error = "not installed"))
  }
  tryCatch(
    {
      loadNamespace(pkg)
      list(version = format(utils::packageVersion(pkg)), error = NULL)
    },
    error = function(e) list(version = NULL, error = conditionMessage(e))
  )
}

probe_cmdstan_version <- function() {
  cmdstanr::cmdstan_version(error_on_NA = FALSE)
}

probe_cmdstan_toolchain <- function() {
  tryCatch(
    {
      cmdstanr::check_cmdstan_toolchain(fix = FALSE, quiet = TRUE)
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

# cmdstanr draws the sampler seed from R's generator; restoring it leaves the
# user's random stream where it was
probe_smoke_test <- function(backend) {
  withr::local_envvar(PKG_CPPFLAGS = NA, PKG_LIBS = NA, USE_CXX17 = NA)
  withr::local_options(try.outFile = nullfile())
  seconds <- system.time(
    fit <- capture_check_conditions(withr::with_preserve_seed(bmm(
      bmmformula(kappa ~ 1, thetat ~ 1),
      data.frame(y = seq(-1, 1, length.out = 50)),
      mixture2p(resp_error = "y"),
      backend = backend, chains = 1, iter = 200, refresh = 0,
      silent = 2, sort_data = FALSE
    )))
  )[["elapsed"]]
  list(
    error = fit$error,
    seconds = seconds,
    output = unique(trimws(c(fit$messages, fit$warnings)))
  )
}

#' @export
print.bmm_setup <- function(x, color = getOption("bmm.color_summary", TRUE), ...) {
  withr::local_options(bmm.color_summary = color)
  os <- c(windows = "Windows", macos = "macOS", linux = "Linux")[[attr(x, "os")]]
  cat(glue("bmm setup check ({os}, R {attr(x, 'r_version')})"), "\n\n", sep = "")
  checks <- format(x$check)
  for (i in seq_len(nrow(x))) {
    print_setup_row(x$status[i], checks[i], x$detail[i])
    if (!is.na(x$fix[i])) {
      cat_wrapped(x$fix[i], "      Fix: ")
    }
  }
  output <- utils::tail(attr(x, "smoke_output"), 10)
  if (length(output) > 0) {
    cat("\nLast messages of the failed smoke test:\n")
    cat(paste0("  ", output), sep = "\n")
  }
  cat("\n")
  cat_wrapped(setup_verdict(x), "")
  invisible(x)
}

# the label is colored after wrapping, so its escape codes do not count
# towards the indent of continuation lines
print_setup_row <- function(status, check, detail) {
  label <- toupper(status)
  prefix <- paste0(label, "  ", check, "  ")
  lines <- strwrap(gsub("\\s+", " ", detail), width = 80, initial = prefix, exdent = nchar(prefix))
  colors <- c(pass = "green", fail = "red", skip = "grey50")
  lines[1] <- sub(label, style(colors[[status]])(label), lines[1], fixed = TRUE)
  cat(lines, sep = "\n")
}

setup_verdict <- function(x) {
  n_fail <- sum(x$status == "fail")
  backend <- attr(x, "backend")
  if (isTRUE(any(x$status[x$check == "Smoke test"] == "pass"))) {
    ready <- glue("bmm can fit models with {backend} on this machine.")
    if (n_fail == 0) {
      return(ready)
    }
    other <- setdiff(names(backend_requirements), backend)
    return(c(ready, glue("The failures above matter only for the {other} backend.")))
  }
  if (n_fail == 0) {
    return("No check failed; run bmm_setup() to compile and sample a test model too.")
  }
  if (n_fail == 1) {
    return("1 check failed. Fix it, restart R, and run bmm_setup() again.")
  }
  glue("{n_fail} checks failed. Fix them from the top, restart R, and run bmm_setup() again.")
}

# texts of failures that bmm_setup() can diagnose; each one was provoked on a
# broken machine or read from the installed sources of brms, cmdstanr, inline
setup_error_patterns <- c(
  "CmdStan path has not been set",
  "Please install the '(cmdstanr|rstan)' package",
  "error occurr?ed during compilation",
  "Compilation ERROR",
  "check_cmdstan_toolchain",
  # the call and the package name read the same in every language R speaks
  "requirePackage\\(package\\).*.(rstan|StanHeaders).",
  "loadNamespace\\(.*.(rstan|StanHeaders|cmdstanr|RcppParallel).",
  # a corrupt shared object is raised by dyn.load(), unwrapped, from loadNamespace()
  "dyn\\.load\\(.*[/\\\\](rstan|StanHeaders|RcppParallel)[/\\\\]libs[/\\\\]",
  # rstan 2.32 reports any failed compilation as an invalid connection here
  'sink\\(type = "output"\\)'
)

is_setup_error <- function(e) {
  text <- paste(paste(deparse(conditionCall(e)), collapse = " "), conditionMessage(e))
  any(vapply(setup_error_patterns, grepl, logical(1), x = text))
}

# re-signals the same condition, so its class and call survive and handlers
# further up see only the longer message
add_setup_hint <- function(e) {
  if (is_setup_error(e)) {
    e$message <- paste0(
      e$message, "\n\nThis error can come from the C++ toolchain, the Stan backend, or the Stan code ",
      "bmm generates. Run bmm_setup() to check the first two; if it finds nothing, please report the error ",
      "at https://github.com/popov-lab/bmm/issues."
    )
    stop(e)
  }
}
