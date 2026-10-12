test_that("warn_deprecated() uses one fixed template", {
  expect_identical(
    collect_warnings(warn_deprecated("fit_model()", "bmm()", since = "1.5.0"))$warnings,
    paste(
      "`fit_model()` is deprecated since bmm 1.5.0 and will be removed in bmm",
      "1.6.0 at the earliest. Use `bmm()` instead."
    )
  )
})

test_that("removed_in defaults to the next minor release after since", {
  removed_in <- function(since) {
    msg <- collect_warnings(warn_deprecated("f()", "g()", since = since))$warnings
    sub(".*removed in bmm ([0-9.]+) at the earliest.*", "\\1", msg)
  }

  expect_identical(removed_in("1.5.0"), "1.6.0")
  expect_identical(removed_in("1.5.2"), "1.6.0")
  expect_identical(removed_in("1.4.0"), "1.5.0")
  expect_identical(removed_in("1.9.0"), "1.10.0")
})

test_that("an explicit removed_in replaces the default", {
  expect_match(
    collect_warnings(warn_deprecated("f()", "g()", since = "1.5.0", removed_in = "2.0.0"))$warnings,
    "removed in bmm 2.0.0 at the earliest",
    fixed = TRUE
  )
})

test_that("details are appended on a new line after the template", {
  expect_identical(
    collect_warnings(
      warn_deprecated("f()", "g()", since = "1.5.0", details = "g() also centres x.")
    )$warnings,
    paste0(
      "`f()` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest. ",
      "Use `g()` instead.\n",
      "g() also centres x."
    )
  )
})

test_that("many deprecated forms in one call give one warning with one statement each", {
  out <- collect_warnings(
    warn_deprecated(
      what = c("a ~ 1", "b ~ 1", "a ~ 1"),
      with = c("A ~ 1", "B ~ 1", "A ~ 1"),
      since = "1.5.0"
    )
  )

  expect_length(out$warnings, 1)
  expect_identical(
    out$warnings,
    paste(
      "`a ~ 1` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest. Use `A ~ 1` instead.",
      "`b ~ 1` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest. Use `B ~ 1` instead.",
      sep = "\n"
    )
  )
})

test_that("braces in the code are printed, not evaluated", {
  expect_identical(
    collect_warnings(
      warn_deprecated("f(x = '{y}')", "g(x = '{y}')", since = "1.5.0", details = "see {z}")
    )$warnings,
    paste0(
      "`f(x = '{y}')` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest. ",
      "Use `g(x = '{y}')` instead.\n",
      "see {z}"
    )
  )
})

test_that("no deprecated form gives no warning, not an empty one", {
  expect_silent(warn_deprecated(character(0), character(0), since = "1.5.0"))
  expect_silent(
    warn_deprecated(character(0), character(0), since = "1.5.0", details = "g() also centres x.")
  )
})

test_that("the warning names no call, so users do not see an internal function", {
  w <- tryCatch(warn_deprecated("f()", "g()", since = "1.5.0"), warning = identity)

  expect_null(conditionCall(w))
})

test_that("a deprecation without replacement ends the message at the removal release", {
  expect_identical(
    collect_warnings(warn_deprecated("f()", NULL, since = "1.5.0"))$warnings,
    "`f()` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest."
  )
  expect_identical(
    collect_warnings(
      warn_deprecated(c("f()", "g()"), NULL, since = "1.5.0", details = "Nothing replaces them.")
    )$warnings,
    paste(
      "`f()` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest.",
      "`g()` is deprecated since bmm 1.5.0 and will be removed in bmm 1.6.0 at the earliest.",
      "Nothing replaces them.",
      sep = "\n"
    )
  )
})

test_that("more than 20 forms list the first 20 and count the rest", {
  what <- sprintf("f%d()", 1:25)
  msg <- collect_warnings(warn_deprecated(what, "g()", since = "1.5.0", details = "See NEWS."))$warnings
  lines <- strsplit(msg, "\n", fixed = TRUE)[[1]]

  expect_length(lines, 22)
  expect_match(lines[1], "^`f1\\(\\)` is deprecated", perl = TRUE)
  expect_match(lines[20], "^`f20\\(\\)` is deprecated", perl = TRUE)
  expect_identical(lines[21], "Not listed: 5 more.")
  expect_identical(lines[22], "See NEWS.")
})

test_that("the warning is raised with R's longest warning length", {
  withr::local_options(warning.length = 1000L)
  length_when_raised <- NULL

  withCallingHandlers(
    warn_deprecated("f()", "g()", since = "1.5.0"),
    warning = function(w) {
      length_when_raised <<- getOption("warning.length")
      invokeRestart("muffleWarning")
    }
  )

  expect_identical(length_when_raised, 8170L)
})

test_that("the warning length option is the caller's again after the call", {
  withr::local_options(warning.length = 1000L)
  suppressWarnings(warn_deprecated("f()", "g()", since = "1.5.0"))

  expect_identical(getOption("warning.length"), 1000L)
})
