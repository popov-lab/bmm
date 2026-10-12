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
