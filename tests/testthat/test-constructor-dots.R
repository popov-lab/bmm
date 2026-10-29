# #552: model constructors report arguments they do not use. Each row is a
# constructor with valid arguments; `foo = 1` is added to it.
constructor_calls <- function() {
  tree <- mpt_tree("study", list(
    C = "cp + (1 - cp) * rp * rp",
    E = "2 * (1 - cp) * rp * (1 - rp)",
    U = "(1 - cp) * (1 - rp) * (1 - rp)"
  ))
  list(
    cswald = quote(cswald("rt", "resp")),
    ddm = quote(ddm("rt", "resp")),
    ezdm = quote(ezdm("m", "v", "u", "n")),
    imm = quote(imm("y", "nt", "d", 2)),
    IMMfull = quote(IMMfull("y", "nt", "d", 2)),
    IMMbsc = quote(IMMbsc("y", "nt", "d", 2)),
    IMMabc = quote(IMMabc("y", "nt", 2)),
    m3 = quote(m3(c("a", "b"), c(1, 2), version = "ss")),
    mixture2p = quote(mixture2p("y")),
    mixture3p = quote(mixture3p("y", nt_features = "x", set_size = 2)),
    mpt = bquote(mpt(.(tree))),
    sdm = quote(sdm("y")),
    sdmSimple = quote(sdmSimple("y")),
    sdt_cdp = quote(sdt_cdp(stimulus = "s", n_new = 2, n_old = 2)),
    sdt_mafc = quote(sdt_mafc("k", "n", m = 3)),
    sdt_ranking = quote(sdt_ranking(c("a", "b", "c"), m = 3)),
    sdt_rating = quote(sdt_rating(c("r1", "r2", "r3", "r4"), "s")),
    sdt_yn = quote(sdt_yn("r", "s", "n"))
  )
}

with_args <- function(call, ...) as.call(c(as.list(call), list(...)))

warnings_of <- function(call) collect_warnings(eval(call))$warnings

# only for calls that stop after the warning: a required argument is missing,
# or the call head is a function object, where stop_missing_args() fails
warnings_before_error <- function(call) collect_warnings(try(eval(call), silent = TRUE))$warnings

unused_message <- function(constructor, name, hint = "") {
  paste0("`", name, "` is not an argument of ", constructor, "() and is ignored.", hint, "\n",
         "See ?", constructor, " for its arguments. From bmm 1.6.0 on this is an error.")
}

two_close <- function(aaaaab, aaaaac, ...) .check_constructor_dots(list(...))

test_that("every constructor names the argument it does not use and itself", {
  calls <- constructor_calls()
  for (name in names(calls)) {
    msg <- grep("is not an argument", warnings_of(with_args(calls[[name]], foo = 1)), value = TRUE)
    expect_identical(msg, paste0(
      "`foo` is not an argument of ", name, "() and is ignored.\n",
      "See ?", name, " for its arguments. From bmm 1.6.0 on this is an error."
    ), label = name)
  }
})

test_that("all unused arguments are named in one warning", {
  msg <- warnings_of(quote(sdt_yn("r", "s", "n", foo = 1, bar = 2)))
  expect_length(msg, 1)
  expect_match(msg, "`foo` is not an argument of sdt_yn()", fixed = TRUE)
  expect_match(msg, "`bar` is not an argument of sdt_yn()", fixed = TRUE)
})

test_that("the names a constructor passes on to its model still work", {
  calls <- constructor_calls()
  passed_on <- list(
    imm = "links", IMMfull = "links", IMMbsc = "links",
    IMMabc = c("links", "nt_distances"), m3 = c("links", "default_priors"),
    mixture2p = "links", mixture3p = "links", mpt = "default_priors",
    sdm = "links", sdmSimple = "links"
  )
  for (name in names(passed_on)) {
    for (arg in passed_on[[name]]) {
      call <- with_args(calls[[name]], NULL)
      names(call)[length(call)] <- arg
      expect_false(any(grepl("is not an argument", warnings_of(call))), label = paste(name, arg))
    }
  }
  # R matches `link` to `links` by its prefix, so it is used
  expect_silent(sdm("y", link = list(kappa = "log")))
})

test_that("a prefix is matched only to an argument no exact name already took", {
  m3_call <- quote(m3(c("a", "b"), c(1, 2), "softmax", "ss", links = list(c = "log")))
  expect_identical(
    warnings_of(with_args(m3_call, l = list(a = "log"))), unused_message("m3", "l")
  )
  expect_length(warnings_of(m3_call), 0)
  expect_length(warnings_of(quote(m3(c("a", "b"), c(1, 2), "softmax", "ss", l = list(a = "log")))), 0)
})

test_that("setsize keeps its deprecation warning and nothing more", {
  setsize_warning <- "The argument 'setsize' is deprecated. Please use 'set_size' instead."
  deprecated_function <- function(name, version) {
    paste0("The function `", name, "()` is deprecated. Please use `imm(version = '", version, "')` instead.")
  }
  expect_identical(
    warnings_of(quote(mixture3p("y", nt_features = "x", setsize = 2))), setsize_warning
  )
  expect_identical(warnings_of(quote(imm("y", "nt", "d", setsize = 2))), setsize_warning)
  expect_identical(
    warnings_of(quote(IMMfull("y", "nt", "d", setsize = 2))),
    c(deprecated_function("IMMfull", "full"), setsize_warning)
  )
  expect_identical(
    warnings_of(quote(IMMbsc("y", "nt", "d", setsize = 2))),
    c(deprecated_function("IMMbsc", "bsc"), setsize_warning)
  )
  expect_identical(
    warnings_of(quote(IMMabc("y", "nt", setsize = 2))),
    c(deprecated_function("IMMabc", "abc"), setsize_warning)
  )
})

test_that("a misspelt argument is pointed to the close one", {
  expect_identical(
    warnings_of(quote(mixture2p("y", lnks = list(kappa = "log")))),
    unused_message("mixture2p", "lnks", " Did you mean `links`?")
  )
  expect_identical(
    warnings_of(quote(mixture2p("y", Links = list(kappa = "log")))),
    unused_message("mixture2p", "Links", " Did you mean `links`?")
  )
  expect_identical(
    warnings_before_error(quote(mixture3p("y", nt_features = "x", set_sze = 2))),
    unused_message("mixture3p", "set_sze", " Did you mean `set_size`?")
  )
})

test_that("a suggestion reaches as far as the length of the name allows", {
  # a transposition costs two edits
  expect_identical(
    warnings_of(quote(sdm("y", ilnks = 1))), unused_message("sdm", "ilnks", " Did you mean `links`?")
  )
  expect_identical(
    warnings_of(quote(sdt_yn("r", "s", "n", idst = 1))),
    unused_message("sdt_yn", "idst", " Did you mean `dist`?")
  )
  expect_identical(
    warnings_of(quote(sdt_yn("r", "s", "n", dst = 1))),
    unused_message("sdt_yn", "dst", " Did you mean `dist`?")
  )
  expect_identical(
    warnings_of(quote(sdm("y", linkzzs = 1))), unused_message("sdm", "linkzzs", " Did you mean `links`?")
  )
  # three letters at two edits, and eight letters at three edits
  expect_identical(warnings_of(quote(sdm("y", lik = 1))), unused_message("sdm", "lik"))
  expect_identical(warnings_of(quote(sdm("y", linkzzzs = 1))), unused_message("sdm", "linkzzzs"))
  # a two-letter name is two edits from the one-letter `m`
  expect_identical(
    warnings_before_error(quote(sdt_mafc("k", "n", xy = 1))), unused_message("sdt_mafc", "xy")
  )
})

test_that("a suggestion never offers an argument the call already gave", {
  expect_identical(
    warnings_of(quote(sdt_yn("r", "s", "n", dist = "normal", distt = "logistic"))),
    unused_message("sdt_yn", "distt")
  )
  expect_identical(
    warnings_of(quote(mixture3p("y", nt_features = "x", set_size = 2, set_sze = 2))),
    unused_message("mixture3p", "set_sze")
  )
  expect_identical(
    warnings_of(quote(sdm("y", links = list(kappa = "log"), lnks = 1))), unused_message("sdm", "lnks")
  )
  # `link` is `links` by prefix
  expect_identical(
    warnings_of(quote(sdm("y", link = list(kappa = "log"), lnks = 1))), unused_message("sdm", "lnks")
  )
})

test_that("two formals within reach of the misspelling give no guess", {
  expect_identical(
    warnings_of(quote(two_close(aaaaad = 1))), unused_message("two_close", "aaaaad")
  )
})

test_that("an argument given twice is reported once", {
  expect_identical(warnings_of(quote(sdm("y", foo = 1, foo = 2))), unused_message("sdm", "foo"))
})

test_that("a call through a function object does not print the function", {
  wording <- paste0("`foo` is not an argument of this model constructor and is ignored.\n",
                    "From bmm 1.6.0 on this is an error.")
  expect_identical(warnings_before_error(quote(do.call(sdm, list("y", foo = 1)))), wording)
  expect_identical(warnings_before_error(quote(Map(sdm, "y", foo = 1))), wording)
  expect_identical(warnings_before_error(quote(lapply("y", sdm, foo = 1))), wording)
  # sdm() takes four arguments by position (#558), so the fifth reaches `...`
  expect_match(
    warnings_before_error(quote(do.call(sdm, list("y", 1, 2, 3, 4)))),
    "This model constructor takes further arguments by name only", fixed = TRUE
  )
  expect_match(
    warnings_of(quote(bmm::sdm("y", foo = 1))),
    "`foo` is not an argument of bmm::sdm() and is ignored.\nSee ?bmm::sdm for",
    fixed = TRUE
  )
})

test_that("unnamed extra arguments are reported", {
  unnamed_message <- function(n) paste0(
    "sdt_yn() takes further arguments by name only; it was given ", n, ".\n",
    "See ?sdt_yn for its arguments. From bmm 1.6.0 on this is an error."
  )
  expect_identical(
    warnings_of(quote(sdt_yn("r", "s", "n", "normal", NULL, 5))), unnamed_message("1 unnamed argument")
  )
  expect_identical(
    warnings_of(quote(sdt_yn("r", "s", "n", "normal", NULL, 5, 6))), unnamed_message("2 unnamed arguments")
  )
})
