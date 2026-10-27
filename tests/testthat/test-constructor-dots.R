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

test_that("setsize keeps its deprecation warning and nothing more", {
  expect_identical(
    warnings_of(quote(mixture3p("y", nt_features = "x", setsize = 2))),
    "The argument 'setsize' is deprecated. Please use 'set_size' instead."
  )
  expect_identical(
    warnings_of(quote(imm("y", "nt", "d", setsize = 2))),
    "The argument 'setsize' is deprecated. Please use 'set_size' instead."
  )
})

test_that("a misspelt argument is pointed to the close one", {
  expect_match(
    warnings_of(quote(mixture2p("y", lnks = list(kappa = "log")))),
    "`lnks` is not an argument of mixture2p() and is ignored. Did you mean `links`?",
    fixed = TRUE
  )
})

test_that("unnamed extra arguments are reported", {
  expect_match(
    warnings_of(quote(sdt_yn("r", "s", "n", "normal", NULL, 5))),
    "sdt_yn() takes further arguments by name only; it was given 1 unnamed argument.",
    fixed = TRUE
  )
})
