# The Stan chunks in inst/stan_chunks reach a compiler only when a model that
# uses them is fitted, and the compile tier loads only the chunks its parity
# and gradient tests need. These checks read every chunk as text and need no
# CmdStan, so an empty file, an unbalanced brace or an include that cannot
# resolve fails on every job instead of at a user's first fit (#563).

chunk_dir <- system.file("stan_chunks", package = "bmm")
chunk_files <- list.files(chunk_dir, pattern = "\\.stan$")

read_chunk <- function(file) read_lines2(file.path(chunk_dir, file))

# string literals go first so that the // of a URL inside reject() does not
# start a comment; the comments go because m3_gaussian_funs.stan and
# sdt_dist_funs.stan write braces in them (prod_{j != k}, F^{-1})
strip_stan_comments <- function(code) {
  gsub("(?s)\"[^\"\n]*\"|/\\*.*?\\*/|//[^\n]*", "", code, perl = TRUE)
}

# wrapped in one more pair, the code balances exactly when that outer brace
# closes at the last character: an unclosed brace leaves it open, a stray
# closing brace closes it early
braces_balance <- function(code) {
  wrapped <- paste0("{", strip_stan_comments(code), "}")
  close <- tryCatch(find_matching_brace(wrapped, 1L), error = function(e) NA_integer_)
  identical(close, nchar(wrapped))
}

test_that("every Stan chunk has code", {
  expect_gt(length(chunk_files), 0)
  has_code <- vapply(chunk_files, function(f) grepl("\\S", strip_stan_comments(read_chunk(f))), logical(1))
  empty <- chunk_files[!has_code]
  expect_identical(empty, character(0), info = paste("no code in", toString(empty)))
})

test_that("the braces of every Stan chunk balance", {
  expect_gt(length(chunk_files), 0)
  balanced <- vapply(chunk_files, function(f) braces_balance(read_chunk(f)), logical(1))
  unbalanced <- chunk_files[!balanced]
  expect_identical(unbalanced, character(0),
                   info = paste("unbalanced braces in", toString(unbalanced)))
})

test_that("every #include in a Stan chunk names a brms chunk", {
  # bmm passes no include path: brms's .stancode() replaces each
  # #include '<file>' with that file from its own chunks folder, matching
  # exactly this form, so a bmm chunk cannot be included and any other
  # spelling reaches stanc unexpanded; brms matches anywhere in the text, so an
  # include commented out with // is still expanded
  brms_chunks <- system.file("chunks", package = "brms")
  expect_true(dir.exists(brms_chunks))
  unresolved <- as.character(unlist(lapply(chunk_files, function(f) {
    code <- read_chunk(f)
    includes <- regmatches(code, gregexpr("#include[^\n]*", code))[[1]]
    target <- sub("^#include '([^']+)'.*", "\\1", includes)
    ok <- grepl("^#include '[^']+'", includes) & file.exists(file.path(brms_chunks, target))
    if (any(!ok)) paste0(f, ": ", includes[!ok])
  })))
  expect_identical(unresolved, character(0), info = paste("unresolved:", toString(unresolved)))
})
