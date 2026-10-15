ranks5 <- paste0("rank", 1:5)

test_that("check_data rejects NA in the set-size column", {
  model <- sdt_ranking(ranks5, m = "set_size")
  dat <- data.frame(id = 1, set_size = NA_integer_,
                    rank1 = 30, rank2 = 15, rank3 = 5, rank4 = 0, rank5 = 0)
  expect_error(check_data(model, dat, bmf(d ~ 1)),
               "Set-size column 'set_size' must not contain NA")
})

test_that("dsdt_ranking rejects NA in counts", {
  counts <- rbind(c(40, 30, 20, 10), c(10, 20, NA, 40))
  expect_error(dsdt_ranking(counts, m = 4, d = c(1.5, 0.2)),
               "counts must not contain NA")
})
