# Test unified SDT API wrappers and trial-level prep helpers

test_that("sdt() dispatches to the expected specialized constructors", {
  expect_s3_class(
    sdt(response = c("r1", "r2", "r3", "r4"), stimulus = "stimulus"),
    "sdt_rating"
  )
  expect_s3_class(
    sdt(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"),
    "sdt_binary"
  )
  expect_s3_class(
    sdt(version = "dpsdt", response = c("r1", "r2", "r3", "r4"),
        stimulus = "stimulus"),
    "sdt_dp"
  )
})

test_that("rsdt() and dsdt() dispatch correctly", {
  dat <- rsdt(version = "rating", n_per_cell = 20, n_subjects = 2,
              dprime = 1, criterion = 0, n_ratings = 4, spacing = 0.5)
  expect_true(all(c("r1", "r2", "r3", "r4") %in% names(dat)))

  dens <- dsdt(version = "binary", n_old = 5, n_trials = 10,
               stimulus = 1, dprime = 1, criterion = 0)
  expect_true(dens > 0)
})

test_that("prepare_sdt_data aggregates binary trial data", {
  dat <- data.frame(
    id = rep(1:2, each = 4),
    stimulus = rep(c(0, 1), each = 2, times = 2),
    response = c(0, 1, 1, 1, 0, 0, 1, 1)
  )

  out <- prepare_sdt_data(
    dat,
    stimulus = "stimulus",
    response = "response",
    id_cols = "id",
    outcome = "binary"
  )

  expect_equal(names(out), c("id", "stimulus", "n_old", "n_trials"))
  expect_equal(out$n_trials, c(2L, 2L, 2L, 2L))
  expect_equal(out$n_old, c(1L, 2L, 0L, 2L))
})

test_that("prepare_sdt_data aggregates rating trial data", {
  dat <- data.frame(
    id = rep(1:2, each = 8),
    stimulus = rep(c(0, 1), each = 4, times = 2),
    response = c(0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1),
    confidence = c(2, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2)
  )

  out <- prepare_sdt_data(
    dat,
    stimulus = "stimulus",
    response = "response",
    confidence = "confidence",
    id_cols = "id",
    outcome = "rating",
    n_ratings = 4
  )

  expect_true(all(c("id", "stimulus", "r1", "r2", "r3", "r4", "nTrials") %in%
                    names(out)))
  expect_equal(out$nTrials, c(4L, 4L, 4L, 4L))
  expect_equal(rowSums(out[, paste0("r", 1:4)]), out$nTrials)
})
