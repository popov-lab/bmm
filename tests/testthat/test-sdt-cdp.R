# Tests for CDP (Continuous Dual-Process) SDT model

############################################################################# !
# CONSTRUCTOR TESTS                                                      ####
############################################################################# !

test_that("sdt_cdp constructor creates correct class chain", {
  model <- sdt_cdp(
    new_response = c("n1", "n2", "n3"),
    old_know = c("k4", "k5", "k6"),
    old_remember = c("r4", "r5", "r6"),
    stimulus = "stimulus"
  )
  expect_s3_class(model, "bmmodel")
  expect_s3_class(model, "sdt")
  expect_s3_class(model, "sdt_cdp")
  expect_equal(class(model), c("bmmodel", "sdt", "sdt_cdp"))
})

test_that("sdt_cdp constructor sets correct parameters for 2-way model", {
  model <- sdt_cdp(
    new_response = c("n1", "n2", "n3"),
    old_know = c("k4", "k5", "k6"),
    old_remember = c("r4", "r5", "r6"),
    stimulus = "stimulus"
  )
  expected_params <- c("dprimef", "dprimer", "criterion", "spacing",
                       "rcrit", "sigmar")
  expect_equal(sort(names(model$parameters)), sort(expected_params))
  expect_equal(model$fixed_parameters$mu, 0)
  expect_equal(model$fixed_parameters$sigmar, 0)
  expect_false(model$other_vars$has_guess)
})

test_that("sdt_cdp constructor sets correct parameters for 3-way model", {
  model <- sdt_cdp(
    new_response = c("n1", "n2", "n3"),
    old_know = c("k4", "k5", "k6"),
    old_remember = c("r4", "r5", "r6"),
    old_guess = c("g4", "g5", "g6"),
    stimulus = "stimulus"
  )
  expected_params <- c("dprimef", "dprimer", "criterion", "spacing",
                       "rcrit", "sigmar", "kcrit")
  expect_equal(sort(names(model$parameters)), sort(expected_params))
  expect_equal(model$fixed_parameters$kcrit, -100)
  expect_true(model$other_vars$has_guess)
})

test_that("sdt_cdp validates response vector lengths", {
  expect_error(
    sdt_cdp(new_response = c("n1", "n2"),
            old_know = c("k3", "k4", "k5"),
            old_remember = c("r3", "r4"),
            stimulus = "stimulus"),
    "same length"
  )
})

test_that("sdt_cdp validates minimum confidence levels", {
  expect_error(
    sdt_cdp(new_response = c("n1"),
            old_know = c("k2"),
            old_remember = c("r2"),
            stimulus = "stimulus"),
    "at least 2"
  )
})

test_that("sdt_cdp supports all noise distributions", {
  for (d in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_cdp(
      new_response = c("n1", "n2"),
      old_know = c("k3", "k4"),
      old_remember = c("r3", "r4"),
      stimulus = "stimulus",
      dist = d
    )
    expect_equal(model$other_vars$dist, d)
  }
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that(".sdt_cdp_category_probs sums to 1 for 2-way model", {
  thresholds <- c(-1.5, -0.5, 0.5, 1.5, 2.5)
  for (stim in c(0, 1)) {
    probs <- bmm:::.sdt_cdp_category_probs(
      thresholds, dprimef = 0.8, dprimer = 1.0, sigmar = 0,
      rcrit = 0.5, kcrit = NULL, stimulus = stim, n_conf = 3L,
      dist = "normal"
    )
    expect_equal(sum(probs), 1, tolerance = 1e-6,
                 label = paste("2-way probs sum, stim =", stim))
    expect_true(all(probs > 0))
    expect_length(probs, 9)  # 3 new + 3 know + 3 remember
  }
})

test_that(".sdt_cdp_category_probs sums to 1 for 3-way model", {
  thresholds <- c(-1.5, -0.5, 0.5, 1.5, 2.5)
  for (stim in c(0, 1)) {
    probs <- bmm:::.sdt_cdp_category_probs(
      thresholds, dprimef = 0.8, dprimer = 1.0, sigmar = 0,
      rcrit = 0.5, kcrit = 0.0, stimulus = stim, n_conf = 3L,
      dist = "normal"
    )
    expect_equal(sum(probs), 1, tolerance = 1e-6,
                 label = paste("3-way probs sum, stim =", stim))
    expect_true(all(probs > 0))
    expect_length(probs, 12)  # 3 new + 3 guess + 3 know + 3 remember
  }
})

test_that("3-way model with kcrit=-100 matches 2-way model", {
  thresholds <- c(-1.5, -0.5, 0.5, 1.5, 2.5)
  probs_2way <- bmm:::.sdt_cdp_category_probs(
    thresholds, dprimef = 0.8, dprimer = 1.0, sigmar = 0,
    rcrit = 0.5, kcrit = NULL, stimulus = 1, n_conf = 3L
  )
  probs_3way <- bmm:::.sdt_cdp_category_probs(
    thresholds, dprimef = 0.8, dprimer = 1.0, sigmar = 0,
    rcrit = 0.5, kcrit = -100, stimulus = 1, n_conf = 3L
  )
  # 3-way: n(3), g(3), k(3), r(3) → g should be ~0, k should match 2-way k
  # 2-way: n(3), k(3), r(3)
  expect_equal(probs_3way[1:3], probs_2way[1:3], tolerance = 1e-4,
               label = "new probs match")
  expect_true(all(probs_3way[4:6] < 1e-4),
              label = "guess probs near zero when kcrit=-100")
  expect_equal(probs_3way[7:9], probs_2way[4:6], tolerance = 1e-4,
               label = "know probs match")
  expect_equal(probs_3way[10:12], probs_2way[7:9], tolerance = 1e-4,
               label = "remember probs match")
})

test_that("rsdt_cdp generates data with correct structure", {
  dat <- rsdt_cdp(n_per_cell = 50, n_subjects = 5,
                  dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = 0.5,
                  rcrit = 0.5, n_ratings = 6)
  expect_equal(nrow(dat), 10)  # 5 subjects x 2 stimulus conditions
  expect_true(all(c("id", "stimulus", "n1", "n2", "n3",
                     "k4", "k5", "k6", "r4", "r5", "r6",
                     "nTrials") %in% names(dat)))
  expect_equal(unique(dat$nTrials), 50L)
})

test_that("rsdt_cdp 3-way generates data with guess columns", {
  dat <- rsdt_cdp(n_per_cell = 50, n_subjects = 3,
                  dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = 0.5,
                  rcrit = 0.5, kcrit = 0.0, n_ratings = 6)
  expect_true(all(c("g4", "g5", "g6") %in% names(dat)))
})

test_that("dsdt_cdp returns valid density", {
  counts <- c(50, 30, 20, 10, 15, 25, 5, 20, 75)
  dens <- dsdt_cdp(counts, stimulus = 1, dprimef = 0.8, dprimer = 1.0,
                   thresholds = c(-1.5, -0.5, 0, 0.5, 1.5),
                   rcrit = 0.5)
  expect_true(is.numeric(dens))
  expect_true(dens >= 0)
})


############################################################################# !
# PIPELINE TESTS (MOCK BACKEND)                                          ####
############################################################################# !

test_that("sdt_cdp 2-way pipeline runs with mock backend", {
  dat <- rsdt_cdp(n_per_cell = 50, n_subjects = 5,
                  dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = 0.5,
                  rcrit = 0.5, n_ratings = 6)

  model <- sdt_cdp(
    new_response = c("n1", "n2", "n3"),
    old_know = c("k4", "k5", "k6"),
    old_remember = c("r4", "r5", "r6"),
    stimulus = "stimulus"
  )

  f <- bmf(dprimef ~ 1, dprimer ~ 1, criterion ~ 1,
           spacing ~ 1, rcrit ~ 1)

  expect_silent(
    bmm(f, dat, model, backend = "mock", mock = 1, rename = FALSE)
  )
})
