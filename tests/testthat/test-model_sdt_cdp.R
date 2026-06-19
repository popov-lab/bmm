test_that("sdt_cdp model can be created", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "sdt")
  expect_s3_class(m, "sdt_cdp")
})

test_that("sdt_cdp has the CDP parameters and no mu", {
  m <- sdt_cdp(judgment = "j", confidence = "c", count = "n", stimulus = "s")
  expect_setequal(
    names(m$parameters),
    c("dprimef", "dprimer", "criterion", "spacing", "rcrit", "sigmar",
      "rho", "kcrit")
  )
  expect_false("mu" %in% names(m$parameters))
  expect_false(isTRUE(m$void_mu))
})

test_that("sdt_cdp fixes sigmar, rho, and kcrit by default", {
  m <- sdt_cdp(judgment = "j", confidence = "c", count = "n", stimulus = "s")
  expect_equal(m$fixed_parameters$sigmar, 0)
  expect_equal(m$fixed_parameters$rho, 0)
  expect_equal(m$fixed_parameters$kcrit, -100)
})

test_that("sdt_cdp has identity links and accepts custom links", {
  m <- sdt_cdp(judgment = "j", confidence = "c", count = "n", stimulus = "s")
  expect_true(all(unlist(m$links) == "identity"))
  m2 <- sdt_cdp(judgment = "j", confidence = "c", count = "n", stimulus = "s",
                links = list(dprimef = "log"))
  expect_equal(m2$links$dprimef, "log")
})

test_that("sdt_cdp only supports dist = 'normal'", {
  expect_error(
    sdt_cdp(judgment = "j", confidence = "c", stimulus = "s", dist = "logistic"),
    "supports only dist"
  )
})

test_that("sdt_cdp requires judgment, confidence, and stimulus", {
  expect_error(sdt_cdp(confidence = "c", stimulus = "s"))
  expect_error(sdt_cdp(judgment = "j", stimulus = "s"))
  expect_error(sdt_cdp(judgment = "j", confidence = "c"))
})

test_that("sdt_cdp init_ranges cover the estimated parameters and exclude mu", {
  m <- sdt_cdp(judgment = "j", confidence = "c", count = "n", stimulus = "s")
  expect_true(all(names(m$parameters) %in% names(m$init_ranges)))
  expect_false("mu" %in% names(m$init_ranges))
})

test_that(".cdp_make_thresholds anchors the old/new boundary at criterion", {
  # symmetric: boundary is the middle threshold
  thr <- .cdp_make_thresholds(0.3, -0.2, n_new = 3, n_old = 3, "parsimonious")
  expect_length(thr, 5)
  expect_equal(thr[3], 0.3)             # boundary between bin 3 and 4
  expect_true(!is.unsorted(thr))
  # asymmetric (1 new + 5 old): boundary is thr[1]
  thr_a <- .cdp_make_thresholds(0.5, -0.1, n_new = 1, n_old = 5, "equidistant")
  expect_length(thr_a, 5)
  expect_equal(thr_a[1], 0.5)
  expect_true(!is.unsorted(thr_a))
})

test_that(".sdt_cdp_category_probs sums to 1 across variants", {
  thr <- .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious")
  # R/K, target and lure
  pt <- .sdt_cdp_category_probs(thr, 0.8, 1.0, 0, 0.5, NULL, 1, 3, 3, "normal", 0)
  pl <- .sdt_cdp_category_probs(thr, 0.8, 1.0, 0, 0.5, NULL, 0, 3, 3, "normal", 0)
  expect_length(pt, 9)
  expect_equal(sum(pt), 1, tolerance = 1e-10)
  expect_equal(sum(pl), 1, tolerance = 1e-10)
  # unequal variance + correlation
  puv <- .sdt_cdp_category_probs(thr, 0.8, 1.0, log(1.5), 0.5, NULL, 1, 3, 3,
                                 "normal", 0.4)
  expect_equal(sum(puv), 1, tolerance = 1e-10)
  # 3-way R/K/G
  p3 <- .sdt_cdp_category_probs(thr, 0.8, 1.0, 0, 0.6, 0.1, 1, 3, 3, "normal", 0)
  expect_length(p3, 12)
  expect_equal(sum(p3), 1, tolerance = 1e-10)
})

test_that(".sdt_cdp_category_probs guess + know recover the 2-way know total", {
  thr <- .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious")
  p2 <- .sdt_cdp_category_probs(thr, 0.8, 1.0, 0, 0.6, NULL, 1, 3, 3, "normal", 0)
  p3 <- .sdt_cdp_category_probs(thr, 0.8, 1.0, 0, 0.6, 0.1, 1, 3, 3, "normal", 0)
  # know(2-way) == guess + know (3-way) per old bin, after renormalising the
  # 2-way vector is identical because the bin totals are unchanged
  know2 <- p2[4:6]
  guess3 <- p3[4:6]
  know3 <- p3[7:9]
  expect_equal(know2, guess3 + know3, tolerance = 1e-8)
})

test_that("sdt_cdp_logmu matches .sdt_cdp_category_probs and preserves shape", {
  thr_name <- "parsimonious"
  args <- list(n_new = 3, n_old = 3, thresh = 1L, has_guess = 0L,
               dprimef = 0.8, dprimer = 1.0, criterion = 0.1, spacing = -0.3,
               rcrit = 0.5, sigmar = 0, rho = 0.2, kcrit = -100, stimulus = 1)
  thr <- .cdp_make_thresholds(args$criterion, args$spacing, 3, 3, thr_name)
  probs <- .sdt_cdp_category_probs(thr, args$dprimef, args$dprimer, args$sigmar,
                                   args$rcrit, NULL, 1, 3, 3, "normal", args$rho)
  for (cat in seq_len(9)) {
    val <- sdt_cdp_logmu(cat, 3, 3, 1L, 0L, args$dprimef, args$dprimer,
                         args$criterion, args$spacing, args$rcrit, args$sigmar,
                         args$rho, args$kcrit, args$stimulus)
    expect_equal(val, log(probs[cat]), tolerance = 1e-10)
  }
  # draws-by-observation shape is preserved
  mat <- matrix(0.8, nrow = 4, ncol = 2)
  out <- sdt_cdp_logmu(1, 3, 3, 1L, 0L, mat, mat, mat, matrix(-0.3, 4, 2),
                       mat, matrix(0, 4, 2), matrix(0, 4, 2),
                       matrix(-100, 4, 2), matrix(1, 4, 2))
  expect_equal(dim(out), c(4, 2))
})

test_that("rsdt_cdp returns tidy long data ready for sdt_cdp", {
  dat <- rsdt_cdp(n_per_cell = 100, n_subjects = 5, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.5, n_new = 3, n_old = 3)
  expect_true(all(c("id", "stimulus", "judgment", "confidence", "count") %in%
                    names(dat)))
  expect_setequal(unique(dat$judgment), c("new", "know", "remember"))
  # each subject x stimulus cell totals n_per_cell
  tot <- tapply(dat$count, list(dat$id, dat$stimulus), sum)
  expect_true(all(tot == 100))
})

test_that("dsdt_cdp computes a valid multinomial density", {
  thr <- .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious")
  d <- dsdt_cdp(counts = c(40, 20, 10, 15, 25, 35, 5, 20, 80), stimulus = 1,
                dprimef = 0.8, dprimer = 1.0, thresholds = thr, rcrit = 0.5,
                n_new = 3, n_old = 3)
  expect_true(is.finite(d))
  expect_true(d > 0)
  dl <- dsdt_cdp(counts = c(40, 20, 10, 15, 25, 35, 5, 20, 80), stimulus = 1,
                 dprimef = 0.8, dprimer = 1.0, thresholds = thr, rcrit = 0.5,
                 n_new = 3, n_old = 3, log = TRUE)
  expect_equal(log(d), dl, tolerance = 1e-10)
})

test_that("check_data pivots tidy long to a wide multinomial matrix", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 200, n_subjects = 4, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.5, n_new = 3, n_old = 3)
  cd <- check_data(m, dat, bmf(dprimef ~ 1))
  expect_equal(attr(cd, "n_new"), 3)
  expect_equal(attr(cd, "n_old"), 3)
  expect_false(attr(cd, "has_guess"))
  expect_equal(ncol(cd$Y), 9)
  expect_equal(nrow(cd), 8)                  # 4 subjects x 2 stimulus classes
  expect_equal(sum(cd$Y), sum(dat$count))    # counts preserved
  expect_true(all(cd$nTrials == 200))
})

test_that("check_data infers asymmetric n_new / n_old", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 300, n_subjects = 3, dprimef = 0.9, dprimer = 0.8,
                  criterion = 0.2, spacing = -0.2, rcrit = 0.9, sigmar = log(1.5),
                  n_new = 1, n_old = 5)
  cd <- check_data(m, dat, bmf(dprimef ~ 1))
  expect_equal(attr(cd, "n_new"), 1)
  expect_equal(attr(cd, "n_old"), 5)
  expect_equal(ncol(cd$Y), 11)
})

test_that("check_data aggregates trial-level data when count is omitted", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 150, n_subjects = 3, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.5, n_new = 3, n_old = 3)
  trials <- dat[rep(seq_len(nrow(dat)), dat$count),
                c("id", "stimulus", "judgment", "confidence")]
  cd <- check_data(m, trials, bmf(dprimef ~ 1))
  expect_equal(nrow(cd), 6)
  expect_equal(sum(cd$Y), sum(dat$count))
})

test_that("check_data detects the Know/Guess split from a 'guess' level", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 200, n_subjects = 3, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.6, kcrit = 0.0,
                  n_new = 3, n_old = 3)
  cd <- check_data(m, dat, bmf(dprimef ~ 1))
  expect_true(attr(cd, "has_guess"))
  expect_equal(ncol(cd$Y), 12)
})

test_that("check_data validates judgment labels, stimulus, and confidence partition", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  base <- data.frame(
    id = 1, stimulus = 1,
    judgment = c("new", "know", "remember"),
    confidence = c(1, 2, 2), count = c(10, 5, 7)
  )
  # bad judgment label
  bad_j <- base
  bad_j$judgment[1] <- "maybe"
  expect_error(check_data(m, bad_j, bmf(dprimef ~ 1)), "not allowed")
  # missing a required judgment type
  miss <- base[base$judgment != "remember", ]
  expect_error(check_data(m, miss, bmf(dprimef ~ 1)), "must contain")
  # bad stimulus coding
  bad_s <- base
  bad_s$stimulus <- 2
  expect_error(check_data(m, bad_s, bmf(dprimef ~ 1)), "0 .* and 1")
})

test_that("bmf2bf builds a multinomial non-linear formula with one mu per category", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 100, n_subjects = 3, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.5, n_new = 3, n_old = 3)
  cd <- check_data(m, dat, bmf(dprimef ~ 1))
  m$other_vars$n_new <- attr(cd, "n_new")
  m$other_vars$n_old <- attr(cd, "n_old")
  m$other_vars$has_guess <- attr(cd, "has_guess")
  bf <- bmf2bf(m, bmf(dprimef ~ 1))
  expect_s3_class(bf, "brmsformula")
  # 9 categories => 1 main + 8 mu non-linear formulas referencing sdt_cdp_logmu
  pforms <- vapply(bf$pforms, function(f) paste(deparse(f), collapse = " "),
                   character(1))
  expect_true(sum(grepl("sdt_cdp_logmu", pforms)) >= 8)
})

test_that("sdt_cdp produces multinomial stancode with the CDP functions", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 100, n_subjects = 4, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.5, n_new = 3, n_old = 3)
  sc <- stancode(bmf(dprimef ~ 1, dprimer ~ 1, criterion ~ 1, spacing ~ 1,
                     rcrit ~ 1), dat, m)
  expect_true(grepl("sdt_cdp_logmu", sc))
  expect_true(grepl("cdp_Phi2", sc))
  expect_true(grepl("owens_t", sc))
  expect_true(grepl("multinomial", sc))
})

test_that("freeing rho and sigmar adds them to the stancode", {
  m <- sdt_cdp(judgment = "judgment", confidence = "confidence",
               count = "count", stimulus = "stimulus")
  dat <- rsdt_cdp(n_per_cell = 100, n_subjects = 4, dprimef = 0.8, dprimer = 1.0,
                  criterion = 0, spacing = -0.3, rcrit = 0.5, n_new = 3, n_old = 3)
  sc <- stancode(bmf(dprimef ~ 1, dprimer ~ 1, criterion ~ 1, spacing ~ 1,
                     rcrit ~ 1, sigmar ~ 1, rho ~ 1), dat, m)
  expect_true(grepl("rho", sc))
  expect_true(grepl("sigmar", sc))
})
