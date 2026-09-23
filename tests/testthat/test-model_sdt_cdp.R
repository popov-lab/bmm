sim_cdp_data <- function(n_subjects = 4, n_trials = 100, n_new = 3, n_old = 3,
                         kcrit = NULL, thresholds = NULL, ...) {
  dat <- expand.grid(id = seq_len(n_subjects), stimulus = c(0L, 1L))
  thresholds <- thresholds %||%
    .cdp_make_thresholds(0, -0.3, n_new, n_old, "parsimonious")
  cbind(dat, rsdt_cdp(nrow(dat), n_trials, dat$stimulus,
                      dfam = 0.8, drec = 1.0, thresholds = thresholds,
                      rcrit = 0.5, kcrit = kcrit, n_new = n_new, ...))
}

test_that("sdt_cdp model can be created", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  expect_s3_class(m, "bmmodel")
  expect_s3_class(m, "sdt")
  expect_s3_class(m, "sdt_cdp")
})

test_that("sdt_cdp has the CDP parameters and no mu", {
  m <- sdt_cdp(stimulus = "s", n_new = 3, n_old = 3)
  expect_setequal(
    names(m$parameters),
    c("dfam", "drec", "criterion", "spacing", "rcrit", "sigmar",
      "rho", "kcrit")
  )
  expect_false("mu" %in% names(m$parameters))
  expect_false(isTRUE(m$void_mu))
})

test_that("sdt_cdp fixes sigmar, rho, and kcrit by default", {
  m <- sdt_cdp(stimulus = "s", n_new = 3, n_old = 3)
  expect_equal(m$fixed_parameters$sigmar, 0)
  expect_equal(m$fixed_parameters$rho, 0)
  expect_equal(m$fixed_parameters$kcrit, -100)
})

test_that("sdt_cdp has identity links and accepts custom links", {
  m <- sdt_cdp(stimulus = "s", n_new = 3, n_old = 3)
  expect_true(all(unlist(m$links) == "identity"))
  m2 <- sdt_cdp(stimulus = "s", n_new = 3, n_old = 3,
                links = list(dfam = "log"))
  expect_equal(m2$links$dfam, "log")
})

test_that("sdt_cdp only supports dist = 'normal'", {
  expect_error(
    sdt_cdp(stimulus = "s", n_new = 3, n_old = 3, dist = "logistic"),
    "supports only dist"
  )
})

test_that("sdt_cdp requires stimulus, n_new, and n_old", {
  expect_error(sdt_cdp(n_new = 3, n_old = 3))
  expect_error(sdt_cdp(stimulus = "s", n_old = 3))
  expect_error(sdt_cdp(stimulus = "s", n_new = 3))
  expect_error(sdt_cdp(stimulus = "s", n_new = 1, n_old = 1), "sum to >= 3")
})

test_that("sdt_cdp init_ranges cover the estimated parameters and exclude mu", {
  m <- sdt_cdp(stimulus = "s", n_new = 3, n_old = 3)
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

test_that(".cdp_make_thresholds log_distance builds cumulative log-distances", {
  # symmetric 3/3: anchor at thr[3]; distances are exp(delta) away from criterion
  d <- log(c(0.8, 0.6, 0.6, 0.8))   # delta1, delta2, delta4, delta5
  thr <- .cdp_make_thresholds(0.5, 0, 3, 3, "log_distance", deltas = d)
  expect_length(thr, 5)
  expect_equal(thr[3], 0.5)
  expect_equal(thr[4] - thr[3], 0.6)   # exp(delta4)
  expect_equal(thr[5] - thr[4], 0.8)   # exp(delta5)
  expect_equal(thr[3] - thr[2], 0.6)   # exp(delta2)
  expect_equal(thr[2] - thr[1], 0.8)   # exp(delta1)
  expect_true(!is.unsorted(thr))
  # asymmetric 1/5: anchor thr[1], all four distances above
  thr_a <- .cdp_make_thresholds(0.2, 0, 1, 5, "log_distance",
                                deltas = log(c(0.6, 0.7, 0.8, 0.9)))
  expect_equal(thr_a[1], 0.2)
  expect_equal(thr_a[2] - thr_a[1], 0.6)
  expect_true(!is.unsorted(thr_a))
  # asymmetric 4/2: anchor thr[4], three distances below + one above
  thr_b <- .cdp_make_thresholds(0.1, 0, 4, 2, "log_distance",
                                deltas = log(c(0.5, 0.6, 0.7, 0.8)))
  expect_equal(thr_b[4], 0.1)
  expect_equal(thr_b[5] - thr_b[4], 0.8)
  expect_true(!is.unsorted(thr_b))
  # wrong delta length errors
  expect_error(.cdp_make_thresholds(0.5, 0, 3, 3, "log_distance",
                                    deltas = log(c(0.8, 0.6))))
})

test_that(".cdp_make_thresholds vectorizes over draws", {
  crit <- c(0.1, -0.2, 0.4)
  sp <- c(-0.3, 0, 0.2)
  thr <- .cdp_make_thresholds(crit, sp, 3, 3, "parsimonious")
  expect_equal(dim(thr), c(3, 5))
  for (i in 1:3) {
    expect_equal(thr[i, ],
                 .cdp_make_thresholds(crit[i], sp[i], 3, 3, "parsimonious"))
  }
  deltas <- rbind(log(c(0.8, 0.6, 0.6, 0.8)), log(c(0.5, 0.7, 0.9, 1.1)))
  thr_ld <- .cdp_make_thresholds(c(0.5, -0.1), 0, 3, 3, "log_distance",
                                 deltas = deltas)
  expect_equal(dim(thr_ld), c(2, 5))
  for (i in 1:2) {
    expect_equal(thr_ld[i, ],
                 .cdp_make_thresholds(deltas = deltas[i, ],
                                      c(0.5, -0.1)[i], 0, 3, 3,
                                      "log_distance"))
  }
})

test_that("sdt_cdp log_distance declares per-distance deltas and no spacing", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3,
               threshold_type = "log_distance")
  expect_true(all(c("delta1", "delta2", "delta4", "delta5") %in%
                    names(m$parameters)))
  expect_false("delta3" %in% names(m$parameters))
  expect_false("spacing" %in% names(m$parameters))
  expect_true(all(c("delta1", "delta2", "delta4", "delta5") %in%
                    names(m$init_ranges)))
})

test_that("sdt_cdp log_distance integrates with the pipeline via mock backend", {
  thr <- .cdp_make_thresholds(0.3, 0, 3, 3, "log_distance",
                              deltas = log(c(0.8, 0.6, 0.6, 0.8)))
  dat <- sim_cdp_data(thresholds = thr)
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3,
               threshold_type = "log_distance")
  expect_silent(
    bmm(bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, rcrit ~ 1,
            delta1 ~ 1, delta2 ~ 1, delta4 ~ 1, delta5 ~ 1),
        dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that(".sdt_cdp_logmu_stan generates a log_distance wrapper with delta args", {
  m <- sdt_cdp(stimulus = "s", n_new = 3, n_old = 3,
               threshold_type = "log_distance")
  sc <- .sdt_cdp_logmu_stan(m)
  expect_true(grepl("real sdt_cdp_logmu", sc))
  expect_true(grepl("real delta1.*real delta2.*real delta4.*real delta5", sc))
  expect_true(grepl("array\\[4\\] real deltas", sc))
  expect_true(grepl("cdp_make_thresholds\\(criterion, spacing, deltas", sc))
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

test_that(".sdt_cdp_category_probs vectorizes over observations", {
  n <- 10
  thr <- .cdp_make_thresholds(rnorm(n, 0, 0.2), rnorm(n, -0.3, 0.1), 3, 3,
                              "parsimonious")
  df <- rnorm(n, 0.8, 0.2)
  dr <- rnorm(n, 1.0, 0.2)
  sg <- rnorm(n, 0.1, 0.1)
  rh <- rnorm(n, 0.2, 0.2)
  rc <- rnorm(n, 0.5, 0.1)
  kc <- rnorm(n, -0.3, 0.1)
  st <- rep_len(c(0L, 1L), n)
  pv <- .sdt_cdp_category_probs(thr, df, dr, sg, rc, kc, st, 3, 3, "normal", rh)
  expect_equal(dim(pv), c(n, 12))
  for (j in seq_len(n)) {
    expect_equal(pv[j, ],
                 .sdt_cdp_category_probs(thr[j, ], df[j], dr[j], sg[j], rc[j],
                                         kc[j], st[j], 3, 3, "normal", rh[j]))
  }
})

test_that("sdt_cdp_logmu matches .sdt_cdp_category_probs and preserves shape", {
  thr_name <- "parsimonious"
  args <- list(n_new = 3, n_old = 3, thresh = 1L, has_guess = 0L,
               dfam = 0.8, drec = 1.0, criterion = 0.1, spacing = -0.3,
               rcrit = 0.5, sigmar = 0, rho = 0.2, kcrit = -100, stimulus = 1)
  thr <- .cdp_make_thresholds(args$criterion, args$spacing, 3, 3, thr_name)
  probs <- .sdt_cdp_category_probs(thr, args$dfam, args$drec, args$sigmar,
                                   args$rcrit, NULL, 1, 3, 3, "normal", args$rho)
  for (cat in seq_len(9)) {
    val <- sdt_cdp_logmu(cat, 3, 3, 1L, 0L, args$dfam, args$drec,
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

test_that("sdt_cdp_logmu vectorized path matches per-element evaluation", {
  nd <- 6
  nobs <- 4
  set.seed(11)
  mk <- function(center, spread) matrix(rnorm(nd * nobs, center, spread), nd, nobs)
  a <- list(dfam = mk(0.8, 0.2), drec = mk(1.0, 0.2),
            criterion = mk(0.1, 0.1), spacing = mk(-0.3, 0.1),
            rcrit = mk(0.5, 0.1), sigmar = mk(0.1, 0.1), rho = mk(0.2, 0.1),
            kcrit = mk(-0.3, 0.1),
            stimulus = matrix(rep(c(0, 0, 1, 1), each = nd), nd, nobs))
  for (cat in c(1L, 5L, 9L, 12L)) {
    vec <- sdt_cdp_logmu(cat, 3, 3, 1L, 1L, a$dfam, a$drec, a$criterion,
                         a$spacing, a$rcrit, a$sigmar, a$rho, a$kcrit,
                         a$stimulus)
    ref <- matrix(0, nd, nobs)
    for (i in seq_len(nd)) {
      for (j in seq_len(nobs)) {
        ref[i, j] <- sdt_cdp_logmu(cat, 3, 3, 1L, 1L, a$dfam[i, j],
                                   a$drec[i, j], a$criterion[i, j],
                                   a$spacing[i, j], a$rcrit[i, j],
                                   a$sigmar[i, j], a$rho[i, j], a$kcrit[i, j],
                                   a$stimulus[i, j])
      }
    }
    expect_equal(vec, ref, tolerance = 1e-12)
  }
})

test_that("rsdt_cdp returns the count columns sdt_cdp expects", {
  dat <- expand.grid(id = 1:5, stimulus = c(0L, 1L))
  thr <- .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious")
  cnts <- rsdt_cdp(nrow(dat), 100, dat$stimulus, 0.8, 1.0, thresholds = thr,
                   rcrit = 0.5, n_new = 3)
  expect_equal(colnames(cnts),
               c("new1", "new2", "new3", "know4", "know5", "know6",
                 "remember4", "remember5", "remember6"))
  expect_true(all(rowSums(cnts) == 100))
  # guess columns appear with a finite kcrit
  cnts_g <- rsdt_cdp(nrow(dat), 100, dat$stimulus, 0.8, 1.0, thresholds = thr,
                     rcrit = 0.5, kcrit = -0.3, n_new = 3)
  expect_true(all(c("guess4", "guess5", "guess6") %in% colnames(cnts_g)))
  # per-observation parameters are honoured
  df_v <- seq(0.2, 1.4, length.out = nrow(dat))
  cnts_v <- rsdt_cdp(nrow(dat), 10000, dat$stimulus, df_v, 1.0, thresholds = thr,
                     rcrit = 0.5, n_new = 3)
  expect_equal(dim(cnts_v), c(nrow(dat), 9))
})

test_that("dsdt_cdp computes a valid multinomial density and vectorizes", {
  thr <- .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious")
  cnt <- c(40, 20, 10, 15, 25, 35, 5, 20, 80)
  d <- dsdt_cdp(counts = cnt, stimulus = 1, dfam = 0.8, drec = 1.0,
                thresholds = thr, rcrit = 0.5, n_new = 3)
  expect_true(is.finite(d))
  expect_true(d > 0)
  dl <- dsdt_cdp(counts = cnt, stimulus = 1, dfam = 0.8, drec = 1.0,
                 thresholds = thr, rcrit = 0.5, n_new = 3, log = TRUE)
  expect_equal(log(d), dl, tolerance = 1e-10)
  # matches dmultinom against the category probabilities
  probs <- .sdt_cdp_category_probs(thr, 0.8, 1.0, 0, 0.5, NULL, 1, 3, 3,
                                   "normal", 0)
  expect_equal(dl, dmultinom(cnt, prob = probs, log = TRUE), tolerance = 1e-10)
  # rows vectorize with recycled parameters
  cnts <- rbind(cnt, c(60, 15, 5, 10, 20, 30, 10, 25, 75))
  dv <- dsdt_cdp(cnts, c(0L, 1L), 0.8, 1.0, thr, 0.5, n_new = 3, log = TRUE)
  expect_length(dv, 2)
  expect_equal(dv[2],
               dsdt_cdp(cnts[2, ], 1L, 0.8, 1.0, thr, 0.5, n_new = 3,
                        log = TRUE))
  # column count is validated
  expect_error(dsdt_cdp(cnt[-1], 1, 0.8, 1.0, thr, 0.5, n_new = 3),
               "columns")
})

test_that("aggregate_sdt_cdp_data pivots long data to the wide count columns", {
  long <- expand.grid(id = 1:4, stimulus = c(0L, 1L),
                      judgment = c("new", "know", "remember"),
                      confidence = 1:6)
  long <- long[(long$judgment == "new") == (long$confidence <= 3), ]
  long$count <- sample(5:20, nrow(long), replace = TRUE)
  wide <- aggregate_sdt_cdp_data(long, "judgment", "confidence", "count")
  expect_equal(nrow(wide), 8)                # 4 subjects x 2 stimulus classes
  expect_true(all(c("new1", "know4", "remember6") %in% colnames(wide)))
  expect_equal(sum(wide[grep("new|know|remember", colnames(wide))]),
               sum(long$count))
  # a response prefix is applied to the count columns
  wide_p <- aggregate_sdt_cdp_data(long, "judgment", "confidence", "count",
                                   response = "cdp")
  expect_true(all(c("cdpnew1", "cdpknow4") %in% colnames(wide_p)))
})

test_that("aggregate_sdt_cdp_data aggregates trial-level data without counts", {
  long <- expand.grid(id = 1:3, stimulus = c(0L, 1L),
                      judgment = c("new", "know", "remember"),
                      confidence = 1:6)
  long <- long[(long$judgment == "new") == (long$confidence <= 3), ]
  long$count <- sample(2:8, nrow(long), replace = TRUE)
  trials <- long[rep(seq_len(nrow(long)), long$count),
                 c("id", "stimulus", "judgment", "confidence")]
  wide <- aggregate_sdt_cdp_data(trials, "judgment", "confidence")
  expect_equal(nrow(wide), 6)
  expect_equal(sum(wide[grep("new|know|remember", colnames(wide))]),
               sum(long$count))
})

test_that("aggregate_sdt_cdp_data emits guess columns and validates input", {
  base <- data.frame(
    id = 1, stimulus = 1,
    judgment = c("new", "know", "remember", "guess"),
    confidence = c(1, 2, 2, 2), count = c(10, 5, 7, 3)
  )
  wide <- aggregate_sdt_cdp_data(base, "judgment", "confidence", "count")
  expect_true("guess2" %in% colnames(wide))
  # bad judgment label
  bad_j <- base
  bad_j$judgment[1] <- "maybe"
  expect_error(aggregate_sdt_cdp_data(bad_j, "judgment", "confidence", "count"),
               "not allowed")
  # missing a required judgment type
  miss <- base[base$judgment != "remember", ]
  expect_error(aggregate_sdt_cdp_data(miss, "judgment", "confidence", "count"),
               "must contain")
  # confidence partition violated: an old judgment on a new level
  bad_c <- base
  bad_c$confidence[2] <- 1
  expect_error(aggregate_sdt_cdp_data(bad_c, "judgment", "confidence", "count"),
               "must occupy levels")
})

test_that("check_data.sdt_cdp binds the count columns into a multinomial matrix", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  dat <- sim_cdp_data(n_trials = 200)
  cd <- check_data(m, dat, bmf(dfam ~ 1))
  expect_false(attr(cd, "has_guess"))
  expect_equal(ncol(cd$Y), 9)
  expect_equal(nrow(cd), 8)                  # 4 subjects x 2 stimulus classes
  expect_true(all(cd$nTrials == 200))
})

test_that("check_data.sdt_cdp handles asymmetric scales and guess columns", {
  m_a <- sdt_cdp(stimulus = "stimulus", n_new = 1, n_old = 5)
  thr_a <- .cdp_make_thresholds(0.2, -0.2, 1, 5, "parsimonious")
  dat_a <- sim_cdp_data(n_subjects = 3, n_trials = 300, n_new = 1, n_old = 5,
                        thresholds = thr_a, sigmar = log(1.5))
  cd_a <- check_data(m_a, dat_a, bmf(dfam ~ 1))
  expect_equal(ncol(cd_a$Y), 11)

  m_g <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  dat_g <- sim_cdp_data(n_subjects = 3, n_trials = 200, kcrit = 0.0)
  cd_g <- check_data(m_g, dat_g, bmf(dfam ~ 1))
  expect_true(attr(cd_g, "has_guess"))
  expect_equal(ncol(cd_g$Y), 12)
})

test_that("check_data.sdt_cdp validates columns, stimulus, and counts", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  dat <- sim_cdp_data()
  # missing response columns point at the aggregation helper
  expect_error(check_data(m, dat[-which(colnames(dat) == "new1")],
                          bmf(dfam ~ 1)),
               "aggregate_sdt_cdp_data")
  # a declared scale that mismatches the columns is a missing-column error
  m_big <- sdt_cdp(stimulus = "stimulus", n_new = 2, n_old = 4)
  expect_error(check_data(m_big, dat, bmf(dfam ~ 1)), "missing in the data")
  # partial guess columns
  dat_pg <- dat
  dat_pg$guess4 <- 1L
  expect_error(check_data(m, dat_pg, bmf(dfam ~ 1)), "all of them or none")
  # bad stimulus coding
  bad_s <- dat
  bad_s$stimulus <- 2
  expect_error(check_data(m, bad_s, bmf(dfam ~ 1)), "0 .* and 1")
  # negative counts
  bad_n <- dat
  bad_n$new1[1] <- -1L
  expect_error(check_data(m, bad_n, bmf(dfam ~ 1)), "non-negative")
})

test_that("sdt_cdp accepts a response prefix for the count columns", {
  dat <- sim_cdp_data()
  resp_cols <- grep("new|know|remember", colnames(dat))
  colnames(dat)[resp_cols] <- paste0("cdp", colnames(dat)[resp_cols])
  m <- sdt_cdp(response = "cdp", stimulus = "stimulus", n_new = 3, n_old = 3)
  cd <- check_data(m, dat, bmf(dfam ~ 1))
  expect_equal(ncol(cd$Y), 9)
  expect_silent(
    bmm(bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, spacing ~ 1, rcrit ~ 1),
        dat, m, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("bmf2bf builds a multinomial non-linear formula with one mu per category", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  dat <- sim_cdp_data(n_subjects = 3)
  cd <- check_data(m, dat, bmf(dfam ~ 1))
  m$other_vars$has_guess <- attr(cd, "has_guess")
  bf <- bmf2bf(m, bmf(dfam ~ 1))
  expect_s3_class(bf, "brmsformula")
  # 9 categories => 1 main + 8 mu non-linear formulas referencing sdt_cdp_logmu
  pforms <- vapply(bf$pforms, function(f) paste(deparse(f), collapse = " "),
                   character(1))
  expect_true(sum(grepl("sdt_cdp_logmu", pforms)) >= 8)
})

test_that("sdt_cdp produces multinomial stancode with the CDP functions", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  dat <- sim_cdp_data()
  sc <- stancode(bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, spacing ~ 1,
                     rcrit ~ 1), dat, m)
  expect_true(grepl("sdt_cdp_logmu", sc))
  expect_true(grepl("cdp_Phi2", sc))
  expect_true(grepl("owens_t", sc))
  expect_true(grepl("multinomial", sc))
})

test_that("freeing rho and sigmar adds them to the stancode", {
  m <- sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3)
  dat <- sim_cdp_data()
  sc <- stancode(bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, spacing ~ 1,
                     rcrit ~ 1, sigmar ~ 1, rho ~ 1), dat, m)
  expect_true(grepl("rho", sc))
  expect_true(grepl("sigmar", sc))
})
