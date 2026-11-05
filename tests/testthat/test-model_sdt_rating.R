# Test Rating SDT model specification and integration

############################################################################# !
# MODEL CONSTRUCTOR TESTS                                                ####
############################################################################# !

test_that("sdt_rating takes the number of categories from the response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(model$other_vars$n_ratings, 4L)
  model <- sdt_rating(c("r1", "r2", "r3", "r4", "r5", "r6"), "stimulus")
  expect_equal(model$other_vars$n_ratings, 6L)
})

test_that("sdt_rating model rejects fewer than 3 response columns", {
  expect_error(sdt_rating(c("r1", "r2"), "stimulus"), "more than 2")
})

test_that("one rating column with n_ratings declares one row per trial", {
  trial <- sdt_rating("rating", "stimulus", n_ratings = 6)
  expect_true(trial$other_vars$trial_level)
  expect_identical(trial$other_vars$n_ratings, 6L)
  expect_match(trial$requirements, "one row per trial")
  expect_identical(.sdt_categories(trial), c("rating1", "rating2", "rating3", "rating4", "rating5", "rating6"))
  cells <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_false(cells$other_vars$trial_level)
  expect_identical(.sdt_categories(cells), c("r1", "r2", "r3", "r4"))
  # the thresholds follow K whichever layout gives it
  expect_identical(
    names(sdt_rating("rating", "stimulus", n_ratings = 5, threshold_type = "log_distance")$parameters),
    names(sdt_rating(paste0("r", 1:5), "stimulus", threshold_type = "log_distance")$parameters)
  )
})

test_that("sdt_rating names the layout it recognised when n_ratings does not fit it", {
  expect_error(sdt_rating("rating", "stimulus"),
               "response names one column, which declares one row per trial")
  expect_error(sdt_rating(c("r1", "r2", "r3"), "stimulus", n_ratings = 3),
               "response names 3 count columns, one per rating category")
  expect_error(sdt_rating("rating", "stimulus", n_ratings = 2), "single whole number above 2")
  expect_error(sdt_rating("rating", "stimulus", n_ratings = 4.5), "single whole number above 2")
  expect_error(sdt_rating("conf_rating", "stimulus", n_ratings = 4), "'conf_rating'")
  expect_error(sdt_rating("rating", "stimulus", n_ratings = 3, threshold_type = "log_ratio"),
               "at least 4 rating categories")
  expect_error(sdt_rating("rating", "stimulus", n_ratings = 5, version = "metad"),
               "The scale has 5 categories")
})

test_that("sdt_rating refuses column names brms cannot turn into parameters", {
  expect_error(sdt_rating(paste0("conf_", 1:4), "stimulus"),
               "must not contain '.' or '_'")
  expect_error(sdt_rating(c("r1", "r.2", "r3"), "stimulus"), "'r.2'")
})

test_that("sdt_rating model has log_distance threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_true(all(c("delta1", "delta2") %in% names(model$parameters)))
  expect_false("delta3" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
})

test_that("sdt_rating stores threshold_type correctly", {
  model_pa <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_equal(model_pa$other_vars$threshold_type, "parsimonious")

  model_ld <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                         threshold_type = "log_distance")
  expect_equal(model_ld$other_vars$threshold_type, "log_distance")
})

test_that("sdt_rating model has log_ratio threshold params", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_ratio")
  expect_true(all(c("delta1", "delta2") %in% names(model$parameters)))
  expect_true("criterion" %in% names(model$parameters))
  expect_false("spacing" %in% names(model$parameters))
  expect_false("delta3" %in% names(model$parameters))
  # the roles differ, so the labels must say which delta is the spread
  expect_match(model$parameters$delta2, "spread")
  expect_match(model$parameters$delta1, "ratio")
})

test_that("sdt_rating model has softmax threshold params", {
  model <- sdt_rating(paste0("r", 1:6), "stimulus", threshold_type = "softmax")
  expect_true("spacing" %in% names(model$parameters))
  # K = 6 -> K - 3 = 3 allocation deltas
  expect_true(all(c("delta1", "delta2", "delta3") %in% names(model$parameters)))
})

test_that("sdt_rating gives every parameter an sd default prior", {
  # Ro and Rn are logit-scale probabilities, so they join d at rate 1, as does
  # criterion, which carries half of d's spread; everything else is a
  # log-scale quantity or the midpoint-relative bias and takes rate 2
  rate1 <- c("d", "criterion", "Ro", "Rn")
  for (v in names(bmm:::.sdt_rating_variants)) {
    for (tt in c("parsimonious", "log_distance", "softmax")) {
      model <- sdt_rating(paste0("r", 1:6), "stimulus",
                          threshold_type = tt, version = v)
      sds <- vapply(model$default_priors, function(p) p$sd %||% NA_character_,
                    character(1))
      info <- paste(v, tt)
      expect_false(anyNA(sds), info = info)
      expect_true(all(sds[intersect(names(sds), rate1)] == "exponential(1)"), info = info)
      expect_true(all(sds[setdiff(names(sds), rate1)] == "exponential(2)"), info = info)
    }
  }
})

test_that("every threshold type places criterion the same way for odd and even K", {
  # even K: criterion is the middle threshold; odd K: the centre of the middle
  # category, with the two thresholds around it half an interval away. The R
  # builder must agree with sdt_place_thresholds_rating() in Stan (parity test
  # below), and the two closed forms (parsimonious, equidistant) with both.
  for (K in 3:8) {
    for (tt in c("parsimonious", "equidistant", "log_distance", "log_ratio", "softmax")) {
      if (tt == "log_ratio" && K < 4) next
      nd <- if (tt == "softmax") max(0L, K - 3L) else K - 2L
      deltas <- if (nd > 0) seq(0.1, 0.3, length.out = nd) else NULL
      thr <- .sdt_make_thresholds(0.3, K, tt, spacing = 0.2, deltas = deltas)
      expect_false(is.unsorted(thr), info = paste(tt, "K =", K))
      if (K %% 2 == 0) {
        expect_equal(thr[K / 2], 0.3, info = paste(tt, "K =", K))
      } else {
        g <- (K - 1) / 2
        expect_equal((thr[g] + thr[g + 1]) / 2, 0.3, info = paste(tt, "K =", K))
      }
    }
  }
  expect_equal(.sdt_make_thresholds(0, 5L, "equidistant", spacing = 0), c(-1.5, -0.5, 0.5, 1.5))
  expect_equal(.sdt_make_thresholds(0, 7L, "equidistant", spacing = 0), seq(-2.5, 2.5))
  expect_equal(.sdt_make_thresholds(0, 6L, "equidistant", spacing = 0), c(-2, -1, 0, 1, 2))
  # with equal intervals the anchored types coincide with equidistant
  for (K in 5:6) {
    expect_equal(.sdt_make_thresholds(0, K, "log_distance", deltas = rep(0, K - 2)),
                 .sdt_make_thresholds(0, K, "equidistant", spacing = 0))
    expect_equal(.sdt_make_thresholds(0, K, "log_ratio", deltas = rep(0, K - 2)),
                 .sdt_make_thresholds(0, K, "equidistant", spacing = 0))
    expect_equal(.sdt_make_thresholds(0, K, "softmax", spacing = 0, deltas = rep(0, K - 3)),
                 .sdt_make_thresholds(0, K, "equidistant", spacing = 0))
  }
})

test_that("log_ratio deltas act as spread and ratios, in interval order", {
  # even K = 6: delta3 is the spread (interval above the criterion threshold),
  # delta2 the ratio of the interval below to it, delta1 and delta4 ratios to
  # the first interval on their side
  gaps <- function(thr) diff(thr)
  expect_equal(gaps(.sdt_make_thresholds(0, 6L, "log_ratio", deltas = c(0, 0, 1, 0))),
               rep(exp(1), 4))
  expect_equal(gaps(.sdt_make_thresholds(0, 6L, "log_ratio", deltas = c(0, 1, 0, 0))),
               c(exp(1), exp(1), 1, 1))
  expect_equal(gaps(.sdt_make_thresholds(0, 6L, "log_ratio", deltas = c(1, 0, 0, 0))),
               c(exp(1), 1, 1, 1))
  # odd K = 5: delta2 is the middle category, delta1 and delta3 ratios to it
  expect_equal(gaps(.sdt_make_thresholds(0, 5L, "log_ratio", deltas = c(0, 1, 0))),
               rep(exp(1), 3))
  expect_equal(gaps(.sdt_make_thresholds(0, 5L, "log_ratio", deltas = c(0, 0, 1))),
               c(1, 1, exp(1)))
  # odd K = 7: delta3 is the middle category, delta2 and delta4 ratios to it,
  # delta1 and delta5 ratios to the first interval on their side
  expect_equal(gaps(.sdt_make_thresholds(0, 7L, "log_ratio", deltas = c(1, 0, 0, 0, 0))),
               c(exp(1), 1, 1, 1, 1))
  expect_equal(gaps(.sdt_make_thresholds(0, 7L, "log_ratio", deltas = c(0, 1, 0, 0, 0))),
               c(exp(1), exp(1), 1, 1, 1))
  expect_equal(gaps(.sdt_make_thresholds(0, 7L, "log_ratio", deltas = c(0, 0, 0, 1, 0))),
               c(1, 1, 1, exp(1), exp(1)))
})

test_that("deltas are named in interval order", {
  for (tt in c("log_distance", "log_ratio")) {
    model <- sdt_rating(paste0("r", 1:5), "stimulus", threshold_type = tt)
    expect_equal(grep("^delta", names(model$parameters), value = TRUE),
                 paste0("delta", 1:3), info = tt)
  }
})

############################################################################# !
# CHECK_DATA TESTS                                                       ####
############################################################################# !

sim_rating <- function(n_subjects, n_trials, d, criterion, n_ratings,
                       spacing = NULL, deltas = NULL,
                       threshold_type = "parsimonious",
                       sdratio = 1, dist = "normal") {
  thr <- bmm:::.sdt_make_thresholds(criterion, n_ratings, threshold_type,
                                    spacing, deltas)
  dat <- expand.grid(id = seq_len(n_subjects), stimulus = c(0L, 1L))
  cbind(dat, as.data.frame(rsdt_rating(nrow(dat), n_trials, dat$stimulus,
                                       d, thr, sdratio = sdratio,
                                       dist = dist)))
}

test_that("sdt_rating check_data validates response columns", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)

  valid_data <- data.frame(
    r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30), r4 = c(5, 55),
    stimulus = c(0L, 1L)
  )
  result <- check_data(model, valid_data, formula)
  # multinomial encoding: one row per observation with a count matrix Y
  expect_equal(nrow(result), nrow(valid_data))
  expect_true(all(c("Y", "nTrials", "stimulus") %in% colnames(result)))
  expect_equal(ncol(result$Y), 4)
  expect_equal(result$nTrials, c(50, 100))
})

test_that("sdt_rating check_data refuses an all-NA column as NA, not a type", {
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  dat <- data.frame(r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30),
                    r4 = NA, stimulus = c(0L, 1L))
  expect_error(check_data(model, dat, bmf(d ~ 1, criterion ~ 1, spacing ~ 1)),
               "must not contain NA counts")
})


############################################################################# !
# DISTRIBUTION FUNCTION TESTS                                            ####
############################################################################# !

test_that("rsdt_rating recycles vectorized parameters per observation", {
  withr::local_seed(445)
  counts <- rsdt_rating(6, c(50, 100), c(0L, 1L), d = rnorm(6, 1.5, 0.2),
                        thresholds = c(-0.5, 0, 0.5), sdratio = c(1, 1.3))
  expect_equal(dim(counts), c(6L, 4L))
  expect_true(all(rowSums(counts) == c(50, 100)))

  thr <- rbind(c(-0.5, 0, 0.5), c(-1, 0, 1))
  counts2 <- rsdt_rating(2, 100, c(0L, 1L), d = 1.5, thresholds = thr)
  expect_equal(dim(counts2), c(2L, 4L))
})

test_that("dsdt_rating matches dmultinom", {
  thr <- c(-0.5, 0, 0.5)
  probs <- bmm:::.sdt_category_probs(thr, 1.5, 1, 1, "normal")
  expect_equal(dsdt_rating(c(10, 20, 30, 40), 1, 1.5, thr),
               dmultinom(c(10, 20, 30, 40), prob = probs), tolerance = 1e-12)
})

test_that(".sdt_make_thresholds vectorized path matches per-draw evaluation", {
  cr <- c(-0.2, 0.1, 0.4)
  sp <- c(0.1, 0.3, -0.1)
  for (tt in c("parsimonious", "equidistant", "log_distance",
               "log_ratio", "softmax")) {
    K <- 5L
    nd <- if (tt == "softmax") K - 3L else K - 2L
    dl <- if (tt %in% c("log_distance", "log_ratio", "softmax")) {
      matrix(seq(-0.3, 0.3, length.out = 3 * nd), 3, nd)
    }
    vec <- bmm:::.sdt_make_thresholds(cr, K, tt, sp, dl)
    ref <- t(vapply(1:3, function(i) {
      bmm:::.sdt_make_thresholds(cr[i], K, tt, sp[i],
                                 if (is.null(dl)) NULL else dl[i, ])
    }, numeric(K - 1L)))
    expect_equal(vec, ref, tolerance = 1e-12, info = tt)
  }
})

test_that(".sdt_category_probs vectorized path matches per-row evaluation", {
  thr <- rbind(c(-0.5, 0, 0.5), c(-1, 0.2, 0.9))
  probs <- bmm:::.sdt_category_probs(thr, c(1.2, 0.8), c(1.3, 1), c(1L, 0L),
                                     "normal")
  ref <- t(vapply(1:2, function(i) {
    bmm:::.sdt_category_probs(thr[i, ], c(1.2, 0.8)[i], c(1.3, 1)[i],
                              c(1L, 0L)[i], "normal")
  }, numeric(4)))
  expect_equal(probs, ref, tolerance = 1e-12)
})

test_that("a bias is the criterion moved by half the separation", {
  # anchor 2 reads the middle threshold from the midpoint between the means,
  # anchor 1 from the noise mean; on the same thresholds both give one model
  d <- matrix(c(1, 1.5, 0.5, 2), nrow = 2)
  sdratio <- 0.3
  bias <- -0.2
  criterion <- bias + d * sqrt((1 + exp(sdratio)^2) / 2) / 2
  for (fun in list(sdt_rating_logmu, sdt_dpsdt_logmu, sdt_metad_logmu)) {
    extra <- intersect(names(formals(fun)), c("Ro", "Rn", "logmratio"))
    extra <- stats::setNames(rep(list(0.4), length(extra)), extra)
    for (k in 1:4) {
      args <- c(list(cat = k, K = 4L, dist = 1L, thresh = 1L, d = d,
                     spacing = 0.3, sdratio = sdratio, stimulus = c(0, 1)),
                extra)
      from_bias <- do.call(fun, c(args, list(anchor = 2L, criterion = bias)))
      from_crit <- do.call(fun, c(args, list(anchor = 1L, criterion = criterion)))
      expect_equal(from_bias, from_crit, tolerance = 1e-12)
    }
  }
})

test_that("sdt_rating_logmu is vectorized and preserves the draw shape", {
  d <- matrix(c(1, 1.5, 0.5, 2), nrow = 2)
  out <- sdt_rating_logmu(2L, 4L, 1L, 1L, 1L, d, criterion = 0.1,
                          spacing = 0.3, sdratio = 0.2, stimulus = c(0, 1))
  expect_equal(dim(out), dim(d))

  thr <- bmm:::.sdt_make_thresholds(0.1, 4L, "parsimonious", 0.3)
  stim <- rep_len(c(0, 1), 4)
  ref <- vapply(seq_along(d), function(j) {
    log(bmm:::.sdt_category_probs(thr, as.vector(d)[j], exp(0.2),
                                  stim[j], "normal")[2])
  }, numeric(1))
  expect_equal(as.vector(out), ref, tolerance = 1e-12)
})

test_that("sdt_rating_logmu passes each draw's deltas to the delta threshold builders", {
  # log_distance is threshold builder 3; the deltas differ by draw, so a
  # delta read from the wrong element moves a threshold
  d <- c(1, 1.5, 0.5)
  delta1 <- c(-0.2, 0.1, 0.3)
  delta2 <- c(0.4, -0.1, 0)
  out <- sdt_rating_logmu(3L, 4L, 1L, 3L, 1L, d, criterion = 0.1, spacing = 0,
                          sdratio = 0, stimulus = 1, delta1, delta2)
  ref <- vapply(seq_along(d), function(j) {
    thr <- bmm:::.sdt_make_thresholds(0.1, 4L, "log_distance", deltas = c(delta1[j], delta2[j]))
    log(bmm:::.sdt_category_probs(thr, d[j], 1, 1, "normal")[3])
  }, numeric(1))
  expect_equal(out, ref, tolerance = 1e-12)
})

test_that("log_ratio requires at least 4 rating categories", {
  expect_error(sdt_rating(c("r1", "r2", "r3"), "stimulus",
                          threshold_type = "log_ratio"),
               "at least 4 rating categories")
})

test_that("K = 3 softmax reduces to a single interval centred on criterion", {
  thr <- bmm:::.sdt_make_thresholds(0.1, 3L, "softmax", spacing = 0.2)
  expect_equal(thr, 0.1 + c(-0.5, 0.5) * exp(0.2), tolerance = 1e-12)
})


############################################################################# !
# FORMULA CONSTRUCTION TESTS                                              ####
############################################################################# !

test_that("bmf2bf maps every rating category to its own logit", {
  # names that do not sort in category order, so a mapping by name or by sorted
  # position would show up as a shifted category index
  cats <- c("surenew", "new", "old", "sureold")
  bf <- bmf2bf(sdt_rating(cats, "stimulus"), bmf(d ~ 1))
  expect_match(deparse(bf$formula, width.cutoff = 500),
               "sdt_rating_logmu(1, 4,", fixed = TRUE)
  for (k in 2:4) {
    pform <- deparse(bf$pforms[[paste0("mu", cats[k])]], width.cutoff = 500)
    expect_match(pform, paste0("sdt_rating_logmu(", k, ", 4,"), fixed = TRUE)
  }
})

test_that("a link on d or criterion is inverted inside the rating formula", {
  # the multinomial family has no link of its own for d and criterion, so a
  # link the user sets reaches the kernel only through the formula
  bf <- bmf2bf(sdt_rating(paste0("r", 1:4), "stimulus",
                          links = list(d = "log", criterion = "softplus")),
               bmf(d ~ 1))
  for (f in c(list(bf$formula), bf$pforms[paste0("mur", 2:4)])) {
    txt <- paste(deparse(f, width.cutoff = 500), collapse = "")
    expect_match(txt, "exp(d)", fixed = TRUE)
    expect_match(txt, "log1p_exp(criterion)", fixed = TRUE)
  }
  bf <- bmf2bf(sdt_rating(paste0("r", 1:4), "stimulus"), bmf(d ~ 1))
  expect_false(grepl("exp(d)", deparse(bf$formula, width.cutoff = 500), fixed = TRUE))
})

# the second category's call in the generated program: K, the noise
# distribution (1, normal), the threshold builder (1 parsimonious, 3
# log_distance, 4 log_ratio, 5 softmax) and the anchor, then the spacing (0
# for the delta builders) and the deltas the builder reads
test_that("the generated rating call passes K, the threshold builder and its deltas", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  counts <- function(K) {
    cells <- expand.grid(id = 1:2, stimulus = c(0L, 1L))
    cbind(cells, as.data.frame(matrix(10L, nrow(cells), K, dimnames = list(NULL, paste0("r", 1:K)))))
  }
  base <- "nlp_d[n] , nlp_criterion[n] ,"
  cases <- list(
    list(3, "parsimonious", bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
         "(2 , 3 , 1 , 1 , 1 ,", "nlp_spacing[n] , nlp_sdratio[n] , C_mur2_1[n]))"),
    list(3, "softmax", bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
         "(2 , 3 , 1 , 5 , 1 ,", "nlp_spacing[n] , nlp_sdratio[n] , C_mur2_1[n]))"),
    list(4, "log_distance", bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1),
         "(2 , 4 , 1 , 3 , 1 ,", "0 , nlp_sdratio[n] , C_mur2_1[n] , nlp_delta1[n] , nlp_delta2[n]))"),
    list(4, "log_ratio", bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1),
         "(2 , 4 , 1 , 4 , 1 ,", "0 , nlp_sdratio[n] , C_mur2_1[n] , nlp_delta1[n] , nlp_delta2[n]))"),
    list(6, "softmax", bmf(d ~ 1, criterion ~ 1, spacing ~ 1, delta1 ~ 1, delta2 ~ 1, delta3 ~ 1),
         "(2 , 6 , 1 , 5 , 1 ,",
         "nlp_spacing[n] , nlp_sdratio[n] , C_mur2_1[n] , nlp_delta1[n] , nlp_delta2[n] , nlp_delta3[n]))"),
    list(6, "parsimonious", bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
         "(2 , 6 , 1 , 1 , 1 ,", "nlp_spacing[n] , nlp_sdratio[n] , C_mur2_1[n]))")
  )
  for (case in cases) {
    model <- sdt_rating(paste0("r", seq_len(case[[1]])), "stimulus", threshold_type = case[[2]])
    code <- stancode(case[[3]], data = counts(case[[1]]), model = model)
    expected <- paste0("sdt_rating_logmu", case[[4]], " ", base, " ", case[[5]])
    expect_match(gsub("\\s+", " ", code), gsub("\\s+", " ", expected), fixed = TRUE,
                 label = paste(case[[1]], case[[2]]))
  }
})

test_that("sdt_rating stancode loads the shared distribution dispatch", {
  # the noise-distribution dispatch lives in sdt_dist_functions.stan; the rating
  # likelihood relies on the log-scale dispatchers from that shared chunk.
  dat <- sim_rating(3, 50, d = 1.5, criterion = 0, n_ratings = 4,
                    spacing = 0.5)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  code <- stancode(formula, data = dat, model = model)
  expect_true(grepl("sdt_log_cumprob", code, fixed = TRUE))
  expect_true(grepl("sdt_log_one_minus_cumprob", code, fixed = TRUE))
  expect_true(grepl("log_diff_exp", code, fixed = TRUE))
})

# The Stan kernels run through a fixed_param generated-quantities program, as
# the ranking parity tests do. sig_figs = 17 round-trips a double. Every row
# builds K thresholds with its own builder and returns one category's log
# probability, so both the builders and the kernel meet their R counterparts.
rating_stan_logmu <- function(grid, deltas, K) {
  sc <- system.file("stan_chunks", package = "bmm")
  funs <- paste(read_lines2(file.path(sc, "sdt_dist_functions.stan")),
                read_lines2(file.path(sc, "sdt_rating_functions.stan")), sep = "\n")
  program <- paste0(
    "functions {\n", funs, "\n}\n",
    "data {\n  int N; int K; array[N] int cat; array[N] int thresh;\n",
    "  vector[N] criterion; vector[N] spacing; array[N, K - 2] real deltas;\n",
    "  vector[N] d; vector[N] sdratio; vector[N] stimulus; array[N] int dist;\n}\n",
    "generated quantities {\n  vector[N] lp;\n  for (i in 1:N) {\n",
    "    vector[K - 1] thr = sdt_make_thresholds_rating(criterion[i], spacing[i],\n",
    "      deltas[i], K, thresh[i]);\n",
    "    lp[i] = sdt_rating_logmu_cat(cat[i], thr, d[i], sdratio[i], stimulus[i], dist[i]);\n",
    "  }\n}\n"
  )
  data <- c(list(N = nrow(grid), K = K, deltas = deltas),
            as.list(grid[c("cat", "thresh", "criterion", "spacing", "d",
                           "sdratio", "stimulus", "dist")]))
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = data, fixed_param = TRUE, chains = 1, iter_sampling = 1,
    iter_warmup = 0, refresh = 0, show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#",
                         check.names = FALSE)
  as.numeric(csv[1, paste0("lp.", seq_len(nrow(grid)))])
}

test_that("the Stan rating kernel matches its R counterpart in both tails", {
  skip_unless_stan_tests()
  skip_on_cran()

  # criterion = 12 puts every threshold far into the upper tail, where the log
  # cdf of both bounds of an interval rounds to 0; -12 does the same below.
  # K = 5 and K = 6 exercise the odd and even placement of the criterion; K = 7
  # is the smallest odd K with log_ratio intervals beyond the first on each side.
  for (K in c(5L, 6L, 7L)) {
    grid <- expand.grid(cat = seq_len(K), thresh = seq_along(bmm:::.sdt_threshold_types),
                        criterion = c(-12, -0.5, 0.3, 12), spacing = c(-0.5, 0.1),
                        d = c(0, 2), sdratio = c(-0.6, 0, 0.9),
                        stimulus = c(0, 1), dist = seq_along(bmm:::.sdt_dists))
    deltas <- matrix(seq(-0.6, 0.6, length.out = K - 2L), nrow(grid), K - 2L, byrow = TRUE)
    stan <- rating_stan_logmu(grid, deltas, K)

    r <- vapply(seq_len(nrow(grid)), function(i) {
      g <- grid[i, ]
      type <- bmm:::.sdt_threshold_types[g$thresh]
      n_deltas <- if (type == "softmax") K - 3L else K - 2L
      thr <- bmm:::.sdt_make_thresholds(
        g$criterion, K, type, g$spacing,
        if (type %in% c("log_distance", "log_ratio", "softmax")) deltas[1, seq_len(n_deltas)]
      )
      bmm:::.sdt_category_log_probs(thr, g$d, exp(g$sdratio), g$stimulus,
                                    names(bmm:::.sdt_dists)[g$dist])[g$cat]
    }, numeric(1))

    expect_true(all(is.finite(stan)), info = paste("K =", K))
    expect_equal(stan, r, tolerance = 1e-10, info = paste("K =", K))
  }
})

test_that("an upper-tail rating category keeps its mass", {
  # under gumbel_min the log cdf rounds to 0 from eta ~ 6.6, so an interval
  # taken as the difference of two log cdfs came out -Inf. Its survival
  # function exp(-exp(x)) gives the interval's mass in closed form.
  lo <- c(4, 7.5, 12)
  hi <- lo + 0.5
  expected <- -exp(lo) + log1m_exp(exp(lo) - exp(hi))
  got <- vapply(seq_along(lo), function(i) {
    bmm:::.sdt_category_log_probs(c(lo[i] - 1, lo[i], hi[i]), 0, 1, 0,
                                  "gumbel_min")[3]
  }, numeric(1))
  expect_equal(got, expected, tolerance = 1e-12)
})


############################################################################# !
# UV-SDT TESTS (sdratio overridable fixed parameter)                     ####
############################################################################# !

test_that("an estimated sdratio gets a prior and a fixed one is the constant 0", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  counts <- expand.grid(id = 1:2, stimulus = c(0L, 1L))
  counts[paste0("r", 1:4)] <- 10L
  for (type in c("parsimonious", "log_distance")) {
    model <- sdt_rating(paste0("r", 1:4), "stimulus", threshold_type = type)
    thresholds <- if (type == "parsimonious") bmf(spacing ~ 1) else bmf(delta1 ~ 1, delta2 ~ 1)
    estimated <- stancode(bmf(d ~ 1, criterion ~ 1, sdratio ~ 1) + thresholds, counts, model)
    fixed <- stancode(bmf(d ~ 1, criterion ~ 1) + thresholds, counts, model)
    expect_match(estimated, "normal_lpdf(b_sdratio[1] | 0, 0.3)", fixed = TRUE, label = type)
    expect_no_match(estimated, "b_sdratio[1] = 0;", fixed = TRUE, label = type)
    expect_match(fixed, "b_sdratio[1] = 0;", fixed = TRUE, label = type)
  }
})


############################################################################# !
# PIPELINE INTEGRATION TESTS (mock backend)                              ####
############################################################################# !

test_that("sdt_rating integrates with the bmm pipeline via mock backend", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  dat <- sim_rating(6, 80, d = 1.2, criterion = 0, n_ratings = 4,
                    spacing = 0.4)
  for (di in c("normal", "logistic", "gumbel_min", "gumbel_max")) {
    model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", dist = di)
    expect_silent(
      bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
          dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
    )
  }
})

test_that("sdt_rating UV-SDT integrates with the bmm pipeline via mock backend", {
  dat <- sim_rating(6, 80, d = 1.2, criterion = 0, n_ratings = 4,
                    spacing = 0.4, sdratio = 1.3)
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus")
  expect_silent(
    bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1, sdratio ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("sdt_rating log_distance integrates with the bmm pipeline via mock backend", {
  # log_distance has variable-arity delta parameters, so this exercises the
  # generated multinomial formula and Stan function for that case.
  dat <- sim_rating(5, 60, d = 1.5, criterion = 0, n_ratings = 4,
                    deltas = c(0.5, 0.5), threshold_type = "log_distance")
  model <- sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus",
                      threshold_type = "log_distance")
  expect_silent(
    bmm(bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1),
        dat, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
})

test_that("the generated Stan logmu function reads either anchor as its R companion does", {
  skip_unless_stan_tests()
  skip_on_cran()

  model <- sdt_rating(paste0("r", 1:4), "stimulus")
  sc <- system.file("stan_chunks", package = "bmm")
  program <- paste0(
    "functions {\n", read_lines2(file.path(sc, "sdt_dist_functions.stan")), "\n",
    read_lines2(file.path(sc, "sdt_rating_functions.stan")), "\n",
    bmm:::.sdt_rating_logmu_stan(model), "\n}\n",
    "data {\n  int N; array[N] int cat; array[N] int anchor; array[N] int dist;\n",
    "  vector[N] d; vector[N] crit; vector[N] sdratio; vector[N] stimulus;\n}\n",
    "generated quantities {\n  vector[N] lp;\n  for (i in 1:N) {\n",
    "    lp[i] = sdt_rating_logmu(cat[i], 4, dist[i], 1, anchor[i], d[i], crit[i],\n",
    "      -0.4, sdratio[i], stimulus[i]);\n",
    "  }\n}\n"
  )
  grid <- expand.grid(cat = 1:4, anchor = 1:2, dist = seq_along(bmm:::.sdt_dists),
                      d = c(0.4, 2), crit = c(-0.5, 0.8), sdratio = c(-0.4, 0, 0.5),
                      stimulus = c(0, 1))
  fit <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(program))$sample(
    data = c(list(N = nrow(grid)), as.list(grid)), fixed_param = TRUE,
    chains = 1, iter_sampling = 1, iter_warmup = 0, refresh = 0,
    show_messages = FALSE, sig_figs = 17
  )
  csv <- utils::read.csv(fit$output_files()[1], comment.char = "#",
                         check.names = FALSE)
  r <- vapply(seq_len(nrow(grid)), function(i) with(grid[i, ], sdt_rating_logmu(
    cat, 4L, dist, 1L, anchor, d, crit, -0.4, sdratio, stimulus)), numeric(1))
  expect_equal(as.numeric(csv[1, paste0("lp.", seq_len(nrow(grid)))]), r,
               tolerance = 1e-10)
})

test_that("a rating formula with bias passes anchor 2 and the bias parameter", {
  model <- sdt_rating(paste0("r", 1:4), "stimulus")
  bf_bias <- bmm:::bmf2bf(model, bmf(d ~ 1, bias ~ 1, spacing ~ 1))
  expect_match(deparse1(bf_bias$formula), "sdt_rating_logmu(1, 4, 1, 1, 2, d, bias,",
               fixed = TRUE)
  bf_crit <- bmm:::bmf2bf(model, bmf(d ~ 1, criterion ~ 1, spacing ~ 1))
  expect_match(deparse1(bf_crit$formula), "sdt_rating_logmu(1, 4, 1, 1, 1, d, criterion,",
               fixed = TRUE)
})


############################################################################# !
# TRIAL-LEVEL DATA                                                       ####
############################################################################# !

# 2 participants x 2 stimulus types x 6 trials, ratings 1..4
rating_trials <- function() {
  data.frame(
    id = rep(1:2, each = 12), stimulus = rep(rep(0:1, each = 6), 2),
    rating = c(1, 2, 1, 3, 2, 1,  3, 4, 4, 2, 4, 3,
               2, 1, 1, 1, 4, 2,  4, 3, 4, 4, 2, 4)
  )
}

test_that("trial-level check_data() names the ratings outside 1..n_ratings", {
  model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  dat <- rating_trials()
  dat$rating[c(2, 5)] <- c(0, 2.5)
  expect_error(check_data(model, dat, bmf(d ~ 1, criterion ~ 1, spacing ~ 1)),
               "the rating column 'rating' must hold one of '1', '2', '3', '4' in every row. Found: '0', '2.5'")
  dat <- rating_trials()
  dat$rating[3] <- NA
  expect_error(check_data(model, dat, bmf(d ~ 1, criterion ~ 1, spacing ~ 1)),
               "'rating' must not contain NA")
})

# match() reads TRUE as category 1; the rank column goes through the same helper
test_that("trial-level check_data() refuses a logical rating or rank column", {
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  dat <- rating_trials()
  dat$rating <- TRUE
  expect_error(check_data(model, dat, formula),
               "The rating column 'rating' is logical (TRUE/FALSE), but it must hold whole numbers from 1 to 4.",
               fixed = TRUE)
  dat$rating <- rep(c(TRUE, FALSE), 12)
  expect_error(check_data(model, dat, formula), "'rating' is logical", fixed = TRUE)
  ranks <- data.frame(id = 1:4, rank = c(TRUE, FALSE, TRUE, TRUE))
  expect_error(check_data(sdt_ranking("rank", m = 3), ranks, bmf(d ~ 1)),
               "The rank column 'rank' is logical (TRUE/FALSE), but it must hold whole numbers from 1 to 3.",
               fixed = TRUE)
})

test_that("trial-level check_data() keeps the rating as the category number", {
  model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  checked <- check_data(model, rating_trials(), bmf(d ~ 1, criterion ~ 1, spacing ~ 1))
  expect_identical(checked$rating, as.integer(rating_trials()$rating))
  expect_false(any(c("Y", "nTrials") %in% names(checked)))
  expect_null(data_trials(model, rating_trials()))
})

test_that("rating trial rows fit through the trial family and cells keep the multinomial", {
  skip_on_cran()
  trials <- rating_trials()
  cells <- data.frame(id = rep(1:2, each = 2), stimulus = rep(0:1, 2),
                      rating1 = c(3, 0, 3, 0), rating2 = c(2, 1, 2, 1),
                      rating3 = c(1, 2, 0, 1), rating4 = c(0, 3, 1, 4))
  formula <- bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1)
  mock <- function(data, model) suppressMessages(
    bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE)
  )
  trial_code <- brms::stancode(mock(trials, sdt_rating("rating", "stimulus", n_ratings = 4)))
  expect_match(trial_code, paste0(
    "target += sdt_rating_trial_standard_parsimonious_criterion_lpmf(Y[n] | mu[n], d[n], ",
    "criterion[n], spacing[n], sdratio[n], vint1[n], 4, 1, 1, 1);"
  ), fixed = TRUE)
  expect_match(trial_code, paste0(
    "return sdt_rating_logmu(y, K, dist_type, thresh_type, anchor, d, criterion, ",
    "spacing, sdratio, stimulus);"
  ), fixed = TRUE)
  expect_no_match(trial_code, "multinomial")
  cell_code <- brms::stancode(mock(cells, sdt_rating(paste0("rating", 1:4), "stimulus")))
  expect_match(cell_code, "target += multinomial_logit2_lpmf(Y[n] | mu[n]);", fixed = TRUE)
  expect_match(cell_code, "murating4[n] = (sdt_rating_logmu(4 , 4 , 1 , 1 , 1 , nlp_d[n] ,", fixed = TRUE)
  expect_no_match(cell_code, "sdt_rating_trial_")
})

# one prep of a mock trial fit, with the dpars set to fixed values per draw
rating_trial_prep <- function(model, formula, data, dpars) {
  fit <- suppressMessages(bmm(formula, data, model, backend = "mock", mock_fit = 1, rename = FALSE))
  sdata <- brms::standata(fit)
  dpars <- lapply(dpars, function(v) matrix(v, length(v), sdata$N))
  structure(list(dpars = dpars, data = sdata, ndraws = nrow(dpars[[1]]), nobs = sdata$N,
                 family = fit$family),
            class = "brmsprep")
}

# equidistant thresholds at criterion 0.3 and spacing log(0.5): -0.2, 0.3, 0.8
test_that("the trial family's log_lik is the category log-probability of the multinomial route", {
  skip_on_cran()
  dat <- data.frame(stimulus = c(0, 1, 1, 0), rating = c(1, 3, 4, 2))
  prep <- rating_trial_prep(
    sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "equidistant"),
    bmf(d ~ 1, criterion ~ 1, spacing ~ 1), dat,
    list(d = c(1.2, 0.5), criterion = c(0.3, 0.3), spacing = log(c(0.5, 0.5)), sdratio = c(0, 0))
  )
  ll <- vapply(1:4, function(i) bmm:::log_lik_sdt_rating_trial(i, prep), numeric(2))
  expect_equal(ll[1, 2], log(pnorm(0.8 - 1.2) - pnorm(0.3 - 1.2)), tolerance = 1e-10)
  expected <- rbind(
    dsdt_rating(diag(4)[dat$rating, ], dat$stimulus, d = 1.2, thresholds = c(-0.2, 0.3, 0.8), log = TRUE),
    dsdt_rating(diag(4)[dat$rating, ], dat$stimulus, d = 0.5, thresholds = c(-0.2, 0.3, 0.8), log = TRUE)
  )
  expect_equal(ll, expected, tolerance = 1e-10)
  epred <- bmm:::posterior_epred_sdt_rating_trial(prep)
  expect_identical(dim(epred), c(2L, 4L, 4L))
  expect_equal(apply(epred, c(1, 2), sum), matrix(1, 2, 4), tolerance = 1e-12)
  expect_equal(ll, vapply(1:4, function(i) log(epred[, i, dat$rating[i]]), numeric(2)),
               tolerance = 1e-10)
  predicted <- vapply(1:4, function(i) bmm:::posterior_predict_sdt_rating_trial(i, prep), numeric(2))
  expect_true(all(predicted %in% 1:4))
})

# 4000 draws of one parameter point: the share of draws predicting each
# category is the category probability, within 0.04 (five Monte Carlo SEs at
# p = 0.5). A reversed or shifted category mapping misses it by over 0.3
test_that("the trial family's posterior_predict draws each category at its probability", {
  skip_on_cran()
  withr::local_seed(667)
  n <- 4000
  dat <- data.frame(stimulus = c(0, 1, 1, 0), rating = c(1, 3, 4, 2))
  prep <- rating_trial_prep(
    sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "equidistant"),
    bmf(d ~ 1, criterion ~ 1, spacing ~ 1), dat,
    list(d = rep(1.2, n), criterion = rep(0.3, n), spacing = rep(log(0.5), n), sdratio = rep(0, n))
  )
  predicted <- vapply(1:4, function(i) bmm:::posterior_predict_sdt_rating_trial(i, prep), numeric(n))
  share <- vapply(1:4, function(k) colMeans(predicted == k), numeric(4))
  probs <- bmm:::posterior_epred_sdt_rating_trial(prep)[1, , ]
  expect_equal(probs[2, ], c(pnorm(-1.4), pnorm(-0.9) - pnorm(-1.4), pnorm(-0.4) - pnorm(-0.9),
                             1 - pnorm(-0.4)), tolerance = 1e-10)
  expect_lt(max(abs(share - probs)), 0.04)
})

# brms can hand a dpar to the family's methods as draws x rows, as one value
# per draw or as one value; each is recycled over the rows
test_that("the trial family's methods take a dpar as one value per draw or one value", {
  skip_on_cran()
  dat <- data.frame(stimulus = c(0, 1, 1, 0), rating = c(1, 3, 4, 2))
  prep <- rating_trial_prep(
    sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "equidistant"),
    bmf(d ~ 1, criterion ~ 1, spacing ~ 1), dat,
    list(d = c(1.2, 0.5), criterion = c(0.3, 0.3), spacing = log(c(0.5, 0.5)), sdratio = c(0, 0))
  )
  short <- prep
  short$dpars$d <- c(1.2, 0.5)
  short$dpars$sdratio <- 0
  expect_equal(bmm:::posterior_epred_sdt_rating_trial(short), bmm:::posterior_epred_sdt_rating_trial(prep))
  expect_equal(vapply(1:4, function(i) bmm:::log_lik_sdt_rating_trial(i, short), numeric(2)),
               vapply(1:4, function(i) bmm:::log_lik_sdt_rating_trial(i, prep), numeric(2)))
})

# a bias is the criterion moved by half the separation: bias -0.3 at d 1.2 is
# criterion 0.3, so both anchors give the same probabilities
test_that("every version builds on trial rows with either anchor, and the bias anchor converts", {
  skip_on_cran()
  dat <- rating_trials()
  for (version in c("standard", "dpsdt", "metad")) {
    extra <- switch(version, standard = list(), dpsdt = list(Ro ~ 1), metad = list(logmratio ~ 1))
    for (anchor in c("criterion", "bias")) {
      formula <- do.call(bmf, c(list(d ~ 1, stats::as.formula(paste(anchor, "~ 1")), spacing ~ 1), extra))
      fit <- suppressMessages(bmm(formula, dat, sdt_rating("rating", "stimulus", n_ratings = 4,
                                                               version = version),
                                  backend = "mock", mock_fit = 1, rename = FALSE))
      name <- paste("sdt_rating_trial", version, "parsimonious", anchor, sep = "_")
      expect_identical(fit$family$name, name)
      expect_match(brms::stancode(fit), paste0("target += ", name, "_lpmf(Y[n] |"), fixed = TRUE,
                   info = name)
    }
  }
  model <- sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "equidistant")
  dpars <- list(d = 1.2, spacing = log(0.5), sdratio = 0)
  by_criterion <- rating_trial_prep(model, bmf(d ~ 1, criterion ~ 1, spacing ~ 1), dat,
                                    c(dpars, criterion = 0.3))
  by_bias <- rating_trial_prep(model, bmf(d ~ 1, bias ~ 1, spacing ~ 1), dat, c(dpars, bias = -0.3))
  expect_equal(bmm:::posterior_epred_sdt_rating_trial(by_bias),
               bmm:::posterior_epred_sdt_rating_trial(by_criterion), tolerance = 1e-12)
})

# the multinomial route substitutes inv_link() into its nlf call; the trial
# family hands the same links to brms, which applies them to the dpars
test_that("a link on d or criterion reaches the trial family's dpars", {
  skip_on_cran()
  fit <- suppressMessages(bmm(
    bmf(d ~ 1, criterion ~ 1, spacing ~ 1), rating_trials(),
    sdt_rating("rating", "stimulus", n_ratings = 4, links = list(d = "softplus", criterion = "probit")),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  expect_identical(fit$family$link_d, "softplus")
  code <- brms::stancode(fit)
  expect_match(code, "d = log1p_exp(d);", fixed = TRUE)
  expect_match(code, "criterion = Phi(criterion);", fixed = TRUE)
})

test_that("threshold deltas stay nlpars on trial rows and reach the family by letter", {
  skip_on_cran()
  fit <- suppressMessages(bmm(
    bmf(d ~ 1, criterion ~ 1, delta1 ~ 1 + (1 | id), delta2 ~ 1),
    rating_trials(), sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "log_distance"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  bterms <- brms::brmsterms(fit$formula)
  expect_setequal(names(bterms$nlpars), c("delta1", "delta2"))
  expect_setequal(names(bterms$dpars), c("mu", "d", "criterion", "sdratio", "deltaa", "deltab"))
  code <- brms::stancode(fit)
  expect_match(code, "deltaa[n] = (nlp_delta1[n]);", fixed = TRUE)
  expect_match(code, "sdratio, stimulus, deltaa, deltab);", fixed = TRUE)
  expect_identical(bmm:::.get_parameter_info(fit, "d")$type, "dpar")
  expect_identical(bmm:::.get_parameter_info(fit, "delta1")$type, "nlpar")
})

# brms centres a dpar's intercept when the formula has predictors; the trial
# family's dpars are uncentred, as the aggregated fit's nlpars are, so the
# defaults and a prior written for the aggregated fit land on the same rows
test_that("a prior written for the aggregated fit reaches the trial fit unchanged", {
  skip_on_cran()
  dat <- rating_trials()
  dat$cond <- factor(rep(c("a", "b"), 12))
  formula <- bmf(d ~ 1 + cond + (1 | id), criterion ~ 1, spacing ~ 1)
  user <- brms::set_prior("normal(2, 0.3)", nlpar = "d") +
    brms::set_prior("normal(0, 9)", nlpar = "d", coef = "condb") +
    brms::set_prior("exponential(4)", nlpar = "d", class = "sd")
  fit <- suppressMessages(bmm(formula, dat, sdt_rating("rating", "stimulus", n_ratings = 4),
                              prior = user, backend = "mock", mock_fit = 1, rename = FALSE))
  rows <- as.data.frame(fit$prior)[nzchar(fit$prior$prior) & fit$prior$dpar == "d",
                                   c("prior", "class", "coef", "lb")]
  expect_identical(sort(rows$prior), c("exponential(4)", "normal(0, 9)", "normal(2, 0.3)"))
  expect_identical(rows$lb[rows$class == "sd"], "0")
  expect_false("Intercept" %in% fit$prior$class)
  code <- brms::stancode(fit)
  expect_match(code, "lprior += normal_lpdf(b_d[1] | 2, 0.3);", fixed = TRUE)
  expect_match(code, "lprior += normal_lpdf(b_d[2] | 0, 9);", fixed = TRUE)
  expect_match(code, "lprior += normal_lpdf(b_criterion[1] | 0.5, 1.5);", fixed = TRUE)
  expect_match(code, "lprior += exponential_lpdf(sd_1 | 4);", fixed = TRUE)
  expect_match(code, "b_sdratio[1] = 0;", fixed = TRUE)
  expect_no_match(code, "Intercept_d")
})

# one-hot rows carry no multinomial coefficient: 6! / (2! 1! 0! 3!) = 60
# orderings of the cell's 6 trials
test_that("rating trial rows sum to the cell log-likelihood less the log multinomial coefficient", {
  ratings <- c(4, 1, 4, 2, 1, 4)
  one_hot <- diag(4)[ratings, ]
  thresholds <- c(-0.4, 0.3, 1.1)
  for (stim in c(0, 1)) {
    expect_equal(
      sum(dsdt_rating(one_hot, stim, d = 1.3, thresholds = thresholds, log = TRUE)),
      dsdt_rating(c(2, 1, 0, 3), stim, d = 1.3, thresholds = thresholds, log = TRUE) - log(60),
      tolerance = 1e-10
    )
    expect_equal(
      sum(dsdt_rating(one_hot, stim, d = 0.6, thresholds = thresholds, sdratio = 1.4,
                      dist = "gumbel_min", log = TRUE)),
      dsdt_rating(c(2, 1, 0, 3), stim, d = 0.6, thresholds = thresholds, sdratio = 1.4,
                  dist = "gumbel_min", log = TRUE) - log(60),
      tolerance = 1e-10
    )
  }
})

# the helpers that read counts take the trial family's rating through one
# bridge, .observed_counts(), and expand its predicted categories to one-hot
test_that("the count-reading helpers read the trial family's ratings", {
  skip_on_cran()
  trials <- rating_trials()
  model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  fit <- suppressMessages(bmm(bmf(d ~ 1, criterion ~ 1, spacing ~ 1), trials, model,
                              backend = "mock", mock_fit = 1, rename = FALSE))
  expect_true(bmm:::.is_sdt_rating_trial_fit(fit))
  one_hot <- diag(4)[trials$rating, ]
  observed <- bmm:::.observed_counts(fit$bmm$model, fit$data)
  expect_equal(unname(observed), one_hot)
  expect_identical(colnames(observed), paste0("rating", 1:4))
  predicted <- matrix(c(1L, 4L), 2, 24)
  epred <- array(0.25, c(2, 24, 4))
  local_mocked_bindings(posterior_predict = function(...) predicted,
                        posterior_epred = function(...) epred, .package = "brms")

  # ppp(): observed and replicated one-hot, expected passed through
  counts <- bmm:::.ppp_counts.default(fit$bmm$model, fit, NULL)
  expect_equal(unname(counts$observed), one_hot)
  expect_identical(dim(counts$replicated), c(2L, 24L, 4L))
  expect_equal(unname(counts$replicated[1, , ]), matrix(rep(c(1, 0, 0, 0), each = 24), 24))
  expect_equal(unname(counts$replicated[2, , ]), matrix(rep(c(0, 0, 0, 1), each = 24), 24))
  expect_identical(counts$expected, epred)

  # roc_observed(): the cumulative rates of the cells the trials aggregate to
  roc <- bmm:::.roc_observed_rating(fit, fit$bmm$model)
  noise <- tabulate(trials$rating[trials$stimulus == 0], 4) / 12
  signal <- tabulate(trials$rating[trials$stimulus == 1], 4) / 12
  expect_equal(roc$FA, c(1, 1 - cumsum(noise)[1:3], 0))
  expect_equal(roc$Hit, c(1, 1 - cumsum(signal)[1:3], 0))

  # pp_check(): the multinomial path, with the predicted ratings as one-hot rows
  plot_data <- bmm:::.pp_check_multinomial(fit, ndraws = 2)$data
  expect_equal(plot_data$proportion[order(plot_data$category)],
               as.vector(colSums(one_hot)) / 24)
  expect_equal(sort(unique(plot_data$pred_median)), c(0, 0.5))
})

test_that("sdt post-processing asks a trial fit for dpars and its deltas for nlpars", {
  skip_on_cran()
  fit <- suppressMessages(bmm(
    bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1), rating_trials(),
    sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "log_distance"),
    backend = "mock", mock_fit = 1, rename = FALSE
  ))
  asked <- list()
  local_mocked_bindings(posterior_linpred = function(object, dpar = NULL, nlpar = NULL, ...) {
    asked[[length(asked) + 1]] <<- c(dpar = dpar %||% "", nlpar = nlpar %||% "")
    matrix(0, 2, 1)
  }, .package = "brms")
  for (par in c("d", "criterion", "delta1")) bmm:::.sdt_par_draws(fit, par, data.frame())
  expect_identical(do.call(rbind, asked),
                   rbind(c(dpar = "d", nlpar = ""), c(dpar = "criterion", nlpar = ""),
                         c(dpar = "", nlpar = "delta1")))
})

test_that("revert_check_data() round-trips trial rows and a 1.5.0 one-hot frame", {
  model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  checked <- check_data(model, rating_trials(), formula)
  reverted <- revert_check_data(model, checked)
  expect_equal(reverted$rating, rating_trials()$rating)
  expect_identical(check_data(model, reverted, formula), checked)
  # a trial fit of 1.5.0 stored the one-hot Y the multinomial read
  one_hot <- checked[setdiff(names(checked), "rating")]
  one_hot$Y <- diag(4)[checked$rating, ]
  one_hot$nTrials <- 1L
  old <- revert_check_data(model, one_hot)
  expect_equal(old$rating, rating_trials()$rating)
  expect_false(any(c("Y", "nTrials") %in% names(old)))
})

# A trial fit made through the multinomial (1.5.0) stored a one-hot Y and
# nTrials; newdata in the user's layout get them built again, which the frame
# decides, while newdata for a trial-family fit go to brms as they are
test_that("user-layout newdata predict from a 1.5.0 one-hot trial fit and from the trial family", {
  skip_on_cran()
  trials <- rating_trials()
  one_hot <- trials[c("id", "stimulus")]
  one_hot[paste0("rating", 1:4)] <- as.data.frame(diag(4)[trials$rating, ])
  formula <- bmf(d ~ 1, criterion ~ 1, spacing ~ 1)
  trial_model <- sdt_rating("rating", "stimulus", n_ratings = 4)
  old <- suppressMessages(bmm(formula, one_hot, sdt_rating(paste0("rating", 1:4), "stimulus"),
                              backend = "mock", mock_fit = 1, rename = FALSE))
  old$bmm$model <- trial_model
  newdata <- trials[c(3, 8, 20), ]

  expect_true(bmm:::lacks_built_columns(old, newdata))
  checked <- bmm:::check_newdata(old, newdata)
  expect_equal(unname(checked$Y), unname(diag(4)[newdata$rating, ]))
  expect_identical(checked$nTrials, rep(1L, 3))
  expect_null(attr(checked, "fit_one_hot"))
  old_data <- brms::standata(old, newdata = newdata)
  expect_equal(unname(old_data$Y), unname(diag(4)[newdata$rating, ]))
  expect_equal(as.vector(old_data$trials), rep(1, 3))

  new <- suppressMessages(bmm(formula, trials, trial_model, backend = "mock", mock_fit = 1, rename = FALSE))
  expect_false(bmm:::lacks_built_columns(new, newdata))
  expect_equal(as.vector(brms::standata(new, newdata = newdata)$Y), newdata$rating)
})

# sdt_rating()'s trial family plots a parameter's own effects, as the
# aggregated fit does, not every predictor of the model
test_that("conditional_effects() asks brms for a trial-family parameter's own effects only", {
  skip_on_cran()
  dat <- rating_trials()
  dat$cond <- factor(rep(c("a", "b"), 12))
  fit <- suppressMessages(bmm(bmf(d ~ 1 + cond, criterion ~ 1, delta1 ~ 1, delta2 ~ 1), dat,
                              sdt_rating("rating", "stimulus", n_ratings = 4, threshold_type = "log_distance"),
                              backend = "mock", mock_fit = 1, rename = FALSE))
  asked <- list()
  local_mocked_bindings(.brms_conditional_effects = function(x, ...) {
    asked[[length(asked) + 1]] <<- list(...)
    structure(list(), names = character(0), class = "brms_conditional_effects")
  })
  for (par in c("d", "criterion", "delta1")) {
    ce <- bmm:::.ce_compute_and_transform(fit, par, bmm:::.get_parameter_info(fit, par), "native")
    expect_length(ce, 0)
  }
  expect_length(asked, 1)
  expect_identical(asked[[1]][c("effects", "dpar")], list(effects = "cond", dpar = "d"))
  bmm:::.ce_compute_and_transform(fit, "criterion", bmm:::.get_parameter_info(fit, "criterion"),
                                  "native", effects = "cond")
  expect_identical(asked[[2]][c("dpar", "effects")], list(dpar = "criterion", effects = "cond"))
})

# The trial family and the one-hot multinomial (aggregated rows of one trial
# each, the program 1.5.0 fits to trial rows) share every parameter and prior,
# so their total log densities agree at any point. log_prob through the
# executable: cmdstanr's init_model_methods() can clash with RcppParallel's TBB.
# About 28 s: it compiles both programs, which no other test builds
test_that("the trial family's log density equals the one-hot multinomial's", {
  skip_unless_stan_tests()
  withr::local_seed(11)
  dir <- withr::local_tempdir()
  trials <- data.frame(id = rep(1:4, each = 12), item = rep(1:12, 4),
                       cond = factor(rep(c("a", "b"), 24)), stimulus = rep(rep(0:1, each = 6), 4),
                       rating = sample.int(6, 48, replace = TRUE))
  one_hot <- trials
  one_hot[paste0("rating", 1:6)] <- as.data.frame(diag(6)[trials$rating, ])
  formula <- bmf(d ~ 1 + cond + (1 | id) + (1 | item), criterion ~ 1 + (1 | id), spacing ~ 1)
  fits <- list(
    trial = suppressMessages(bmm(formula, trials, sdt_rating("rating", "stimulus", n_ratings = 6),
                                 backend = "mock", mock_fit = 1, rename = FALSE)),
    one_hot = suppressMessages(bmm(formula, one_hot, sdt_rating(paste0("rating", 1:6), "stimulus"),
                                   backend = "mock", mock_fit = 1, rename = FALSE))
  )
  sdata <- brms::standata(fits$trial)
  values <- list(b_d = array(c(1.1, 0.4)), b_criterion = array(0.3), b_spacing = array(-0.5),
                 sd_1 = array(0.4), z_1 = matrix(stats::rnorm(4), 1),
                 sd_2 = array(0.2), z_2 = matrix(stats::rnorm(12), 1),
                 sd_3 = array(0.3), z_3 = matrix(stats::rnorm(4), 1))
  log_prob <- function(name) {
    mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(brms::stancode(fits[[name]]), dir),
                                   quiet = TRUE)
    files <- file.path(dir, paste0(name, c("_data.json", "_params.json", "_out.csv")))
    data <- lapply(unclass(brms::standata(fits[[name]])), function(x) {
      if (is.null(dim(x)) && length(x) > 1) as.array(x) else unclass(x)
    })
    cmdstanr::write_stan_json(data, files[1])
    cmdstanr::write_stan_json(values, files[2])
    status <- system2(mod$exe_file(), c("log_prob", paste0("constrained_params=", files[2]),
                                        "jacobian=1", "data", paste0("file=", files[1]),
                                        "output", paste0("file=", files[3]), "sig_figs=18"),
                      stdout = FALSE, stderr = FALSE)
    expect_identical(status, 0L)
    utils::read.csv(files[3], comment.char = "#")$lp[1]
  }
  expect_identical(c(sdata$K_d, sdata$N_1, sdata$N_2, sdata$N_3), c(2L, 4L, 12L, 4L))
  expect_equal(log_prob("trial"), log_prob("one_hot"), tolerance = 1e-8)
})
