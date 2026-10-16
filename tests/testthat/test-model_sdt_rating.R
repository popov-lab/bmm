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
  expect_error(sdt_rating("rating", "stimulus"), "more than 2")
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
  # the noise-distribution dispatch lives in sdt_dist_funs.stan; the rating
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
  funs <- paste(read_lines2(file.path(sc, "sdt_dist_funs.stan")),
                read_lines2(file.path(sc, "sdt_rating_funs.stan")), sep = "\n")
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
    "functions {\n", read_lines2(file.path(sc, "sdt_dist_funs.stan")), "\n",
    read_lines2(file.path(sc, "sdt_rating_funs.stan")), "\n",
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
