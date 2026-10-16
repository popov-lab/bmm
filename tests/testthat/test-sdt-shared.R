# Test SDT shared infrastructure

############################################################################# !
# CONTRACTS EVERY SDT FORMAT SHARES                                      ####
############################################################################# !

# One row per format. ranking is two rows because dist = "normal" adds
# sdratio and takes every link away; the rating versions only add rows where
# they change the answer.
sdt_formats <- list(
  sdt_yn = list(
    make = function(...) sdt_yn("n_old", "stimulus", "n_trials", ...),
    class = c("bmmodel", "sdt", "sdt_yn"),
    params = c("d", "criterion", "bias", "sdratio"),
    fixed = list(sdratio = 0),
    links = list(d = "identity", criterion = "identity", bias = "identity",
                 sdratio = "log"),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    custom = list(d = "log"),
    settable = c("d", "criterion", "bias"),
    # sdratio is fixed at 0 and read on whatever link the model carries; only
    # log maps that 0 to the equal-variance ratio of 1
    refused = list(sdratio = c("identity", "softplus", "logit", "sqrt")),
    # links reach Stan through custom_family(), so sqrt is a valid link here
    unknown_link = "nonsense"
  ),
  sdt_rating = list(
    make = function(...) sdt_rating(c("r1", "r2", "r3", "r4"), "stimulus", ...),
    class = c("bmmodel", "sdt", "sdt_rating", "sdt_rating_standard"),
    params = c("d", "criterion", "bias", "spacing", "sdratio"),
    fixed = list(sdratio = 0),
    links = list(d = "identity", criterion = "identity", bias = "identity",
                 spacing = "identity", sdratio = "identity"),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    custom = list(d = "log"),
    settable = c("d", "criterion", "bias"),
    # sdratio and the thresholds are read through exp() in Stan, and sdratio is
    # additionally fixed at 0, so neither survives a change of link
    refused = list(sdratio = "log", spacing = "log"),
    unknown_link = "sqrt"
  ),
  sdt_mafc = list(
    make = function(...) sdt_mafc("n_correct", "n_trials", m = 4, ...),
    class = c("bmmodel", "sdt", "sdt_mafc"),
    params = "d",
    fixed = list(),
    links = list(d = "identity"),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    custom = list(d = "log"),
    settable = "d",
    refused = NULL,
    unknown_link = "nonsense"
  ),
  sdt_ranking = list(
    make = function(...) sdt_ranking(paste0("rank", 1:4), m = 4, ...),
    class = c("bmmodel", "sdt", "sdt_ranking"),
    params = "d",
    fixed = list(),
    links = list(d = "identity"),
    dists = c("gumbel_min", "normal"),
    custom = list(d = "log"),
    settable = "d",
    refused = NULL,
    unknown_link = "sqrt"
  ),
  sdt_ranking_normal = list(
    make = function(...) sdt_ranking(paste0("rank", 1:4), m = 4,
                                     dist = "normal", ...),
    class = c("bmmodel", "sdt", "sdt_ranking"),
    params = c("d", "sdratio"),
    fixed = list(sdratio = 0),
    links = list(d = "identity", sdratio = "identity"),
    dists = NULL,
    custom = NULL,
    settable = character(0),
    # Stan reads exp(sdratio), and the fixed 0 means equal variance only on the
    # identity link; d is relinkable only for gumbel_min, because the Gaussian
    # quadrature fails at the d a log link reaches
    refused = list(sdratio = c("log", "softplus", "logit"),
                   d = c("log", "softplus", "logit")),
    unknown_link = NULL
  ),
  sdt_cdp = list(
    make = function(...) sdt_cdp(stimulus = "s", n_new = 3, n_old = 3, ...),
    class = c("bmmodel", "sdt", "sdt_cdp"),
    params = c("dfam", "drec", "criterion", "spacing", "rcrit", "sigmar",
               "rho", "kcrit"),
    fixed = list(sigmar = 0, rho = 0, kcrit = -100),
    links = list(dfam = "identity", drec = "identity", criterion = "identity",
                 spacing = "identity", rcrit = "identity", sigmar = "identity",
                 rho = "identity", kcrit = "identity"),
    dists = "normal",
    custom = list(dfam = "log"),
    settable = c("dfam", "drec", "criterion", "rcrit"),
    # sigmar, rho and kcrit are fixed, and spacing is read through exp()
    refused = list(sigmar = "log", rho = "log", kcrit = "log", spacing = "log"),
    unknown_link = "sqrt"
  )
)

sdt_link_rows <- c(sdt_formats, list(
  sdt_rating_dpsdt = list(
    make = function(...) sdt_rating(paste0("r", 1:4), "stimulus",
                                    version = "dpsdt", ...),
    settable = c("d", "criterion", "bias"),
    refused = list(sdratio = "log", spacing = "log", Ro = "log", Rn = "log"),
    unknown_link = "sqrt"
  ),
  sdt_rating_metad = list(
    make = function(...) sdt_rating(paste0("r", 1:4), "stimulus",
                                    version = "metad", ...),
    settable = c("d", "criterion", "bias"),
    refused = list(sdratio = "log", spacing = "log", logmratio = "log"),
    unknown_link = "sqrt"
  )
))

# check_data() rows: a valid two-row data set, the count column the type
# checks probe (named in the message as `label`), and the columns each format
# reads besides the counts
sdt_data_rows <- list(
  sdt_yn = list(
    make = sdt_formats$sdt_yn$make,
    formula = bmf(d ~ 1, criterion ~ 1),
    data = data.frame(n_old = c(30, 40), stimulus = c(0L, 1L),
                      n_trials = c(50, 50)),
    count = "n_old", label = "Response variable 'n_old'",
    trials = "n_trials", stim = "stimulus",
    drop = c("stimulus", "n_trials"), missing = "missing in the data"
  ),
  sdt_rating = list(
    make = sdt_formats$sdt_rating$make,
    formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1),
    data = data.frame(r1 = c(10, 5), r2 = c(20, 10), r3 = c(15, 30),
                      r4 = c(5, 55), stimulus = c(0L, 1L)),
    count = "r1", label = "Response column 'r1'",
    trials = NULL, stim = "stimulus",
    drop = c("r3", "r4"), missing = "missing in the data"
  ),
  sdt_mafc = list(
    make = sdt_formats$sdt_mafc$make,
    formula = bmf(d ~ 1),
    data = data.frame(n_correct = c(70, 65), n_trials = c(100, 100)),
    count = "n_correct", label = "Response variable 'n_correct'",
    trials = "n_trials", stim = NULL,
    drop = "n_correct", missing = "missing in the data"
  ),
  sdt_ranking = list(
    make = sdt_formats$sdt_ranking$make,
    formula = bmf(d ~ 1),
    data = data.frame(rank1 = c(40, 30), rank2 = c(30, 30), rank3 = c(20, 20),
                      rank4 = c(10, 20)),
    count = "rank1", label = "Response column 'rank1'",
    trials = NULL, stim = NULL,
    drop = c("rank3", "rank4"), missing = "missing in the data"
  ),
  sdt_cdp = list(
    make = sdt_formats$sdt_cdp$make,
    formula = bmf(dfam ~ 1),
    data = data.frame(s = c(0L, 1L), new1 = c(30, 5), new2 = c(20, 5),
                      new3 = c(10, 10), know4 = c(10, 15), know5 = c(10, 15),
                      know6 = c(5, 10), remember4 = c(5, 10),
                      remember5 = c(5, 10), remember6 = c(5, 20)),
    count = "know4", label = "Response column 'know4'",
    trials = NULL, stim = "s",
    # missing count columns point at the aggregation helper
    drop = "new1",
    missing = "missing in the data. Use aggregate_sdt_cdp_data\\(\\)"
  )
)

# d/r rows: `d_vec` holds two-observation calls whose per-observation
# arguments have length 2 and are checked against one call per row; the
# `bad` lists hold argument overrides and the message each one must raise
sdt_cdp_thr <- .cdp_make_thresholds(0, -0.3, 3, 3, "parsimonious")
sdt_cdp_counts <- c(40, 20, 10, 15, 25, 35, 5, 20, 80)
sdt_dr_rows <- list(
  sdt_yn = list(
    r = rsdt_yn,
    r_args = list(n = 10, n_trials = 50, stimulus = rep(c(0L, 1L), 5),
                  d = 1.5, criterion = 0),
    cols = NULL,
    r_bad = list(list(args = list(n = 2, stimulus = c(0L, 2L)), msg = "0.*1")),
    d = dsdt_yn,
    d_args = list(n_old = 80, n_trials = 100, stimulus = 1, d = 1.5,
                  criterion = 0.2),
    d_vec = list(
      list(n_old = c(30, 80), n_trials = c(100, 100), stimulus = c(0, 1),
           d = 1.5, criterion = 0.2),
      list(n_old = c(30, 80), n_trials = c(100, 100), stimulus = c(0, 1),
           d = 1.5, criterion = 0.2, sdratio = c(1.0, 1.3))
    ),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    d_bad = list(
      list(args = list(n_old = -1), msg = "non-negative"),
      list(args = list(n_old = 150), msg = "must not exceed"),
      list(args = list(stimulus = 2), msg = "0.*1"),
      list(args = list(sdratio = 0), msg = "sdratio must be positive"),
      list(args = list(dist = "cauchy"), msg = "should be one of"),
      list(args = list(n_old = NA), msg = "n_old must not contain NA"),
      list(args = list(n_trials = NA), msg = "n_trials must not contain NA")
    )
  ),
  sdt_mafc = list(
    r = rsdt_mafc,
    r_args = list(n = 5, n_trials = 50, m = 4, d = 1.5),
    cols = NULL,
    r_bad = list(list(args = list(n = 2, m = 1), msg = "m must be")),
    d = dsdt_mafc,
    d_args = list(n_correct = 80, n_trials = 100, d = 1.5, m = 4),
    d_vec = list(
      list(n_correct = c(60, 80), n_trials = c(100, 100), d = 1.5, m = 4),
      list(n_correct = c(60, 80), n_trials = 100, d = c(1, 2), m = c(4, 4))
    ),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    d_bad = list(
      list(args = list(m = 1), msg = "m must be"),
      list(args = list(n_correct = 120), msg = "must not exceed"),
      list(args = list(m = NA), msg = "m must not contain NA"),
      list(args = list(n_correct = NA), msg = "n_correct must not contain NA"),
      list(args = list(n_trials = NA), msg = "n_trials must not contain NA")
    )
  ),
  sdt_rating = list(
    r = rsdt_rating,
    r_args = list(n = 10, n_trials = 50, stimulus = rep(c(0L, 1L), 5),
                  d = 1.5, thresholds = c(-0.5, 0, 0.5)),
    cols = paste0("r", 1:4),
    r_bad = list(
      list(args = list(n = 2, stimulus = c(0L, 2L)),
           msg = "0 \\(noise\\) or 1 \\(signal\\)"),
      list(args = list(n = 2, stimulus = c(0L, 1L), dist = "cauchy"),
           msg = "should be one of")
    ),
    d = dsdt_rating,
    d_args = list(counts = c(10, 20, 30, 40), stimulus = 1, d = 1.5,
                  thresholds = c(-0.5, 0, 0.5)),
    d_vec = list(
      list(counts = rbind(c(10, 20, 30, 40), c(40, 30, 20, 10)),
           stimulus = c(1L, 0L), d = 1.5, thresholds = c(-0.5, 0, 0.5))
    ),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    d_bad = list(
      list(args = list(counts = c(-1, 20, 30, 40)), msg = "non-negative"),
      list(args = list(counts = c(10, 20, 30)), msg = "K - 1"),
      list(args = list(counts = c(NA, 20, 30, 40)), msg = "must not contain NA")
    )
  ),
  sdt_ranking = list(
    r = rsdt_ranking,
    r_args = list(n = 5, n_trials = 100, m = 4, d = 1.5),
    cols = paste0("rank", 1:4),
    r_bad = list(list(args = list(n = 2, m = 1), msg = "m must be")),
    d = dsdt_ranking,
    d_args = list(counts = c(40, 30, 20, 10), d = 1.5, m = 4),
    d_vec = list(
      list(counts = rbind(c(40, 30, 20, 10), c(10, 20, 30, 40)), m = 4,
           d = c(1.5, 0.2))
    ),
    dists = c("gumbel_min", "normal"),
    d_bad = list(
      list(args = list(counts = rbind(c(30, 15, 5, 6), c(30, 12, 6, 2)),
                       m = c(3, 4), d = 1.2),
           msg = "beyond the row's set size"),
      list(args = list(m = NA), msg = "m must not contain NA")
    )
  ),
  sdt_cdp = list(
    r = rsdt_cdp,
    r_args = list(n = 10, n_trials = 100, stimulus = rep(c(0L, 1L), each = 5),
                  dfam = 0.8, drec = 1.0, thresholds = sdt_cdp_thr,
                  rcrit = 0.5, n_new = 3),
    cols = c("new1", "new2", "new3", "know4", "know5", "know6",
             "remember4", "remember5", "remember6"),
    r_bad = list(),
    d = dsdt_cdp,
    d_args = list(counts = sdt_cdp_counts, stimulus = 1, dfam = 0.8,
                  drec = 1.0, thresholds = sdt_cdp_thr, rcrit = 0.5,
                  n_new = 3),
    d_vec = list(
      list(counts = rbind(sdt_cdp_counts, c(60, 15, 5, 10, 20, 30, 10, 25, 75)),
           stimulus = c(0L, 1L), dfam = 0.8, drec = 1.0,
           thresholds = sdt_cdp_thr, rcrit = 0.5, n_new = 3)
    ),
    dists = "normal",
    d_bad = list(
      list(args = list(counts = sdt_cdp_counts[-1]), msg = "columns"),
      list(args = list(counts = replace(sdt_cdp_counts, 2, NA)),
           msg = "must not contain NA"),
      list(args = list(counts = replace(sdt_cdp_counts, 2, -1)),
           msg = "non-negative"),
      list(args = list(stimulus = 2), msg = "stimulus must be 0")
    )
  ),
  sdt_dpsdt = list(
    r = rsdt_dpsdt,
    r_args = list(n = 8, n_trials = 50, stimulus = rep(c(0L, 1L), 4), d = 1.4,
                  thresholds = .sdt_make_thresholds(0, 4L, "parsimonious", 0.3),
                  Ro = 0.3, Rn = 0.1),
    cols = paste0("r", 1:4),
    r_bad = list(),
    d = dsdt_dpsdt,
    d_args = list(counts = c(2, 8, 20, 70), stimulus = 1L, d = 1.5,
                  thresholds = c(-0.5, 0, 0.5), Ro = 0.3, Rn = 0),
    d_vec = list(
      list(counts = rbind(c(2, 8, 20, 70), c(70, 20, 8, 2)),
           stimulus = c(1L, 0L), d = 1.5, thresholds = c(-0.5, 0, 0.5),
           Ro = c(0.3, 0.2), Rn = c(0, 0.1))
    ),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    d_bad = list(list(args = list(Ro = 1.5), msg = "Ro must be a probability"))
  ),
  sdt_metad = list(
    r = rsdt_metad,
    r_args = list(n = 8, n_trials = 50, stimulus = rep(c(0L, 1L), 4), d = 1.4,
                  thresholds = .sdt_make_thresholds(0, 4L, "parsimonious", 0.3),
                  metad = 1.0),
    cols = paste0("r", 1:4),
    r_bad = list(),
    d = dsdt_metad,
    d_args = list(counts = c(5, 15, 25, 55), stimulus = 1L, d = 1.5,
                  thresholds = c(-0.5, 0, 0.5), metad = 1.0),
    d_vec = list(),
    dists = c("normal", "logistic", "gumbel_min", "gumbel_max"),
    d_bad = list()
  )
)

# the i-th observation of a d_vec call: count matrices give their row, and
# length-2 arguments their element
sdt_dr_obs <- function(args, i) {
  lapply(args, function(a) {
    if (is.matrix(a)) a[i, ] else if (length(a) == 2) a[i] else a
  })
}

test_that("every SDT constructor builds the class, parameters and links it documents", {
  for (name in names(sdt_formats)) {
    row <- sdt_formats[[name]]
    expect_silent(model <- row$make())
    expect_equal(class(model), row$class, info = name)
    expect_setequal(names(model$parameters), row$params)
    expect_false(isTRUE(model$void_mu), info = name)
    expect_setequal(names(model$default_priors), row$params)
    expect_equal(model$fixed_parameters, row$fixed, info = name)
    expect_equal(model$links, row$links, info = name)
    for (dist in row$dists) {
      expect_silent(model <- row$make(dist = dist))
      expect_equal(model$other_vars$dist, dist, info = paste(name, dist))
    }
    if (!is.null(row$custom)) {
      expect_equal(row$make(links = row$custom)$links,
                   utils::modifyList(row$links, row$custom), info = name)
    }
  }
})

test_that("every SDT model supplies an init range for each of its parameters", {
  # create_initfun() looks up init_ranges per parameter; a missing entry yields
  # NA inits, and the flexible rating threshold types reject brms' default
  # random init, so each of them gets its own row
  thresholds <- c("parsimonious", "equidistant", "log_distance", "log_ratio",
                  "softmax")
  models <- c(
    lapply(sdt_link_rows, function(row) row$make()),
    lapply(stats::setNames(nm = paste0("sdt_rating_", thresholds)), function(tt) {
      sdt_rating(paste0("r", 1:6), "stimulus",
                 threshold_type = sub("sdt_rating_", "", tt))
    })
  )
  for (name in names(models)) {
    model <- models[[name]]
    expect_true(all(names(model$parameters) %in% names(model$init_ranges)),
                info = name)
    expect_false("mu" %in% names(model$init_ranges), info = name)
    for (par in names(model$init_ranges)) {
      range <- model$init_ranges[[par]]
      expect_length(range, 2)
      expect_true(range[1] < range[2], info = paste(name, par))
    }
  }
})

test_that("every SDT constructor refuses the links it cannot apply, before and after construction", {
  for (name in names(sdt_link_rows)) {
    row <- sdt_link_rows[[name]]
    model <- row$make()
    # check_links() returns early without this attribute, which is how sdt_yn
    # once opted out of every check on a link assigned after construction
    expect_equal(attr(model, "links_default"), model$links, info = name)
    expect_equal(settable_links(model), row$settable, info = name)
    for (par in names(row$refused)) {
      refusal <- paste0("link of '", par, "' cannot be changed")
      for (link in row$refused[[par]]) {
        expect_error(row$make(links = stats::setNames(list(link), par)),
                     refusal, info = paste(name, par, link))
      }
      relinked <- row$make()
      relinked$links[[par]] <- row$refused[[par]][1]
      expect_error(check_links(relinked), refusal, info = paste(name, par))
    }
    # naming the link a parameter already has asks for no change
    for (par in names(model$links)) {
      expect_silent(row$make(links = stats::setNames(list(model$links[[par]]), par)))
    }
    expect_error(row$make(links = list(sensitivity = "log")),
                 "Unrecognized link target", info = name)
    relinked <- row$make()
    relinked$links$typo <- "log"
    expect_error(check_links(relinked), "Unrecognized link target", info = name)
    for (par in row$settable) {
      expect_error(row$make(links = stats::setNames(list(row$unknown_link), par)),
                   "Unknown link function", info = paste(name, par))
      relinked <- row$make()
      relinked$links[[par]] <- row$unknown_link
      expect_error(check_links(relinked), "Unknown link function",
                   info = paste(name, par))
    }
  }
})

test_that("every SDT check_data() refuses malformed counts and stimulus columns", {
  for (name in names(sdt_data_rows)) {
    row <- sdt_data_rows[[name]]
    model <- row$make()
    check <- function(data) check_data(model, data, row$formula)
    with_col <- function(col, value) {
      data <- row$data
      data[[col]] <- value
      data
    }
    # integer-valued doubles are counts too
    expect_silent(check(row$data))
    expect_error(check(row$data[setdiff(names(row$data), row$drop)]),
                 row$missing, info = name)
    for (conv in list(as.character, as.factor)) {
      expect_error(check(with_col(row$count, conv(row$data[[row$count]]))),
                   paste(row$label, "must be numeric"), info = name)
    }
    expect_warning(check(with_col(row$count, row$data[[row$count]] + 0.5)),
                   paste(row$label, "should contain integer counts"))
    expect_error(check(with_col(row$count, c(-1, row$data[[row$count]][-1]))),
                 paste(row$label, "must contain non-negative counts"),
                 info = name)

    if (!is.null(row$trials)) {
      trials <- row$data[[row$trials]]
      for (conv in list(as.character, as.factor)) {
        expect_error(check(with_col(row$trials, conv(trials))),
                     paste0("Variable '", row$trials, "' must be numeric"),
                     info = name)
      }
      expect_error(
        check(with_col(row$count, c(trials[1] + 10, row$data[[row$count]][-1]))),
        paste0("Response counts in '", row$count, "' must not exceed '",
               row$trials, "'"),
        info = name
      )
      zero <- with_col(row$count, c(0, 0))
      zero[[row$trials]] <- c(-10, trials[-1])
      expect_error(check(zero),
                   paste0("Variable '", row$trials, "' must contain positive values"),
                   info = name)
    }

    if (!is.null(row$stim)) {
      stim <- paste0("Stimulus variable '", row$stim, "'")
      for (coding in list(c(1, 2), c(0, 2))) {
        expect_error(check(with_col(row$stim, coding)),
                     paste(stim, "must be coded as 0 \\(noise\\) and 1 \\(signal\\)"),
                     info = name)
      }
      # one message for every problem told the user nothing about theirs
      expect_error(check(with_col(row$stim, c(NA_integer_, 1L))),
                   "but 1 of 2 values are NA", info = name)
      # a constant stimulus was accepted once, and d and criterion were then
      # "estimated" from one rate
      for (value in 0:1) {
        expect_error(check(with_col(row$stim, value)),
                     paste(stim, "is", value, "in every row"), info = name)
      }
    }
  }
})

test_that("every SDT d/r pair simulates the documented shape and refuses bad input", {
  withr::local_seed(2101)
  for (name in names(sdt_dr_rows)) {
    row <- sdt_dr_rows[[name]]
    call_r <- function(...) do.call(row$r, utils::modifyList(row$r_args, list(...)))
    call_d <- function(...) do.call(row$d, utils::modifyList(row$d_args, list(...)))

    sim <- call_r()
    n <- row$r_args$n
    if (is.null(row$cols)) {
      expect_length(sim, n)
      expect_true(all(sim >= 0 & sim <= row$r_args$n_trials), info = name)
    } else {
      expect_true(is.matrix(sim), info = name)
      expect_equal(dim(sim), c(n, length(row$cols)), info = name)
      expect_equal(colnames(sim), row$cols, info = name)
      expect_true(all(rowSums(sim) == row$r_args$n_trials), info = name)
    }
    expect_error(call_r(n = c(2, 3)), "single positive integer", info = name)
    for (bad in row$r_bad) {
      expect_error(do.call(call_r, bad$args), bad$msg, info = name)
    }

    dens <- call_d()
    expect_true(dens > 0 && dens <= 1, info = name)
    expect_equal(call_d(log = TRUE), log(dens), tolerance = 1e-10, info = name)
    for (dist in row$dists) {
      expect_true(call_d(dist = dist) > 0, info = paste(name, dist))
    }
    for (args in row$d_vec) {
      dens <- do.call(row$d, args)
      expect_length(dens, 2)
      expect_true(all(dens > 0), info = name)
      expect_equal(do.call(row$d, c(args, log = TRUE)), log(dens),
                   tolerance = 1e-10, info = name)
      for (i in 1:2) {
        expect_equal(dens[i], do.call(row$d, sdt_dr_obs(args, i)),
                     info = paste(name, "observation", i))
      }
    }
    for (bad in row$d_bad) {
      expect_error(do.call(call_d, bad$args), bad$msg, info = name)
    }
  }
})

############################################################################# !
# CDF HELPER TESTS                                                       ####
############################################################################# !

test_that(".sdt_cdf returns correct values for known inputs", {
  # Normal: Phi(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "normal"), 0.5)
  # Logistic: inv_logit(0) = 0.5
  expect_equal(bmm:::.sdt_cdf(0, "logistic"), 0.5)
  # Gumbel_min (smallest extreme value) at 0: 1 - exp(-exp(0)) = 1 - exp(-1)
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_min"), 1 - exp(-1), tolerance = 1e-10)
  # Gumbel_max (largest extreme value) at 0: exp(-exp(0)) = exp(-1)
  expect_equal(bmm:::.sdt_cdf(0, "gumbel_max"), exp(-1), tolerance = 1e-10)
})

test_that("gumbel labels follow the extreme-value convention", {
  # the min-type is the mirror image of the max-type about zero
  eta <- c(-2, -0.5, 0, 0.5, 2)
  expect_equal(bmm:::.sdt_cdf(eta, "gumbel_min"),
               1 - bmm:::.sdt_cdf(-eta, "gumbel_max"))
  # gumbel_max is the log-log / evd::pgumbel parameterisation
  expect_equal(bmm:::.sdt_cdf(eta, "gumbel_max"), exp(-exp(-eta)))
})

test_that("every dist argument offers exactly the registry's distributions", {
  # the registry order defines the dist_type integer passed to Stan, so a
  # signature that drifts out of step with it becomes an off-by-one
  fns <- list(sdt_yn, dsdt_yn, rsdt_yn, sdt_d, sdt_criterion,
              sdt_mafc, dsdt_mafc, rsdt_mafc,
              sdt_rating, dsdt_rating, rsdt_rating,
              dsdt_dpsdt, rsdt_dpsdt, dsdt_metad, rsdt_metad)
  for (f in fns) {
    expect_equal(eval(formals(f)$dist), names(bmm:::.sdt_dists))
  }
})

test_that("the noise distribution reaches Stan as the id the Stan switch reads", {
  # sdt_dist_funs.stan dispatches on 1 normal, 2 gumbel_min, 3 gumbel_max,
  # 4 logistic; the code is the same for every distribution, so the id in the
  # data is the only thing that selects it
  ids <- c(normal = 1L, gumbel_min = 2L, gumbel_max = 3L, logistic = 4L)
  yn_data <- data.frame(n_old = c(10, 40), stimulus = c(0L, 1L), n_trials = 50)
  mafc_data <- data.frame(n_correct = c(70, 65), n_trials = 100)
  for (dist in names(ids)) {
    yn <- standata(bmf(d ~ 1, criterion ~ 1), yn_data,
                   sdt_yn("n_old", "stimulus", "n_trials", dist = dist))
    mafc <- standata(bmf(d ~ 1), mafc_data,
                     sdt_mafc("n_correct", "n_trials", m = 4, dist = dist))
    expect_equal(as.vector(yn$vint2), rep(ids[[dist]], 2), label = paste("sdt_yn", dist))
    expect_equal(as.vector(mafc$vint2), rep(ids[[dist]], 2), label = paste("sdt_mafc", dist))
  }
})

test_that("a predictor and a random effect on d reach every SDT model's program", {
  skip_on_cran() # over 0.4 s in CRAN mode (ADR 0004, R10)
  cells <- data.frame(id = rep(1:2, each = 2), stimulus = c(0L, 1L, 0L, 1L),
                      condition = c("A", "A", "B", "B"))
  cases <- list(
    sdt_yn = list(sdt_yn("n_old", "stimulus", "n_trials"),
                  cbind(cells, n_old = c(10, 40, 15, 35), n_trials = 50), bmf(criterion ~ 1)),
    sdt_mafc = list(sdt_mafc("n_correct", "n_trials", m = 4),
                    cbind(cells, n_correct = c(70, 65, 80, 60), n_trials = 100), NULL),
    sdt_ranking = list(sdt_ranking(paste0("rank", 1:4), m = 4),
                       cbind(cells, rank1 = 20L, rank2 = 15L, rank3 = 10L, rank4 = 5L), NULL),
    sdt_rating = list(sdt_rating(paste0("r", 1:4), "stimulus"),
                      cbind(cells, r1 = 20L, r2 = 15L, r3 = 10L, r4 = 5L),
                      bmf(criterion ~ 1, spacing ~ 1))
  )
  for (name in names(cases)) {
    case <- cases[[name]]
    formula <- bmf(d ~ condition + (1 | id))
    if (!is.null(case[[3]])) formula <- formula + case[[3]]
    code <- stancode(formula, data = case[[2]], model = case[[1]])
    # Kc_d and Xc_d where d is a dpar with a centred design, K_d and X_d where
    # it is an nlpar
    for (term in c("vector\\[Kc?_d\\] b_d;", "Xc?_d \\* b_d", "r_1_d_1\\[J_1\\[n\\]\\] \\* Z_1_d_1\\[n\\]")) {
      expect_match(code, term, label = paste(name, term))
    }
  }
})

test_that("a two-category rating reduces to the binary likelihood", {
  # both models must implement the same decision rule -- evidence above the
  # criterion is an "old" response -- which for the asymmetric extreme-value
  # distributions is only true if neither evaluates its cdf on the mirrored axis
  for (di in names(bmm:::.sdt_dists)) {
    for (stim in c(0L, 1L)) {
      expect_equal(
        bmm:::.sdt_category_probs(-0.2, 1.5, 1.3, stim, di)[2],
        dsdt_yn(1, 1, stim, d = 1.5, criterion = -0.2,
                sdratio = 1.3, dist = di),
        info = paste(di, "stimulus", stim)
      )
    }
  }
})

test_that("the m-AFC probabilities agree with the registry's own cdf", {
  # P(correct) is the probability the signal variate beats m - 1 distractors,
  # int f(x - d) F(x)^(m - 1) dx, so a closed form is only right if it agrees
  # with the cdf the rest of the model uses. Taking the density off that same
  # cdf by central difference ties each branch of .mafc_pc_r() to the registry
  # entry it claims: relabelling gumbel_min and gumbel_max leaves every block
  # in test-model_sdt_mafc.R green, because its closed-form tests restate the
  # formula the implementation already uses, and fails here. gumbel_min and
  # gumbel_max coincide at m = 2, which is why the grid goes past it.
  # 1e-7 is the implementation's side: the 40-point Gauss-Hermite normal branch
  # sits 4.5e-8 from the integral at m = 8, d = 0, which is the rule's own
  # error and not the central difference's (h = 1e-6 and an exact density both
  # give 4.5e-8); every other cell stays below 7e-9.
  pc_integral <- function(d, m, dist) {
    cdf <- bmm:::.sdt_dists[[dist]]$cdf
    dens <- function(x) (cdf(x + 1e-5) - cdf(x - 1e-5)) / 2e-5
    stats::integrate(function(x) dens(x - d) * cdf(x)^(m - 1),
                     -Inf, Inf, rel.tol = 1e-10)$value
  }
  for (dist in names(bmm:::.sdt_dists)) {
    for (m in c(2L, 3L, 5L, 8L)) {
      for (d in c(0, 0.8, 2)) {
        expect_equal(bmm:::.mafc_pc_r(d, m, dist), pc_integral(d, m, dist),
                     tolerance = 1e-7,
                     info = paste(dist, "m =", m, "d =", d))
      }
    }
  }
})

test_that("the default sdratio prior brackets the group-level estimates", {
  # two group-level SD ratios are on record: Mickes et al. (2007) Table 1
  # averages sd(lure)/sd(target) = 0.79 over their 13 retained subjects, i.e.
  # 1.26 signal over noise, and the broeder_schuetz_2009_e3 posterior gives
  # 1.46 [1.25, 1.71]. The prior has to reach past the widest ratio those
  # support and still keep its mass off ratios above 2, which recognition does
  # not produce -- one-sided, "covers the empirical range" cannot fail.
  models <- list(sdt_yn("n_old", "stimulus", "n_trials"),
                 sdt_ranking(paste0("rank", 1:4), m = 4, dist = "normal"))
  for (model in models) {
    main <- model$default_priors$sdratio$main
    sd <- as.numeric(sub("normal\\(0, ([0-9.]+)\\)", "\\1", main))
    expect_true(is.finite(sd), info = main)
    expect_gt(exp(qnorm(0.975, 0, sd)), 1.71)
    expect_gt(diff(pnorm(log(c(0.5, 2)), 0, sd)), 0.95)
  }
})

test_that("sdratio priors agree across the SDT models", {
  yn <- sdt_yn("n_old", "stimulus", "n_trials")$default_priors$sdratio
  rk <- sdt_ranking(paste0("rank", 1:4), m = 4,
                    dist = "normal")$default_priors$sdratio
  rt <- sdt_rating(paste0("r", 1:4), "stimulus")$default_priors$sdratio
  expect_identical(yn, rk)
  expect_identical(yn, rt)
})

test_that("quantile functions invert their cdfs", {
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  for (d in names(bmm:::.sdt_dists)) {
    expect_equal(bmm:::.sdt_cdf(bmm:::.sdt_dists[[d]]$qf(p), d), p,
                 tolerance = 1e-10, info = d)
  }
})


############################################################################# !
# LOG-SCALE CDF TESTS                                                    ####
############################################################################# !

test_that("lcdf and lccdf agree with the probability scale where it is exact", {
  # only the central range: further out the naive log(cdf) is the inaccurate
  # side of the comparison, which is the whole reason the log-scale pair exists
  eta <- seq(-2, 2, by = 0.5)
  for (d in names(bmm:::.sdt_dists)) {
    p <- bmm:::.sdt_cdf(eta, d)
    expect_equal(bmm:::.sdt_log_cdf(eta, d), log(p), tolerance = 1e-10, info = d)
    expect_equal(bmm:::.sdt_log_ccdf(eta, d), log(1 - p),
                 tolerance = 1e-10, info = d)
  }
})

test_that("lcdf and lccdf stay finite where the probability scale underflows", {
  # the naive path returns log(0) = -Inf here; Stan and the log-scale path do not
  for (d in names(bmm:::.sdt_dists)) {
    expect_true(is.finite(bmm:::.sdt_log_cdf(-40, d)), info = d)
    expect_true(is.finite(bmm:::.sdt_log_ccdf(40, d)), info = d)
  }
  expect_equal(log(bmm:::.sdt_cdf(-40, "normal")), -Inf)
})

test_that(".sdt_cdf is vectorized over eta", {
  eta <- c(-1, 0, 1)
  result <- bmm:::.sdt_cdf(eta, "normal")
  expect_length(result, 3)
  expect_equal(result, pnorm(eta))
})


############################################################################# !
# PDF HELPER TESTS                                                       ####
############################################################################# !

test_that(".sdt_dists pdf is the derivative of its cdf and integrates to 1", {
  x <- seq(-4, 4, by = 0.1)
  h <- 1e-5
  for (dist in names(bmm:::.sdt_dists)) {
    cdf <- bmm:::.sdt_dists[[dist]]$cdf
    pdf <- bmm:::.sdt_dists[[dist]]$pdf
    fd  <- (cdf(x + h) - cdf(x - h)) / (2 * h)
    expect_equal(pdf(x), fd, tolerance = 1e-5, info = dist)
    expect_equal(integrate(pdf, -Inf, Inf)$value, 1, tolerance = 1e-6,
                 info = dist)
  }
})


############################################################################# !
# DPRIME AND CRITERION TESTS                                              ####
############################################################################# !

test_that("sdt_d computes correct values for normal distribution", {
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_d(hr, far, dist = "normal"), 2, tolerance = 1e-10)
})

test_that("sdt_criterion and sdt_bias compute correct values for normal distribution", {
  hr <- pnorm(1)
  far <- pnorm(-1)
  expect_equal(sdt_criterion(hr, far, dist = "normal"), 1, tolerance = 1e-10)
  expect_equal(sdt_bias(hr, far, dist = "normal"), 0, tolerance = 1e-10)
})

test_that("sdt_d validates input", {
  expect_error(sdt_d(0, 0.5), "between 0 and 1")
  expect_error(sdt_d(1, 0.5), "between 0 and 1")
  expect_error(sdt_d(0.5, 0), "between 0 and 1")
})


############################################################################# !
# SHARED FORMULA CHECKS                                                    ####
############################################################################# !

sdt_stim_data <- function() {
  data.frame(
    n_old = c(10, 40, 15, 35),
    stimulus = c(0L, 1L, 0L, 1L),
    n_trials = c(50, 50, 50, 50),
    condition = c("A", "A", "B", "B"),
    id = c(1, 1, 2, 2)
  )
}

# check_model() must run first, or add_missing_parameters() overwrites a formula
# the user gave a fixed parameter: on a raw model `sdratio ~ stimulus` comes back
# as `sdratio ~ 1` and a direct check_formula() call tests nothing
sdt_check_formula <- function(formula, data = sdt_stim_data(),
                              model = sdt_yn("n_old", "stimulus",
                                             "n_trials")) {
  model <- check_model(model, data, formula)
  check_formula(model, data, formula)
}

test_that("a stimulus term on a parameter formula is refused", {
  # the likelihood already consumes stimulus to set each row's role, so the term
  # is confounded with the parameter's own intercept -- profiled flat to 6e-11
  expect_error(
    sdt_check_formula(bmf(d ~ stimulus, criterion ~ 1)),
    "'d'.*uses the stimulus variable 'stimulus'"
  )
  expect_error(
    sdt_check_formula(bmf(d ~ 1, criterion ~ stimulus)),
    "'criterion'.*uses the stimulus variable"
  )
  expect_error(
    sdt_check_formula(bmf(d ~ 1, criterion ~ 1, sdratio ~ stimulus)),
    "'sdratio'.*uses the stimulus variable"
  )
  # both offenders are named, not just the first
  expect_error(
    sdt_check_formula(bmf(d ~ stimulus, criterion ~ stimulus)),
    "'d', 'criterion'"
  )
})

test_that("a stimulus grouping factor is refused with the same message", {
  # rhs_vars() reports random-effect grouping variables, and a random intercept
  # per stimulus level is the same confound as a fixed effect
  expect_error(
    sdt_check_formula(bmf(d ~ 1 + (1 | stimulus), criterion ~ 1)),
    "as a predictor or as a grouping factor"
  )
})

test_that("the stimulus check leaves legitimate predictors alone", {
  expect_no_error(
    sdt_check_formula(bmf(d ~ condition, criterion ~ 1 + (1 | id)))
  )
  # neither the confound check nor the all-intercept sdratio heuristic applies
  # to a design with a real predictor
  expect_no_warning(
    sdt_check_formula(bmf(d ~ condition, criterion ~ 1 + (1 | id)))
  )
  # and it fires through the full pipeline, not only on a direct call
  expect_error(
    bmm(bmf(d ~ stimulus, criterion ~ 1), sdt_stim_data(),
        sdt_yn("n_old", "stimulus", "n_trials"),
        backend = "mock", mock_fit = 1, rename = FALSE),
    "uses the stimulus variable"
  )
})

test_that("a format without a stimulus other-var skips the confound check", {
  # sdt_mafc and sdt_ranking have no stimulus other-var, and check_formula.sdt
  # must fall through to NextMethod() instead of failing on a NULL stim_var
  model <- structure(
    list(name = "stub", parameters = list(par1 = ""),
         resp_vars = list(response = "n_old"),
         other_vars = list(stimulus = NULL)),
    class = c("bmmodel", "sdt", "stub")
  )
  formula <- bmf(par1 ~ 1)
  expect_identical(check_formula(model, sdt_stim_data(), formula), formula)
})


############################################################################# !
# VECTORIZATION TESTS                                                      ####
############################################################################# !

test_that(".sdt_eta is vectorized over all arguments", {
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -0.2)
  expect_equal(eta[2], 1.5 - 0.2)

  eta <- bmm:::.sdt_eta(c(1.0, 2.0), 0.2, c(0, 1))
  expect_length(eta, 2)
  expect_equal(eta[1], -0.2)
  expect_equal(eta[2], 2.0 - 0.2)

  # `d` is d_a, so the separation in noise units is d * sqrt((1 + r^2) / 2)
  s <- sqrt((1 + c(1.0, 1.3)^2) / 2)
  eta <- bmm:::.sdt_eta(1.5, 0.2, c(0, 1), sdratio = c(1.0, 1.3))
  expect_length(eta, 2)
  expect_equal(eta[1], -0.2)
  expect_equal(eta[2], (1.5 * s[2] - 0.2) / 1.3)
})

test_that(".sdt_eta divides only the signal trials by sdratio", {
  # the noise trial carries no scale division, so its eta must match an
  # equal-variance call with the separation widened to d * s
  s <- sqrt((1 + 1.5^2) / 2)
  eta_uv <- bmm:::.sdt_eta(1.5, 0.2, 0, sdratio = 1.5)
  eta_ev <- bmm:::.sdt_eta(1.5 * s, 0.2, 0, sdratio = 1)
  expect_equal(eta_uv, eta_ev)
})
