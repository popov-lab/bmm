# Shared by test-update.R and test-reliability.R.

# one mock-fittable (model, formula, data) per model, each chosen so that the
# stored model frame is the hard case: m3 with a category that has zero options
# in some rows and with generated option columns, and the two non-target models
# with a formula that leaves brms no reason to keep the set_size column;
# mixture3p_set_size covers the opposite case, where the formula names set_size
# directly and brms keeps the real column instead
cdp_data <- function(guess, prefix = "rk") {
  cols <- .sdt_cdp_response_cols(2, 2, guess, prefix)
  counts <- matrix(c(30L, 20L, 5L, 7L, 9L, 11L, 4L, 14L), 12, 8, byrow = TRUE)
  cbind(
    data.frame(stimulus = rep(0:1, 6), id = factor(rep(1:6, each = 2))),
    stats::setNames(as.data.frame(counts[, seq_along(cols)]), cols)
  )
}

stored_frame_cases <- function() {
  rt_data <- data.frame(
    rt = rep(c(0.6, 0.8, 1.1, 0.7), 5),
    response = rep(c(1, 0), 10),
    id = factor(rep(1:5, each = 4))
  )
  rt_formula <- bmf(drift ~ 1, bound ~ 1, ndt ~ 1)
  ez_data <- data.frame(
    mean_rt = rep(c(0.5, 0.6), 10), var_rt = rep(c(0.02, 0.03), 10),
    n_upper = rep(c(60, 70), 10), n_trials = 100, id = factor(rep(1:10, each = 2))
  )
  nt_features <- paste0("col_nt", 1:7)
  # two participants still carry every set size from 1 to 8, so they exercise
  # the same LureIdx columns at a tenth of the mock-fitting cost
  lin_2017 <- oberauer_lin_2017[oberauer_lin_2017$ID %in% 1:2, ]
  m3_cats <- c("corr", "other", "dist", "npl")
  m3_formula <- bmf(
    corr ~ b + a + c, other ~ b + a, dist ~ b + d, npl ~ b,
    c ~ 1, a ~ 1, d ~ 1
  )
  m3_links <- list(c = "log", a = "log", d = "log")
  mafc_data <- data.frame(
    n_correct = c(80, 55, 78, 60, 85, 52, 81, 58), n_trials = 100,
    n_afc = rep(c(2, 4), 4), cond = factor(rep(c("a", "b"), each = 4))
  )
  ranking_data <- meyer_grant_jakob_2025[as.integer(meyer_grant_jakob_2025$id) <= 4, ]
  rating_data <- data.frame(
    stimulus = rep(c(0L, 1L), 4), id = factor(rep(1:4, each = 2)),
    r1 = c(30, 8, 26, 10, 33, 6, 28, 9), r2 = c(25, 12, 27, 14, 22, 11, 24, 13),
    r3 = c(20, 15, 21, 16, 19, 17, 22, 15), r4 = c(15, 25, 16, 22, 17, 28, 14, 26),
    r5 = c(10, 40, 10, 38, 9, 38, 12, 37)
  )

  list(
    cswald = list(
      model = cswald("rt", "response", version = "simple"),
      formula = rt_formula, data = rt_data
    ),
    ddm = list(model = ddm("rt", "response"), formula = rt_formula, data = rt_data),
    ezdm = list(
      model = ezdm("mean_rt", "var_rt", "n_upper", "n_trials", version = "3par"),
      formula = rt_formula, data = ez_data
    ),
    imm = list(
      model = imm("dev_rad",
        nt_features = nt_features,
        nt_distances = paste0("dist_nt", 1:7), set_size = "set_size"
      ),
      formula = bmf(c ~ 1, a ~ 1, s ~ 1, kappa ~ 1), data = lin_2017
    ),
    m3 = list(
      model = m3(
        resp_cats = m3_cats, num_options = paste0("n_", m3_cats),
        choice_rule = "simple", links = m3_links
      ),
      formula = m3_formula, data = oberauer_lewandowsky_2019_e1
    ),
    m3_num_options = list(
      model = m3(
        resp_cats = m3_cats,
        num_options = c(opt_corr = 1, opt_other = 4, opt_dist = 5, opt_npl = 5),
        choice_rule = "simple", links = m3_links
      ),
      formula = m3_formula, data = oberauer_lewandowsky_2019_e1
    ),
    m3_num_options_by_category = list(
      model = m3(
        resp_cats = m3_cats,
        num_options = c(other = 4, corr = 1, npl = 5, dist = 5),
        choice_rule = "simple", links = m3_links
      ),
      formula = m3_formula, data = oberauer_lewandowsky_2019_e1
    ),
    m3_char_num_options_by_category = list(
      model = m3(
        resp_cats = m3_cats,
        num_options = c(other = "n_other", corr = "n_corr", npl = "n_npl", dist = "n_dist"),
        choice_rule = "simple", links = m3_links
      ),
      formula = m3_formula, data = oberauer_lewandowsky_2019_e1
    ),
    mixture2p = list(
      model = mixture2p("dev_rad"), formula = bmf(kappa ~ 1, thetat ~ 1),
      data = lin_2017
    ),
    mixture3p = list(
      model = mixture3p("dev_rad", nt_features = nt_features, set_size = "set_size"),
      formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 1), data = lin_2017
    ),
    mixture3p_set_size = list(
      model = mixture3p("dev_rad", nt_features = nt_features, set_size = "set_size"),
      formula = bmf(kappa ~ 1, thetat ~ 1, thetant ~ 0 + set_size), data = lin_2017
    ),
    sdm = list(
      model = sdm("dev_rad"), formula = bmf(c ~ 1, kappa ~ 1),
      data = lin_2017
    ),
    sdt_yn = list(
      model = sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"),
      formula = bmf(d ~ 1, criterion ~ 1, sdratio ~ 1), data = broeder_schuetz_2009_e3
    ),
    sdt_mafc = list(
      model = sdt_mafc("n_correct", "n_trials", m = 4),
      formula = bmf(d ~ 1 + cond), data = mafc_data
    ),
    sdt_mafc_m_column = list(
      model = sdt_mafc("n_correct", "n_trials", m = "n_afc"),
      formula = bmf(d ~ 1 + cond), data = mafc_data
    ),
    sdt_mafc_m_predictor = list(
      model = sdt_mafc("n_correct", "n_trials", m = "n_afc"),
      formula = bmf(d ~ 1 + n_afc), data = mafc_data
    ),
    sdt_ranking = list(
      model = sdt_ranking(paste0("rank", 1:5), m = 5),
      formula = bmf(d ~ 1), data = ranking_data[ranking_data$set_size == 5, ]
    ),
    sdt_ranking_m_column = list(
      model = sdt_ranking(paste0("rank", 1:5), m = "set_size", dist = "normal"),
      formula = bmf(d ~ 1, sdratio ~ 1), data = ranking_data
    ),
    sdt_ranking_m_predictor = list(
      model = sdt_ranking(paste0("rank", 1:5), m = "set_size"),
      formula = bmf(d ~ 1 + set_size), data = ranking_data
    ),
    sdt_rating = list(
      model = sdt_rating(paste0("r", 1:5), "stimulus"),
      formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, sdratio ~ 1),
      data = rating_data
    ),
    sdt_rating_deltas = list(
      model = sdt_rating(paste0("r", 1:5), "stimulus", threshold_type = "log_distance"),
      formula = bmf(d ~ 1, criterion ~ 1, delta1 ~ 1, delta2 ~ 1, delta3 ~ 1),
      data = rating_data
    ),
    sdt_rating_dpsdt = list(
      model = sdt_rating(paste0("r", 1:5), "stimulus", version = "dpsdt"),
      formula = bmf(d ~ 1 + (1 | id), criterion ~ 1, spacing ~ 1, Ro ~ 1),
      data = rating_data
    ),
    # a column prefix and the Know/Guess split, which check_data() infers from
    # the guess columns alone
    sdt_cdp = list(
      model = sdt_cdp("rk", "stimulus", n_new = 2, n_old = 2),
      formula = bmf(dfam ~ 1 + (1 | id), drec ~ 1, criterion ~ 1, spacing ~ 1,
                    rcrit ~ 1, kcrit ~ 1),
      data = cdp_data(guess = TRUE)
    ),
    sdt_cdp_deltas = list(
      model = sdt_cdp(stimulus = "stimulus", n_new = 2, n_old = 2,
                      threshold_type = "log_distance"),
      formula = bmf(dfam ~ 1, drec ~ 1, criterion ~ 1, rcrit ~ 1, sigmar ~ 1,
                    delta1 ~ 1, delta2 ~ 1),
      data = cdp_data(guess = FALSE, prefix = "")
    ),
    # meta-d' needs an even number of categories
    sdt_rating_metad = list(
      model = sdt_rating(paste0("r", 1:4), "stimulus", version = "metad"),
      formula = bmf(d ~ 1, criterion ~ 1, spacing ~ 1, logmratio ~ 1 + (1 | id)),
      data = rating_data[setdiff(names(rating_data), "r5")]
    )
  )
}

stored_frame_fit <- function(case) {
  # the toy rt_data has a 50% error rate, which cswald "simple" warns about
  suppressWarnings(suppressMessages(
    bmm(case$formula, case$data, case$model,
      backend = "mock", mock_fit = 1, rename = FALSE
    )
  ))
}
