# One small fit per model, shared by make-references.R (which stores their
# fixed effects in references.rds), test-references.R (which compares a refit
# with them) and test-conditional_effects-models.R. A case changes the
# reference: after editing one, rerun
#   Rscript tests/local/make-references.R <model>
# Every formula gives at least one parameter a two-level predictor, so
# conditional_effects() has an effect to show for each model.

reference_cases <- function() {
  two_conds <- function(x) factor(x, levels = c("A", "B"))
  first_ids <- function(data, id, n) {
    data[data[[id]] %in% utils::head(sort(unique(data[[id]])), n), ]
  }
  tree_target <- mpt_tree("target", old ~ D + (1 - D) * g, new ~ (1 - D) * (1 - g))
  tree_lure <- mpt_tree("lure", old ~ (1 - D) * g, new ~ D + (1 - D) * (1 - g))

  ol17 <- first_ids(oberauer_lin_2017, "ID", 3)
  ol17 <- ol17[ol17$set_size %in% c(2, 4), ]
  ol17$set_size <- droplevels(ol17$set_size)
  zl08 <- first_ids(zhang_luck_2008, "subID", 4)
  zl08 <- zl08[zl08$setsize %in% c(2, 6), ]
  zl08$setsize <- droplevels(zl08$setsize)
  ole19 <- first_ids(oberauer_lewandowsky_2019_e1, "ID", 20)
  ole19 <- ole19[ole19$cond %in% c("old reordered", "old same"), ]
  ole19$cond <- droplevels(ole19$cond)

  list(
    cswald = withr::with_seed(101, {
      cond <- two_conds(rep(c("A", "B"), each = 150))
      dat <- rcswald(300, drift = ifelse(cond == "A", 2, 3), bound = 1.5, ndt = 0.3)
      list(
        model = cswald(rt = "rt", response = "response"),
        formula = bmf(drift ~ cond, bound ~ 1, ndt ~ 1),
        data = cbind(dat, cond = cond), seed = 101
      )
    }),
    ddm = withr::with_seed(102, {
      cond <- two_conds(rep(c("A", "B"), each = 150))
      dat <- rddm(300, drift = ifelse(cond == "A", 1, 2), bound = 1.2, ndt = 0.3)
      list(
        model = ddm(rt = "rt", response = "response"),
        formula = bmf(drift ~ cond, bound ~ 1, ndt ~ 1),
        data = cbind(dat, cond = cond), seed = 102
      )
    }),
    ezcdm = withr::with_seed(115, {
      dat <- data.frame(id = rep(1:10, 2), cond = two_conds(rep(c("A", "B"), each = 10)))
      dat <- cbind(dat, rezcdm(20, 100, drift = ifelse(dat$cond == "A", 1.5, 2.5),
                               bound = 1.5, ndt = 0.3))
      list(
        model = ezcdm(mean_angle = "mean_angle", var_angle = "var_angle",
                      mean_rt = "mean_rt", var_rt = "var_rt", n_trials = "n_trials"),
        formula = bmf(drift ~ cond, bound ~ 1, ndt ~ 1),
        data = dat, seed = 115
      )
    }),
    ezdm = withr::with_seed(103, {
      dat <- data.frame(id = rep(1:10, 2), cond = two_conds(rep(c("A", "B"), each = 10)))
      dat <- cbind(dat, rezdm(20, 100, drift = ifelse(dat$cond == "A", 1, 2),
                              bound = 1.5, ndt = 0.3))
      list(
        model = ezdm(mean_rt = "mean_rt", var_rt = "var_rt",
                     n_upper = "n_upper", n_trials = "n_trials"),
        formula = bmf(drift ~ cond, bound ~ 1, ndt ~ 1),
        data = dat, seed = 103
      )
    }),
    imm = list(
      model = imm(resp_error = "dev_rad", nt_features = paste0("col_nt", 1:3),
                  nt_distances = paste0("dist_nt", 1:3), set_size = "set_size"),
      formula = bmf(c ~ 0 + set_size, a ~ 1, s ~ 1, kappa ~ 1),
      data = ol17, seed = 104
    ),
    m3 = list(
      model = m3(resp_cats = c("corr", "other", "npl"),
                 num_options = c("n_corr", "n_other", "n_npl"), version = "ss"),
      formula = bmf(c ~ 0 + cond, a ~ 1),
      data = ole19, seed = 105
    ),
    mixture2p = list(
      model = mixture2p(resp_error = "response_error"),
      formula = bmf(kappa ~ 0 + setsize, thetat ~ 0 + setsize),
      data = zl08, seed = 106
    ),
    mixture3p = list(
      model = mixture3p(resp_error = "response_error", nt_features = "col_lure",
                        set_size = "setsize", regex = TRUE),
      formula = bmf(kappa ~ 1, thetat ~ 0 + setsize, thetant ~ 0 + setsize),
      data = zl08, seed = 107
    ),
    mpt = withr::with_seed(108, {
      dat <- expand.grid(id = 1:5, item_type = c("target", "lure"),
                         cond = two_conds(c("A", "B")), stringsAsFactors = FALSE)
      D <- ifelse(dat$cond == "A", 0.4, 0.7)
      p_old <- ifelse(dat$item_type == "target", D + (1 - D) * 0.4, (1 - D) * 0.4)
      dat$old <- stats::rbinom(nrow(dat), 50, p_old)
      dat$new <- 50L - dat$old
      list(
        model = mpt(list(tree_target, tree_lure), tree_id = "item_type"),
        formula = bmf(D ~ cond, g ~ 1),
        data = dat, seed = 108
      )
    }),
    sdm = list(
      model = sdm(resp_error = "dev_rad"),
      formula = bmf(c ~ 0 + set_size, kappa ~ 0 + set_size),
      data = ol17, seed = 109
    ),
    sdt_cdp = withr::with_seed(110, {
      dat <- expand.grid(id = 1:6, stimulus = c(0L, 1L), cond = two_conds(c("A", "B")))
      dat <- cbind(dat, rsdt_cdp(nrow(dat), 200, dat$stimulus,
                                 dfam = ifelse(dat$cond == "A", 0.5, 1.2), drec = 1,
                                 thresholds = c(-1.1, -0.5, 0, 0.6, 1.3),
                                 rcrit = 0.7, n_new = 3))
      list(
        model = sdt_cdp(stimulus = "stimulus", n_new = 3, n_old = 3),
        formula = bmf(dfam ~ cond, drec ~ 1, criterion ~ 1, spacing ~ 1, rcrit ~ 1),
        data = dat, seed = 110
      )
    }),
    sdt_mafc = withr::with_seed(111, {
      dat <- data.frame(cond = two_conds(rep(c("A", "B"), each = 10)), n_trials = 100L)
      dat$n_correct <- rsdt_mafc(nrow(dat), dat$n_trials, m = 4,
                                 d = ifelse(dat$cond == "A", 1, 2))
      list(
        model = sdt_mafc(response = "n_correct", n_trials = "n_trials", m = 4),
        formula = bmf(d ~ cond),
        data = dat, seed = 111
      )
    }),
    sdt_ranking = withr::with_seed(112, {
      dat <- data.frame(cond = two_conds(rep(c("A", "B"), each = 10)))
      dat <- cbind(dat, rsdt_ranking(nrow(dat), 100, m = 4,
                                     d = ifelse(dat$cond == "A", 0.8, 1.6)))
      list(
        model = sdt_ranking(response = paste0("rank", 1:4), m = 4),
        formula = bmf(d ~ cond),
        data = dat, seed = 112
      )
    }),
    sdt_rating = withr::with_seed(113, {
      dat <- expand.grid(id = 1:6, stimulus = c(0L, 1L), cond = two_conds(c("A", "B")))
      dat <- cbind(dat, rsdt_rating(nrow(dat), 150, dat$stimulus,
                                    d = ifelse(dat$cond == "A", 1, 2),
                                    thresholds = c(-0.5, 0, 0.5)))
      list(
        model = sdt_rating(response = paste0("r", 1:4), stimulus = "stimulus"),
        formula = bmf(d ~ cond, criterion ~ 1, spacing ~ 1),
        data = dat, seed = 113
      )
    }),
    sdt_yn = withr::with_seed(114, {
      dat <- expand.grid(id = 1:10, stimulus = c(0L, 1L), cond = two_conds(c("A", "B")))
      dat$n_trials <- 100L
      dat$n_old <- rsdt_yn(nrow(dat), dat$n_trials, dat$stimulus,
                           d = ifelse(dat$cond == "A", 1, 2), criterion = 0.5)
      list(
        model = sdt_yn(response = "n_old", stimulus = "stimulus", n_trials = "n_trials"),
        formula = bmf(d ~ cond, criterion ~ 1),
        data = dat, seed = 114
      )
    })
  )
}

# the parameters of a case whose formula has a predictor
reference_ce_pars <- function(case) {
  pars <- intersect(names(case$formula), names(case$model$parameters))
  pars[vapply(pars, function(p) length(rhs_vars(case$formula[[p]])) > 0, logical(1))]
}

# Fits live in an option, not in this file's environment: testthat sources the
# helpers again for every test_file() call, and the two test files are meant
# to share one fit per model within an R session.
reference_fit <- function(name) {
  cache <- getOption("bmm.local.reference_fits")
  if (is.null(cache)) {
    cache <- new.env()
    options(bmm.local.reference_fits = cache)
  }
  if (is.null(cache[[name]])) {
    case <- reference_cases()[[name]]
    cache[[name]] <- bmm(
      case$formula, case$data, case$model,
      backend = "cmdstanr", chains = 2, iter = 1000, warmup = 500,
      refresh = 0, silent = 2, seed = case$seed
    )
  }
  cache[[name]]
}
