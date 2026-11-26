# Recovery plan for sdt_ranking() with the gumbel_min default on one row per
# trial (ADR 0006). The generators run the observer forward and share no code
# with bmm: on each trial the target's evidence is d plus a smallest-extreme-value
# draw, each of the m - 1 lures' evidence is such a draw, the items are ranked
# from the largest evidence (rank 1 = most likely target), and the response is
# the target's rank. The draw is the inverse CDF of F(x) = 1 - exp(-exp(x)),
# x = log(-log(1 - u)) with u uniform. bmm computes the rank probabilities in
# closed form from gamma-function ratios, so an error in dsdt_ranking(),
# rsdt_ranking() or the Stan likelihood does not cancel.
sdt_ranking_m <- 4L

gumbel_min_draw <- function(n) {
  log(-log(1 - stats::runif(n)))
}

# the target's rank (1 to m); one d per trial
sdt_ranking_respond <- function(d, m = sdt_ranking_m) {
  n <- length(d)
  target <- d + gumbel_min_draw(n)
  lures <- matrix(gumbel_min_draw(n * (m - 1L)), nrow = n)
  1L + as.integer(rowSums(lures > target))
}

# one subject's trials; no design variable (each trial ranks one target among
# m - 1 lures, and m is constant)
sdt_ranking_generator <- function(pars, n_trials, model) {
  out <- data.frame(sdt_ranking_respond(rep(pars$d, n_trials), model$other_vars$m))
  names(out) <- model$resp_vars$response
  out
}

# the calibration's data sets: one prior draw under the fit's names, the
# ranks drawn again from the subject values of that draw
sdt_ranking_sbc_generator <- function(draws, data) {
  d <- draws[["b_d_Intercept"]] +
    unname(draws[sprintf("r_id__d[%s,Intercept]", as.character(data$id))])
  data$rank <- sdt_ranking_respond(d)
  data
}

bmmtools::validation_plan(
  bmm::sdt_ranking(response = "rank", m = sdt_ranking_m, dist = "gumbel_min"),
  # d = 1.6 sits 0.6 prior SDs above the prior centre (1); m = 4 so that the
  # rank carries more than the correct / incorrect bit that m = 2 would
  pars = c(d = 1.6),
  # below the prior's mean SD (exponential(1) has mean 1)
  sds = c(d = 0.4),
  # bmmtools' default grid
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(30L, 100L, 300L)),
  generator = sdt_ranking_generator,
  sbc = list(
    generator = sdt_ranking_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
