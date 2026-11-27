# Recovery plan for mixture3p() in its version "simple": a response error comes
# from the target (probability e^thetat / Z), from one of the lures (each lure
# e^(thetant - log(set_size - 1)) / Z, so the lures share e^thetant / Z equally)
# or from guessing (e^0 / Z), with Z the sum of those weights; a set size of 1
# has no lure. thetat and thetant are on the softmax link, so they are log
# weights against the guessing weight 0.
#
# How this differs from bmm: each trial's source (target, one lure, guess) is
# drawn from those weights, then the response from a von Mises sampled by
# rejection (accept x with probability exp(kappa * (cos(x) - 1))) around that
# source's location, or from the uniform for a guess. bmm's rmixture3p(),
# dmixture3p() and the Stan likelihood are not called, nor any von Mises
# routine of bmm, stats or circular.
#
# Design: the set size varies from 1 to 4 over the trials (an equal share
# each), with three lure columns whose locations are uniform on the circle and
# are only read up to the trial's set size. A constant set size could not tell
# a likelihood that shares thetant by set_size - 1 from one that shares it by
# another count, and set size 1 is the trial without lures; both only show up
# when the set size varies. Trial numbers are higher than bmmtools' standard
# grid because a trial holds three sources, and 30 trials leave thetant to
# about 3 trials per subject that report a lure.
m3p_rvm <- function(kappa) {
  out <- numeric(length(kappa))
  open <- seq_along(out)
  while (length(open) > 0L) {
    x <- stats::runif(length(open), -pi, pi)
    keep <- stats::runif(length(open)) < exp(kappa[open] * (cos(x) - 1))
    out[open[keep]] <- x[keep]
    open <- open[!keep]
  }
  out
}

# the trials' lure columns and set sizes; a function of the trial count, which
# bmmtools calls once per subject
m3p_design <- function(n_trials) {
  set_size <- sample(rep_len(1:4, n_trials))
  lures <- matrix(stats::runif(n_trials * 3L, -pi, pi), n_trials, 3L)
  out <- data.frame(set_size = set_size, lures)
  names(out)[2:4] <- paste0("nt", 1:3)
  out
}

# One response error per trial. thetat, thetant and kappa are vectors with one
# value per trial; target at 0, lures in the columns nt1 to nt3
m3p_respond <- function(thetat, thetant, kappa, set_size, lures) {
  n <- length(kappa)
  centre <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    n_lure <- set_size[t] - 1L
    weight <- c(exp(thetat[t]), rep(exp(thetant[t] - log(n_lure)), n_lure), 1)
    pick <- sample.int(length(weight), 1L, prob = weight)
    centre[t] <- if (pick == 1L) 0 else if (pick <= n_lure + 1L) lures[t, pick - 1L] else NA_real_
  }
  memory <- !is.na(centre)
  y <- stats::runif(n, -pi, pi)
  y[memory] <- centre[memory] + m3p_rvm(kappa[memory])
  (y + pi) %% (2 * pi) - pi
}

# bmmtools hands pars over with kappa natural and thetat, thetant as log weights
m3p_generator <- function(pars, n_trials, model, trial_design) {
  n <- n_trials
  out <- data.frame(m3p_respond(
    rep(pars$thetat, n), rep(pars$thetant, n), rep(pars$kappa, n),
    trial_design$set_size, as.matrix(trial_design[paste0("nt", 1:3)])
  ))
  names(out) <- model$resp_vars$resp_error
  out
}

# the calibration's data sets: one prior draw under the fit's names (kappa log,
# thetat and thetant log weights) and the design data, whose set sizes and lure
# columns are kept; the responses are drawn again from the subject values
m3p_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    unname(draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))])
  }
  data$y <- m3p_respond(
    subject("thetat"), subject("thetant"), exp(subject("kappa")),
    data$set_size, as.matrix(data[paste0("nt", 1:3)])
  )
  data
}

bmmtools::validation_plan(
  bmm::mixture3p(
    resp_error = "y", nt_features = paste0("nt", 1:3), set_size = "set_size"
  ),
  # log kappa 3 (kappa 20), 1 SD above the centre of the default prior
  # normal(2, 1). thetat 1.2 sits 1.2 logistic units above the prior centre 0:
  # at set size 4 (one target, three lures, guess) the target weight is
  # e^1.2 = 3.3 against the lures' 0.8 in total and the guess 1, so a target
  # response has probability 0.65, lure 0.16, guess 0.19. thetant -0.2, 0.2
  # below its prior centre 0, keeps lure responses frequent enough to count
  pars = c(kappa = 3, thetat = 1.2, thetant = -0.2),
  sds = c(kappa = 0.3, thetat = 0.5, thetant = 0.4),
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(60L, 150L, 300L)),
  trial_design = m3p_design,
  generator = m3p_generator,
  sbc = list(
    prior = NULL, n_trials = 100L,
    generator = m3p_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
