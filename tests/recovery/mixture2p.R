# Recovery plan for mixture2p() in its default version "simple" (Zhang & Luck,
# 2008): a response error comes from memory (a von Mises centred on the target
# with concentration kappa) with probability thetat, otherwise from guessing
# (uniform on the circle).
#
# How this differs from bmm: the generators draw the source of each trial first
# (memory with probability thetat, guess otherwise) and then the response, with
# a von Mises sampled by rejection from the uniform proposal, accepting x with
# probability exp(kappa * (cos(x) - 1)). bmm's rmixture2p() and its density
# (dmixture2p(), the Stan likelihood) are not called, and no von Mises routine
# of bmm, stats or circular is used, so a bug in bmm's density or sampler
# cannot cancel against the generator.

# von Mises(0, kappa) draws, kappa a vector with one value per draw
m2p_rvm <- function(kappa) {
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

# one response error per trial; thetat and kappa are vectors (one subject's
# value repeated, or one value per trial in the calibration)
m2p_respond <- function(thetat, kappa) {
  n <- length(kappa)
  memory <- stats::runif(n) < thetat
  y <- stats::runif(n, -pi, pi)
  y[memory] <- m2p_rvm(kappa[memory])
  y
}

# one subject's trials, for the recovery grid; thetat arrives as a
# probability and kappa on its natural scale
m2p_generator <- function(pars, n_trials, model) {
  out <- data.frame(m2p_respond(rep(pars$thetat, n_trials), rep(pars$kappa, n_trials)))
  names(out) <- model$resp_vars$resp_error
  out
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names (kappa on the log scale, thetat on the logit scale) and the
# design data; the responses are drawn again from the subject values of that draw
m2p_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    unname(draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))])
  }
  data$y <- m2p_respond(stats::plogis(subject("thetat")), exp(subject("kappa")))
  data
}

bmmtools::validation_plan(
  bmm::mixture2p(resp_error = "y"),
  # kappa on the log scale: e^3 = 20 (SD of the memory error 0.22 rad), 1 SD
  # above the centre of the default prior (normal(2, 1), kappa 7.4); thetat at
  # 0.8 (logit 1.39), 1.4 SD above the centre of logistic(0, 1) (0.5), a
  # typical set-size-3 memory rate
  pars = c(kappa = 3, thetat = stats::qlogis(0.8)),
  # between-person SDs on the link scale, 0.3 on log kappa and 0.6 on the
  # logit of thetat, moderate against the exponential(1) default prior
  sds = c(kappa = 0.3, thetat = 0.6),
  # bmmtools' standard grid: two weights in the mixture plus kappa are
  # identified from 30 to 300 trials per subject
  generator = m2p_generator,
  # prior = NULL puts bmm's own default priors back: bmmtools 0.3.1's entry
  # for mixture2p narrows them for an rmixture2p() overflow that this plan's
  # generator cannot have, and it names kappa and thetat as nlpars, which are
  # dpars in bmm 1.5.0, so that prior would not reach the fit
  sbc = list(
    prior = NULL, generator = m2p_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
