# Recovery plan for psychometric_function() with the default "normal" sigmoid
# on its default layout: one 0/1 response per trial, guess and lapse estimated
# (a yes/no detection task).
#
# The generator is the observer written out. On a trial the observer lapses
# with probability `lapse` and then answers "yes" at the guess rate; otherwise
# the stimulus is detected with probability F(x), a cumulative normal with
# mean `midpoint` and SD width / (2 * qnorm(0.95)), because bmm defines
# `width` as the distance from F = 0.05 to F = 0.95. A non-lapse trial that
# is not detected is also a guess. The sum of those branches is
# guess + (1 - guess) * (1 - lapse) * F(x), the model's psi, but bmm's closed
# form and its Stan likelihood are never used here.
psy_sigma <- function(width) width / (2 * stats::qnorm(0.95))

psy_respond <- function(intensity, midpoint, width, guess, lapse) {
  n <- length(intensity)
  lapsed <- stats::runif(n) < lapse
  detected <- stats::runif(n) < stats::pnorm(intensity, midpoint, psy_sigma(width))
  guessed <- stats::runif(n) < guess
  as.integer(ifelse(lapsed, guessed, detected | guessed))
}

# Five levels, each an equal share of the trials (30, 100 and 300 are
# multiples of 5). Two levels sit on the plateaus, so that guess (at 0) and
# lapse (at 1) are identified apart from the midpoint; three sit in the rise,
# at about F = 0.1, 0.5 and 0.9 for the truth below, because the slope is
# identified by how fast the rate climbs across them, not by the plateaus.
# Levels with one trial in five at each end cost trials in the rise, which is
# why there are not more of them.
psy_levels <- c(0, 0.2, 0.3, 0.4, 1)

psy_generator <- function(pars, n_trials, model, trial_design) {
  data.frame(y = psy_respond(
    trial_design$intensity, pars$midpoint, pars$width, pars$guess, pars$lapse
  ))
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names (midpoint identity, width log, guess and lapse logit) and the
# design data; the responses are drawn again from the subject values
psy_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]
  }
  data$y <- psy_respond(
    data$intensity, unname(subject("midpoint")), exp(unname(subject("width"))),
    stats::plogis(unname(subject("guess"))), stats::plogis(unname(subject("lapse")))
  )
  data
}

bmmtools::validation_plan(
  bmm::psychometric_function(response = "y", intensity = "intensity"),
  # on the link scale. The default priors follow the intensity range (0 to 1
  # here): midpoint normal(0.5, 0.5), log width normal(log 0.5, 1), guess and
  # lapse normal(-3, 1.3). Midpoint 0.3 is 0.4 SD below the prior centre;
  # width 0.25 (log -1.39) is 0.7 SD below it, a steep observer; guess 0.15
  # (logit -1.73) and lapse 0.10 (logit -2.20) are 1.0 and 0.6 SD above the
  # prior centre, both well above the typical 5% so they have to be recovered
  pars = c(midpoint = 0.3, width = log(0.25), guess = stats::qlogis(0.15), lapse = stats::qlogis(0.10)),
  # between-person SDs: midpoint 0.04 of a unit range keeps every subject's
  # rise between the levels at 0 and 1; the rest are moderate on their scales
  sds = c(midpoint = 0.04, width = 0.2, guess = 0.4, lapse = 0.4),
  trial_design = function(n_trials) {
    data.frame(intensity = rep(psy_levels, length.out = n_trials))
  },
  generator = psy_generator,
  sbc = list(
    generator = psy_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
