# Recovery plan for tcc_cd(). The generators run the observer forward and
# share no code with bmm: each probe lies a whole number of degrees from the
# target (half of them at the target, the rest 1 to 180 degrees away), its
# similarity is read from the vector by that number, its familiarity is a
# normal draw with mean d * f and SD 1, and the observer answers "change" (1)
# when the familiarity falls short of the criterion. bmm computes
# P(change) = pnorm(criterion - d * f) instead, so an error in dtcc_cd(),
# rtcc_cd() or the Stan likelihood does not cancel.
tcc_cd_similarity <- (exp(-(0:180) / 20) - exp(-9)) / (1 - exp(-9))

# 1 = "change" for probes `degrees` away, one d and criterion per trial
tcc_cd_respond <- function(degrees, d, criterion) {
  familiarity <- stats::rnorm(length(degrees), mean = d * tcc_cd_similarity[degrees + 1], sd = 1)
  as.integer(familiarity < criterion)
}

# one subject's trials, for the recovery grid
tcc_cd_generator <- function(pars, n_trials, model) {
  at_target <- stats::runif(n_trials) < 0.5
  degrees <- ifelse(at_target, 0L, sample.int(180L, n_trials, replace = TRUE))
  sign <- sample(c(-1, 1), n_trials, replace = TRUE)
  out <- data.frame(tcc_cd_respond(degrees, pars$d, pars$criterion), sign * degrees * pi / 180)
  names(out) <- c(model$resp_vars$response, model$other_vars$distance)
  out
}

# the calibration's data sets: bmmtools hands over one prior draw under the
# fit's names and the design data (from tcc_cd_generator), and the responses
# are drawn again from the subject values of that draw
tcc_cd_sbc_generator <- function(draws, data) {
  subject <- function(par) draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]
  d <- exp(draws[["b_d_Intercept"]] + subject("d"))
  criterion <- draws[["b_criterion_Intercept"]] + subject("criterion")
  data$change <- tcc_cd_respond(round(abs(data$dist) * 180 / pi), d, criterion)
  data
}

bmmtools::validation_plan(
  bmm::tcc_cd(response = "change", distance = "dist", similarity = tcc_cd_similarity),
  # d on the log scale: e^0.8 = 2.2, away from the centre of the default
  # prior (0.5); the criterion at 1, near the unbiased observer (d / 2 = 1.1)
  pars = c(d = 0.8, criterion = 1),
  sds = c(d = 0.3, criterion = 0.4),
  grid = expand.grid(n_subjects = c(20, 50), n_trials = c(100, 300)),
  generator = tcc_cd_generator,
  # bmmtools 0.3.0 has no default calibration generator for tcc_cd. Cores per
  # fit, one per chain, are given here, because SBC's default reads
  # future::availableCores(), whose deprecation warning (future 1.82) fails
  # on the call bmmtools builds with rlang::exec(). bmmtools' default of 2
  # chains of 500 at adapt_delta 0.8 left 4 of 100 fits with Rhat above 1.05
  # and 30 with divergences
  sbc = list(
    generator = tcc_cd_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
