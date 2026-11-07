# Recovery plan for sdm_cd(), criterion fixed at 0. The generator runs the
# observer forward and shares no code with bmm: it draws the remembered
# feature from the SDM density by its own rejection sampler, takes the
# normalising constant from integrate() (bmm uses a Gauss-Legendre rule), and
# answers "change" when Lin and Oberauer's (2022) log-likelihood ratio of
# "change" against "same" at that feature exceeds the criterion, bounded to
# the trial's range of evidence as the model bounds it. It reads
# pars$criterion, so a free-criterion arm needs only a model formula and pars.
sdm_cd_generator <- function(pars, n_trials, model) {
  peak <- pars$c * sqrt(pars$kappa / (2 * pi))
  # the activation minus its peak, the log density up to log Z - peak
  excess <- function(u) peak * expm1(pars$kappa * (cos(u) - 1))
  draw <- function(n) {
    out <- numeric(0)
    while (length(out) < n) {
      x <- stats::runif(10 * n, -pi, pi)
      out <- c(out, x[log(stats::runif(10 * n)) < excess(x)])
    }
    out[seq_len(n)]
  }
  width <- 1 / sqrt(pars$kappa * max(peak, 1))
  breaks <- sort(unique(c(0, pmin(pi, width * c(1, 4, 16, 64)), pi)))
  log_z <- peak + log(2 * sum(vapply(seq_len(length(breaks) - 1), function(i) {
    stats::integrate(function(u) exp(excess(u)), breaks[i], breaks[i + 1],
                     rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value
  }, numeric(1))))

  target <- stats::runif(n_trials, -pi, pi)
  change <- stats::runif(n_trials) < 0.5
  probe <- ifelse(change, stats::runif(n_trials, -pi, pi), target)
  x <- target + pars$mu + draw(n_trials)

  # log of 2 pi times the observer's density at distance d from the probe
  log_2pi_f <- function(d) log(2 * pi) + peak + excess(d) - log_z
  llr <- -log_2pi_f(x - probe)
  llr_max <- -log_2pi_f(pi)
  llr_min <- -log_2pi_f(0)
  criterion <- pars$criterion
  bounded <- if (criterion > 0) {
    llr_max * tanh(criterion / llr_max)
  } else if (criterion < 0) {
    -llr_min * tanh(criterion / -llr_min)
  } else {
    0
  }

  out <- data.frame(as.integer(llr > bounded), probe, target)
  names(out) <- unlist(model$resp_vars[c("response", "probe", "target")], use.names = FALSE)
  out
}

bmmtools::validation_plan(
  bmm::sdm_cd(response = "change", probe = "probe", target = "target"),
  # log scale: c e^1.2 = 3.3 and kappa e^2.5 = 12.2, away from the centres of
  # the default priors (log c 2, log kappa 1.75)
  pars = c(c = 1.2, kappa = 2.5),
  sds = c(c = 0.3, kappa = 0.3),
  grid = expand.grid(n_subjects = c(20, 40), n_trials = c(200, 500)),
  generator = sdm_cd_generator,
  seed = 2026
)
