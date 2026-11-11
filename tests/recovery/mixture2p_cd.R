# Recovery plan for mixture2p_cd(), version "simple" with the criterion fixed
# at 0. The generator runs the observer forward and shares no code with bmm:
# it draws the remembered feature (its own von Mises sampler), computes Lin and
# Oberauer's (2022) log-likelihood ratio of "change" against "same" at that
# feature, and answers "change" when it exceeds the criterion, bounded to the
# trial's range of evidence as the model bounds it. It handles every version
# of the model, so the plan can be rebuilt for "slot" or "slot_averaging" by
# changing the model and pars below.
mixture2p_cd_generator <- function(pars, n_trials, model) {
  rvm <- function(n, kappa) {
    out <- numeric(0)
    while (length(out) < n) {
      x <- stats::runif(n, -pi, pi)
      out <- c(out, x[stats::runif(n) < exp(kappa * (cos(x) - 1))])
    }
    out[seq_len(n)]
  }
  # the concentration whose Fisher information kappa I1 / I0 is `slots`
  # times that of one slot of precision kappa
  slot_kappa <- function(kappa, slots) {
    fisher <- function(k) k * besselI(k, 1, TRUE) / besselI(k, 0, TRUE)
    if (slots <= 1) return(kappa)
    stats::uniroot(function(k) fisher(k) - slots * fisher(kappa), c(kappa, slots * kappa + 1),
                   tol = 1e-12)$root
  }
  vm_2pi <- function(cos_d, kappa) exp(kappa * (cos_d - 1)) / besselI(kappa, 0, TRUE)

  version <- model$version
  target <- stats::runif(n_trials, -pi, pi)
  change <- stats::runif(n_trials) < 0.5
  probe <- ifelse(change, stats::runif(n_trials, -pi, pi), target)
  set_size <- if (version == "simple") rep(1, n_trials) else sample(c(2, 3, 4, 6), n_trials, TRUE)

  # the observer's precision states: kappa_lo with weight 1 - w_hi, kappa_hi
  # with w_hi; p_s is its prior that the probed item is stored
  if (version == "slot_averaging") {
    share <- pars$K / set_size
    whole <- floor(share)
    extra <- share - whole
    slots <- whole + (stats::runif(n_trials) < extra)
    kappa_lo <- vapply(whole, function(s) slot_kappa(pars$kappa, s), 0)
    kappa_hi <- vapply(whole + 1, function(s) slot_kappa(pars$kappa, s), 0)
    stored <- slots >= 1
    kappa_x <- ifelse(slots == whole, kappa_lo, kappa_hi)
    p_s <- ifelse(whole >= 1, 1, extra)
    w_hi <- if (identical(model$knowledge, "rich")) as.numeric(slots > whole) else extra
    w_hi[whole < 1] <- 1
  } else {
    p_s <- if (version == "slot") pmin(1, pars$K / set_size) else rep(pars$thetat, n_trials)
    stored <- stats::runif(n_trials) < p_s
    kappa_lo <- kappa_hi <- kappa_x <- rep(pars$kappa, n_trials)
    w_hi <- rep(0, n_trials)
  }
  x <- stats::runif(n_trials, -pi, pi)
  x[stored] <- target[stored] + pars$mu + vapply(kappa_x[stored], function(k) rvm(1, k), 0)

  density_2pi <- function(cos_d) {
    hi <- w_hi * vm_2pi(cos_d, kappa_hi)
    lo <- ifelse(w_hi < 1, (1 - w_hi) * vm_2pi(cos_d, kappa_lo), 0)
    hi + lo
  }
  llr <- -log(1 - p_s + p_s * density_2pi(cos(x - probe)))
  llr_max <- -log(1 - p_s + p_s * density_2pi(-1))
  llr_min <- -log(1 - p_s + p_s * density_2pi(1))
  criterion <- pars$criterion
  bounded <- if (criterion > 0) {
    llr_max * tanh(criterion / llr_max)
  } else if (criterion < 0) {
    -llr_min * tanh(criterion / -llr_min)
  } else {
    0
  }

  # the model reads the probe relative to the target, in [-pi, pi)
  out <- data.frame(as.integer(llr > bounded), (probe - target + pi) %% (2 * pi) - pi)
  names(out) <- unlist(model$resp_vars[c("response", "probe")], use.names = FALSE)
  if (version != "simple") out[[model$other_vars$set_size]] <- set_size
  out
}

bmmtools::validation_plan(
  bmm::mixture2p_cd(response = "change", probe = "probe"),
  # link scale: kappa e^1.4 = 4.1 and P(memory) 0.73, both away from the
  # centres of the default priors (log kappa 2, logit thetat 0)
  pars = c(kappa = 1.4, thetat = 1),
  sds = c(kappa = 0.3, thetat = 0.6),
  grid = expand.grid(n_subjects = c(20, 40), n_trials = c(200, 500)),
  generator = mixture2p_cd_generator,
  seed = 2026
)
