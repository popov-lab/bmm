# Recovery plan for imm() in its version "full" (Oberauer & Lin, 2017):
# parameters kappa, a, c and s are estimated; mu is fixed at 0 and the
# background activation b at 1 (log b = 0), which is what makes the others
# identified. A trial of set size k has the target and k - 1 lures with
# distances d_i to the cue. The activation of the target is a + c, that of lure
# i is a + c * exp(-s * d_i) and the background is b; a response is the
# normalised choice among those sources (Luce), then a von Mises around the
# chosen item or, for the background, the uniform.
#
# How this differs from bmm: each trial's source is drawn from those weights
# and then the response from a von Mises sampled by rejection (accept x with
# probability exp(kappa * (cos(x) - 1))); bmm's rimm(), dimm() and the Stan
# likelihood are not called, nor any von Mises routine of bmm, stats or
# circular. The weights are taken from the activation equations in bmm's IMM
# article (A = b * A_b + a * A_a + c * A_c), with b the log-0 background.
#
# Design: set sizes 2 to 5 (an equal share each; four lure columns), lure
# locations uniform on the circle and lure distances uniform on (0, pi). The
# spread of the distances is what separates a (flat in distance) from c (falling
# with distance at rate s), and the several set sizes put different numbers of
# lures against the one target. Trial numbers are higher than bmmtools'
# standard grid because four parameters are read from the lure responses
# (about a quarter of the trials).
imm_rvm <- function(kappa) {
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

imm_design <- function(n_trials) {
  set_size <- sample(rep_len(2:5, n_trials))
  lures <- matrix(stats::runif(n_trials * 4L, -pi, pi), n_trials, 4L)
  dist <- matrix(stats::runif(n_trials * 4L, 0, pi), n_trials, 4L)
  out <- data.frame(set_size = set_size, lures, dist)
  names(out)[2:9] <- c(paste0("nt", 1:4), paste0("dist", 1:4))
  out
}

# One response error per trial; kappa, a, c and s are vectors with one value
# per trial, the background is 1
imm_respond <- function(kappa, a, c, s, set_size, lures, dist) {
  n <- length(kappa)
  centre <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    n_lure <- set_size[t] - 1L
    weight <- c(a[t] + c[t], a[t] + c[t] * exp(-s[t] * dist[t, seq_len(n_lure)]), 1)
    pick <- sample.int(length(weight), 1L, prob = weight)
    centre[t] <- if (pick == 1L) 0 else if (pick <= n_lure + 1L) lures[t, pick - 1L] else NA_real_
  }
  memory <- !is.na(centre)
  y <- stats::runif(n, -pi, pi)
  y[memory] <- centre[memory] + imm_rvm(kappa[memory])
  (y + pi) %% (2 * pi) - pi
}

# bmmtools hands pars over on the natural scale (b = 1, mu = 0 are among them)
imm_generator <- function(pars, n_trials, model, trial_design) {
  n <- n_trials
  out <- data.frame(imm_respond(
    rep(pars$kappa, n), rep(pars$a, n), rep(pars$c, n), rep(pars$s, n),
    trial_design$set_size, as.matrix(trial_design[paste0("nt", 1:4)]),
    as.matrix(trial_design[paste0("dist", 1:4)])
  ))
  names(out) <- model$resp_vars$resp_error
  out
}

# the calibration's data sets: one prior draw under the fit's names (every
# parameter on its log scale) and the design data, whose set sizes, lure
# locations and distances are kept; the responses are drawn again
imm_sbc_generator <- function(draws, data) {
  subject <- function(par) {
    exp(unname(draws[[sprintf("b_%s_Intercept", par)]] +
      draws[sprintf("r_id__%s[%s,Intercept]", par, as.character(data$id))]))
  }
  data$y <- imm_respond(
    subject("kappa"), subject("a"), subject("c"), subject("s"), data$set_size,
    as.matrix(data[paste0("nt", 1:4)]), as.matrix(data[paste0("dist", 1:4)])
  )
  data
}

bmmtools::validation_plan(
  bmm::imm(
    resp_error = "y", nt_features = paste0("nt", 1:4),
    nt_distances = paste0("dist", 1:4), set_size = "set_size", version = "full"
  ),
  # all on the log scale, against default priors that are normal(0, 1) for a,
  # c and s and normal(2, 1) for kappa: kappa 3 (kappa 20) is 1 SD above the
  # centre; c 1.2 (c 3.3) 1.2 SD above, a -1.2 (a 0.30) 1.2 SD below, and s 0.7
  # (s 2) 0.7 SD above. At set size 4 (target a + c = 3.6, three lures about
  # 0.3 + 0.5 each, background 1) a response comes from the target with
  # probability 0.5, from a lure 0.34 and from the background 0.14
  pars = c(kappa = 3, a = -1.2, c = 1.2, s = 0.7),
  sds = c(kappa = 0.3, a = 0.3, c = 0.3, s = 0.3),
  grid = expand.grid(n_subjects = c(20L, 50L), n_trials = c(100L, 300L, 600L)),
  trial_design = imm_design,
  generator = imm_generator,
  sbc = list(
    prior = NULL, n_trials = 150L,
    generator = imm_sbc_generator, cores_per_fit = 4L,
    sampler = list(chains = 4L, iter = 1000L, warmup = 500L, control = list(adapt_delta = 0.95))
  ),
  seed = 2026
)
