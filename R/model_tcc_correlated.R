############################################################################# !
# TCC WITH CORRELATED AND MOTOR NOISE: THE TABLE FROM R                  ####
############################################################################# !

# The correlated-noise TCC (Schurgin et al., 2020): signals m = d f + z with
# z ~ N(0, C), C circulant over the n options, the response the argmax of m,
# then optional motor noise round(x + e), e ~ N(0, sd^2) in option steps. Its
# response probabilities have no closed form, so the likelihood table that
# tcc_lpmf() reads is a Monte Carlo estimate, built here once per fit and
# passed to Stan as data (table = "r"). Options are indexed 0, ..., n - 1 with
# the target at 0, so a rotation of a draw is a cyclic shift of its rows.
#
# The estimator conditions on the drift (#566, correlated-noise spike v2).
# Dropping the constant mode, which cannot move an argmax, and with f_R = f
# centred, z = Delta f_R + z_perp, Delta = b'z with b = C^+ f_R / kappa^2 and
# kappa^2 = f_R' C^+ f_R, so Delta ~ N(0, 1 / kappa^2) independent of z_perp.
# Given z_perp, the winner at s = d + Delta is the upper envelope of
# z_perp,x + s f_R,x, made of segments [lo, hi) with one winner each, and an
# option's segment adds Phi(kappa (hi - d)) - Phi(kappa (lo - d)) to its P at
# d. The endpoints are binned on an s grid, and one matrix product with
# Phi(kappa (s_b - d)) gives the table at every d at once.

# the estimator's version enters the table seed, so that a changed engine
# never reuses a table built by an older one
.tcc_corr_version <- "drift-1"

# draws are taken in groups: one draw of z and its rotations, each rotation an
# exact draw again because C is circulant. Six rotations cost about the same
# per certified trial as none but cut the peak memory from about 1 GB to
# 270 MB at 24,000 envelopes per chunk (spike v2, Q1-Q2)
.tcc_corr_rotations <- 6L

# A value by distance 0, ..., n / 2 as the first row of a circulant matrix
.tcc_circ_row <- function(by_distance) {
  n <- 2 * (length(by_distance) - 1)
  lag <- 0:(n - 1)
  by_distance[pmin(lag, n - lag) + 1]
}

# Everything the build needs from the correlation and the similarity by
# distance. Eigenvalues below 1e-10 of the largest count as null modes. The
# similarity's part on null modes is dropped (the estimator conditions on
# f_R alone), and `dropped` is reported for tcc() to refuse a correlation
# that would lose more than 1e-6 of the similarity's energy
.tcc_corr_setup <- function(correlation, similarity) {
  n <- 2 * (length(similarity) - 1)
  lam <- Re(stats::fft(.tcc_circ_row(correlation)))
  keep <- lam > 1e-10 * max(lam)
  keep[1] <- FALSE
  f <- .tcc_circ_row(similarity)
  f_c <- f - mean(f)
  f_R <- if (sum(!keep) > 1) Re(stats::fft(stats::fft(f_c) * keep, inverse = TRUE)) / n else f_c
  inv_lam <- ifelse(keep, 1 / pmax(lam, 1e-300), 0)
  Cp_f <- Re(stats::fft(stats::fft(f_R) * inv_lam, inverse = TRUE)) / n
  kappa2 <- sum(f_R * Cp_f)
  list(
    n = n, f_R = f_R, b = Cp_f / kappa2, kappa = sqrt(kappa2),
    sqrt_lam = sqrt(pmax(lam, 0)) * keep,
    dropped = sum((f_c - f_R)^2) / sum(f_c^2)
  )
}

# m draws of z (n x m) without the constant mode, by the circulant square root
.tcc_corr_draw <- function(m, setup) {
  G <- matrix(stats::rnorm(setup$n * m), setup$n, m)
  Re(stats::mvfft(stats::mvfft(G) * setup$sqrt_lam, inverse = TRUE)) / setup$n
}

# z (n x m) -> z_perp (m x n), the draws rotated by `shift` options first
.tcc_corr_split <- function(Z, setup, shift = 0) {
  if (shift) {
    Z <- Z[((seq_len(setup$n) - 1 + shift) %% setup$n) + 1, , drop = FALSE]
  }
  t(Z) - outer(drop(crossprod(setup$b, Z)), setup$f_R)
}

# The upper envelope of W[i, k] + t f_k over t >= 0, for all rows at once: from
# the winner at t = 0, each step moves to the option that overtakes it first,
# until the winner is an option with the largest f, which wins for ever
.tcc_corr_walk <- function(W, f) {
  N <- nrow(W)
  fmax <- max(f)
  j <- max.col(W, ties.method = "first")
  t_cur <- numeric(N)
  segs <- list()
  active <- which(f[j] < fmax)
  while (length(active)) {
    Wa <- W[active, , drop = FALSE]
    ja <- j[active]
    den <- outer(-f[ja], f, "+")
    Tm <- (Wa[cbind(seq_along(active), ja)] - Wa) / den
    Tm[den <= 0] <- Inf
    k <- max.col(-Tm, ties.method = "first")
    tk <- pmax(Tm[cbind(seq_along(active), k)], t_cur[active])
    segs[[length(segs) + 1]] <- list(draw = active, opt = ja, lo = t_cur[active], hi = tk)
    j[active] <- k
    t_cur[active] <- tk
    active <- active[f[k] < fmax]
  }
  segs[[length(segs) + 1]] <- list(draw = seq_len(N), opt = j, lo = t_cur, hi = rep(Inf, N))
  lapply(c(draw = "draw", opt = "opt", lo = "lo", hi = "hi"),
         function(v) unlist(lapply(segs, `[[`, v)))
}

# the segments over all s of z_perp + s f_R: one walk up along +f_R for s >= 0
# and one along -f_R for s <= 0
.tcc_corr_envelope <- function(Rho, f_R) {
  up <- .tcc_corr_walk(Rho, f_R)
  down <- .tcc_corr_walk(Rho, -f_R)
  list(draw = c(up$draw, down$draw), opt = c(up$opt, down$opt),
       lo = c(up$lo, -down$hi), hi = c(up$hi, -down$lo))
}

# A step of 0.01 keeps the binning error below 1.2e-4 relative down to P =
# 1e-7, 10 to 100 times below the Monte Carlo SE; 8 / kappa of padding puts
# Phi at 0 or 1 to 6e-16 beyond the grid (spike v2, Q1)
.tcc_corr_sgrid <- function(d_values, kappa, ds = 0.01, pad = 8) {
  lo <- min(d_values) - pad / kappa
  nb <- ceiling((max(d_values) + pad / kappa - lo) / ds) + 1
  list(lo = lo, ds = ds, nb = nb, s = lo + (seq_len(nb) - 1) * ds)
}

.tcc_corr_wsum <- function(idx, w, len) {
  out <- numeric(len)
  if (length(idx)) {
    s <- rowsum(w, idx, reorder = FALSE)
    out[as.integer(rownames(s))] <- s[, 1]
  }
  out
}

# Adds a batch of segments: + at hi, - at lo, split linearly between the two
# nearest s bins. An endpoint above the s grid has Phi = 1 at every d of the
# table, one below it Phi = 0
.tcc_corr_accumulate <- function(acc, env, sg) {
  n <- length(acc$above)
  s <- c(env$hi, env$lo)
  o <- c(env$opt, env$opt)
  sgn <- rep(c(1, -1), each = length(env$opt))
  pos <- (s - sg$lo) / sg$ds
  top <- pos > sg$nb - 1
  acc$above <- acc$above + .tcc_corr_wsum(o[top], sgn[top], n)
  keep <- pos >= 0 & !top
  pos <- pos[keep]
  o <- o[keep]
  sgn <- sgn[keep]
  b0 <- pmin(floor(pos), sg$nb - 2)
  wt <- pos - b0
  acc$bins <- acc$bins +
    .tcc_corr_wsum(c(b0 * n + o, (b0 + 1) * n + o), c(sgn * (1 - wt), sgn * wt), n * sg$nb)
  acc
}

# P of each option (n x length(d_values)) from `draws` envelopes, rounded up to
# whole groups of rotations, in chunks of 24,000 envelopes. Draws from the
# current RNG stream (.tcc_noise_table() sets the seed). For the certification
# it returns the segments of every chunk that starts before group
# `keep_groups`: whole chunks, so at least `keep_groups` groups and up to a
# chunk more (12,000 at the defaults, from 3 chunks of 4,000 groups)
.tcc_corr_build <- function(setup, d_values, draws, keep_groups = 1e4, chunk = 24000) {
  n <- setup$n
  n_rot <- .tcc_corr_rotations
  shifts <- floor((seq_len(n_rot) - 1) * n / n_rot)
  sg <- .tcc_corr_sgrid(d_values, setup$kappa)
  groups <- ceiling(draws / n_rot)
  acc <- list(bins = numeric(n * sg$nb), above = numeric(n))
  kept <- list()
  done <- 0
  while (done < groups) {
    m <- min(chunk / n_rot, groups - done)
    Z <- .tcc_corr_draw(m, setup)
    for (shift in shifts) {
      env <- .tcc_corr_envelope(.tcc_corr_split(Z, setup, shift), setup$f_R)
      acc <- .tcc_corr_accumulate(acc, env, sg)
      if (done < keep_groups) {
        env$draw <- env$draw + done
        kept[[length(kept) + 1]] <- env
      }
    }
    done <- done + m
  }
  envelopes <- groups * n_rot
  H <- stats::pnorm(outer(sg$s, d_values, "-") * setup$kappa)
  P <- (matrix(acc$bins, n) %*% H + acc$above) / envelopes
  kept <- lapply(c(draw = "draw", opt = "opt", lo = "lo", hi = "hi"),
                 function(v) unlist(lapply(kept, `[[`, v)))
  nlist(P, envelopes, kept)
}

# Motor noise as the generator has it: the chosen position plus e ~ N(0,
# sd^2), rounded to the nearest option, so K(k) = Phi((k + 1/2) / sd) -
# Phi((k - 1/2) / sd), out to 6 SD. The authors' fitting code samples
# dnorm(-6:6, 0, 2) instead, which differs by up to 0.15 in log P at large d
# and is 0 beyond 6 steps (spike v2, Q5)
.tcc_motor_kernel <- function(sd_steps) {
  k <- seq(-ceiling(6 * sd_steps), ceiling(6 * sd_steps))
  w <- stats::pnorm((k + 0.5) / sd_steps) - stats::pnorm((k - 0.5) / sd_steps)
  w / sum(w)
}

# Circular convolution over the options (rows) of log P by direct summation:
# every term is positive, so a P of 1e-300 keeps its relative accuracy, where
# an FFT would bury it under 1e-17 of round-off; in logs, so that an exact
# table at a large d_max does not underflow. A cell that is -Inf at every tap
# stays -Inf
.tcc_convolve_log <- function(L, kernel) {
  n <- nrow(L)
  half <- (length(kernel) - 1) / 2
  shifted <- function(s) L[((seq_len(n) - 1 - s) %% n) + 1, , drop = FALSE]
  # one shifted copy at a time: a list of all of them takes 4 MB per tap at
  # d_max = 30, over 1 GB for a motor SD of 20 options
  top <- L
  for (s in -half:half) top <- pmax(top, shifted(s))
  safe <- ifelse(is.finite(top), top, 0)
  total <- 0
  for (s in -half:half) total <- total + kernel[s + half + 1] * exp(shifted(s) - safe)
  ifelse(is.finite(top), safe + log(total), -Inf)
}

# log P of each option (n x G) -> one option at each distance 0, ..., n / 2,
# the two mirror options pooled (their P is equal for a symmetric C and f, and
# pooling halves the Monte Carlo variance)
.tcc_by_distance <- function(L) {
  n <- nrow(L)
  k <- 0:(n / 2)
  a <- L[k + 1, , drop = FALSE]
  b <- L[((n - k) %% n) + 1, , drop = FALSE]
  top <- pmax(a, b)
  ifelse(is.finite(top), top + log((exp(a - top) + exp(b - top)) / 2), -Inf)
}

# The table's seed is a function of everything that defines the table, so the
# same model always gets the same table, whatever the user's RNG state
.tcc_corr_seed <- function(key) {
  strtoi(substr(rlang::hash(key), 1, 7), 16L)
}

.tcc_motor_steps <- function(other_vars) {
  motor_sd <- other_vars$motor_sd %||% 0
  if (identical(other_vars$units, "bins")) motor_sd else motor_sd * other_vars$n_options / (2 * pi)
}

# The table that tcc_lpmf() reads when it comes from R (table = "r"): rows the
# d grid, columns the distances, log P. Without correlation or motor noise it
# is .tcc_table(). Built once per R process and inputs: update(), loo() refits
# and kfold() in that process reuse it, and so does postprocess_brm(), which
# reads `info`; parallel future workers each build their own
.tcc_noise_table <- function(model, grid) {
  ov <- model$other_vars
  similarity <- .tcc_similarity(model)
  motor <- .tcc_motor_steps(ov)
  if (is.null(ov$correlation) && motor == 0) {
    return(list(tab = .tcc_table(similarity, grid), info = NULL))
  }
  key <- list(as.numeric(ov$correlation), as.numeric(similarity), ov$correlation_draws,
              motor, grid, .tcc_corr_version)
  hash <- rlang::hash(key)
  cache <- .tcc_cache$noise %||% list()
  if (is.null(cache[[hash]])) {
    built <- .tcc_build_noise_table(ov$correlation, similarity, ov$correlation_draws, motor,
                                    grid, .tcc_corr_seed(key))
    # tables are 2 MB at d_max = 30; a few cover the fits of one analysis
    .tcc_cache$noise <- utils::tail(c(cache, stats::setNames(list(built), hash)), 4)
    return(built)
  }
  cache[[hash]]
}

# `seed` is not a user's seed: .tcc_noise_table() derives it from the hash of
# the table's inputs, so that one model always gets one table
.tcc_build_noise_table <- function(correlation, similarity, draws, motor, grid, seed) {
  d_values <- grid$lo + (seq_len(grid$rows) - 2) * grid$step
  info <- NULL
  if (is.null(correlation)) {
    n <- 2 * (length(similarity) - 1)
    L <- vapply(d_values, .tcc_logp, numeric(length(similarity)), similarity = similarity)
    L <- L[pmin(0:(n - 1), n - 0:(n - 1)) + 1, , drop = FALSE]
  } else {
    setup <- .tcc_corr_setup(correlation, similarity)
    # the engine draws under its own RNG kinds, so the table does not depend
    # on the user's RNGkind(), and the user's stream is left as it was
    build <- withr::with_seed(
      seed, .tcc_corr_build(setup, d_values, draws),
      .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
      .rng_sample_kind = "Rejection"
    )
    # the bins carry signed sums, so a P that is 0 can come out a few ulp below
    L <- log(pmax(build$P, 0))
  }
  kernel <- NULL
  if (motor > 0) {
    kernel <- .tcc_motor_kernel(motor)
    L <- .tcc_convolve_log(L, kernel)
  }
  L <- .tcc_by_distance(L)
  if (!is.null(correlation)) {
    # A cell no envelope reached is not impossible, only rarer than the draws
    # can show. Flooring it at half a draw keeps the likelihood finite; only
    # weak correlations leave such cells inside the grid (spike v2, Q2)
    zero <- !is.finite(L)
    L[zero] <- log(0.5 / build$envelopes)
    info <- c(
      list(draws = build$envelopes, seed = seed, kappa = setup$kappa,
           zeros = sum(zero),
           zeros_from = if (any(zero)) min(d_values[col(L)[zero]]) else NA_real_),
      .tcc_corr_certify(build$kept, setup, L, d_values, grid, kernel, build$envelopes)
    )
  }
  tab <- t(L)
  tab[grid$rows, ] <- pmin(tab[grid$rows, ], tab[grid$rows - 2, ])
  nlist(tab, info)
}

# How many trials one estimate may pool before the table's own Monte Carlo
# error matters. To first order the shift a table P-hat induces in an
# estimate of d, in units of its posterior SD sqrt(1 / (n_pool I(d))), is
# sqrt(n_pool) delta(d) / sqrt(I(d)) with delta(d) = -sum_k mult_k s_k(d)
# P-hat_k(d), minus the score's mean under the table, which is linear in
# P-hat. Each group of draws therefore contributes psi = -sum over its
# segments of g_opt(d) c_seg(d), with g the score of the option's distance
# convolved back through the motor kernel, and SD(delta) = SD(psi) /
# sqrt(rotations x envelopes). The certified n_pool is the largest at which
# that SD stays below 0.1 posterior SD at every d from 0.3 to 12 (or d_max):
# a property of the build method, which the table in hand meets on average.
# The scores come from the table itself. At 1e6 draws with the authors' C
# this gave 12,996 trials, and the spread of four independent 4e6 builds was
# 0.98 times the prediction (spike v2, Q2)
.tcc_corr_certify <- function(kept, setup, L, d_values, grid, kernel, envelopes) {
  d_hi <- min(12, max(d_values) - 2 * grid$step)
  at <- unique(round((seq(min(0.3, d_hi / 2), d_hi, length.out = 25) - grid$lo) / grid$step) + 2)
  at <- at[at > 1 & at < length(d_values)]
  D <- exp(L)
  S <- (L[, at + 1, drop = FALSE] - L[, at - 1, drop = FALSE]) / (2 * grid$step)
  info_per_trial <- colSums(.tcc_mult(nrow(L)) * D[, at, drop = FALSE] * S^2)
  n <- setup$n
  score <- S[pmin(0:(n - 1), n - 0:(n - 1)) + 1, , drop = FALSE]
  if (!is.null(kernel)) {
    score <- Reduce(`+`, lapply(seq_along(kernel), function(i) {
      kernel[i] * score[((seq_len(n) - 1 + i - (length(kernel) + 1) / 2) %% n) + 1, , drop = FALSE]
    }))
  }
  psi <- matrix(0, max(kept$draw), length(at))
  for (idx in split(seq_along(kept$opt), ceiling(seq_along(kept$opt) / 1e5))) {
    contrib <- stats::pnorm(setup$kappa * outer(kept$hi[idx], d_values[at], "-")) -
      stats::pnorm(setup$kappa * outer(kept$lo[idx], d_values[at], "-"))
    sums <- rowsum(-score[kept$opt[idx], , drop = FALSE] * contrib, kept$draw[idx])
    rows <- as.integer(rownames(sums))
    psi[rows, ] <- psi[rows, ] + sums
  }
  rho <- matrixStats::colSds(psi) / sqrt(.tcc_corr_rotations * envelopes) / sqrt(info_per_trial)
  # with every group kept, colSums(psi) / envelopes is the table's first-order
  # shift against its own scores, -sum_k mult_k s_k P-hat_k: about 0, and a
  # check of the bookkeeping above
  list(n_pool = floor((0.1 / max(rho))^2), rho = rho, rho_d = d_values[at],
       shift = colSums(psi) / envelopes)
}

############################################################################# !
# AFTER THE FIT: CERTIFICATION, SUMMARY AND PRINT                        ####
############################################################################# !

# A Monte Carlo table is certified for n_pool trials per estimate (see
# .tcc_corr_certify()); n_pool grows in proportion to the draws, which the
# warning scales to cover the data. All trials are counted, as a
# population-level d pools every one of them. The suggested draws are rounded
# up, and the 10 % margin covers the n_pool of the refit's own seed, which
# differed by 4 % between two seeds at 1e6 draws (13,104 and 12,564): with 5 %
# rounded to the nearest value, a refit could warn again
.tcc_check_certified <- function(info, nobs, d_max) {
  needed <- info$draws * nobs / info$n_pool * 1.1
  unit <- 10^(floor(log10(needed)) - 1)
  needed <- ceiling(needed / unit) * unit
  warnif(
    nobs > info$n_pool,
    "The likelihood table of the correlated noise is a Monte Carlo estimate from \\
    {format(info$draws, big.mark = ',', scientific = FALSE)} draws, certified for \\
    {format(info$n_pool, big.mark = ',', scientific = FALSE)} trials per estimate: an estimate that pools \\
    more trials can move by more than 0.1 posterior SD with the table's seed. The \\
    data hold {format(nobs, big.mark = ',', scientific = FALSE)} trials. If an estimate pools more than \\
    {format(info$n_pool, big.mark = ',', scientific = FALSE)} of them, refit with correlation_draws = \\
    {format(needed, big.mark = ',', scientific = FALSE)} (one million draws took \\
    about 70 s to build on one core of a laptop under other load; see the \\
    section Correlated noise of ?tcc)."
  )
  warnif(
    isTRUE(info$zeros_from <= d_max),
    "No Monte Carlo draw of the correlated noise reached {info$zeros} cells of the \\
    likelihood table, the first at d = {signif(info$zeros_from, 3)}, which is below \\
    d_max = {d_max}. Their probability was set to half a draw, which makes the likelihood \\
    there a floor, not the model. This happens when the correlation is weak; check \\
    that the posterior of d stays well below {signif(info$zeros_from, 3)}, or raise \\
    correlation_draws."
  )
}

.tcc_print_noise <- function(other_vars) {
  parts <- c(
    if (!is.null(other_vars$correlation)) {
      glue("correlated by distance (likelihood table from \\
           {format(other_vars$correlation_draws, big.mark = ',', scientific = FALSE)} \\
           Monte Carlo draws)")
    },
    if (!is.null(other_vars$motor_sd)) {
      glue("motor noise with SD {signif(other_vars$motor_sd, 3)} \\
           {if (identical(other_vars$units, 'bins')) 'option steps' else 'radians'}")
    }
  )
  if (length(parts)) {
    cat(glue("Noise:      {paste(parts, collapse = '; ')}."), "\n")
  }
}

# One line for each correlated fit: the table is a Monte Carlo estimate, and
# its seed and certification make the fit reproducible and its limits visible
#' @export
summary_notes.tcc <- function(model, x) {
  info <- model$other_vars$table_info
  if (is.null(info)) {
    return(NULL)
  }
  paste0(
    "Note: the likelihood of the correlated noise reads a Monte Carlo table from ",
    format(info$draws, big.mark = ",", scientific = FALSE), " draws (seed ", info$seed,
    "),\n      certified for ", format(info$n_pool, big.mark = ",", scientific = FALSE), " trials per estimate."
  )
}
