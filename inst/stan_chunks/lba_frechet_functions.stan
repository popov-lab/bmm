real lba_frechet_log_F(real x, real v, real s) {
  return -pow(x / s, -v);
}

real lba_frechet_log_M(real t, real v, real b, real A, real s) {
  array[16] real nodes = {
    -0.9894009349916499, -0.9445750230732326, -0.8656312023878318,
    -0.7554044083550030, -0.6178762444026438, -0.4580167776572274,
    -0.2816035507792589, -0.0950125098376374,  0.0950125098376374,
     0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
     0.7554044083550030,  0.8656312023878318,  0.9445750230732326,
     0.9894009349916499
  };
  array[16] real weights = {
    0.0271524594117541, 0.0622535239386479, 0.0951585116824928,
    0.1246289712555339, 0.1495959888165767, 0.1691565193950025,
    0.1826034150449236, 0.1894506104550685, 0.1894506104550685,
    0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
    0.1246289712555339, 0.0951585116824928, 0.0622535239386479,
    0.0271524594117541
  };
  real lo = (b - A) / t;
  real hi = b / t;
  real mid = 0.5 * (hi + lo);
  real half_range = 0.5 * (hi - lo);
  array[16] real log_terms;
  for (j in 1:16) {
    real u = mid + half_range * nodes[j];
    real log_z = log(u / s);
    real log_integrand = log(v) - (v * log_z) - exp(-v * log_z);
    log_terms[j] = log(weights[j]) + log_integrand;
  }
  return log(half_range) + log_sum_exp(log_terms);
}

// Single-accumulator log-PDF for LBA with Frechet drift
real lba_frechet_single_lpdf(real t, real v, real b, real A, real s) {
  return lba_frechet_log_M(t, v, b, A, s) - log(A);
}

// Single-accumulator log-survival for LBA with Frechet drift
real lba_frechet_single_lccdf(real t, real v, real b, real A, real s) {
  real log_M = lba_frechet_log_M(t, v, b, A, s);
  real hi = b / t;
  real lo = (b - A) / t;
  real log_u = log_diff_exp(
    log(b) + lba_frechet_log_F(hi, v, s),
    log(b - A) + lba_frechet_log_F(lo, v, s)
  );
  // survival numerator (u - t*M) in log space: log_diff_exp is exact where
  // log_u > log(t) + log_M; the cancellation tail (survival numerically 0) floors
  // to the same penalty as lba_log_clip (log(1e-300)).
  real log_tM = log(t) + log_M;
  return (log_u > log_tM ? log_diff_exp(log_u, log_tM) : log(1e-300)) - log(A);
}
