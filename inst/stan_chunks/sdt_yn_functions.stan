// yes/no SDT likelihood over all rows at once (loop = FALSE), for trial rows and
// aggregated counts alike; needs sdt_dist_functions.stan
// d is d_a, sdratio the SD ratio on the natural scale, stimulus 0 = new and 1 = old
// criterion is measured from the noise mean, in noise-SD units
// dist_type is one code for all rows: check_data() writes the same value in each

int sdt_n_between(array[] int k, int lo, int hi) {
  int out = 0;
  for (n in 1:size(k)) out += k[n] >= lo && k[n] <= hi;
  return out;
}

array[] int sdt_which_between(array[] int k, int lo, int hi) {
  array[sdt_n_between(k, lo, hi)] int out;
  int j = 1;
  for (n in 1:size(k)) {
    if (k[n] >= lo && k[n] <= hi) {
      out[j] = n;
      j += 1;
    }
  }
  return out;
}

// sum over rows of k log P(old) (old = 1) or k log P(new) (old = 0), where
// k counts that row's responses on the side; rows without any are skipped, so
// a probability that underflows to 0 where its count is 0 adds nothing. Rows
// with one response, all of them on trial rows, are summed in one call
real sdt_yn_side(array[] int k, vector x, int dist_type, int old) {
  array[sdt_n_between(k, 1, 1)] int single = sdt_which_between(k, 1, 1);
  array[sdt_n_between(k, 2, 2147483647)] int many = sdt_which_between(k, 2, 2147483647);
  if (old == 1) {
    return sdt_log_one_minus_cumprob_sum(x[single], dist_type) +
      dot_product(to_vector(k[many]), sdt_log_one_minus_cumprob(x[many], dist_type));
  }
  return sdt_log_cumprob_sum(x[single], dist_type) +
    dot_product(to_vector(k[many]), sdt_log_cumprob(x[many], dist_type));
}

// x is -eta, the criterion in units of each row's evidence distribution.
// P(old) is the survivor at x, which differs from F(-x) for the asymmetric
// Gumbels. lchoose(n, 0) and lchoose(n, n) are exactly 0, so only rows with
// responses on both sides add it
real sdt_yn_binomial(array[] int y, vector x, int dist_type, array[] int trials) {
  int N = size(y);
  array[N] int n_new;
  real out = 0;
  for (n in 1:N) {
    n_new[n] = trials[n] - y[n];
    if (y[n] > 0 && n_new[n] > 0) out += lchoose(trials[n], y[n]);
  }
  return out + sdt_yn_side(y, x, dist_type, 1) + sdt_yn_side(n_new, x, dist_type, 0);
}

// equal_variance = 1 says sdratio is the constant 1 on every row (it is fixed at
// its default); sdt_rms_scale(1) and the division by 1 are exact, so skipping
// them gives the same bits for less work
real sdt_yn_lpmf(array[] int y, vector mu, vector d, vector criterion,
                 vector sdratio, array[] int stimulus, int dist_type,
                 array[] int trials, int equal_variance) {
  array[sdt_n_between(stimulus, 1, 1)] int signal = sdt_which_between(stimulus, 1, 1);
  vector[size(y)] x = criterion;
  if (equal_variance == 1) {
    x[signal] = criterion[signal] - d[signal];
  } else {
    x[signal] = (criterion[signal] - d[signal] .* sdt_rms_scale(sdratio[signal])) ./
                sdratio[signal];
  }
  return sdt_yn_binomial(y, x, dist_type, trials);
}

// the same likelihood with the criterion measured from the midpoint between the means
real sdt_yn_bias_lpmf(array[] int y, vector mu, vector d, vector bias,
                      vector sdratio, array[] int stimulus, int dist_type,
                      array[] int trials, int equal_variance) {
  array[sdt_n_between(stimulus, 1, 1)] int signal = sdt_which_between(stimulus, 1, 1);
  if (equal_variance == 1) {
    vector[size(y)] x = bias + d / 2.0;
    x[signal] = x[signal] - d[signal];
    return sdt_yn_binomial(y, x, dist_type, trials);
  }
  vector[size(y)] separation = d .* sdt_rms_scale(sdratio);
  vector[size(y)] x = bias + separation / 2.0;
  x[signal] = (x[signal] - separation[signal]) ./ sdratio[signal];
  return sdt_yn_binomial(y, x, dist_type, trials);
}
