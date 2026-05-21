// Single-accumulator log-PDF for LBA with normal drift
// Exact formulas follow the R reference implementation in distributions.R
real lba_normal_single_lpdf(real t, real v, real b, real A, real s) {
  real hi = b / t;
  real lo = (b - A) / t;
  real z_hi = (hi - v) / s;
  real z_lo = (lo - v) / s;
  real log_F_hi = std_normal_lcdf(z_hi | );
  real log_F_lo = std_normal_lcdf(z_lo | );
  real Phi_hi = exp(log_F_hi);
  real Phi_lo = exp(log_F_lo);
  real phi_hi = exp(std_normal_lpdf(z_hi | ));
  real phi_lo = exp(std_normal_lpdf(z_lo | ));
  real M = v * (Phi_hi - Phi_lo) + s * (phi_lo - phi_hi);
  real log_denom = std_normal_lcdf(v / s | );
  return lba_log_positive(M) - log(A) - log_denom;
}

// Single-accumulator log-survival for LBA with normal drift
// Implements posdrift (truncated normal): S_corrected = 1 - F_standard/Phi(v/s)
//   = (surv_num/A - Phi(-v/s)) / Phi(v/s)
real lba_normal_single_lccdf(real t, real v, real b, real A, real s) {
  real hi = b / t;
  real lo = (b - A) / t;
  real z_hi = (hi - v) / s;
  real z_lo = (lo - v) / s;
  real log_F_hi = std_normal_lcdf(z_hi | );
  real log_F_lo = std_normal_lcdf(z_lo | );
  real Phi_hi = exp(log_F_hi);
  real Phi_lo = exp(log_F_lo);
  real phi_hi = exp(std_normal_lpdf(z_hi | ));
  real phi_lo = exp(std_normal_lpdf(z_lo | ));
  real M = v * (Phi_hi - Phi_lo) + s * (phi_lo - phi_hi);
  real surv_num = (b * Phi_hi) - ((b - A) * Phi_lo) - (t * M);
  real log_denom = std_normal_lcdf(v / s | );
  real corrected_surv_num = surv_num - A * exp(std_normal_lcdf(-v / s | ));
  return lba_log_positive(corrected_surv_num) - log(A) - log_denom;
}
