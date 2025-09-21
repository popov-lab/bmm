// Specify additional hyperbolic functions
real tanh_over_x(real x) {
  if (abs(x) < 1e-6) {
    real x2 = x * x;
    return 1 - x2/3 + 2*x2*x2/15 - 17*x2*x2*x2/315;
  }
  return tanh(x) / x;
}

real phi(real x) {
  if (abs(x) < 1e-6) {
    real x2 = x * x;
    return 1 + x2/3 - x2*x2/45 + 2*x2*x2*x2/945;
  }
  return x * cosh(x) / fmax(1e-12, sinh(x));
}

real psi(real x) {
  if (abs(x) < 1e-6) {
    real x2 = x * x;
    return 1 - x2/3 + x2*x2/15 - 2*x2*x2*x2/189;
  }
  real sh = sinh(x);
  return (x * x) / fmax(1e-24, sh * sh);
}

// ---------- kappas in 0–a convention ----------
real kz_0a(real v, real a, real s) { return (v * a) / square(s); }
real kx_0a(real v, real z, real s) { return (v * z) / square(s); }

// ---------- choice probability for upper (stable) ----------
real ddm_p_upper(real v, real a, real z, real s) {
  real kz = kz_0a(v, a, s);
  if (abs(kz) < 1e-8) {
    real zr = (z / a);
    return zr + kz * zr * (1 - zr);
  } else {
    real kx = kx_0a(v, z, s);
    real num = -expm1(-2 * kx);  // 1 - exp(-2 kx)
    real den = -expm1(-2 * kz);  // 1 - exp(-2 kz)
    return num / den;
  }
}

real ddm_mdt_upper(real v, real a, real z, real s) {
  real kz = kz_0a(v, a, s);
  real kx = kx_0a(v, z, s);
  return square(s) / square(v) * ( phi(kz) - phi(kx) );
}
real ddm_mdt_lower(real v, real a, real z, real s) {
  real kz = kz_0a(v, a, s);
  real kx = kx_0a(v, z, s);
  return square(s) / square(v) * ( phi(kz) - phi(kz - kx) );
}
real ddm_vrt_upper(real v, real a, real z, real s) {
  real kz = kz_0a(v, a, s);
  real kx = kx_0a(v, z, s);
  real factor = square(square(s) / square(v));
  return factor * ( psi(kz) + phi(kz) - psi(kx) - phi(kx) );
}
real ddm_vrt_lower(real v, real a, real z, real s) {
  real kz = kz_0a(v, a, s);
  real kx = kx_0a(v, z, s);
  real factor = square(square(s) / square(v));
  return factor * ( psi(kz) + phi(kz) - psi(kz - kx) - phi(kz - kx) );
}

real ezDMupper_lpdf (real mrt, real drift, real bound, real ndt, real zr, real s, real vrt, int n_upper, int trials) {
  // compute additional parameters
  real z  = zr * bound;
  real k_z = (drift * bound) / square(s);   // = v a / s^2
  real k_x = (drift * z) / square(s);   // = v z / s^2

  // Calculate moments based on DM parameters
  real MDT_upper = (square(s) / square(drift)) * ( phi(k_z) - phi(k_x) );
  real VRT_upper = square(square(s) / square(drift)) * ( psi(k_z) + phi(k_z) - psi(k_x) - phi(k_x) );

  if (abs(k_z) < 1e-8) {
    real p_upper = zr + k_z * zr * (1 - zr);
    return binomial_lpmf(n_upper | trials, p_upper) + normal_lpdf(mrt | ndt + MDT_upper, sqrt(VRT_upper/n_upper)) + gamma_lpdf(vrt | ((n_upper-1)/2.0), ((n_upper-1)/(2*VRT_upper)));
  } else {
    real num = -expm1(-2 * k_x);   // 1 - exp(-2 k_x)
    real den = -expm1(-2 * k_z);   // 1 - exp(-2 k_z)
    real p_upper  = num / den;
    return binomial_lpmf(n_upper | trials, p_upper) + normal_lpdf(mrt | ndt + MDT_upper, sqrt(VRT_upper/n_upper)) + gamma_lpdf(vrt | ((n_upper-1)/2.0), ((n_upper-1)/(2*VRT_upper)));
  }
}

real ezDMlower_lpdf (real mrt, real drift, real bound, real ndt, real zr, real s, real vrt, int n_upper, int trials) {
  // compute additional parameters
  real z  = zr * bound;
  real k_z = (drift * bound) / square(s);   // = v a / s^2
  real k_x = (drift * z) / square(s);   // = v z / s^2

  // Calculate moments based on DM parameters
  real MDT_lower = (square(s) / square(drift)) * ( phi(k_z) - phi(k_z - k_x) );
  real VRT_lower = square(square(s) / square(drift)) * ( psi(k_z) + phi(k_z) - psi(k_z - k_x) - phi(k_z - k_x) );

  if (abs(k_z) < 1e-8) {
    real p_upper = zr + k_z * zr * (1 - zr);
    return binomial_lpmf(n_upper | trials, 1 - p_upper) + normal_lpdf(mrt | ndt + MDT_lower, sqrt(VRT_lower/(trials - n_upper))) + gamma_lpdf(vrt | (((trials - n_upper)-1)/2.0), (((trials - n_upper)-1)/(2*VRT_lower)));
  } else {
    real p_upper  = -expm1(-2 * k_x) / -expm1(-2 * k_z);
    return binomial_lpmf(n_upper | trials, 1 - p_upper) + normal_lpdf(mrt | ndt + MDT_lower, sqrt(VRT_lower/(trials - n_upper))) + gamma_lpdf(vrt | (((trials - n_upper)-1)/2.0), (((trials - n_upper)-1)/(2*VRT_lower)));
  }
}

// Specify likelihood for ezDM
real ezdm_4par_lpdf (real mrt, real mu, real drift, real bound, real ndt, real zr, real s, real vrt, int hits, int trials, int response) {
  if (response == 1) {
    return ezDMupper_lpdf(mrt | drift, bound, ndt, zr, s, vrt, hits, trials);
  } else {
    return ezDMlower_lpdf(mrt | drift, bound, ndt, zr, s, vrt, hits, trials);
  }
}
