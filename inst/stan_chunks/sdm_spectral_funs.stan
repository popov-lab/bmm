  real sdm_spectral_kernel_from_cos(real cos_diff, real kappa) {
    return exp(kappa * (cos_diff - 1)) * sqrt(kappa) * inv_sqrt(2 * pi());
  }

  real sdm_spectral_activation_obs(real y, real mu, real kappa,
                                   array[] real item_angles,
                                   array[] real item_weights,
                                   int J) {
    real out = 0;

    for (j in 1:J) {
      out += item_weights[j] * sdm_spectral_kernel_from_cos(
        cos(y - (mu + item_angles[j])),
        kappa
      );
    }

    return out;
  }

  real sdm_spectral_weight_sum(array[] real item_weights, int J) {
    real out = 0;

    for (j in 1:J) {
      out += item_weights[j];
    }

    return out;
  }

  real sdm_spectral_weight_max(array[] real item_weights, int J) {
    real out = item_weights[1];

    for (j in 2:J) {
      out = fmax(out, item_weights[j]);
    }

    return out;
  }

  real sdm_spectral_circular_distance(real x) {
    real y = abs(x);
    return fmin(y, 2 * pi() - y);
  }

  real sdm_spectral_resultant_length(array[] real item_angles,
                                     array[] real item_weights,
                                     int J) {
    real cx = 0;
    real sx = 0;

    for (j in 1:J) {
      cx += item_weights[j] * cos(item_angles[j]);
      sx += item_weights[j] * sin(item_angles[j]);
    }

    return sqrt(square(cx) + square(sx));
  }

  real sdm_spectral_min_separation(array[] real item_angles, int J) {
    real out = pi();

    if (J < 2) {
      return out;
    }

    for (j in 1:(J - 1)) {
      for (k in (j + 1):J) {
        out = fmin(out, sdm_spectral_circular_distance(item_angles[j] - item_angles[k]));
      }
    }

    return out;
  }

  int sdm_spectral_bump_L(int L) {
    if (L <= 32) {
      return 64;
    } else if (L <= 64) {
      return 128;
    } else if (L <= 128) {
      return 256;
    } else {
      return 512;
    }
  }

  int sdm_get_L_general(real kappa, int J, real w_sum, real w_max,
                        real R, real delta_min) {
    real rho = R / fmax(w_sum, 1e-12);
    real logA = log1p(sqrt(kappa / (2 * pi())) * w_sum);
    real closeness = -log(fmax(delta_min, 1e-3));
    real score = 2.6520 +
      0.9349 * log(kappa) +
      0.1815 * logA +
      0.4803 * log1p(w_max) +
      0.2386 * rho +
      0.02224 * closeness -
      0.06635 * log(J);
    int L;

    if (score < 5.5) {
      L = 32;
    } else if (score < 6.5) {
      L = 64;
    } else if (score < 7.5) {
      L = 128;
    } else if (score < 8.5) {
      L = 256;
    } else {
      L = 512;
    }

    if (
      kappa > 10 ||
      logA > 5 ||
      (delta_min < 0.25 && rho > 0.8)
    ) {
      L = sdm_spectral_bump_L(L);
    }

    if (
      kappa > 13 ||
      (delta_min < 0.15 && logA > 4.2) ||
      (w_max > 150 && rho > 0.75)
    ) {
      L = 512;
    }

    return L;
  }

  real sdm_spectral_logz(real mu, real kappa,
                         array[] real item_angles,
                         array[] real item_weights,
                         int J,
                         int L,
                         vector grid_cos,
                         vector grid_sin) {
    int stride = num_elements(grid_cos) %/% L;
    vector[L] act;

    for (m in 1:L) {
      int idx = 1 + (m - 1) * stride;
      real out = 0;

      for (j in 1:J) {
        real shift = mu + item_angles[j];
        real cos_shift = cos(shift);
        real sin_shift = sin(shift);
        real cos_diff = grid_cos[idx] * cos_shift + grid_sin[idx] * sin_shift;

        out += item_weights[j] * sdm_spectral_kernel_from_cos(cos_diff, kappa);
      }

      act[m] = out;
    }

    return log_sum_exp(act) + log(2 * pi() / L);
  }
