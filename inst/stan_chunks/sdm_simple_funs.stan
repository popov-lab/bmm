  #include 'fun_tan_half.stan'

  real sdm_simple_lpdf(vector y, vector mu, vector c, vector kappa) {
    int N = num_elements(y);
    real out = 0;

    for (n in 1:N) {
      array[1] real item_angles;
      array[1] real item_weights;

      item_angles[1] = 0;
      item_weights[1] = exp(c[n]);
      out += sdm_spectral_activation_obs(
        y[n],
        mu[n],
        kappa[n],
        item_angles,
        item_weights,
        1
      );
    }

    return out;
  }
