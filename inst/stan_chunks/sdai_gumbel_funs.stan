  // Gumbel-min CDF: P(X <= x) = 1 - exp(-exp(x - mu))
  real gumbelmin_cdf_sdai(real x, real mu) {
    return 1 - exp(-exp(x - mu));
  }

  // Internal helper for source confusion probability
  real sdai_pinternal1(real a, real b, real g1) {
    return exp(-exp(-g1 + a) + b);
  }

  // Internal helper for source confusion probability (two-source)
  real sdai_pinternal2(real a, real b, real g1, real g2) {
    return exp(-exp(-g1 + a) - exp(-g2 + a) + b);
  }

  // Source confusion probability in interval [l, u] for sources at g1, g2
  real sdai_pfun(real l, real u, real g1, real g2) {
    real denominator = exp(g1) + exp(g2);
    return (sdai_pinternal1(l, g1, g1) - sdai_pinternal1(u, g1, g1) +
            sdai_pinternal1(l, g2, g1) - sdai_pinternal1(u, g2, g1) +
            sdai_pinternal2(u, g2, g1, g2) - sdai_pinternal2(l, g2, g1, g2)) /
           denominator;
  }

  // Discrimination probability (both sources at location 0)
  real sdai_p_hit_gumbel(real l, real u) {
    return gumbelmin_cdf_sdai(u, 0) * gumbelmin_cdf_sdai(u, 0) -
           gumbelmin_cdf_sdai(l, 0) * gumbelmin_cdf_sdai(l, 0);
  }
