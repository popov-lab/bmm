// M(t) (the truncated first moment) and the survival numerator are analytically
// non-negative; they can dip below 0 only through floating-point cancellation at
// extreme parameter values. Clip at a small positive floor before taking the log:
// legitimate positive values pass through exactly (a softplus would distort small
// positives), while cancellation noise is floored to a large finite penalty rather
// than -inf, keeping HMC gradients well-behaved during warmup.
real lba_log_clip(real x) {
  return log(fmax(x, 1e-300));
}
