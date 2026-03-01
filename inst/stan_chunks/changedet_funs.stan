  // Shared utility functions for change detection models
  // Based on Lin & Oberauer (2022, Cognitive Psychology)

  // Trapezoidal quadrature integration helper
  // Integrates f(x) over a uniform grid with spacing dx
  real trapz(vector f, real dx) {
    int N = size(f);
    return (sum(f) - 0.5 * (f[1] + f[N])) * dx;
  }
