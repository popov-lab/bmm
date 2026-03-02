  // precompute quadrature grid for SDM change detection
  vector[51] SDM_CD_GRID;
  for (i in 1:51) {
    SDM_CD_GRID[i] = -pi() + (i - 1) * 2 * pi() / 50;
  }
