  int SDM_GRID_SIZE = 256;
  vector[SDM_GRID_SIZE] SDM_GRID;
  real SDM_GRID_STEP = 2 * pi() / SDM_GRID_SIZE;
  for (m in 1:SDM_GRID_SIZE) {
    SDM_GRID[m] = -pi() + (m - 0.5) * SDM_GRID_STEP;
  }
