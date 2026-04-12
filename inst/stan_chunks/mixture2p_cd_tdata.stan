  vector[101] CD_GRID;
  real CD_DX = 2 * pi() / 100;
  for (i in 1:101) {
    CD_GRID[i] = -pi() + (i - 1) * CD_DX;
  }
