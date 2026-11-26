// log density of a response time and decision under the 4-parameter diffusion model
real ddm_lpdf (real rt, real mu, real drift, real bound, real ndt, real zr, int dec) {
  if (dec == 1) {
    return wiener_lpdf (rt | bound, ndt, zr, drift);
  } else {
    return wiener_lpdf (rt | bound, ndt, 1 - zr, -drift);
  }
}
