/**
# My functions */

void vof_advection2 (scalar f2, double dt2, face vector uf2, int i) {
  face vector uf_stored = uf;
  double dt_stored = dt;
  uf = uf2;
  dt = dt2;
  vof_advection ({f2}, i);
  dt = dt_stored;
  uf = uf_stored;
}
