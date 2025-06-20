// geometry_utility.cpp

#include "geometry_utility.hpp"
#include <cmath>

// Interpolation function matching Fortran's interp_p

// Modernized interpolation using Vec3 I/O and internal logic
void interp_p(const Vec3& p1, const Vec3& p2, const Vec3& p3,
              const Vec3& p4, const Vec3& p5, const Vec3& p6,
              double al, double be, double ga,
              double xi, double eta,
              Vec3& x,
              Vec3& dx_dxi, Vec3& dx_deta,
              Vec3& n,
              double& hxi, double& het, double& hs,
              int Ichoose) {
  
  // Intermediate shape parameters
  double alc = 1.0 - al;
  double bec = 1.0 - be;
  double gac = 1.0 - ga;

  double alalc = al * alc;
  double bebec = be * bec;
  double gagac = ga * gac;

  // Shape functions
  double ph2 = xi * (xi - al + eta * (al - ga) / gac) / alc;
  double ph3 = eta * (eta - be + xi * (be + ga - 1.0) / ga) / bec;
  double ph4 = xi * (1.0 - xi - eta) / alalc;
  double ph5 = xi * eta / gagac;
  double ph6 = eta * (1.0 - xi - eta) / bebec;
  double ph1 = 1.0 - ph2 - ph3 - ph4 - ph5 - ph6;

  // Position
  x = ph1 * p1 + ph2 * p2 + ph3 * p3 + ph4 * p4 + ph5 * p5 + ph6 * p6;

  if (Ichoose <= 1) return;

  // Derivatives w.r.t xi
  double dph2 = (2.0 * xi - al + eta * (al - ga) / gac) / alc;
  double dph3 = eta * (be + ga - 1.0) / (ga * bec);
  double dph4 = (1.0 - 2.0 * xi - eta) / alalc;
  double dph5 = eta / gagac;
  double dph6 = -eta / bebec;
  double dph1 = -dph2 - dph3 - dph4 - dph5 - dph6;

  dx_dxi = dph1 * p1 + dph2 * p2 + dph3 * p3 + dph4 * p4 + dph5 * p5 + dph6 * p6;

  // Derivatives w.r.t eta
  double pph2 = xi * (al - ga) / (alc * gac);
  double pph3 = (2.0 * eta - be + xi * (be + ga - 1.0) / ga) / bec;
  double pph4 = -xi / alalc;
  double pph5 = xi / gagac;
  double pph6 = (1.0 - xi - 2.0 * eta) / bebec;
  double pph1 = -pph2 - pph3 - pph4 - pph5 - pph6;

  dx_deta = pph1 * p1 + pph2 * p2 + pph3 * p3 + pph4 * p4 + pph5 * p5 + pph6 * p6;

  // Norms of tangent vectors
  hxi = dx_dxi.norm();
  het = dx_deta.norm();

  // Cross product (area-weighted normal)
  Vec3 raw_n = dx_dxi.cross(dx_deta);
  hs = raw_n.norm();

  n = (hs > 1e-14) ? raw_n / hs : Vec3(0.0, 0.0, 0.0);
}
