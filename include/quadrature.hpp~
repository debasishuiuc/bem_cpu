// quadrature.hpp
#pragma once

#include <vector>
#include <array>

struct QuadratureData {
  std::vector<double> zq, wq;      // Gauss-Legendre
  std::vector<double> xiq, etq;    // Triangle barycentric coords
  std::vector<double> wwq;         // Triangle weights
};

QuadratureData gauss_legendre(int N);
QuadratureData gauss_triangle(int N);
