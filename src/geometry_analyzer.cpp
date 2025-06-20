
// geometry_analyzer.cpp

#include "geometry_analyzer.hpp"
#include "geometry_utility.hpp"
#include "quadrature.hpp"
#include <iostream>
#include <cstdlib>
#include <algorithm>

void GeometryAnalyzer::compute(const Mesh& mesh, int quad_order) {
  const int num_nodes = mesh.p.size();
  const int num_elements = mesh.n.size();

  totalArea = 0.0;
  totalVolume = 0.0;
  surfaceCentroid = Vec3(0.0, 0.0, 0.0);
  Mmat.setZero();

  QuadratureData quad = gauss_triangle(quad_order);
  int mint = quad.xiq.size();

  // Resize output containers
  nodeNormals.assign(num_nodes, Vec3(0.0, 0.0, 0.0));
  std::vector<int> nodeContribCount(num_nodes, 0);
  elementCurvature.assign(num_elements, 0.0);
  elementNormals.assign(num_elements, Vec3(0.0, 0.0, 0.0));
  elementCentroids.assign(num_elements, Vec3(0.0, 0.0, 0.0));
  elementArea.assign(num_elements, 0.0);

  for (int k = 0; k < num_elements; ++k) {
    std::vector<Vec3> nodes(6);
    for (int i = 0; i < 6; ++i) {
      nodes[i] = mesh.p[mesh.n[k][i]];
    }

    double d42 = (nodes[3] - nodes[1]).norm();
    double d41 = (nodes[3] - nodes[0]).norm();
    double d63 = (nodes[5] - nodes[2]).norm();
    double d61 = (nodes[5] - nodes[0]).norm();
    double d52 = (nodes[4] - nodes[1]).norm();
    double d53 = (nodes[4] - nodes[2]).norm();

    double al = 1.0 / (1.0 + d42 / d41);
    double be = 1.0 / (1.0 + d63 / d61);
    double ga = 1.0 / (1.0 + d52 / d53);

    double arel = 0.0;
    Vec3 position_moment(0.0, 0.0, 0.0);

    for (int i = 0; i < mint; ++i) {
      Vec3 x, dx_dxi, dx_deta, n;
      double hxi, het, hs;

      interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
               al, be, ga, quad.xiq[i], quad.etq[i],
               x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

      double cf = 0.5 * hs * quad.wwq[i];
      arel += cf;
      position_moment += cf * x;
      totalVolume += x.dot(n) * cf;
    }

    totalArea += arel;
    surfaceCentroid += position_moment;
    elementArea[k] = arel;
    elementCentroids[k] = position_moment / arel;

    std::array<double, 6> xxi = {0.0, 1.0, 0.0, al, ga, 0.0};
    std::array<double, 6> eet = {0.0, 0.0, 1.0, 0.0, 1.0 - ga, be};
    std::array<Vec3, 6> dxdxi, dxdeta, vnorm;

    for (int i = 0; i < 6; ++i) {
      Vec3 x, dx_dxi, dx_deta, n;
      double hxi, het, hs;

      interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
               al, be, ga, xxi[i], eet[i],
               x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

      dxdxi[i] = dx_dxi;
      dxdeta[i] = dx_deta;
      vnorm[i] = n;

      int m = mesh.n[k][i];
      nodeNormals[m] += n;
      nodeContribCount[m]++;
    }

    Vec3 crv(0.0, 0.0, 0.0);
    auto cross_dx = [&](int i) -> Vec3 { return vnorm[i].cross(dxdxi[i]); };
    auto cross_de = [&](int i) -> Vec3 { return vnorm[i].cross(dxdeta[i]); };

    crv += al * cross_dx(0) + cross_dx(3) + (1.0 - al) * cross_dx(1);
    crv -= (1.0 - ga) * cross_dx(1) + cross_dx(4) + ga * cross_dx(2);
    crv += (1.0 - ga) * cross_de(1) + cross_de(4) + ga * cross_de(2);
    crv -= be * cross_de(0) + cross_de(5) + (1.0 - be) * cross_de(2);

    Vec3 x_c, dx_dxi, dx_deta, n_centroid;
    double hxi, het, hs;
    interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
             al, be, ga, 1.0 / 3.0, 1.0 / 3.0,
             x_c, dx_dxi, dx_deta, n_centroid,
             hxi, het, hs, 2);

    elementNormals[k] = n_centroid;
    elementCurvature[k] = 0.25 * crv.dot(n_centroid) / arel;
  }

  // Normalize node normals
  for (int i = 0; i < num_nodes; ++i) {
    if (nodeContribCount[i] > 0) {
      nodeNormals[i] /= static_cast<double>(nodeContribCount[i]);
      if (nodeNormals[i].norm() > 1e-12) {
        nodeNormals[i] = nodeNormals[i].normalized();
      }
    }
  }

  totalVolume /= 3.0;
  surfaceCentroid /= totalArea;

  // compute node curvature from element curvature

  nodeCurvature.assign(nodeNormals.size(), 0.0);      // Initialize node curvature vector
  std::vector<double> nodeAreaSum(nodeNormals.size(), 0.0);  // To accumulate area weights per node

  for (int elem = 0; elem < mesh.Nelm; ++elem) {
    double curv = elementCurvature[elem];
    double area = elementArea[elem]; // element area computed previously

    // Loop over all nodes of this element (assuming 6-node triangle)
    for (int localNodeIdx = 0; localNodeIdx < 6; ++localNodeIdx) {
      int globalNodeIdx = mesh.n[elem][localNodeIdx];
      nodeCurvature[globalNodeIdx] += curv * area;
      nodeAreaSum[globalNodeIdx] += area;
    }
  }

  // Normalize to get weighted average
  for (int i = 0; i < nodeCurvature.size(); ++i) {
    if (nodeAreaSum[i] > 1e-14) {
      nodeCurvature[i] /= nodeAreaSum[i];
    } else {
      nodeCurvature[i] = 0.0;  // Or handle isolated nodes differently
    }
  }

  // Compute Mmat (moment of inertia)
  for (int k = 0; k < num_elements; ++k) {
    std::vector<Vec3> nodes(6);
    for (int i = 0; i < 6; ++i) {
      nodes[i] = mesh.p[mesh.n[k][i]];
    }

    double d42 = (nodes[3] - nodes[1]).norm();
    double d41 = (nodes[3] - nodes[0]).norm();
    double d63 = (nodes[5] - nodes[2]).norm();
    double d61 = (nodes[5] - nodes[0]).norm();
    double d52 = (nodes[4] - nodes[1]).norm();
    double d53 = (nodes[4] - nodes[2]).norm();

    double al = 1.0 / (1.0 + d42 / d41);
    double be = 1.0 / (1.0 + d63 / d61);
    double ga = 1.0 / (1.0 + d52 / d53);

    for (int i = 0; i < mint; ++i) {
      Vec3 x, dx_dxi, dx_deta, n;
      double hxi, het, hs;
      interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
               al, be, ga, quad.xiq[i], quad.etq[i],
               x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

      double cf = 0.5 * hs * quad.wwq[i];
      Vec3 xhat = x - surfaceCentroid;
      double r2 = xhat.dot(xhat);

      for (int m = 0; m < 3; ++m) {
        for (int n = 0; n < 3; ++n) {
          double delta = (m == n) ? 1.0 : 0.0;
          Mmat(m, n) += (r2 * delta - xhat[m] * xhat[n]) * cf;
        }
      }
    }
  }

  Mmat = Mmat.lu().inverse();
}

void GeometryAnalyzer::print_summary(const Mesh& mesh, int num_print) const {
  std::cout << "Total Surface Area: " << totalArea << "\n";
  std::cout << "Total Volume: " << totalVolume << "\n";
  std::cout << "Centroid: (" << surfaceCentroid[0] << ", "
            << surfaceCentroid[1] << ", " << surfaceCentroid[2] << ")\n";
  std::cout << "Inverse of Mmat:\n" << Mmat << "\n";

  std::cout << "\nFirst " << num_print << " Node Normals (averaged):\n";
  for (int i = 0; i < std::min(num_print, (int)mesh.p.size()); ++i) {
    const auto& n = nodeNormals[i];
    const auto& p = mesh.p[i];
    std::cout << "Node " << i << ": Pos (" << p[0] << ", " << p[1] << ", " << p[2]
              << "), Normal (" << n[0] << ", " << n[1] << ", " << n[2] << ")\n";
  }

  std::cout << "\nFirst " << num_print << " Element Normals and Centroids:\n";
  for (int i = 0; i < std::min(num_print, (int)elementNormals.size()); ++i) {
    const auto& c = elementCentroids[i];
    const auto& n = elementNormals[i];
    std::cout << "Element " << i << ": Centroid (" << c[0] << ", " << c[1] << ", " << c[2]
              << "), Normal (" << n[0] << ", " << n[1] << ", " << n[2] << ")\n";
  }

  std::cout << "\nFirst " << num_print << " Element Curvatures:\n";
  for (int i = 0; i < std::min(num_print, (int)elementCurvature.size()); ++i) {
    std::cout << "Element " << i << ": Curvature = " << elementCurvature[i] << "\n";
  }
}
