// geometry.cpp

#include "geometry.hpp"
#include "geometry_utility.hpp"
#include "quadrature.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>


// void check_element_orientations(const Mesh& mesh) {
//   std::cout << "\n=== Element Orientation Check ===\n";

//   for (size_t k = 0; k < mesh.n.size(); ++k) {
//     const auto& n0 = mesh.p[mesh.n[k][0]]; // A
//     const auto& n1 = mesh.p[mesh.n[k][1]]; // B
//     const auto& n2 = mesh.p[mesh.n[k][2]]; // C

//     // Edge vectors
//     double v1[3] = { n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2] };
//     double v2[3] = { n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2] };

//     // Cross product v1 x v2
//     double nx = v1[1]*v2[2] - v1[2]*v2[1];
//     double ny = v1[2]*v2[0] - v1[0]*v2[2];
//     double nz = v1[0]*v2[1] - v1[1]*v2[0];

//     // Centroid
//     double cx = (n0[0] + n1[0] + n2[0]) / 3.0;
//     double cy = (n0[1] + n1[1] + n2[1]) / 3.0;
//     double cz = (n0[2] + n1[2] + n2[2]) / 3.0;

//     // Dot product between normal and centroid
//     double dot = nx * cx + ny * cy + nz * cz;

//     int orientation = (dot > 0.0) ? 1 : 0; // 1 = CCW, 0 = CW

//     std::cout << "Element " << k << ": Orientation = " << orientation << "\n";
//   }
// }

void compute_geometry(const Mesh& mesh, int quad_order) {
  const int num_elements = mesh.n.size();

  double area = 0.0, vlm = 0.0;
  Vec3 centroid{0.0, 0.0, 0.0};
  
  QuadratureData quad = gauss_triangle(quad_order);
  int mint = quad.xiq.size();

  std::vector<Vec3> vna(mesh.p.size(), Vec3{0.0, 0.0, 0.0});
  std::vector<int> itally(mesh.p.size(), 0);
  std::vector<double> crvmel(num_elements, 0.0);
  std::vector<Vec3> vnc(num_elements, Vec3{0.0, 0.0, 0.0});
  std::vector<Vec3> xcoll(num_elements, Vec3{0.0, 0.0, 0.0});

  for (int k = 0; k < num_elements; ++k) {
    std::vector<Vec3> nodes(6);
    for (int i = 0; i < 6; ++i) {
      int nid = mesh.n[k][i];
      nodes[i] = mesh.p[nid];
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
    Vec3 position_moment{0.0, 0.0, 0.0};
    
    for (int i = 0; i < mint; ++i) {
      double xi = quad.xiq[i];
      double eta = quad.etq[i];
      double w = quad.wwq[i];

      Vec3 x, dx_dxi, dx_deta, n;
      double hxi, het, hs;

      interp_p(nodes[0], nodes[1], nodes[2],
               nodes[3], nodes[4], nodes[5],
               al, be, ga, xi, eta,
               x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

      double cf = 0.5 * hs * w;

      arel += cf;
      position_moment += cf * x;
      vlm += x.dot(n) * cf;
    }

    area += arel;
    centroid += position_moment; 

    std::array<double, 6> xxi = {0.0, 1.0, 0.0, al, ga, 0.0};
    std::array<double, 6> eet = {0.0, 0.0, 1.0, 0.0, 1.0 - ga, be};

    std::array<Vec3, 6> dxdxi, dxdeta, vnorm;

    for (int i = 0; i < 6; ++i) {
      double xi = xxi[i];
      double eta = eet[i];

      Vec3 x, dx_dxi, dx_deta, n;
      double hxi, het, hs;

      interp_p(nodes[0], nodes[1], nodes[2],
               nodes[3], nodes[4], nodes[5],
               al, be, ga, xi, eta,
               x, dx_dxi, dx_deta, n,
               hxi, het, hs, 2);

      dxdxi[i] = dx_dxi;
      dxdeta[i] = dx_deta;
      vnorm[i] = n;

      int m = mesh.n[k][i];
      vna[m] += n;
      itally[m]++;
    }

    Vec3 crv(0.0, 0.0, 0.0);

    auto cross_dx = [&](int i) -> Vec3 {
      return vnorm[i].cross(dxdxi[i]);
    };

    Vec3 b1 = cross_dx(0), b2 = cross_dx(3), b3 = cross_dx(1);
    crv += al * b1 + b2 + (1.0 - al) * b3;

    b1 = cross_dx(1); b2 = cross_dx(4); b3 = cross_dx(2);
    crv -= (1.0 - ga) * b1 + b2 + ga * b3;

    auto cross_de = [&](int i) -> Vec3 {
      return vnorm[i].cross(dxdeta[i]);
    };

    b1 = cross_de(1); b2 = cross_de(4); b3 = cross_de(2);
    crv += (1.0 - ga) * b1 + b2 + ga * b3;

    b1 = cross_de(0); b2 = cross_de(5); b3 = cross_de(2);
    crv -= be * b1 + b2 + (1.0 - be) * b3;

    Vec3 x_c, dx_dxi, dx_deta, n_centroid;
    double hxi, het, hs;

    interp_p(nodes[0], nodes[1], nodes[2],
             nodes[3], nodes[4], nodes[5],
             al, be, ga, 1.0 / 3.0, 1.0 / 3.0,
             x_c, dx_dxi, dx_deta, n_centroid,
             hxi, het, hs, 2);

    crvmel[k] = 0.25 * crv.dot(n_centroid) / arel;
    vnc[k] = n_centroid;
    xcoll[k] = x_c;
  }

  for (int i = 0; i < mesh.p.size(); ++i) {
    if (itally[i] > 0) {
      vna[i] /= static_cast<double>(itally[i]);
      if (vna[i].norm() > 1e-12) {
        vna[i] = vna[i].normalized();
      }
    }
  }

  vlm /= 3.0;
  centroid /= area;
  
  std::cout << "Total Surface Area: " << area << "\n";
  std::cout << "Total Volume: " << vlm << "\n";
  std::cout << "Centroid: (" << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << ")\n";
  
  Eigen::Matrix3d Mmat = Eigen::Matrix3d::Zero();

  for (int k = 0; k < num_elements; ++k) {
    std::vector<Vec3> nodes(6);
    for (int i = 0; i < 6; ++i) {
      int nid = mesh.n[k][i];
      nodes[i] = mesh.p[nid];
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

    for (int i = 0; i < quad.xiq.size(); ++i) {
      double xi = quad.xiq[i];
      double eta = quad.etq[i];
      double w = quad.wwq[i];

      Vec3 x_vec, dx_dxi, dx_deta, n_vec;
      double hxi, het, hs;

      interp_p(nodes[0], nodes[1], nodes[2],
               nodes[3], nodes[4], nodes[5],
               al, be, ga, xi, eta,
               x_vec, dx_dxi, dx_deta, n_vec,
               hxi, het, hs, 2);

      double cf = 0.5 * hs * w;
      Vec3 xhat = x_vec - centroid;
      double r2 = xhat.dot(xhat);

      Mmat(0, 0) += (r2 - xhat[0] * xhat[0]) * cf;
      Mmat(0, 1) += (-xhat[0] * xhat[1]) * cf;
      Mmat(0, 2) += (-xhat[0] * xhat[2]) * cf;

      Mmat(1, 0) += (-xhat[1] * xhat[0]) * cf;
      Mmat(1, 1) += (r2 - xhat[1] * xhat[1]) * cf;
      Mmat(1, 2) += (-xhat[1] * xhat[2]) * cf;

      Mmat(2, 0) += (-xhat[2] * xhat[0]) * cf;
      Mmat(2, 1) += (-xhat[2] * xhat[1]) * cf;
      Mmat(2, 2) += (r2 - xhat[2] * xhat[2]) * cf;
    }
  }

  Eigen::Matrix3d Mmat_inv = Mmat.lu().inverse();
  Mmat = Mmat_inv;

  std::cout << "Inverse of Mmat:\n" << Mmat << "\n";

  std::cout << "\nFirst 5 Node Normals (averaged):\n";
  for (int i = 0; i < std::min(5, static_cast<int>(mesh.p.size())); ++i) {
    if (itally[i] > 0) {
      std::cout << "Node " << i << ": Pos (" << mesh.p[i][0] << ", "
                << mesh.p[i][1] << ", " << mesh.p[i][2] << "), Normal ("
                << vna[i][0] << ", " << vna[i][1] << ", " << vna[i][2] << ")\n";
    } else {
      std::cout << "Node " << i << ": No contributing elements.\n";
    }
  }

  std::cout << "\nFirst 5 Element Normals and Centroids:\n";
  for (int i = 0; i < std::min(5, static_cast<int>(vnc.size())); ++i) {
    std::cout << "Element " << i << ": Centroid (" << xcoll[i][0] << ", "
              << xcoll[i][1] << ", " << xcoll[i][2] << "), Normal ("
              << vnc[i][0] << ", " << vnc[i][1] << ", " << vnc[i][2] << ")\n";
  }

  std::cout << "\nFirst 5 Element Curvatures:\n";
  for (int i = 0; i < std::min(5, static_cast<int>(crvmel.size())); ++i) {
    std::cout << "Element " << i << ": Curvature = " << crvmel[i] << "\n";
  }

  // check_element_orientations(mesh);
}
