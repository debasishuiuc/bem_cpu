
// geometry_analyzer.cpp

#include "geometry_analyzer.hpp"
#include "geometry_utility.hpp"
#include "quadrature.hpp"
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <omp.h>



void GeometryAnalyzer::compute(const Mesh& mesh, int quad_order) {
  const int num_nodes = mesh.p.size();
  const int num_elements = mesh.n.size();

  totalArea = 0.0;
  totalVolume = 0.0;
  surfaceCentroid = Vec3(0.0, 0.0, 0.0);
  Mmat.setZero();

  QuadratureData quad = gauss_triangle(quad_order);
  int mint = quad.xiq.size();

  initialize_containers(mesh, mint);

  compute_element_geometry(mesh, quad, mint);

  average_node_normals();

  compute_node_curvature(mesh);

  compute_moment_matrix(mesh, quad, mint);
}






void GeometryAnalyzer::initialize_containers(const Mesh& mesh, int mint) {
  int num_nodes = mesh.p.size();
  int num_elements = mesh.n.size();

  nodeNormals.assign(num_nodes, Vec3(0.0, 0.0, 0.0));
  nodeContribCount.assign(num_nodes, 0);  // ✅ Add this line
  elementCurvature.assign(num_elements, 0.0);
  elementNormals.assign(num_elements, Vec3(0.0, 0.0, 0.0));
  elementCentroids.assign(num_elements, Vec3(0.0, 0.0, 0.0));
  elementArea.assign(num_elements, 0.0);
}



void GeometryAnalyzer::compute_element_geometry(const Mesh& mesh, const QuadratureData& quad, int mint) {
  int num_nodes = mesh.p.size();
  int num_elements = mesh.n.size();

  Vec3 positionMomentSum(0.0, 0.0, 0.0);
  std::vector<Vec3> globalNodeNormals(num_nodes, Vec3(0.0, 0.0, 0.0));

#pragma omp parallel
  {
    Vec3 localPositionMoment(0.0, 0.0, 0.0);
    double localArea = 0.0, localVolume = 0.0;
    std::vector<Vec3> threadNodeNormals(num_nodes, Vec3(0.0, 0.0, 0.0));
    std::vector<int> threadContribCount(num_nodes, 0);

#pragma omp for schedule(static)
    for (int k = 0; k < num_elements; ++k) {
      std::vector<Vec3> nodes(6);
      for (int i = 0; i < 6; ++i)
        nodes[i] = mesh.p[mesh.n[k][i]];

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
        localVolume += x.dot(n) * cf;
      }

      localArea += arel;
      localPositionMoment += position_moment;

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
        threadNodeNormals[m] += n;
        threadContribCount[m]++;
      }

      Vec3 crv(0.0, 0.0, 0.0);
      auto cross_dx = [&](int i) { return vnorm[i].cross(dxdxi[i]); };
      auto cross_de = [&](int i) { return vnorm[i].cross(dxdeta[i]); };

      crv += al * cross_dx(0) + cross_dx(3) + (1.0 - al) * cross_dx(1);
      crv -= (1.0 - ga) * cross_dx(1) + cross_dx(4) + ga * cross_dx(2);
      crv += (1.0 - ga) * cross_de(1) + cross_de(4) + ga * cross_de(2);
      crv -= be * cross_de(0) + cross_de(5) + (1.0 - be) * cross_de(2);

      Vec3 x_c, dx_dxi, dx_deta, n_centroid;
      double hxi, het, hs;
      interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
               al, be, ga, 1.0 / 3.0, 1.0 / 3.0,
               x_c, dx_dxi, dx_deta, n_centroid, hxi, het, hs, 2);

      elementNormals[k] = n_centroid;
      elementCurvature[k] = 0.25 * crv.dot(n_centroid) / arel;
    }

#pragma omp critical
    {
      totalArea += localArea;
      totalVolume += localVolume;
      positionMomentSum += localPositionMoment;

      for (int i = 0; i < num_nodes; ++i) {
        globalNodeNormals[i] += threadNodeNormals[i];
        nodeContribCount[i] += threadContribCount[i];
      }
    }
  }

  surfaceCentroid = positionMomentSum / totalArea;
  totalVolume /= 3.0;
  nodeNormals = globalNodeNormals;
}







void GeometryAnalyzer::average_node_normals() {
  int num_nodes = nodeNormals.size();
#pragma omp parallel for
  for (int i = 0; i < num_nodes; ++i) {
    if (nodeContribCount[i] > 0) {
      nodeNormals[i] /= static_cast<double>(nodeContribCount[i]);
      if (nodeNormals[i].norm() > 1e-12)
        nodeNormals[i] = nodeNormals[i].normalized();
    }
  }
}


void GeometryAnalyzer::compute_node_curvature(const Mesh& mesh) {
  int num_nodes = mesh.p.size();
  int num_elements = mesh.n.size();

  std::vector<double> nodeCurvatureLocal(num_nodes, 0.0);
  std::vector<double> nodeAreaSumLocal(num_nodes, 0.0);

#pragma omp parallel
  {
    std::vector<double> curvPrivate(num_nodes, 0.0);
    std::vector<double> areaPrivate(num_nodes, 0.0);

#pragma omp for
    for (int elem = 0; elem < num_elements; ++elem) {
      double curv = elementCurvature[elem];
      double area = elementArea[elem];

      for (int local = 0; local < 6; ++local) {
        int global = mesh.n[elem][local];
        curvPrivate[global] += curv * area;
        areaPrivate[global] += area;
      }
    }

#pragma omp critical
    {
      for (int i = 0; i < num_nodes; ++i) {
        nodeCurvatureLocal[i] += curvPrivate[i];
        nodeAreaSumLocal[i] += areaPrivate[i];
      }
    }
  }

  nodeCurvature = nodeCurvatureLocal;

#pragma omp parallel for
  for (int i = 0; i < num_nodes; ++i) {
    if (nodeAreaSumLocal[i] > 1e-14)
      nodeCurvature[i] /= nodeAreaSumLocal[i];
    else
      nodeCurvature[i] = 0.0;
  }
}





void GeometryAnalyzer::compute_moment_matrix(const Mesh& mesh, const QuadratureData& quad, int mint) {
  int num_elements = mesh.n.size();
  Eigen::Matrix3d Mmat_accum = Eigen::Matrix3d::Zero();

#pragma omp parallel
  {
    Eigen::Matrix3d Mmat_local = Eigen::Matrix3d::Zero();

#pragma omp for
    for (int k = 0; k < num_elements; ++k) {
      std::vector<Vec3> nodes(6);
      for (int i = 0; i < 6; ++i)
        nodes[i] = mesh.p[mesh.n[k][i]];

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
            Mmat_local(m, n) += (r2 * delta - xhat[m] * xhat[n]) * cf;
          }
        }
      }
    }

#pragma omp critical
    {
      Mmat_accum += Mmat_local;
    }
  }

  Mmat = Mmat_accum.lu().inverse();
}


// void GeometryAnalyzer::compute(const Mesh& mesh, int quad_order) {
//   const int num_nodes = mesh.p.size();
//   const int num_elements = mesh.n.size();

//   totalArea = 0.0;
//   totalVolume = 0.0;
//   surfaceCentroid = Vec3(0.0, 0.0, 0.0);
//   Mmat.setZero();

//   QuadratureData quad = gauss_triangle(quad_order);
//   int mint = quad.xiq.size();

//   // Resize output containers
//   nodeNormals.assign(num_nodes, Vec3(0.0, 0.0, 0.0));
//   std::vector<int> nodeContribCount(num_nodes, 0);
//   elementCurvature.assign(num_elements, 0.0);
//   elementNormals.assign(num_elements, Vec3(0.0, 0.0, 0.0));
//   elementCentroids.assign(num_elements, Vec3(0.0, 0.0, 0.0));
//   elementArea.assign(num_elements, 0.0);

//   Vec3 positionMomentSum(0.0, 0.0, 0.0);

//   std::vector<Vec3> globalNodeNormals(num_nodes, Vec3(0.0, 0.0, 0.0));
//   std::vector<int> globalContribCount(num_nodes, 0);

//   // =========================
//   // Parallel Element Loop
//   // =========================
// #pragma omp parallel
//   {
//     Vec3 localPositionMoment(0.0, 0.0, 0.0);
//     double localArea = 0.0;
//     double localVolume = 0.0;

//     std::vector<Vec3> threadNodeNormals(num_nodes, Vec3(0.0, 0.0, 0.0));
//     std::vector<int> threadContribCount(num_nodes, 0);

// #pragma omp for schedule(static)
//     for (int k = 0; k < num_elements; ++k) {
//       std::vector<Vec3> nodes(6);
//       for (int i = 0; i < 6; ++i)
//         nodes[i] = mesh.p[mesh.n[k][i]];

//       double d42 = (nodes[3] - nodes[1]).norm();
//       double d41 = (nodes[3] - nodes[0]).norm();
//       double d63 = (nodes[5] - nodes[2]).norm();
//       double d61 = (nodes[5] - nodes[0]).norm();
//       double d52 = (nodes[4] - nodes[1]).norm();
//       double d53 = (nodes[4] - nodes[2]).norm();

//       double al = 1.0 / (1.0 + d42 / d41);
//       double be = 1.0 / (1.0 + d63 / d61);
//       double ga = 1.0 / (1.0 + d52 / d53);

//       double arel = 0.0;
//       Vec3 position_moment(0.0, 0.0, 0.0);

//       for (int i = 0; i < mint; ++i) {
//         Vec3 x, dx_dxi, dx_deta, n;
//         double hxi, het, hs;

//         interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
//                  al, be, ga, quad.xiq[i], quad.etq[i],
//                  x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

//         double cf = 0.5 * hs * quad.wwq[i];
//         arel += cf;
//         position_moment += cf * x;
//         localVolume += x.dot(n) * cf;
//       }

//       localArea += arel;
//       localPositionMoment += position_moment;

//       elementArea[k] = arel;
//       elementCentroids[k] = position_moment / arel;

//       std::array<double, 6> xxi = {0.0, 1.0, 0.0, al, ga, 0.0};
//       std::array<double, 6> eet = {0.0, 0.0, 1.0, 0.0, 1.0 - ga, be};
//       std::array<Vec3, 6> dxdxi, dxdeta, vnorm;

//       for (int i = 0; i < 6; ++i) {
//         Vec3 x, dx_dxi, dx_deta, n;
//         double hxi, het, hs;

//         interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
//                  al, be, ga, xxi[i], eet[i],
//                  x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

//         dxdxi[i] = dx_dxi;
//         dxdeta[i] = dx_deta;
//         vnorm[i] = n;

//         int m = mesh.n[k][i];
//         threadNodeNormals[m] += n;
//         threadContribCount[m]++;
//       }

//       Vec3 crv(0.0, 0.0, 0.0);
//       auto cross_dx = [&](int i) { return vnorm[i].cross(dxdxi[i]); };
//       auto cross_de = [&](int i) { return vnorm[i].cross(dxdeta[i]); };

//       crv += al * cross_dx(0) + cross_dx(3) + (1.0 - al) * cross_dx(1);
//       crv -= (1.0 - ga) * cross_dx(1) + cross_dx(4) + ga * cross_dx(2);
//       crv += (1.0 - ga) * cross_de(1) + cross_de(4) + ga * cross_de(2);
//       crv -= be * cross_de(0) + cross_de(5) + (1.0 - be) * cross_de(2);

//       Vec3 x_c, dx_dxi, dx_deta, n_centroid;
//       double hxi, het, hs;
//       interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
//                al, be, ga, 1.0 / 3.0, 1.0 / 3.0,
//                x_c, dx_dxi, dx_deta, n_centroid, hxi, het, hs, 2);

//       elementNormals[k] = n_centroid;
//       elementCurvature[k] = 0.25 * crv.dot(n_centroid) / arel;
//     }

// #pragma omp critical
//     {
//       totalArea += localArea;
//       totalVolume += localVolume;
//       positionMomentSum += localPositionMoment;

//       for (int i = 0; i < num_nodes; ++i) {
//         globalNodeNormals[i] += threadNodeNormals[i];
//         globalContribCount[i] += threadContribCount[i];
//       }
//     }
//   }

//   surfaceCentroid = positionMomentSum / totalArea;
//   totalVolume /= 3.0;
//   nodeNormals = globalNodeNormals;
//   nodeContribCount = globalContribCount;

//   // =========================
//   // Normalize Node Normals
//   // =========================
// #pragma omp parallel for
//   for (int i = 0; i < num_nodes; ++i) {
//     if (nodeContribCount[i] > 0) {
//       nodeNormals[i] /= static_cast<double>(nodeContribCount[i]);
//       if (nodeNormals[i].norm() > 1e-12)
//         nodeNormals[i] = nodeNormals[i].normalized();
//     }
//   }

//   // =========================
//   // Compute Node Curvature from Elements
//   // =========================
//   std::vector<double> nodeCurvatureLocal(num_nodes, 0.0);
//   std::vector<double> nodeAreaSumLocal(num_nodes, 0.0);

// #pragma omp parallel
//   {
//     std::vector<double> curvPrivate(num_nodes, 0.0);
//     std::vector<double> areaPrivate(num_nodes, 0.0);

// #pragma omp for
//     for (int elem = 0; elem < num_elements; ++elem) {
//       double curv = elementCurvature[elem];
//       double area = elementArea[elem];

//       for (int local = 0; local < 6; ++local) {
//         int global = mesh.n[elem][local];
//         curvPrivate[global] += curv * area;
//         areaPrivate[global] += area;
//       }
//     }

// #pragma omp critical
//     {
//       for (int i = 0; i < num_nodes; ++i) {
//         nodeCurvatureLocal[i] += curvPrivate[i];
//         nodeAreaSumLocal[i] += areaPrivate[i];
//       }
//     }
//   }

//   nodeCurvature = nodeCurvatureLocal;

//   // Normalize node curvature
// #pragma omp parallel for
//   for (int i = 0; i < num_nodes; ++i) {
//     if (nodeAreaSumLocal[i] > 1e-14)
//       nodeCurvature[i] /= nodeAreaSumLocal[i];
//     else
//       nodeCurvature[i] = 0.0;
//   }

//   // =========================
//   // Compute Mmat (Moment Matrix)
//   // =========================
//   Eigen::Matrix3d Mmat_accum = Eigen::Matrix3d::Zero();

// #pragma omp parallel
//   {
//     Eigen::Matrix3d Mmat_local = Eigen::Matrix3d::Zero();

// #pragma omp for
//     for (int k = 0; k < num_elements; ++k) {
//       std::vector<Vec3> nodes(6);
//       for (int i = 0; i < 6; ++i)
//         nodes[i] = mesh.p[mesh.n[k][i]];

//       double d42 = (nodes[3] - nodes[1]).norm();
//       double d41 = (nodes[3] - nodes[0]).norm();
//       double d63 = (nodes[5] - nodes[2]).norm();
//       double d61 = (nodes[5] - nodes[0]).norm();
//       double d52 = (nodes[4] - nodes[1]).norm();
//       double d53 = (nodes[4] - nodes[2]).norm();

//       double al = 1.0 / (1.0 + d42 / d41);
//       double be = 1.0 / (1.0 + d63 / d61);
//       double ga = 1.0 / (1.0 + d52 / d53);

//       for (int i = 0; i < mint; ++i) {
//         Vec3 x, dx_dxi, dx_deta, n;
//         double hxi, het, hs;
//         interp_p(nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5],
//                  al, be, ga, quad.xiq[i], quad.etq[i],
//                  x, dx_dxi, dx_deta, n, hxi, het, hs, 2);

//         double cf = 0.5 * hs * quad.wwq[i];
//         Vec3 xhat = x - surfaceCentroid;
//         double r2 = xhat.dot(xhat);

//         for (int m = 0; m < 3; ++m) {
//           for (int n = 0; n < 3; ++n) {
//             double delta = (m == n) ? 1.0 : 0.0;
//             Mmat_local(m, n) += (r2 * delta - xhat[m] * xhat[n]) * cf;
//           }
//         }
//       }
//     }

// #pragma omp critical
//     {
//       Mmat_accum += Mmat_local;
//     }
//   }

//   Mmat = Mmat_accum.lu().inverse();
// }


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



void GeometryAnalyzer::check_element_quality(const Mesh& mesh, int quad_order) const {
  QuadratureData quad = gauss_triangle(quad_order);
  int mint = quad.xiq.size();
  int num_elements = mesh.n.size();

  double global_min_area = std::numeric_limits<double>::max();
  double global_max_area = 0.0;
  double global_min_aspect = std::numeric_limits<double>::max();
  double global_max_aspect = 0.0;
  int global_flipped = 0;

#pragma omp parallel
  {
    double min_area = std::numeric_limits<double>::max();
    double max_area = 0.0;
    double min_aspect = std::numeric_limits<double>::max();
    double max_aspect = 0.0;
    int flipped_count = 0;

#pragma omp for
    for (int k = 0; k < num_elements; ++k) {
      const auto& tri = mesh.n[k];
      const Vec3& A = mesh.p[tri[0]];
      const Vec3& B = mesh.p[tri[1]];
      const Vec3& C = mesh.p[tri[2]];

      // Aspect ratio via side lengths and area
      double a = (B - A).norm();
      double b = (C - B).norm();
      double c = (A - C).norm();
      double s = 0.5 * (a + b + c);
      double tri_area = std::sqrt(std::max(s * (s - a) * (s - b) * (s - c), 0.0)); // Heron’s formula

      double h = (tri_area > 1e-14) ? (2.0 * tri_area / std::max({a, b, c})) : 1e-14;
      double aspect = std::max({a, b, c}) / h;

      min_aspect = std::min(min_aspect, aspect);
      max_aspect = std::max(max_aspect, aspect);

      // Area
      double area = elementArea[k];
      min_area = std::min(min_area, area);
      max_area = std::max(max_area, area);

      // Jacobian positivity check
      std::vector<Vec3> nodes(6);
      for (int i = 0; i < 6; ++i)
        nodes[i] = mesh.p[mesh.n[k][i]];

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
        if (hs <= 0.0)
          flipped_count++;
      }
    }

#pragma omp critical
    {
      global_min_area = std::min(global_min_area, min_area);
      global_max_area = std::max(global_max_area, max_area);
      global_min_aspect = std::min(global_min_aspect, min_aspect);
      global_max_aspect = std::max(global_max_aspect, max_aspect);
      global_flipped += flipped_count;
    }
  }

  // Print summary
  std::cout << "\n=== Element Quality Diagnostics ===\n";
  std::cout << "Min Element Area: " << global_min_area << "\n";
  std::cout << "Max Element Area: " << global_max_area << "\n";
  std::cout << "Min Aspect Ratio: " << global_min_aspect << "\n";
  std::cout << "Max Aspect Ratio: " << global_max_aspect << "\n";
  if (global_flipped > 0)
    std::cout << "⚠️  " << global_flipped << " elements have non-positive Jacobians (inverted/degenerate).\n";
  else
    std::cout << "All Jacobians are positive. ✅\n";
}
