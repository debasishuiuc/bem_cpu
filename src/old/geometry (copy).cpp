// geometry.cpp

#include "geometry.hpp"
#include "quadrature.hpp"
#include "mesh_utils.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <omp.h>
#include <set>

using Eigen::Vector3d;
using Eigen::Matrix3d;

// Helper function to compute Euclidean distance between two Vec3
inline double dist(const Vec3& a, const Vec3& b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) +
                   (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]));
}

// Diagnostic function
void check_mesh_integrity(const Mesh& mesh, double near_duplicate_tol, double degenerate_tol) {
  const auto& points = mesh.p;
  const auto& elements = mesh.n;

  std::cout << "=== Mesh Integrity Diagnostics ===\n";

  // Hashable voxel key
  struct VoxelKey {
    int x, y, z;
    bool operator==(const VoxelKey& other) const {
      return x == other.x && y == other.y && z == other.z;
    }
  };

  struct VoxelKeyHash {
    std::size_t operator()(const VoxelKey& key) const {
      return std::hash<int>()(key.x) ^ std::hash<int>()(key.y << 1) ^ std::hash<int>()(key.z << 2);
    }
  };

  int near_duplicates = 0;
  double voxel_size = near_duplicate_tol;

  // Map from voxel -> list of point indices
  std::unordered_map<VoxelKey, std::vector<int>, VoxelKeyHash> voxel_map;

  // Insert points into voxel grid
  for (int i = 0; i < (int)points.size(); ++i) {
    const Vec3& p = points[i];
    VoxelKey key = {
      static_cast<int>(std::floor(p[0] / voxel_size)),
      static_cast<int>(std::floor(p[1] / voxel_size)),
      static_cast<int>(std::floor(p[2] / voxel_size))
    };
    voxel_map[key].push_back(i);
  }

  // Now search each voxel and its 26 neighbors
  for (const auto& [key, indices] : voxel_map) {
    for (int i : indices) {
      const Vec3& pi = points[i];

      for (int dx = -1; dx <= 1; ++dx)
	for (int dy = -1; dy <= 1; ++dy)
	  for (int dz = -1; dz <= 1; ++dz) {
	    VoxelKey neighbor = {key.x + dx, key.y + dy, key.z + dz};
	    if (!voxel_map.count(neighbor)) continue;

	    for (int j : voxel_map[neighbor]) {
	      if (j <= i) continue;  // avoid double count
	      const Vec3& pj = points[j];
	      if (dist(pi, pj) < near_duplicate_tol) {
		++near_duplicates;
	      }
	    }
	  }
    }
  }

  std::cout << "Near-duplicate points (within " << near_duplicate_tol << "): " << near_duplicates << "\n";


  // Check 2: Degenerate triangles
  int degenerate_count = 0;
  double min_area = std::numeric_limits<double>::max();

#pragma omp parallel for reduction(+:degenerate_count) reduction(min:min_area)
  for (int i = 0; i < (int)elements.size(); ++i) {
    const auto& tri = elements[i];
    const Vec3& A = points[tri[0]];
    const Vec3& B = points[tri[1]];
    const Vec3& C = points[tri[2]];

    // Compute triangle area using Heron's formula or cross product
    Vec3 AB = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
    Vec3 AC = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};
    Vec3 cross = {
      AB[1]*AC[2] - AB[2]*AC[1],
      AB[2]*AC[0] - AB[0]*AC[2],
      AB[0]*AC[1] - AB[1]*AC[0]
    };
    double area = 0.5 * std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);

    if (area < degenerate_tol)
      ++degenerate_count;

    if (area < min_area)
      min_area = area;
  }

  std::cout << "Degenerate triangles (area < " << degenerate_tol << "): " << degenerate_count << "\n";
  std::cout << "Minimum triangle area in mesh: " << min_area << "\n";

  // Check 3: Minimum edge length
  double min_edge = std::numeric_limits<double>::max();
#pragma omp parallel for reduction(min:min_edge)
  for (int i = 0; i < (int)elements.size(); ++i) {
    const auto& tri = elements[i];
    const double l1 = dist(points[tri[0]], points[tri[1]]);
    const double l2 = dist(points[tri[1]], points[tri[2]]);
    const double l3 = dist(points[tri[2]], points[tri[0]]);
    min_edge = std::min({min_edge, l1, l2, l3});
  }

  std::cout << "Minimum edge length: " << min_edge << "\n";
  std::cout << "Suggested deduplication tolerance: ~" << min_edge / 2.0 << "\n";
  std::cout << "====================================\n";
}

void fix_triangle_orientation(Mesh& mesh) {
  size_t flipped_total = 0;

#pragma omp parallel
  {
    size_t flipped_local = 0;

#pragma omp for
    for (size_t i = 0; i < mesh.n.size(); ++i) {
      const auto& elem = mesh.n[i];

      // Sanity check: ensure the 6-node triangle doesn't have repeated nodes
      std::set<int> unique(elem.begin(), elem.end());
      if (unique.size() < 6) {
        std::cerr << "‼️ Degenerate element detected before flip: element " << i << "\n";
        for (int j = 0; j < 6; ++j) {
          const Vec3& pt = mesh.p[elem[j]];
          std::cerr << "  [" << j << "] " << pt.x << " " << pt.y << " " << pt.z << "\n";
        }
        continue;
      }

      const Vec3& v1_ = mesh.p[elem[0]];
      const Vec3& v2_ = mesh.p[elem[2]];
      const Vec3& v3_ = mesh.p[elem[4]];
      Eigen::Vector3d v1(v1_.x, v1_.y, v1_.z);
      Eigen::Vector3d v2(v2_.x, v2_.y, v2_.z);
      Eigen::Vector3d v3(v3_.x, v3_.y, v3_.z);

      Eigen::Vector3d normal = (v2 - v1).cross(v3 - v1);
      Eigen::Vector3d center = (v1 + v2 + v3) / 3.0;

      if (normal.dot(center) < 0.0) {
        std::array<int, 6> flipped_elem = elem;
        std::swap(flipped_elem[0], flipped_elem[2]);
        std::swap(flipped_elem[1], flipped_elem[5]);
        mesh.n[i] = flipped_elem;

        // Double-check if it became degenerate after flipping
        std::set<int> unique2(flipped_elem.begin(), flipped_elem.end());
        if (unique2.size() < 6) {
          std::cerr << "‼️ Degenerate element detected after flip: element " << i << "\n";
          for (int j = 0; j < 6; ++j) {
            const Vec3& pt = mesh.p[flipped_elem[j]];
            std::cerr << "  [" << j << "] " << pt.x << " " << pt.y << " " << pt.z << "\n";
          }
          continue;  // Skip counting as flipped
        }

        ++flipped_local;
      }
    }

#pragma omp atomic
    flipped_total += flipped_local;
  }

  std::cout << "🔄 Fixed orientation of " << flipped_total << " elements (if any)." << std::endl;
}



void precompute_shape_functions(int Nq, double al, double be, double ga, ShapeFunctionData& sfd) {
  QuadratureData quad = gauss_triangle(Nq);
  const auto& xiq = quad.xiq;
  const auto& etq = quad.etq;
  int mint = xiq.size();

  assert(sfd.phi.size() == mint);
  assert(sfd.dphidxi.size() == mint);
  assert(sfd.dphideta.size() == mint);


  double alc = 1.0 - al;
  double bec = 1.0 - be;
  double gac = 1.0 - ga;
  double alalc = al * alc;
  double bebec = be * bec;
  double gagac = ga * gac;

  for (int q = 0; q < mint; ++q) {
    double xi = xiq[q], eta = etq[q];
    auto& phi = sfd.phi[q];
    auto& dphidxi = sfd.dphidxi[q];
    auto& dphideta = sfd.dphideta[q];

    phi[1] = xi * (xi - al + eta * (al - ga) / gac) / alc;
    phi[2] = eta * (eta - be + xi * (be + ga - 1.0) / ga) / bec;
    phi[3] = xi * (1.0 - xi - eta) / alalc;
    phi[4] = xi * eta / gagac;
    phi[5] = eta * (1.0 - xi - eta) / bebec;
    phi[0] = 1.0 - phi[1] - phi[2] - phi[3] - phi[4] - phi[5];

    dphidxi[1] = (2.0 * xi - al + eta * (al - ga) / gac) / alc;
    dphidxi[2] = eta * (be + ga - 1.0) / (ga * bec);
    dphidxi[3] = (1.0 - 2.0 * xi - eta) / alalc;
    dphidxi[4] = eta / gagac;
    dphidxi[5] = -eta / bebec;
    dphidxi[0] = -dphidxi[1] - dphidxi[2] - dphidxi[3] - dphidxi[4] - dphidxi[5];

    dphideta[1] = xi * (al - ga) / (alc * gac);
    dphideta[2] = (2.0 * eta - be + xi * (be + ga - 1.0) / ga) / bec;
    dphideta[3] = -xi / alalc;
    dphideta[4] = xi / gagac;
    dphideta[5] = (1.0 - xi - 2.0 * eta) / bebec;
    dphideta[0] = -dphideta[1] - dphideta[2] - dphideta[3] - dphideta[4] - dphideta[5];
  }
}

inline Vector3d curl(const Vector3d& a, const Vector3d& b) {
  return {
    a.y() * b.z() - a.z() * b.y(),
    a.z() * b.x() - a.x() * b.z(),
    a.x() * b.y() - a.y() * b.x()
  };
}

double compute_element_curvature(
				 const std::array<Vector3d, 6>& DxDxi,
				 const std::array<Vector3d, 6>& DxDet,
				 const std::array<Vector3d, 6>& vn,
				 const Vector3d& vn_centroid,
				 double al, double be, double ga,
				 double alc, double bec, double gac,
				 double elem_area)
{
  Vector3d crv = al  * curl(vn[0], DxDxi[0]) +
    curl(vn[3], DxDxi[3]) +
    alc * curl(vn[1], DxDxi[1]);

  crv -= gac * curl(vn[1], DxDxi[1]) +
    curl(vn[4], DxDxi[4]) +
    ga  * curl(vn[2], DxDxi[2]);

  crv += gac * curl(vn[1], DxDet[1]) +
    curl(vn[4], DxDet[4]) +
    ga  * curl(vn[2], DxDet[2]);

  crv -= be  * curl(vn[0], DxDet[0]) +
    curl(vn[5], DxDet[5]) +
    bec * curl(vn[2], DxDet[2]);

  return 0.25 * crv.dot(vn_centroid) / elem_area;
}

void compute_geometry(const Mesh& mesh, Geometry& geom, int Nq) {
  const int Npts = mesh.Npts, Nelm = mesh.Nelm;
  const auto& p = mesh.p;
  const auto& n = mesh.n;

  QuadratureData quad = gauss_triangle(Nq);
  const auto& xiq = quad.xiq;
  const auto& etq = quad.etq;
  const auto& wq = quad.wwq;
  int mint = xiq.size();

  geom.elemArea.assign(Nelm, 0.0);
  geom.elemCentroid.assign(Nelm, Vector3d::Zero());
  geom.elemNormals.assign(Nelm, Vector3d::Zero());
  geom.nodeNormals.assign(Npts, Vector3d::Zero());
  geom.elemCurvature.assign(Nelm, 0.0);

  std::vector<Vector3d> vna(Npts, Vector3d::Zero());
  std::vector<int> tally(Npts, 0);
  std::vector<double> alp(Nelm), bet(Nelm), gam(Nelm);
  std::vector<ShapeFunctionData> sfd_all(Nelm);
  for (int k = 0; k < Nelm; ++k) {
    sfd_all[k].phi.resize(mint);
    sfd_all[k].dphidxi.resize(mint);
    sfd_all[k].dphideta.resize(mint);
  }

  double total_area = 0.0, total_volume = 0.0;
  Vector3d total_centroid = Vector3d::Zero();

  // Preallocate thread-local accumulators
  int nthreads = omp_get_max_threads();
  std::vector<std::vector<Vector3d>> vna_thread(nthreads, std::vector<Vector3d>(Npts, Vector3d::Zero()));
  std::vector<std::vector<int>> tally_thread(nthreads, std::vector<int>(Npts, 0));
  std::vector<double> area_thread(nthreads, 0.0);
  std::vector<double> volume_thread(nthreads, 0.0);
  std::vector<Vector3d> centroid_thread(nthreads, Vector3d::Zero());

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    auto& vna_private = vna_thread[tid];
    auto& tally_private = tally_thread[tid];
    double& area_private = area_thread[tid];
    double& volume_private = volume_thread[tid];
    Vector3d& centroid_private = centroid_thread[tid];

#pragma omp for
    for (int k = 0; k < Nelm; ++k) {
      const auto& ni = n[k];
      std::array<Vector3d, 6> X;
      get_element_coords_eigen(mesh, k, X);

      // check for elements on sphere or not
      for (int i = 0; i < 6; ++i) {
	double r = X[i].norm();
	if (std::abs(r - 1.0) > 1e-6) {
	  std::cout << "⚠️ Non-unit node at element " << k << ", node " << i << ", |X| = " << r << std::endl;
	}
      }
      
      double d42 = (X[3] - X[1]).norm();
      double d41 = (X[3] - X[0]).norm();
      double d63 = (X[5] - X[2]).norm();
      double d61 = (X[5] - X[0]).norm();
      double d52 = (X[4] - X[1]).norm();
      double d53 = (X[4] - X[2]).norm();

      double al = 1.0 / (1.0 + d42 / d41);
      double be = 1.0 / (1.0 + d63 / d61);
      double ga = 1.0 / (1.0 + d52 / d53);
      double alc = 1.0 - al, bec = 1.0 - be, gac = 1.0 - ga;

      alp[k] = al; bet[k] = be; gam[k] = ga;

      ShapeFunctionData& sfd = sfd_all[k];
      precompute_shape_functions(Nq, al, be, ga, sfd);

      double ak = 0.0;
      Vector3d xmom = Vector3d::Zero();

      for (int q = 0; q < mint; ++q) {
        Vector3d x = Vector3d::Zero(), DxDxi = Vector3d::Zero(), DxDet = Vector3d::Zero();
        for (int i = 0; i < 6; ++i) {
          x += sfd.phi[q][i] * X[i];
          DxDxi += sfd.dphidxi[q][i] * X[i];
          DxDet += sfd.dphideta[q][i] * X[i];
        }
        Vector3d vn = DxDxi.cross(DxDet);
        double hs = vn.norm();
        vn.normalize();

        double cf = hs * wq[q];
        ak += cf;
        xmom += cf * x;
        volume_private += x.dot(vn) * cf;
      }

      ak *= 0.5;
      xmom *= 0.5;
      geom.elemArea[k] = ak;
      geom.elemCentroid[k] = xmom / ak;
      area_private += ak;
      centroid_private += xmom;

      std::array<Vector3d, 6> DxDxi_eval, DxDet_eval, vn_eval;

      for (int i = 0; i < 6; ++i) {
        Vector3d DxDxi = Vector3d::Zero(), DxDet = Vector3d::Zero();
        for (int j = 0; j < 6; ++j) {
          DxDxi += sfd.dphidxi[0][j] * X[j];
          DxDet += sfd.dphideta[0][j] * X[j];
        }
        Vector3d vn = DxDxi.cross(DxDet).normalized();
        DxDxi_eval[i] = DxDxi;
        DxDet_eval[i] = DxDet;
        vn_eval[i] = vn;
        vna_private[ni[i]] += vn;
        tally_private[ni[i]]++;
      }

      Vector3d vn_centroid = -(DxDxi_eval[0].cross(DxDet_eval[0])).normalized();
      geom.elemNormals[k] = vn_centroid;

      geom.elemCurvature[k] = compute_element_curvature(
        DxDxi_eval, DxDet_eval, vn_eval, vn_centroid,
        al, be, ga, alc, bec, gac, ak);
    }
  }

  // Manual reduction after parallel region
  for (int t = 0; t < nthreads; ++t) {
    total_area += area_thread[t];
    total_volume += volume_thread[t];
    total_centroid += centroid_thread[t];
    for (int i = 0; i < Npts; ++i) {
      vna[i] += vna_thread[t][i];
      tally[i] += tally_thread[t][i];
    }
  }

  geom.area = total_area;
  geom.volume = 0.5 * total_volume / 3.0;
  geom.centroid = total_centroid / total_area;

#pragma omp parallel for
  for (int i = 0; i < Npts; ++i) {
    if (tally[i] > 0) {
      geom.nodeNormals[i] = -vna[i] / tally[i];
      geom.nodeNormals[i].normalize();
    }
  }

  Matrix3d M = Matrix3d::Zero();

  std::vector<Matrix3d> M_thread(nthreads, Matrix3d::Zero());

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    Matrix3d& M_private = M_thread[tid];

#pragma omp for
    for (int k = 0; k < Nelm; ++k) {
      std::array<Vector3d, 6> X;
      get_element_coords_eigen(mesh, k, X);
      const ShapeFunctionData& sfd = sfd_all[k];
      double al = alp[k], be = bet[k], ga = gam[k];

      for (int q = 0; q < mint; ++q) {
	Vector3d x = Vector3d::Zero(), DxDxi = Vector3d::Zero(), DxDet = Vector3d::Zero();
	for (int i = 0; i < 6; ++i) {
	  x += sfd.phi[q][i] * X[i];
	  DxDxi += sfd.dphidxi[q][i] * X[i];
	  DxDet += sfd.dphideta[q][i] * X[i];
	}

	Vector3d vn = DxDxi.cross(DxDet);
	double hs = vn.norm();
	double cf = 0.5 * hs * wq[q];
	Vector3d d = x - geom.centroid;
	double r2 = d.squaredNorm();

	M_private(0, 0) += (r2 - d[0]*d[0]) * cf;
	M_private(0, 1) += (-d[0]*d[1]) * cf;
	M_private(0, 2) += (-d[0]*d[2]) * cf;
	M_private(1, 0) += (-d[1]*d[0]) * cf;
	M_private(1, 1) += (r2 - d[1]*d[1]) * cf;
	M_private(1, 2) += (-d[1]*d[2]) * cf;
	M_private(2, 0) += (-d[2]*d[0]) * cf;
	M_private(2, 1) += (-d[2]*d[1]) * cf;
	M_private(2, 2) += (r2 - d[2]*d[2]) * cf;
      }
    }
  }

  // Reduce M_thread into M
  for (int t = 0; t < nthreads; ++t) {
    M += M_thread[t];
  }

  geom.Mmat = M.inverse();

}
