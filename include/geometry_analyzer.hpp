// geometry_analyzer.hpp

#pragma once

#include "vec3.hpp"
#include "mesh_utils.hpp"
#include "quadrature.hpp"
#include <vector>
#include <Eigen/Dense>

class GeometryAnalyzer {
public:
  // === Outputs ===
  std::vector<Vec3> elementNormals;     // [1] One normal per element (outward-facing)
  std::vector<double> elementCurvature; // [2] Mean curvature per element
  std::vector<Vec3> nodeNormals;        // [3] Averaged normal vector per node
  std::vector<double> elementArea;      // [4] Surface area per element
  std::vector<Vec3> elementCentroids;   // [5] Centroid of each element
  std::vector<double> nodeCurvature;    // [6] Node-wise interpolated curvature
  double totalArea = 0.0;               // [7] Total surface area
  double totalVolume = 0.0;             // [8] Enclosed volume (signed)
  Vec3 surfaceCentroid = Vec3(0.0, 0.0, 0.0); // [9] Surface centroid (area-weighted)
  Eigen::Matrix3d Mmat = Eigen::Matrix3d::Zero(); // [10] Moment of inertia matrix (default to zero)

  // === Constructor ===
  GeometryAnalyzer() = default;

  // === Main Entry Point ===
  void compute(const Mesh& mesh, int quad_order);

  // === Diagnostics and Debugging ===
  void print_summary(const Mesh& mesh, int num_print = 5) const;
  void check_element_quality(const Mesh& mesh, int quad_order) const;

private:
  // === Internal Helpers ===
  void initialize_containers(const Mesh& mesh, int mint);
  void compute_element_geometry(const Mesh& mesh, const QuadratureData& quad, int mint);
  void average_node_normals();
  void compute_node_curvature(const Mesh& mesh);
  void compute_moment_matrix(const Mesh& mesh, const QuadratureData& quad, int mint);

  std::vector<int> nodeContribCount;  // Helper: number of contributions per node
};

