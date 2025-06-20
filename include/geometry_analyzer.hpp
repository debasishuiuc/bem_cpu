// geometry_analyzer.hpp

#pragma once

#include "vec3.hpp"
#include "mesh_utils.hpp"
#include <vector>
#include <Eigen/Dense>

class GeometryAnalyzer {
public:
  // Member variables to hold geometry outputs
  std::vector<Vec3> elementNormals;     // 1. Normal vector per element
  std::vector<double> elementCurvature; // 2. Curvature per element
  std::vector<Vec3> nodeNormals;        // 3. Averaged normal vector per node
  std::vector<double> elementArea;      // 4. Area per element
  std::vector<Vec3> elementCentroids;   // 5. Centroid per element
  std::vector<double> nodeCurvature;    // 6. Curvature interpolated on nodes
  double totalArea = 0.0;               // 7. Total surface area
  double totalVolume = 0.0;             // 8. Enclosed volume
  Vec3 surfaceCentroid = Vec3();        // 9. Centroid of the surface
  Eigen::Matrix3d Mmat = Eigen::Matrix3d::Zero(); // 10. Moment of inertia matrix

  // Constructor
  GeometryAnalyzer() = default;

  // Compute method: main entry point
  void compute(const Mesh& mesh, int quad_order);

  // Optional: Print summary
  void print_summary(const Mesh& mesh, int num_print = 5) const;

  void check_element_quality(const Mesh& mesh, int quad_order) const;
};
