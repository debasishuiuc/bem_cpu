// geometry.hpp

#pragma once
#include "mesh_utils.hpp"
#include <vector>
#include <Eigen/Dense>

struct Geometry {
  double area = 0.0;
  double volume = 0.0;
  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  std::vector<double> elemArea;
  std::vector<Eigen::Vector3d> elemCentroid;
  std::vector<Eigen::Vector3d> elemNormals;   // vnc
  std::vector<Eigen::Vector3d> nodeNormals;   // vna
  std::vector<double> elemCurvature;          // crvmel
  Eigen::Matrix3d Mmat = Eigen::Matrix3d::Zero();
};


void compute_geometry(const Mesh& mesh, int quad_order);
