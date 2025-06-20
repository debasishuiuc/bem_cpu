// mesh_utils.hpp

#pragma once

#include "trgl6_icos.hpp"
#include "vec3.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>

// Hash and equality for Vec3
struct Vec3Hash {
  std::size_t operator()(const Vec3& v) const {
    auto h1 = std::hash<double>{}(v[0]);
    auto h2 = std::hash<double>{}(v[1]);
    auto h3 = std::hash<double>{}(v[2]);
    return h1 ^ (h2 << 1) ^ (h3 << 2);
  }
};

struct Vec3Equal {
  bool operator()(const Vec3& a, const Vec3& b) const {
    return std::fabs(a[0] - b[0]) < 1e-12 &&
      std::fabs(a[1] - b[1]) < 1e-12 &&
      std::fabs(a[2] - b[2]) < 1e-12;
  }
};

// Define Edge as a pair of ints
using Edge = std::pair<int, int>;

// Hash for Edge (needed for CPU mesh operations and CUDA kernels)
struct pair_hash {
  std::size_t operator()(const std::pair<int, int>& p) const noexcept {
    return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
  }
};

// Collect unique edges from 3-node triangles
void collect_edges(const std::vector<std::array<int, 3>>& triangles,
                   std::vector<Edge>& edges);

// Summary print
inline void print_mesh_summary(const Mesh& mesh, const std::string& label) {
  std::cout << "Mesh generation complete (" << label << ").\n";
  std::cout << "Subdivision level (Ndiv): " << mesh.Ndiv << '\n';
  std::cout << "Points: " << mesh.Npts << ", Elements: " << mesh.Nelm << '\n';
}


bool check_mesh_integrity(const Mesh& mesh);
  
bool check_triangle_orientation_strict(const Mesh& mesh, double threshold = 1e-12);
