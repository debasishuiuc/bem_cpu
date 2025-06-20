// mesh_utils.cpp

#include "mesh_utils.hpp"
#include "vec3.hpp"
#include <omp.h>
#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <parallel/algorithm>
#include <limits>

// Vec3 comparison helpers
bool compare_vec3(const Vec3& a, const Vec3& b) {
  for (int i = 0; i < 3; ++i) {
    if (fabs(a[i] - b[i]) > 1e-12) return a[i] < b[i];
  }
  return false;
}

bool equal_vec3(const Vec3& a, const Vec3& b) {
  return fabs(a[0] - b[0]) < 1e-12 &&
    fabs(a[1] - b[1]) < 1e-12 &&
    fabs(a[2] - b[2]) < 1e-12;
}

// Collect unique edges from triangle list
void collect_edges(const std::vector<std::array<int, 3>>& triangles,
                   std::vector<Edge>& edges) {
  std::unordered_map<Edge, int, pair_hash> edge_set;
  for (const auto& tri : triangles) {
    Edge e1 = std::minmax(tri[0], tri[1]);
    Edge e2 = std::minmax(tri[1], tri[2]);
    Edge e3 = std::minmax(tri[2], tri[0]);
    edge_set[e1] = 1;
    edge_set[e2] = 1;
    edge_set[e3] = 1;
  }

  edges.clear();
  for (const auto& kv : edge_set) {
    edges.push_back(kv.first);
  }
}

// Euclidean distance between Vec3s
inline double dist(const Vec3& a, const Vec3& b) {
  return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) +
                   (a[1] - b[1]) * (a[1] - b[1]) +
                   (a[2] - b[2]) * (a[2] - b[2]));
}
