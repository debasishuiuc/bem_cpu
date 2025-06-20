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
#include <set>

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

bool check_mesh_integrity(const Mesh& mesh) {
  bool all_good = true;

#pragma omp parallel for
  for (int k = 0; k < mesh.n.size(); ++k) {
    std::set<int> node_indices;
    bool duplicate = false;

    for (int i = 0; i < 6; ++i) {
      int idx = mesh.n[k][i];
      if (node_indices.count(idx)) {
#pragma omp critical
        {
          std::cerr << "Duplicate node index " << idx
                    << " in element " << k << ": [";
          for (int j = 0; j < 6; ++j) std::cerr << mesh.n[k][j] << " ";
          std::cerr << "]\n";
        }
        duplicate = true;
        break;
      }
      node_indices.insert(idx);
    }

    if (duplicate) {
#pragma omp critical
      all_good = false;
    }
  }

  return all_good;
}


bool check_triangle_orientation_strict(const Mesh& mesh, double threshold) {
  bool all_ccw = true;

#pragma omp parallel
  {
    bool local_ok = true;

#pragma omp for nowait
    for (int k = 0; k < mesh.n.size(); ++k) {
      const auto& tri = mesh.n[k];

      // Use nodes 0,1,2 to compute orientation
      Vec3 A = mesh.p[tri[0]];
      Vec3 B = mesh.p[tri[1]];
      Vec3 C = mesh.p[tri[2]];

      Vec3 normal = (B - A).cross(C - A).normalized();  // CCW normal
      Vec3 centroid = (A + B + C) / 3.0;
      Vec3 radial = centroid.normalized();              // Outward vector

      double dot = normal.dot(radial);

      if (dot <= threshold) {
#pragma omp critical
        {
          std::cerr << "⚠️ Triangle " << k << " is NOT counter-clockwise (dot = " << dot << ").\n";
        }
        local_ok = false;
      }
    }

#pragma omp critical
    all_ccw = all_ccw && local_ok;
  }

  if (all_ccw) {
    std::cout << "✅ All triangles are correctly oriented (CCW).\n";
  } else {
    std::cerr << "❌ Some triangles are not CCW. Please fix orientation errors.\n";
  }

  return all_ccw;
}
