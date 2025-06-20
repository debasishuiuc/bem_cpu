// deduplicate.cpp

#include "mesh_utils.hpp"
#include "vec3.hpp"
#include <omp.h>
#include <algorithm>
#include <unordered_map>

// Vec3 comparison helpers (can optionally be factored into a common utility)
extern bool compare_vec3(const Vec3& a, const Vec3& b);
extern bool equal_vec3(const Vec3& a, const Vec3& b);

// Deduplicate point list with mapping
void deduplicate_points(const std::vector<Vec3>& all_points,
                        std::vector<Vec3>& unique_points,
                        std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal>& point_to_id) {
  unique_points = all_points;
  std::sort(unique_points.begin(), unique_points.end(), compare_vec3);

  unique_points.erase(std::unique(unique_points.begin(), unique_points.end(), equal_vec3),
                      unique_points.end());

#pragma omp parallel
  {
    std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal> local_map;

#pragma omp for nowait
    for (int i = 0; i < static_cast<int>(all_points.size()); ++i) {
      const Vec3& p = all_points[i];
      auto it = std::lower_bound(unique_points.begin(), unique_points.end(), p, compare_vec3);
      if (it != unique_points.end() && equal_vec3(*it, p)) {
        int idx = std::distance(unique_points.begin(), it);
        local_map[p] = idx;
      }
    }

#pragma omp critical
    for (const auto& kv : local_map)
      point_to_id[kv.first] = kv.second;
  }
}

// Deduplicate and remap mesh node indices
void deduplicate_points(Mesh& mesh) {
  std::vector<Vec3> unique_points;
  std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal> point_to_id;

  deduplicate_points(mesh.p, unique_points, point_to_id);

  for (auto& tri : mesh.n) {
    for (int& node : tri) {
      node = point_to_id[mesh.p[node]];
    }
  }

  mesh.p = std::move(unique_points);
  mesh.Npts = mesh.p.size();
}
