// deduplicate.cpp

// #include "mesh_utils.hpp"
// #include "vec3.hpp"
// #include <omp.h>
// #include <algorithm>
// #include <unordered_map>


// void deduplicate_points(Mesh& mesh) {
//   const auto& all_points = mesh.p;

//   // Step 1: Deduplicate points
//   std::vector<Vec3> unique_points = all_points;
//   std::sort(unique_points.begin(), unique_points.end(), compare_vec3);
//   unique_points.erase(
// 		      std::unique(unique_points.begin(), unique_points.end(), equal_vec3),
// 		      unique_points.end());

//   // Step 2: Build mapping from old point to new index
//   std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal> point_to_id;
//   for (int i = 0; i < static_cast<int>(unique_points.size()); ++i) {
//     point_to_id[unique_points[i]] = i;
//   }

//   // Step 3: Create remap vector (old index → new index)
//   std::vector<int> remap_indices(mesh.p.size(), -1);
// #pragma omp parallel for
//   for (int i = 0; i < static_cast<int>(mesh.p.size()); ++i) {
//     auto it = point_to_id.find(mesh.p[i]);
//     if (it != point_to_id.end()) {
//       remap_indices[i] = it->second;
//     } else {
// #pragma omp critical
//       std::cerr << "ERROR: Could not find point in deduplication map.\n";
//     }
//   }

//   // Step 4: Remap element connectivity
// #pragma omp parallel for
//   for (int k = 0; k < static_cast<int>(mesh.n.size()); ++k) {
//     for (int i = 0; i < 6; ++i) {
//       mesh.n[k][i] = remap_indices[mesh.n[k][i]];
//     }
//   }

//   // Step 5: Update mesh
//   mesh.p = std::move(unique_points);
//   mesh.Npts = mesh.p.size();
// }


#include "mesh_utils.hpp"
#include "vec3.hpp"
#include "trgl6_icos.hpp"  // for Mesh
#include <omp.h>
#include <algorithm>
#include <unordered_map>

void deduplicate_points(Mesh& mesh) {
  const auto& inputCoords = mesh.nodeCoords;

  // Step 1: Deduplicate coordinates
  std::vector<Vec3> uniqueCoords = inputCoords;
  std::sort(uniqueCoords.begin(), uniqueCoords.end(), compare_vec3);
  uniqueCoords.erase(
    std::unique(uniqueCoords.begin(), uniqueCoords.end(), equal_vec3),
    uniqueCoords.end()
  );

  // Step 2: Map unique coordinate to new index
  std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal> coordToIndex;
  for (int i = 0; i < static_cast<int>(uniqueCoords.size()); ++i) {
    coordToIndex[uniqueCoords[i]] = i;
  }

  // Step 3: Remap old indices to new indices
  std::vector<int> remap(mesh.nodeCoords.size(), -1);
#pragma omp parallel for
  for (int i = 0; i < static_cast<int>(mesh.nodeCoords.size()); ++i) {
    auto it = coordToIndex.find(mesh.nodeCoords[i]);
    if (it != coordToIndex.end()) {
      remap[i] = it->second;
    } else {
#pragma omp critical
      std::cerr << "ERROR: Could not find coordinate in deduplication map.\n";
    }
  }

  // Step 4: Apply remap to element connectivity
#pragma omp parallel for
  for (int elem = 0; elem < static_cast<int>(mesh.elementNodes.size()); ++elem) {
    for (int i = 0; i < 6; ++i) {
      mesh.elementNodes[elem][i] = remap[mesh.elementNodes[elem][i]];
    }
  }

  // Step 5: Update mesh
  mesh.nodeCoords = std::move(uniqueCoords);
  mesh.numNodes = mesh.nodeCoords.size();
}
