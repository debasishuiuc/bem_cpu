// // deduplicate.hpp

// #pragma once

// #include "trgl6_icos.hpp"
// #include "vec3.hpp"
// #include <unordered_map>
// #include <vector>

// // Deduplicate a point list
// void deduplicate_points(const std::vector<Vec3>& all_points,
//                         std::vector<Vec3>& unique_points,
//                         std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal>& point_to_id);

// // Deduplicate and remap a mesh
// void deduplicate_points(Mesh& mesh);


#pragma once

#include "trgl6_icos.hpp"
#include "vec3.hpp"
#include <unordered_map>
#include <vector>

// Deduplicate a list of 3D points using a hash-based tolerance map.
// Produces a unique point list and a mapping from original point to new index.
void deduplicate_points(const std::vector<Vec3>& all_points,
                        std::vector<Vec3>& unique_points,
                        std::unordered_map<Vec3, int, Vec3Hash, Vec3Equal>& point_to_id);

// Deduplicate the point list in a Mesh object and remap its triangle connectivity.
void deduplicate_points(Mesh& mesh);
