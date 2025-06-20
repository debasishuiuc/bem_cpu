// trgl6_icos.cpp

// #include "trgl6_icos.hpp"
// #include "mesh_utils.hpp"
// #include "connectivity.hpp"
// #include "deduplicate.hpp"
// #include <unordered_map>
// #include <omp.h>
// #include <vector>
// #include <cmath>
// #include <iostream>
// #include <utility>
// #include <algorithm>
// #include <chrono>

// void trgl6_icos(Mesh& mesh, int Ndiv, int CCW) {
//   using Clock = std::chrono::high_resolution_clock;

//   mesh.Ndiv = Ndiv;

//   std::vector<Vec3> points = {
//     {-0.5257311121191336, 0.0, 0.8506508083520399},
//     { 0.5257311121191336, 0.0, 0.8506508083520399},
//     {-0.5257311121191336, 0.0, -0.8506508083520399},
//     { 0.5257311121191336, 0.0, -0.8506508083520399},
//     { 0.0, 0.8506508083520399, 0.5257311121191336},
//     { 0.0, 0.8506508083520399, -0.5257311121191336},
//     { 0.0, -0.8506508083520399, 0.5257311121191336},
//     { 0.0, -0.8506508083520399, -0.5257311121191336},
//     { 0.8506508083520399, 0.5257311121191336, 0.0},
//     {-0.8506508083520399, 0.5257311121191336, 0.0},
//     { 0.8506508083520399, -0.5257311121191336, 0.0},
//     {-0.8506508083520399, -0.5257311121191336, 0.0}
//   };
//   for (auto& p : points) p = p.normalized();

//   std::vector<std::array<int,3>> triangles;

//   if (CCW) {
//     // Counter-clockwise ordering (for outward normals)
//     triangles = {
//       {0,1,4},  {0,4,9},  {9,4,5},  {4,8,5},  {4,1,8},
//       {8,1,10}, {8,10,3}, {5,8,3},  {5,3,2},  {2,3,7},
//       {7,3,10}, {7,10,6}, {7,6,11}, {11,6,0}, {0,6,1},
//       {6,10,1}, {9,11,0}, {9,2,11}, {9,5,2},  {7,11,2}
//     };
//   } else {
//     // Clockwise ordering (for inward normals)
//     triangles = {
//       {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1}, {8,10,1}, {8,3,10},
//       {5,3,8}, {5,2,3}, {2,7,3}, {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6},
//       {0,1,6}, {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
//     };
//   }

//   std::unordered_map<std::pair<int,int>, int, pair_hash> edge_to_mid;

//   for (int level = 0; level < Ndiv; ++level) {
//     edge_to_mid.clear();

//     // Step 1: Collect unique edges
//     std::vector<std::pair<int,int>> edges;
//     std::unordered_map<std::pair<int,int>, int, pair_hash> edge_map;

//     for (const auto& tri : triangles) {
//       for (int i = 0; i < 3; ++i) {
//         int a = tri[i];
//         int b = tri[(i + 1) % 3];
//         if (a > b) std::swap(a, b);
//         auto e = std::make_pair(a, b);
//         if (edge_map.find(e) == edge_map.end()) {
//           edge_map[e] = edges.size();
//           edges.push_back(e);
//         }
//       }
//     }

//     auto t0 = Clock::now();

//     // Step 2: Compute midpoints
//     int n_edges = edges.size();
//     std::vector<Vec3> midpoints(n_edges);
// #pragma omp parallel for
//     for (int i = 0; i < n_edges; ++i) {
//       int a = edges[i].first;
//       int b = edges[i].second;
//       midpoints[i] = ((points[a] + points[b]) * 0.5).normalized();
//     }

//     auto t1 = Clock::now();

//     // Step 3: Insert midpoints
//     std::vector<int> edge_to_point(n_edges);
//     int start_idx = points.size();
//     points.resize(start_idx + n_edges);
// #pragma omp parallel for
//     for (int i = 0; i < n_edges; ++i) {
//       points[start_idx + i] = midpoints[i];
//       edge_to_point[i] = start_idx + i;
//     }

//     for (const auto& [key, idx] : edge_map)
//       edge_to_mid[key] = edge_to_point[idx];

//     auto t2 = Clock::now();

//     // Step 4: Triangle update
//     int n_old = triangles.size();
//     std::vector<std::array<int,3>> newtris(4 * n_old);
// #pragma omp parallel for
//     for (int i = 0; i < n_old; ++i) {
//       const auto& tri = triangles[i];
//       int a = tri[0], b = tri[1], c = tri[2];
//       int ab = edge_to_mid[{std::min(a,b), std::max(a,b)}];
//       int bc = edge_to_mid[{std::min(b,c), std::max(b,c)}];
//       int ca = edge_to_mid[{std::min(c,a), std::max(c,a)}];

//       newtris[4*i + 0] = {a, ab, ca};
//       newtris[4*i + 1] = {ab, b, bc};
//       newtris[4*i + 2] = {ca, bc, c};
//       newtris[4*i + 3] = {ab, bc, ca};
//     }

//     auto t3 = Clock::now();
//     triangles = std::move(newtris);

//     std::cout << "Level " << level << ": Midpoint compute = "
//               << std::chrono::duration<double>(t1 - t0).count() << " s, "
//               << "Insertion = "
//               << std::chrono::duration<double>(t2 - t1).count() << " s, "
//               << "Triangle update = "
//               << std::chrono::duration<double>(t3 - t2).count() << " s\n";
//   }

//   // Final edge-to-midpoint mapping
//   edge_to_mid.clear();
//   std::unordered_map<std::pair<int,int>, int, pair_hash> final_edge_map;
//   std::vector<std::pair<int,int>> final_edges;

//   for (const auto& tri : triangles) {
//     for (int i = 0; i < 3; ++i) {
//       int a = tri[i];
//       int b = tri[(i + 1) % 3];
//       if (a > b) std::swap(a, b);
//       auto e = std::make_pair(a, b);
//       if (final_edge_map.find(e) == final_edge_map.end()) {
//         final_edge_map[e] = final_edges.size();
//         final_edges.push_back(e);
//       }
//     }
//   }

//   int n_final_edges = final_edges.size();
//   std::vector<Vec3> final_midpoints(n_final_edges);
// #pragma omp parallel for
//   for (int i = 0; i < n_final_edges; ++i) {
//     int a = final_edges[i].first;
//     int b = final_edges[i].second;
//     final_midpoints[i] = ((points[a] + points[b]) * 0.5).normalized();
//   }

//   int final_start_idx = points.size();
//   points.resize(final_start_idx + n_final_edges);
// #pragma omp parallel for
//   for (int i = 0; i < n_final_edges; ++i) {
//     points[final_start_idx + i] = final_midpoints[i];
//   }

//   for (const auto& [key, idx] : final_edge_map)
//     edge_to_mid[key] = final_start_idx + idx;

//   // Construct 6-node triangle elements
//   mesh.p = points;
//   for (const auto& tri : triangles) {
//     int a = tri[0], b = tri[1], c = tri[2];
//     int ab = edge_to_mid[{std::min(a,b), std::max(a,b)}];
//     int bc = edge_to_mid[{std::min(b,c), std::max(b,c)}];
//     int ca = edge_to_mid[{std::min(c,a), std::max(c,a)}];

//     mesh.n.push_back({a, b, c, ab, bc, ca});
//   }

//   mesh.Npts = mesh.p.size();
//   mesh.Nelm = mesh.n.size();

//   auto t_dedup_start = Clock::now();
//   std::cout << "Before deduplicate: " << mesh.p.size() << " points\n";
//   deduplicate_points(mesh);
//   std::cout << "After deduplicate: " << mesh.p.size() << " points\n";
//   auto t_dedup_end = Clock::now();

//   auto t_conn_start = Clock::now();
//   compute_ne(mesh);
//   compute_nbe(mesh.n, mesh.nbe);
//   auto t_conn_end = Clock::now();

//   std::cout << "Deduplication time = "
//             << std::chrono::duration<double>(t_dedup_end - t_dedup_start).count() << " s\n";
//   std::cout << "Connectivity time = "
//             << std::chrono::duration<double>(t_conn_end - t_conn_start).count() << " s\n";

//   print_mesh_summary(mesh, "OpenMP");
// }


#include "trgl6_icos.hpp"
#include "mesh_utils.hpp"
#include "connectivity.hpp"
#include "deduplicate.hpp"
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <omp.h>

// --- Hash for edges (already in mesh_utils.hpp) ---
// struct pair_hash { ... }  // already defined in mesh_utils.hpp

using Clock = std::chrono::high_resolution_clock;

// --- Step 1: Initialize icosahedron ---
void initialize_icosahedron(std::vector<Vec3>& points,
                            std::vector<std::array<int, 3>>& triangles,
                            bool ccw) {
  points = {
    {-0.5257311121191336, 0.0, 0.8506508083520399},
    { 0.5257311121191336, 0.0, 0.8506508083520399},
    {-0.5257311121191336, 0.0, -0.8506508083520399},
    { 0.5257311121191336, 0.0, -0.8506508083520399},
    { 0.0, 0.8506508083520399, 0.5257311121191336},
    { 0.0, 0.8506508083520399, -0.5257311121191336},
    { 0.0, -0.8506508083520399, 0.5257311121191336},
    { 0.0, -0.8506508083520399, -0.5257311121191336},
    { 0.8506508083520399, 0.5257311121191336, 0.0},
    {-0.8506508083520399, 0.5257311121191336, 0.0},
    { 0.8506508083520399, -0.5257311121191336, 0.0},
    {-0.8506508083520399, -0.5257311121191336, 0.0}
  };

  for (auto& p : points) p = p.normalized();

  if (ccw) {
    triangles = {
      {0,1,4},  {0,4,9},  {9,4,5},  {4,8,5},  {4,1,8},
      {8,1,10}, {8,10,3}, {5,8,3},  {5,3,2},  {2,3,7},
      {7,3,10}, {7,10,6}, {7,6,11}, {11,6,0}, {0,6,1},
      {6,10,1}, {9,11,0}, {9,2,11}, {9,5,2},  {7,11,2}
    };
  } else {
    triangles = {
      {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1}, {8,10,1}, {8,3,10},
      {5,3,8}, {5,2,3}, {2,7,3}, {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6},
      {0,1,6}, {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
    };
  }
}

// --- Step 2: Refine triangles ---
void refine_triangles(std::vector<Vec3>& points,
                      std::vector<std::array<int, 3>>& triangles,
                      int Ndiv) {
  for (int level = 0; level < Ndiv; ++level) {
    std::unordered_map<std::pair<int, int>, int, pair_hash> edge_to_mid;

    std::vector<std::pair<int, int>> edges;
    std::unordered_map<std::pair<int, int>, int, pair_hash> edge_map;

    for (const auto& tri : triangles) {
      for (int i = 0; i < 3; ++i) {
        int a = tri[i], b = tri[(i + 1) % 3];
        if (a > b) std::swap(a, b);
        auto e = std::make_pair(a, b);
        if (edge_map.find(e) == edge_map.end()) {
          edge_map[e] = edges.size();
          edges.push_back(e);
        }
      }
    }

    auto t0 = Clock::now();

    int n_edges = edges.size();
    std::vector<Vec3> midpoints(n_edges);
#pragma omp parallel for
    for (int i = 0; i < n_edges; ++i)
      midpoints[i] = ((points[edges[i].first] + points[edges[i].second]) * 0.5).normalized();

    auto t1 = Clock::now();

    std::vector<int> edge_to_point(n_edges);
    int start_idx = points.size();
    points.resize(start_idx + n_edges);
#pragma omp parallel for
    for (int i = 0; i < n_edges; ++i) {
      points[start_idx + i] = midpoints[i];
      edge_to_point[i] = start_idx + i;
    }

    for (const auto& [key, idx] : edge_map)
      edge_to_mid[key] = edge_to_point[idx];

    auto t2 = Clock::now();

    std::vector<std::array<int, 3>> newtris(4 * triangles.size());
#pragma omp parallel for
    for (int i = 0; i < (int)triangles.size(); ++i) {
      const auto& tri = triangles[i];
      int a = tri[0], b = tri[1], c = tri[2];
      int ab = edge_to_mid[{std::min(a, b), std::max(a, b)}];
      int bc = edge_to_mid[{std::min(b, c), std::max(b, c)}];
      int ca = edge_to_mid[{std::min(c, a), std::max(c, a)}];

      newtris[4 * i + 0] = {a, ab, ca};
      newtris[4 * i + 1] = {ab, b, bc};
      newtris[4 * i + 2] = {ca, bc, c};
      newtris[4 * i + 3] = {ab, bc, ca};
    }

    auto t3 = Clock::now();
    triangles = std::move(newtris);

    std::cout << "Level " << level << ": Midpoint compute = "
              << std::chrono::duration<double>(t1 - t0).count() << " s, "
              << "Insertion = "
              << std::chrono::duration<double>(t2 - t1).count() << " s, "
              << "Triangle update = "
              << std::chrono::duration<double>(t3 - t2).count() << " s\n";
  }
}

// --- Step 3: Final midpoint insertion and 6-node construction ---
void construct_elements(Mesh& mesh,
                        std::vector<Vec3>& points,
                        const std::vector<std::array<int, 3>>& triangles) {
  std::unordered_map<std::pair<int,int>, int, pair_hash> final_edge_map;
  std::vector<std::pair<int,int>> final_edges;

  for (const auto& tri : triangles) {
    for (int i = 0; i < 3; ++i) {
      int a = tri[i], b = tri[(i + 1) % 3];
      if (a > b) std::swap(a, b);
      auto e = std::make_pair(a, b);
      if (final_edge_map.find(e) == final_edge_map.end()) {
        final_edge_map[e] = final_edges.size();
        final_edges.push_back(e);
      }
    }
  }

  int n_final_edges = final_edges.size();
  std::vector<Vec3> final_midpoints(n_final_edges);
#pragma omp parallel for
  for (int i = 0; i < n_final_edges; ++i) {
    int a = final_edges[i].first;
    int b = final_edges[i].second;
    final_midpoints[i] = ((points[a] + points[b]) * 0.5).normalized();
  }

  int start_idx = points.size();
  points.resize(start_idx + n_final_edges);
#pragma omp parallel for
  for (int i = 0; i < n_final_edges; ++i) {
    points[start_idx + i] = final_midpoints[i];
  }

  std::unordered_map<std::pair<int, int>, int, pair_hash> edge_to_mid;
  for (const auto& [key, idx] : final_edge_map)
    edge_to_mid[key] = start_idx + idx;

  mesh.p = points;
  mesh.n.clear();

  for (const auto& tri : triangles) {
    int a = tri[0], b = tri[1], c = tri[2];
    int ab = edge_to_mid[{std::min(a, b), std::max(a, b)}];
    int bc = edge_to_mid[{std::min(b, c), std::max(b, c)}];
    int ca = edge_to_mid[{std::min(c, a), std::max(c, a)}];
    mesh.n.push_back({a, b, c, ab, bc, ca});
  }

  mesh.Npts = mesh.p.size();
  mesh.Nelm = mesh.n.size();
}

// --- Top-level mesh generation entry point ---
void trgl6_icos(Mesh& mesh, int Ndiv, int CCW) {
  mesh.Ndiv = Ndiv;

  std::vector<Vec3> points;
  std::vector<std::array<int, 3>> triangles;

  initialize_icosahedron(points, triangles, CCW);
  refine_triangles(points, triangles, Ndiv);
  construct_elements(mesh, points, triangles);

  auto t_dedup_start = Clock::now();
  std::cout << "Before deduplicate: " << mesh.p.size() << " points\n";
  deduplicate_points(mesh);
  std::cout << "After deduplicate: " << mesh.p.size() << " points\n";
  auto t_dedup_end = Clock::now();

  auto t_conn_start = Clock::now();
  compute_ne(mesh);
  compute_nbe(mesh.n, mesh.nbe);
  auto t_conn_end = Clock::now();

  std::cout << "Deduplication time = "
            << std::chrono::duration<double>(t_dedup_end - t_dedup_start).count() << " s\n";
  std::cout << "Connectivity time = "
            << std::chrono::duration<double>(t_conn_end - t_conn_start).count() << " s\n";

  print_mesh_summary(mesh, "OpenMP");
}
