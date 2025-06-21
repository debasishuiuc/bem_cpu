// trgl6_icos.cpp

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

using Clock = std::chrono::high_resolution_clock;

// --- Step 1: Initialize icosahedron ---
void initialize_icosahedron(PointList& nodeCoords,
                            std::vector<std::array<int, 3>>& triangleFaces,
                            bool useCCWOrientation) {
  nodeCoords = {
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

  for (auto& p : nodeCoords) p = p.normalized();

  if (useCCWOrientation) {
    triangleFaces = {
      {0,1,4},  {0,4,9},  {9,4,5},  {4,8,5},  {4,1,8},
      {8,1,10}, {8,10,3}, {5,8,3},  {5,3,2},  {2,3,7},
      {7,3,10}, {7,10,6}, {7,6,11}, {11,6,0}, {0,6,1},
      {6,10,1}, {9,11,0}, {9,2,11}, {9,5,2},  {7,11,2}
    };
  } else {
    triangleFaces = {
      {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1}, {8,10,1}, {8,3,10},
      {5,3,8}, {5,2,3}, {2,7,3}, {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6},
      {0,1,6}, {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
    };
  }
}

// --- Step 2: Refine triangles ---
void refine_triangles(PointList& nodeCoords,
                      std::vector<std::array<int, 3>>& triangleFaces,
                      int subdivisionLevel) {
  for (int level = 0; level < subdivisionLevel; ++level) {
    std::unordered_map<std::pair<int, int>, int, pair_hash> edgeToMidpointIndex;
    std::vector<std::pair<int, int>> uniqueEdges;
    std::unordered_map<std::pair<int, int>, int, pair_hash> edgeMap;

    for (const auto& tri : triangleFaces) {
      for (int i = 0; i < 3; ++i) {
        int a = tri[i], b = tri[(i + 1) % 3];
        if (a > b) std::swap(a, b);
        auto edge = std::make_pair(a, b);
        if (edgeMap.find(edge) == edgeMap.end()) {
          edgeMap[edge] = uniqueEdges.size();
          uniqueEdges.push_back(edge);
        }
      }
    }

    auto t0 = Clock::now();

    int numEdges = uniqueEdges.size();
    std::vector<Vec3> midpoints(numEdges);
#pragma omp parallel for
    for (int i = 0; i < numEdges; ++i) {
      const auto& [a, b] = uniqueEdges[i];
      midpoints[i] = ((nodeCoords[a] + nodeCoords[b]) * 0.5).normalized();
    }

    auto t1 = Clock::now();

    std::vector<int> newNodeIndices(numEdges);
    int newStartIndex = nodeCoords.size();
    nodeCoords.resize(newStartIndex + numEdges);
#pragma omp parallel for
    for (int i = 0; i < numEdges; ++i) {
      nodeCoords[newStartIndex + i] = midpoints[i];
      newNodeIndices[i] = newStartIndex + i;
    }

    for (const auto& [edge, idx] : edgeMap)
      edgeToMidpointIndex[edge] = newNodeIndices[idx];

    auto t2 = Clock::now();

    std::vector<std::array<int, 3>> refinedFaces(4 * triangleFaces.size());
#pragma omp parallel for
    for (int i = 0; i < (int)triangleFaces.size(); ++i) {
      const auto& tri = triangleFaces[i];
      int a = tri[0], b = tri[1], c = tri[2];
      int ab = edgeToMidpointIndex[{std::min(a, b), std::max(a, b)}];
      int bc = edgeToMidpointIndex[{std::min(b, c), std::max(b, c)}];
      int ca = edgeToMidpointIndex[{std::min(c, a), std::max(c, a)}];

      refinedFaces[4 * i + 0] = {a, ab, ca};
      refinedFaces[4 * i + 1] = {ab, b, bc};
      refinedFaces[4 * i + 2] = {ca, bc, c};
      refinedFaces[4 * i + 3] = {ab, bc, ca};
    }

    auto t3 = Clock::now();
    triangleFaces = std::move(refinedFaces);

    std::cout << "Level " << level << ": Midpoint compute = "
              << std::chrono::duration<double>(t1 - t0).count() << " s, "
              << "Insertion = "
              << std::chrono::duration<double>(t2 - t1).count() << " s, "
              << "Triangle update = "
              << std::chrono::duration<double>(t3 - t2).count() << " s\n";
  }
}

// --- Step 3: Construct 6-node triangular elements ---
void construct_elements(Mesh& mesh,
                        PointList& nodeCoords,
                        const std::vector<std::array<int, 3>>& triangleFaces) {
  std::unordered_map<std::pair<int, int>, int, pair_hash> edgeMap;
  std::vector<std::pair<int, int>> edges;

  for (const auto& tri : triangleFaces) {
    for (int i = 0; i < 3; ++i) {
      int a = tri[i], b = tri[(i + 1) % 3];
      if (a > b) std::swap(a, b);
      auto edge = std::make_pair(a, b);
      if (edgeMap.find(edge) == edgeMap.end()) {
        edgeMap[edge] = edges.size();
        edges.push_back(edge);
      }
    }
  }

  int numEdges = edges.size();
  std::vector<Vec3> midpoints(numEdges);
#pragma omp parallel for
  for (int i = 0; i < numEdges; ++i) {
    const auto& [a, b] = edges[i];
    midpoints[i] = ((nodeCoords[a] + nodeCoords[b]) * 0.5).normalized();
  }

  int newStartIndex = nodeCoords.size();
  nodeCoords.resize(newStartIndex + numEdges);
#pragma omp parallel for
  for (int i = 0; i < numEdges; ++i) {
    nodeCoords[newStartIndex + i] = midpoints[i];
  }

  std::unordered_map<std::pair<int, int>, int, pair_hash> edgeToMidpointIndex;
  for (const auto& [edge, idx] : edgeMap)
    edgeToMidpointIndex[edge] = newStartIndex + idx;

  mesh.nodeCoords = nodeCoords;
  mesh.elementNodes.clear();

  for (const auto& tri : triangleFaces) {
    int a = tri[0], b = tri[1], c = tri[2];
    int ab = edgeToMidpointIndex[{std::min(a, b), std::max(a, b)}];
    int bc = edgeToMidpointIndex[{std::min(b, c), std::max(b, c)}];
    int ca = edgeToMidpointIndex[{std::min(c, a), std::max(c, a)}];
    mesh.elementNodes.push_back({a, b, c, ab, bc, ca});
  }

  mesh.numNodes = mesh.nodeCoords.size();
  mesh.numElements = mesh.elementNodes.size();
}

// --- Top-level mesh generation entry point ---
void trgl6_icos(Mesh& mesh, int subdivisionLevel, int useCCWOrientation) {
  mesh.subdivisionLevel = subdivisionLevel;

  PointList nodeCoords;
  std::vector<std::array<int, 3>> triangleFaces;

  initialize_icosahedron(nodeCoords, triangleFaces, useCCWOrientation);
  refine_triangles(nodeCoords, triangleFaces, subdivisionLevel);
  construct_elements(mesh, nodeCoords, triangleFaces);

  auto t_dedup_start = Clock::now();
  std::cout << "Before deduplicate: " << mesh.nodeCoords.size() << " points\n";
  deduplicate_points(mesh);
  std::cout << "After deduplicate: " << mesh.nodeCoords.size() << " points\n";
  auto t_dedup_end = Clock::now();

  auto t_conn_start = Clock::now();
  compute_node_adjacency_for_mesh(mesh);
  compute_element_neighbors(mesh.elementNodes, mesh.elementNeighbors);
  auto t_conn_end = Clock::now();

  std::cout << "Deduplication time = "
            << std::chrono::duration<double>(t_dedup_end - t_dedup_start).count() << " s\n";
  std::cout << "Connectivity time = "
            << std::chrono::duration<double>(t_conn_end - t_conn_start).count() << " s\n";

  print_mesh_summary(mesh, "OpenMP");
}
