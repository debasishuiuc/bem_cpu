// trgl6_icos.hpp

#pragma once
#include <vector>
#include <array>
#include <cmath>
#include "vec3.hpp"

using real = double;

// Mesh data structures
using PointList = std::vector<Vec3>;
using Element6 = std::array<int, 6>;  // 3 vertices + 3 mid-edge nodes
using ElementList = std::vector<Element6>;
using NodeToElementConnectivity = std::vector<std::vector<int>>; // node → elements
using ElementToElementConnectivity = std::vector<std::array<int, 3>>; // element → neighbors

struct Mesh {
  int subdivisionLevel = 0;       // Level of icosahedron refinement
  int numNodes = 0;               // Total number of nodes
  int numElements = 0;            // Total number of surface elements

  PointList nodeCoords;           // Global coordinates of nodes
  ElementList elementNodes;       // Connectivity: global node indices per 6-node triangle

  NodeToElementConnectivity nodeToElements;  // ne(i,j): elements connected to node i
  ElementToElementConnectivity elementNeighbors; // nbe(i,j): neighbors of element i across edge j

  // Geometry (computed post-mesh generation)
  std::vector<Vec3> elementNormals;  // vnc: normals at element centroids
  std::vector<double> elementCurvature;  // crvmel: mean curvature per element
  std::vector<Vec3> nodeNormals;     // vna: normals at each node
};

// Function to generate triangulated spherical mesh
void trgl6_icos(Mesh& mesh, int subdivisionLevel, int useCCWOrientation = 1);
