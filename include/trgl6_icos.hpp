// trgl6_icos.hpp

#pragma once
#include <vector>
#include <array>
#include <cmath>
#include "vec3.hpp"

typedef double real;

// Mesh data structures
typedef std::vector<Vec3> PointList;
typedef std::array<int, 6> Element6;
typedef std::vector<Element6> ElementList;
typedef std::vector<std::vector<int>> Connectivity;

struct Mesh {
  int Ndiv;         // Number of subdivisions
  int Npts;         // Number of points
  int Nelm;         // Number of elements
  PointList p;      // Coordinates of nodes
  ElementList n;    // 6-node triangle connectivity

  
  std::vector<std::vector<int>> ne; // Generate connectivity table: ne(i,j) for elements touching global node i
  std::vector<std::array<int, 3>> nbe; // Element-to-neighbor connectivity

  // Geometry data (computed later)
  std::vector<std::array<double, 3>> vnc;  // Element normals at centroids
  std::vector<double> crvmel;              // Mean curvature per element
  std::vector<std::array<double, 3>> vna;  // Node normals
};

// Function to generate mesh (defined in trgl6_icos.cpp)
void trgl6_icos(Mesh& mesh, int Ndiv, int CCW = 1);
