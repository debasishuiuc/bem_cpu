// write_vtk.cpp

#include "write_vtk.hpp"
#include <fstream>
#include <iostream>

void write_vtk(const Mesh& mesh, const GeometryAnalyzer& analyzer, const std::string& filename) {
  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Could not open file " << filename << " for writing.\n";
    return;
  }

  file << "# vtk DataFile Version 3.0\nSphere Mesh with Geometry\nASCII\n";
  file << "DATASET POLYDATA\n";

  // Write points
  file << "POINTS " << mesh.numNodes << " float\n";
  for (const auto& pt : mesh.nodeCoords)
    file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";

  // Write triangle connectivity (only the first 3 nodes for VTK)
  file << "POLYGONS " << mesh.numElements << " " << mesh.numElements * 4 << "\n";
  for (const auto& tri : mesh.elementNodes)
    file << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";

  // Cell data: normals and curvature
  file << "CELL_DATA " << mesh.numElements << "\n";

  file << "VECTORS elem_normals float\n";
  for (const auto& v : analyzer.elementNormals)
    file << v[0] << " " << v[1] << " " << v[2] << "\n";

  file << "SCALARS curvature float 1\nLOOKUP_TABLE default\n";
  for (const auto& k : analyzer.elementCurvature)
    file << k << "\n";

  file << "POINT_DATA " << mesh.numNodes << "\n";

  file << "SCALARS node_curvature float 1\nLOOKUP_TABLE default\n";
  for (const auto& val : analyzer.nodeCurvature)
    file << val << "\n";

  file << "VECTORS node_normals float\n";
  for (const auto& n : analyzer.nodeNormals)
    file << n[0] << " " << n[1] << " " << n[2] << "\n";

  std::cout << "VTK file written successfully.\n";
}
