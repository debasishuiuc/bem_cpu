// write_txt.cpp

#include "write_txt.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iomanip>

void write_txt(const Mesh& mesh, const GeometryAnalyzer& analyzer, const std::string& outdir) {
  std::ofstream pfile(outdir + "node_coords.txt");
  std::ofstream nfile(outdir + "element_nodes.txt");
  std::ofstream nefile(outdir + "node_to_elements.txt");
  std::ofstream nbefile(outdir + "element_neighbors.txt");
  std::ofstream vncfile(outdir + "element_normals.txt");
  std::ofstream crvfile(outdir + "element_curvature.txt");
  std::ofstream vnafile(outdir + "node_normals.txt");
  std::ofstream nodecurvfile(outdir + "node_curvature.txt");

  if (!pfile || !nfile || !nefile || !nbefile || !vncfile || !crvfile || !vnafile || !nodecurvfile) {
    std::cerr << "Error opening one or more output .txt files in " << outdir << "\n";
    return;
  }

  // Node coordinates
  for (const auto& pt : mesh.nodeCoords)
    pfile << pt[0] << " " << pt[1] << " " << pt[2] << "\n";

  // 6-node triangle connectivity
  for (const auto& tri : mesh.elementNodes)
    nfile << tri[0] << " " << tri[1] << " " << tri[2] << " "
          << tri[3] << " " << tri[4] << " " << tri[5] << "\n";

  // Connectivity: node to elements (ne)
  for (const auto& row : mesh.nodeToElements) {
    // Print the count of connected elements (row[0]) first
    nefile << row[0] << " ";

    // Print the element indices (starting from row[1] onward)
    for (size_t j = 1; j <= row[0]; ++j)  // Start from 1 to skip the count
      nefile << row[j] << " ";

    nefile << "\n";
  }


  // Connectivity: element neighbors (nbe)
  for (const auto& nbs : mesh.elementNeighbors)
    nbefile << nbs[0] << " " << nbs[1] << " " << nbs[2] << "\n";

  // Element normals (vnc)
  for (const auto& n : analyzer.elementNormals)
    vncfile << n[0] << " " << n[1] << " " << n[2] << "\n";

  // Element curvature (crvmel)
  for (double c : analyzer.elementCurvature)
    crvfile << std::fixed << std::setprecision(12) << c << "\n";

  // Node normals (vna)
  for (const auto& n : analyzer.nodeNormals)
    vnafile << n[0] << " " << n[1] << " " << n[2] << "\n";

  // Node curvature (node_curvature)
  for (double val : analyzer.nodeCurvature)
    nodecurvfile << std::fixed << std::setprecision(12) << val << "\n";

  std::cout << "TXT files written to: " << outdir << "\n";
}
