// main.cpp

#include "trgl6_icos.hpp"   
#include "geometry_analyzer.hpp"      // CPU-only geometry computations
#include "connectivity.hpp"  // compute_ne, compute_nbe, collect_edges
#include "deduplicate.hpp"   // deduplicate_points
#include "mesh_utils.hpp"    // Vec3 helpers, mesh integrity checks, etc.
#include "write_txt.hpp"
#include "write_vtk.hpp"
#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <omp.h>

using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

int main(int argc, char* argv[]) {
  int Ndiv = -1;
  bool writeFiles = false;
  int num_threads = 8;
  std::string log_path = "";

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--ndiv" && i + 1 < argc) {
      Ndiv = std::stoi(argv[++i]);
    } else if (arg == "--threads" && i + 1 < argc) {
      num_threads = std::stoi(argv[++i]);
    } else if (arg == "--write") {
      writeFiles = true;
    } else if (arg == "--log" && i + 1 < argc) {
      log_path = argv[++i];
    } else {
      std::cerr << "Unknown or incomplete flag: " << arg << std::endl;
      std::cerr << "Usage: ./meshgen --ndiv N [--write] [--threads N] [--log path]" << std::endl;
      return 1;
    }
  }

  if (Ndiv < 0) {
    std::cerr << "Error: Subdivision level (--ndiv) must be specified.\n";
    return 1;
  }

  omp_set_num_threads(num_threads);
  std::cout << "Running mesh generation with Ndiv = " << Ndiv
            << ", using " << num_threads << " threads...\n";

  Mesh mesh;

  auto t_mesh_start = Clock::now();
  trgl6_icos(mesh, Ndiv, 1);
  auto t_mesh_end = Clock::now();
  Duration mesh_time = t_mesh_end - t_mesh_start;

  auto t_geom_start = Clock::now();
  GeometryAnalyzer analyzer;
  analyzer.compute(mesh, 6);
  
  auto t_geom_end = Clock::now();
  Duration geom_time = t_geom_end - t_geom_start;

  // Print decorated summary
  std::cout << "\n=== Mesh Statistics ===\n";
  std::cout << "Subdivision level (Ndiv): " << Ndiv << "\n";
  std::cout << "Points: " << mesh.numNodes << ", Elements: " << mesh.numElements << "\n";

  std::cout << "\n=== Timing ===\n";
  std::cout << "Mesh generation took " << mesh_time.count() << " seconds.\n";
  std::cout << "Compute Mesh took " << geom_time.count() << " seconds.\n";

  // Print geometry summary
  std::cout << "\n=== Geometry Summary ===\n";
  analyzer.print_summary(mesh, 5); // limit print to first 5 nodes/elements

  analyzer.check_element_quality(mesh, 6);
  
  // Mesh integrity check
  std::cout << "\n=== Mesh Integrity Check ===\n";
  if (!check_mesh_integrity(mesh)) {
    std::cerr << "Mesh integrity check FAILED.\n";
    return 1;
  } else {
    std::cout << "Mesh integrity check passed.\n";
  }

  // Triangle orientation check
  std::cout << "\n=== Triangle Orientation Check ===\n";
  if (!check_triangle_orientation_strict(mesh)) {
    std::cerr << "Orientation check FAILED: Some triangles are not CCW.\n";
    return 1;
  } else {
    std::cout << "All triangles are properly counter-clockwise (CCW).\n";
  }
  
  if (!log_path.empty()) {
    std::filesystem::path log_dir = std::filesystem::path(log_path).parent_path();
    if (!log_dir.empty() && !std::filesystem::exists(log_dir)) {
      std::filesystem::create_directories(log_dir);
    }

    std::ofstream logfile(log_path, std::ios::app);
    if (logfile.is_open()) {
      logfile << "Ndiv=" << Ndiv;
      logfile << ", Mode=OpenMP"
              << ", Threads=" << num_threads;
      logfile << ", MeshTime=" << std::scientific << mesh_time.count()
              << ", GeometryTime=" << geom_time.count()
              << ", Npts=" << mesh.numNodes
              << ", Nelm=" << mesh.numElements
              << ", SurfaceArea=" << std::fixed << analyzer.totalArea
              << ", Volume=" << analyzer.totalVolume
              << ", Centroid=[" << std::scientific
              << analyzer.surfaceCentroid[0] << " "
              << analyzer.surfaceCentroid[1] << " "
              << analyzer.surfaceCentroid[2] << "]"
              << "\n";

      logfile.close();
      std::cout << "📝 Log written to: " << log_path << "\n";
    } else {
      std::cerr << "Warning: could not open " << log_path << " for writing.\n";
    }
  }

  std::cout << "Ndiv = " << Ndiv << ", Threads = " << num_threads
            << ", Mesh Time = " << std::scientific << mesh_time.count()
            << " s, Geometry Time = " << geom_time.count() << " s\n";

  if (writeFiles) {
    std::string outdir = "../output/";
    if (!std::filesystem::exists(outdir)) {
      if (std::filesystem::create_directory(outdir)) {
        std::cout << "Created output directory at " << outdir << "\n";
      } else {
        std::cerr << "Failed to create output directory: " << outdir << "\n";
        return 1;
      }
    }

    compute_node_adjacency_for_mesh(mesh);
    write_txt(mesh, analyzer, outdir);

    std::cout << "Writing VTK to: " << outdir + "mesh.vtk" << std::endl;
    std::cout << "Number of elements = " << mesh.numElements
              << ", Number of vnc = " << analyzer.elementNormals.size()
              << ", Number of crvmel = " << analyzer.elementCurvature.size() << std::endl;

    write_vtk(mesh, analyzer, outdir + "mesh.vtk");
    
    std::cout << "Mesh written to output/ directory.\n";
  }

  return 0;
}
