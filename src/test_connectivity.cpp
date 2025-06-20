// test_connectivity.cpp

#include <gtest/gtest.h>
#include "mesh_utils.hpp"
#include "connectivity_tests.hpp"
#include "connectivity.hpp"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST(MeshConnectivity, NBE_Symmetric) {
  Mesh mesh;
  trgl6_icos(mesh, 1, 1);
  compute_nbe(mesh.n, mesh.nbe);
  EXPECT_TRUE(check_neighbor_symmetry(mesh.nbe));
}

TEST(MeshConnectivity, NBE_NoBoundary) {
  Mesh mesh;
  trgl6_icos(mesh, 1, 1);
  compute_nbe(mesh.n, mesh.nbe);
  EXPECT_TRUE(check_no_boundary_edges(mesh.nbe));
}


TEST(MeshConnectivity, NodeValenceReasonable) {
  Mesh mesh;
  trgl6_icos(mesh, 1, 1);     // generate mesh
  mesh.Npts = mesh.p.size(); // <-- add this line
  compute_ne(mesh);          // now ne_fallback will be populated
  EXPECT_TRUE(check_node_valence(mesh.ne_fallback));
}
