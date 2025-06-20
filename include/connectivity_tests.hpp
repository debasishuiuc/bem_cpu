// tests.hpp

#pragma once
#include "mesh_utils.hpp"

bool check_neighbor_symmetry(const std::vector<std::array<int,3>>& nbe);
bool check_no_boundary_edges(const std::vector<std::array<int,3>>& nbe);
bool check_node_valence(const std::vector<std::vector<int>>& ne_fallback, int min_expected = 4, int max_expected = 8);
void run_connectivity_tests(const Mesh& mesh);
