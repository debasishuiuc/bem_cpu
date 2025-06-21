// connectivity.hpp

#pragma once

#include "trgl6_icos.hpp"
#include <vector>
#include <array>

// Compute element-to-element connectivity (nbe)
// Given the 6-node triangle list 'n', fills 'nbe' with neighboring element indices.
void compute_element_neighbors(const std::vector<std::array<int, 6>>& n,
                 std::vector<std::array<int, 3>>& nbe);

// Compute node-to-element connectivity (ne) using explicit nodes and points
void compute_node_adjacency(const std::vector<std::array<int, 6>>& n,
                const std::vector<Vec3>& p,
                std::vector<std::vector<int>>& ne);

// Overload that works directly with the Mesh object (fills mesh.ne)
void compute_node_adjacency_for_mesh(Mesh& mesh);
