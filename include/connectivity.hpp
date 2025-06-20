// connectivity.hpp

#pragma once

#include "trgl6_icos.hpp"
#include <vector>
#include <array>

void compute_nbe(const std::vector<std::array<int, 6>>&, std::vector<std::array<int, 3>>&);

// Accepts Vec3 points
void compute_ne(const std::vector<std::array<int, 6>>& n,
                const std::vector<Vec3>& p,
                std::vector<std::vector<int>>& ne);

void compute_ne(Mesh& mesh);
