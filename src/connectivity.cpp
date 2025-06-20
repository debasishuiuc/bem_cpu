// connectivity.cpp

#include "mesh_utils.hpp"
#include <omp.h>
#include <vector>
#include <array>
#include <unordered_map>
#include <algorithm>
#include <parallel/algorithm>
#include <iostream>
#include <set>

// Data structure to hold edge info for nbe
struct EdgeInfo {
  int v1, v2;
  int tri;
  int local_edge;

  bool operator<(const EdgeInfo& other) const {
    if (v1 != other.v1) return v1 < other.v1;
    return v2 < other.v2;
  }
};

// Computes nbe: element-to-element adjacency
void compute_nbe(const std::vector<std::array<int, 6>>& elements,
                 std::vector<std::array<int, 3>>& nbe) {
  int Nelm = elements.size();
  std::vector<EdgeInfo> edges(3 * Nelm);

#pragma omp parallel for
  for (int ei = 0; ei < Nelm; ++ei) {
    const auto& elem = elements[ei];
    std::array<std::pair<int, int>, 3> edge_verts = {
      std::make_pair(std::min(elem[0], elem[1]), std::max(elem[0], elem[1])),
      std::make_pair(std::min(elem[1], elem[2]), std::max(elem[1], elem[2])),
      std::make_pair(std::min(elem[2], elem[0]), std::max(elem[2], elem[0]))
    };
    for (int j = 0; j < 3; ++j) {
      edges[3 * ei + j] = {edge_verts[j].first, edge_verts[j].second, ei, j};
    }
  }

  __gnu_parallel::sort(edges.begin(), edges.end());
  nbe.resize(Nelm, {-1, -1, -1});

  for (size_t i = 0; i + 1 < edges.size(); ++i) {
    const auto& e1 = edges[i];
    const auto& e2 = edges[i + 1];
    if (e1.v1 == e2.v1 && e1.v2 == e2.v2) {
      nbe[e1.tri][e1.local_edge] = e2.tri;
      nbe[e2.tri][e2.local_edge] = e1.tri;
      ++i;
    }
  }
}

// Computes ne(i,j): number of elements and their indices for each node
void compute_ne(
    const std::vector<std::array<int, 6>>& n,
    const std::vector<Vec3>& p,
    std::vector<std::vector<int>>& ne)
{
  int Npts = p.size();
  int Nelm = n.size();
  int nthreads = omp_get_max_threads();

  // Per-thread local buffers: [thread][node] = list of element indices
  std::vector<std::vector<std::vector<int>>> thread_ne(nthreads,
      std::vector<std::vector<int>>(Npts));

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
#pragma omp for
    for (int i = 0; i < Nelm; ++i) {
      const auto& tri = n[i];
      for (int j = 0; j < 6; ++j) {
        int node = tri[j];
        if (node < 0 || node >= Npts) {
          std::cerr << "Invalid node index " << node << " in element " << i << std::endl;
          std::abort();
        }
        thread_ne[tid][node].push_back(i);
      }
    }
  }

  // Merge all threads’ buffers into final ne
  ne.assign(Npts, std::vector<int>());
  for (int i = 0; i < Npts; ++i) {
    size_t total_size = 0;
    for (int t = 0; t < nthreads; ++t)
      total_size += thread_ne[t][i].size();

    ne[i].reserve(total_size + 1);  // +1 for count prefix
    for (int t = 0; t < nthreads; ++t)
      ne[i].insert(ne[i].end(),
                   thread_ne[t][i].begin(),
                   thread_ne[t][i].end());
  }

  // Prepend count as first entry in each row
#pragma omp parallel for
  for (int i = 0; i < Npts; ++i) {
    std::vector<int> row;
    row.reserve(ne[i].size() + 1);
    row.push_back(static_cast<int>(ne[i].size()));
    row.insert(row.end(), ne[i].begin(), ne[i].end());
    ne[i] = std::move(row);
  }
}


// Top-level wrapper for Mesh
void compute_ne(Mesh& mesh) {
  if (mesh.Npts == 0) mesh.Npts = mesh.p.size();
  compute_ne(mesh.n, mesh.p, mesh.ne);
}
