// connectivity.cpp

#include "mesh_utils.hpp"
#include "trgl6_icos.hpp"  // for Mesh definition
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

// Computes element-to-element adjacency (nbe)
void compute_element_neighbors(const ElementList& elementNodes,
			       ElementToElementConnectivity& elementNeighbors) {
  int numElements = elementNodes.size();
  std::vector<EdgeInfo> edges(3 * numElements);

#pragma omp parallel for
  for (int ei = 0; ei < numElements; ++ei) {
    const auto& elem = elementNodes[ei];
    std::array<std::pair<int, int>, 3> edgeVerts = {
      std::make_pair(std::min(elem[0], elem[1]), std::max(elem[0], elem[1])),
      std::make_pair(std::min(elem[1], elem[2]), std::max(elem[1], elem[2])),
      std::make_pair(std::min(elem[2], elem[0]), std::max(elem[2], elem[0]))
    };
    for (int j = 0; j < 3; ++j) {
      edges[3 * ei + j] = {edgeVerts[j].first, edgeVerts[j].second, ei, j};
    }
  }

  __gnu_parallel::sort(edges.begin(), edges.end());
  elementNeighbors.resize(numElements, {-1, -1, -1});

  for (size_t i = 0; i + 1 < edges.size(); ++i) {
    const auto& e1 = edges[i];
    const auto& e2 = edges[i + 1];
    if (e1.v1 == e2.v1 && e1.v2 == e2.v2) {
      elementNeighbors[e1.tri][e1.local_edge] = e2.tri;
      elementNeighbors[e2.tri][e2.local_edge] = e1.tri;
      ++i;
    }
  }
}

// Computes node-to-element adjacency (ne)
void compute_node_adjacency(const ElementList& elementNodes,
			    const PointList& nodeCoords,
			    NodeToElementConnectivity& nodeToElements) {
  int numNodes = nodeCoords.size();
  int numElements = elementNodes.size();
  int nthreads = omp_get_max_threads();

  // Per-thread buffers: thread_ne[thread_id][node] = list of element indices
  std::vector<std::vector<std::vector<int>>> thread_ne(nthreads,
						       std::vector<std::vector<int>>(numNodes));

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
#pragma omp for
    for (int i = 0; i < numElements; ++i) {
      const auto& elem = elementNodes[i];
      for (int j = 0; j < 6; ++j) {
        int node = elem[j];
        if (node < 0 || node >= numNodes) {
          std::cerr << "Invalid node index " << node << " in element " << i << std::endl;
          std::abort();
        }
        thread_ne[tid][node].push_back(i);
      }
    }
  }

  // Merge all thread-local buffers into final nodeToElements
  nodeToElements.assign(numNodes, std::vector<int>());
  for (int i = 0; i < numNodes; ++i) {
    size_t total_size = 0;
    for (int t = 0; t < nthreads; ++t)
      total_size += thread_ne[t][i].size();

    nodeToElements[i].reserve(total_size + 1);  // +1 for prefix count
    for (int t = 0; t < nthreads; ++t)
      nodeToElements[i].insert(nodeToElements[i].end(),
                               thread_ne[t][i].begin(),
                               thread_ne[t][i].end());
  }

  // Prefix each row with the count of connected elements
#pragma omp parallel for
  for (int i = 0; i < numNodes; ++i) {
    std::vector<int> row;
    row.reserve(nodeToElements[i].size() + 1);
    row.push_back(static_cast<int>(nodeToElements[i].size()));
    row.insert(row.end(), nodeToElements[i].begin(), nodeToElements[i].end());
    nodeToElements[i] = std::move(row);
  }
}

// Wrapper for Mesh
void compute_node_adjacency_for_mesh(Mesh& mesh) {
  if (mesh.numNodes == 0)
    mesh.numNodes = mesh.nodeCoords.size();
  compute_node_adjacency(mesh.elementNodes, mesh.nodeCoords, mesh.nodeToElements);
}
