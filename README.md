# BEM CPU: OpenMP-Based Boundary Element Method Mesh Generator

This repository implements a fully CPU-based boundary element method (BEM) mesh generator and geometry analyzer using OpenMP for parallelization. The original CUDA-based code has been entirely removed for simplicity and portability.

---

## ✅ Current Status

- Pure C++ implementation with OpenMP parallelism.
- No GPU or CUDA dependencies.
- Generates 6-node curved triangular meshes on a unit sphere.
- Computes:
  - Element-to-element connectivity (`nbe`)
  - Node-to-element connectivity (`ne`)
  - Accurate surface area, volume, and centroid
  - Element and node normals
  - Mean curvature per element
- Supports text and VTK output for visualization in ParaView.
- Fully validated for the base icosahedron (Ndiv = 0) against legacy outputs.

---

## 📁 Project Structure

```
bem_cpu/
├── bin/             # Executables
├── build/           # CMake build artifacts
├── include/         # C++ header files
├── main/            # Main entry point (main.cpp)
├── output/          # Mesh data output: .txt and .vtk
├── src/             # C++ source files
├── tests/           # Benchmark scripts and timing analysis
├── CMakeLists.txt   # Build configuration
├── Dockerfile       # Optional Docker support
├── README.md        # This file
└── run_docker.sh    # Docker helper script
```

---

## ⚙️ Features

### Mesh Generation
- Recursive subdivision of an icosahedron.
- Six-node curved triangle elements (`trgl6_icos.cpp`).
- Midpoint insertion and deduplication (`deduplicate.cpp`).

### Connectivity
- Element-to-element (`nbe`) and node-to-element (`ne`) connectivity.
- Verified against reference Fortran code.

### Geometry Analysis
- Computation of:
  - Total surface area
  - Volume enclosed
  - Centroid
  - Element and node normals
  - Element curvature
  - Inverse moment matrix (`Mmat⁻¹`)

### Output
- `.txt` files: `p.txt`, `n.txt`, `ne.txt`, `nbe.txt`, `crvmel.txt`, etc.
- `.vtk` file for visualization in ParaView.

---

## 🧪 Tests and Benchmarks

The `tests/` folder contains:
- `plot_speedup.py`: script to visualize OpenMP speedup.
- `timing.py`: analyze performance across thread counts.
- Legacy scripts referencing GPU benchmarks (can be cleaned up).

---

## 🛠️ Build Instructions

```bash
mkdir build
cd build
cmake ..
make -j
./bin/meshgen --ndiv 0 --write
```

---

## 🔐 GitHub SSH Setup

SSH access to GitHub has been configured and the repository has been successfully pushed to:

📍 `git@github.com:debasishuiuc/bem_cpu.git`

---

## 🔄 Future Work

- Add formal test cases using `CTest`.
- Benchmark OpenMP performance for higher subdivisions.
- Explore reintroducing GPU parallelism on a separate branch if needed.

---

## 📞 Contact

Maintained by [debasishuiuc](https://github.com/debasishuiuc).
