# BEM CPU vs GPU Mesh Generator

This project implements and benchmarks a recursive 6-node triangle mesh generator for the surface of a unit sphere using:

* 🧠 **CPU version**: Parallelized with OpenMP
* 🚀 **GPU version**: Written in CUDA C for acceleration

Benchmarks compare CPU and GPU performance across different mesh subdivision levels (`Ndiv`).

---

## 📁 Project Structure

```
bem_cpu_gpu/
├── bin/                   # Executables
├── build/                 # CMake build directory
├── include/               # Header files
├── main/                  # main.cpp and driver files
├── output/                # Mesh output files
├── src/                   # Source code (CPU and CUDA)
├── tests/                 # Benchmarking scripts + logs
├── CMakeLists.txt         # CMake config
├── Dockerfile             # (Optional) container support
├── benchmark_cpu_gpu.sh   # Benchmark automation
├── plot_speedup.py        # Benchmark analysis & plotting
└── README.md              # This file
```

---

## ⚙️ Compilation

### CPU Version

```bash
cmake -DUSE_CUDA=OFF -B build -S .
cmake --build build --parallel
```

### GPU Version

```bash
cmake -DUSE_CUDA=ON -B build -S .
cmake --build build --parallel
```

The executable will be at `./bin/meshgen`.

---

## 🚀 Running the Mesh Generator

```bash
./bin/meshgen --ndiv <N> --threads <num_threads> [--write]
```

* `--ndiv <N>`: subdivision level
* `--threads <N>`: number of CPU threads
* `--write`: writes output to `output/` as `.txt` and `.vtk`

---

## 📊 Benchmarking

Run the benchmarking script from the **project root**:

```bash
python tests/plot_speedup.py
```

This:

* Compiles CPU and GPU versions
* Runs each with `Ndiv = 1` to `6`
* Stores logs in `tests/`
* Generates plots:

  * `timing_comparison.png`
  * `speedup.png`

---

## 📆 Dependencies

* CMake ≥ 3.10
* CUDA Toolkit (for GPU build)
* Python ≥ 3.8 with:

  * `matplotlib`
  * `numpy`

---

## 📈 Sample Output

![timing\_comparison.png](tests/timing_comparison.png)
![speedup.png](tests/speedup.png)

---

## 📝 License

MIT License (you may modify this as needed).