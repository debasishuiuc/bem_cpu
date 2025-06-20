## plot_speedup.py

import subprocess
import matplotlib.pyplot as plt
import numpy as np
import re
from collections import defaultdict
from tabulate import tabulate
import os
import glob

# === Paths ===
this_dir = os.getcwd()  # Use current directory for interactive environments
log_dir = this_dir
bash_script = os.path.join(this_dir, "benchmark_cpu_gpu.sh")

# === Ask user for benchmarking range (must match shell script input) ===
print("Provide the SAME input ranges you plan to use for benchmarking:")
ndiv_start = int(input("Ndiv start: "))
ndiv_end   = int(input("Ndiv end (inclusive): "))
ndiv_step  = int(input("Ndiv step: "))
thread_start = int(input("Thread start: "))
thread_end   = int(input("Thread end (inclusive): "))
thread_step  = int(input("Thread step: "))

ndiv_range = range(ndiv_start, ndiv_end + 1, ndiv_step)
thread_range = range(thread_start, thread_end + 1, thread_step)

# === Determine required log files ===
gpu_logs = [
    os.path.join(log_dir, f"log_gpu_ndiv{ndiv:02d}.txt")
    for ndiv in ndiv_range
]
cpu_logs = [
    os.path.join(log_dir, f"log_cpu_ndiv{ndiv:02d}_threads{threads:02d}.txt")
    for ndiv in ndiv_range
    for threads in thread_range
]
missing_logs = [f for f in (gpu_logs + cpu_logs) if not os.path.exists(f)]

# === Run benchmark if needed ===
if missing_logs:
    print("\nSome logs are missing. Running benchmark...")
    subprocess.run([
        "bash", bash_script,
        str(ndiv_start), str(ndiv_end), str(ndiv_step),
        str(thread_start), str(thread_end), str(thread_step)
    ], check=True)
else:
    print("\n✅ All required logs found. Skipping benchmark run.")

# === Parse logs ===
cpu_data = defaultdict(
    lambda: defaultdict(
        lambda: {'mesh': [], 'geom': [], 'total': [], 'area': [], 'volume': []}
    )
)
gpu_data = defaultdict(lambda: {'mesh': [], 'geom': [], 'total': []})

pattern_cpu = re.compile(
    r"Ndiv=(\d+), Mode=OpenMP, Threads=(\d+), MeshTime=([\d\.eE+-]+), GeometryTime=([\d\.eE+-]+), "
    r"Npts=(\d+), Nelm=(\d+), SurfaceArea=([\d\.eE+-]+), Volume=([\d\.eE+-]+), Centroid=\[([^\]]+)\]"
)

pattern_gpu = re.compile(
    r"Ndiv=(\d+), Mode=CUDA, MeshTime=([\d\.eE+-]+), GeometryTime=([\d\.eE+-]+)"
)

print("\n=== DEBUG: Parsing GPU log files ===")
for filepath in glob.glob(os.path.join(log_dir, "log_gpu_ndiv*.txt")):
    print(f"Reading {filepath}")
    with open(filepath) as f:
        for line in f:
            print(f"[GPU line] {line.strip()}")
            match = pattern_gpu.search(line.strip())
            if match:
                ndiv = int(match.group(1))
                mesh = float(match.group(2))
                geom = float(match.group(3))
                gpu_data[ndiv]['mesh'].append(mesh)
                gpu_data[ndiv]['geom'].append(geom)
                gpu_data[ndiv]['total'].append(mesh + geom)

print("\n=== DEBUG: Parsing CPU log files ===")
for filepath in glob.glob(os.path.join(log_dir, "log_cpu_ndiv*_threads*.txt")):
    print(f"Reading {filepath}")
    with open(filepath) as f:
        for line in f:
            print(f"[CPU line] {line.strip()}")
            match = pattern_cpu.search(line.strip())
            if match:
                ndiv = int(match.group(1))
                threads = int(match.group(2))
                mesh = float(match.group(3))
                geom = float(match.group(4))
                area = float(match.group(7))
                vol = float(match.group(8))
                cpu_data[threads][ndiv]['mesh'].append(mesh)
                cpu_data[threads][ndiv]['geom'].append(geom)
                cpu_data[threads][ndiv]['total'].append(mesh + geom)
                cpu_data[threads][ndiv]['area'].append(area)
                cpu_data[threads][ndiv]['volume'].append(vol)

# === Timing Summary Table ===
print("\n=== Timing Summary (Grouped by Ndiv) ===\n")
ndivs_all = sorted(set().union(*(cpu_data[t].keys() for t in cpu_data)))

for ndiv in ndivs_all:
    rows = []
    gpu_time = np.mean(gpu_data.get(ndiv, {}).get('total', [np.nan]))
    for threads in sorted(cpu_data.keys()):
        cpu_times = cpu_data[threads].get(ndiv, {}).get('total', [])
        if cpu_times:
            avg_cpu = np.mean(cpu_times)
            speedup = avg_cpu / gpu_time if not np.isnan(gpu_time) and gpu_time > 0 else float('nan')
            rows.append([threads, f"{avg_cpu:.5e}", f"{gpu_time:.5e}", f"{speedup:.2f}"])
    if rows:
        print(f"\nNdiv = {ndiv}")
        print(tabulate(rows, headers=["Threads", "CPU Time (s)", "GPU Time (s)", "Speedup"]))

# === Plot line plots ===
timing_categories = ['mesh', 'geom', 'total']
titles = {
    'mesh': "Mesh Generation Time vs Ndiv",
    'geom': "Geometry Computation Time vs Ndiv",
    'total': "Total Time vs Ndiv"
}
ylabels = {
    'mesh': "Mesh Time (s)",
    'geom': "Geometry Time (s)",
    'total': "Total Time (s)"
}
filenames = {
    'mesh': "timing_mesh_vs_ndiv.png",
    'geom': "timing_geom_vs_ndiv.png",
    'total': "timing_total_vs_ndiv.png"
}

for cat in timing_categories:
    plt.figure()
    for threads in sorted(cpu_data.keys()):
        ndivs = sorted(cpu_data[threads].keys())
        avg_times = [np.mean(cpu_data[threads][n][cat]) for n in ndivs]
        plt.plot(ndivs, avg_times, marker='o', label=f"OpenMP {threads} threads")

    gpu_ndivs = sorted(gpu_data.keys())
    gpu_avg = [np.mean(gpu_data[n][cat]) for n in gpu_ndivs]
    if gpu_avg:
        plt.plot(gpu_ndivs, gpu_avg, marker='s', linestyle='--', color='black', label="GPU (CUDA)")

    plt.xlabel("Ndiv")
    plt.ylabel(ylabels[cat])
    plt.title(titles[cat])
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    filepath = os.path.join(this_dir, filenames[cat])
    plt.savefig(filepath)
    print(f"✅ Saved: {filepath}")
    plt.close()

# === Bar plots ===
for cat in timing_categories:
    ndivs_sorted = sorted(set(gpu_data.keys()).union(*(cpu_data[t].keys() for t in cpu_data)))
    x = np.arange(len(ndivs_sorted))
    width = 0.35

    fig, ax = plt.subplots()

    max_threads = max(cpu_data.keys())
    cpu_avg = [np.mean(cpu_data[max_threads][n][cat]) if n in cpu_data[max_threads] else np.nan for n in ndivs_sorted]
    gpu_avg = [np.mean(gpu_data[n][cat]) if n in gpu_data else np.nan for n in ndivs_sorted]

    ax.bar(x - width/2, cpu_avg, width, label=f"CPU ({max_threads} threads)")
    ax.bar(x + width/2, gpu_avg, width, label="GPU (CUDA)")

    ax.set_xlabel("Ndiv")
    ax.set_ylabel(ylabels[cat])
    ax.set_title(f"{titles[cat]} (Bar Plot)")
    ax.set_xticks(x)
    ax.set_xticklabels(ndivs_sorted)
    ax.legend()
    ax.grid(True, axis='y')

    fig.tight_layout()
    path = os.path.join(this_dir, f"barplot_{cat}_cpu_vs_gpu.png")
    plt.savefig(path)
    print(f"✅ Saved: {path}")
    plt.close()

# === Plot Surface Area and Volume for max threads ===
max_threads = max(cpu_data.keys())
ndivs_area_vol = sorted(cpu_data[max_threads].keys())

areas = [np.mean(cpu_data[max_threads][n]['area']) for n in ndivs_area_vol]
volumes = [np.mean(cpu_data[max_threads][n]['volume']) for n in ndivs_area_vol]

# Surface Area Plot
plt.figure()
plt.plot(ndivs_area_vol, areas, marker='o', color='green')
plt.xlabel("Ndiv")
plt.ylabel("Surface Area")
plt.title(f"Surface Area vs Ndiv (CPU {max_threads} threads)")
plt.grid(True)
plt.tight_layout()
f_area = os.path.join(this_dir, "surface_area_vs_ndiv.png")
plt.savefig(f_area)
print(f"✅ Saved: {f_area}")
plt.close()

# Volume Plot
plt.figure()
plt.plot(ndivs_area_vol, volumes, marker='o', color='blue')
plt.xlabel("Ndiv")
plt.ylabel("Volume")
plt.title(f"Volume vs Ndiv (CPU {max_threads} threads)")
plt.grid(True)
plt.tight_layout()
f_vol = os.path.join(this_dir, "volume_vs_ndiv.png")
plt.savefig(f_vol)
print(f"✅ Saved: {f_vol}")
plt.close()

