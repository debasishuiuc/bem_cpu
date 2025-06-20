import subprocess
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import os
import re
from tabulate import tabulate

# === Paths ===
this_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(this_dir, ".."))
build_dir = os.path.join(root_dir, "build")
executable = os.path.join(root_dir, "bin", "meshgen")
logfile = os.path.join(this_dir, "log_threads.txt")
max_threads = 12

# === User input ===
max_ndiv = int(input("Enter maximum subdivision level (Ndiv): "))

# === Step 1: Recompile with OpenMP (CUDA OFF) ===
print("🔧 Compiling mesh generator with OpenMP (CUDA OFF)...")
if os.path.exists(build_dir):
    subprocess.run(["rm", "-rf", build_dir], check=True)

subprocess.run(["cmake", "-DUSE_CUDA=OFF", "-B", build_dir, "-S", root_dir], check=True)
subprocess.run(["cmake", "--build", build_dir, "--parallel"], check=True)

# === Optional: clear old log ===
if os.path.exists(logfile):
    os.remove(logfile)

# === Step 2: Run benchmark ===
print(f"\n🚀 Benchmarking: {executable}")
for ndiv in range(1, max_ndiv + 1):
    print(f"\n>>> Ndiv = {ndiv}")
    for t in range(1, max_threads + 1):
        print(f"  Running with {t} thread(s)...")
        result = subprocess.run(
            [executable, "--ndiv", str(ndiv), "--threads", str(t)],
            capture_output=True, text=True
        )
        # Append output to log file
        with open(logfile, "a") as f:
            f.write(result.stdout)
            f.write(result.stderr)

# === Step 3: Parse log file ===
timing_data = defaultdict(lambda: defaultdict(list))
pattern = re.compile(r"Ndiv\s*=\s*(\d+)[^\n]*?Threads\s*=\s*(\d+)[^\n]*?Time\s*=\s*([\d\.eE+-]+)")

with open(logfile) as f:
    for line in f:
        match = pattern.search(line)
        if match:
            ndiv = int(match.group(1))
            threads = int(match.group(2))
            time = float(match.group(3))
            timing_data[ndiv][threads].append(time)

# === Step 4: Plot results ===
plt.figure()
print("\n=== Benchmark Summary ===")
for ndiv in sorted(timing_data.keys()):
    threads = sorted(timing_data[ndiv].keys())
    avg_times = [np.mean(timing_data[ndiv][t]) for t in threads]
    counts = [len(timing_data[ndiv][t]) for t in threads]

    # Print table for this Ndiv
    table_data = [(t, f"{avg:.5f}", c) for t, avg, c in zip(threads, avg_times, counts)]
    print(f"\nNdiv = {ndiv}")
    print(tabulate(table_data, headers=["Threads", "Avg Time (s)", "#Runs"], tablefmt="pretty"))

    # Plot line
    plt.plot(threads, avg_times, marker='o', label=f"Ndiv = {ndiv}")

# Finalize plot
plt.xlabel("Number of Threads")
plt.ylabel("Average Time (s)")
plt.title("Mesh Generation Time vs Number of Threads")
plt.legend()
plt.grid(True)
plt.tight_layout()
plot_path = os.path.join(this_dir, "scaling_all_ndiv.png")
plt.savefig(plot_path)
print(f"\n✅ Saved plot: {plot_path}")
plt.close()
