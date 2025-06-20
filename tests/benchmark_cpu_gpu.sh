#!/bin/bash

# === Always run from project root ===
cd "$(dirname "$0")/.."

# === Prompt user ===
if [ $# -eq 6 ]; then
    NDIV_START=$1
    NDIV_END=$2
    NDIV_STEP=$3
    THREAD_START=$4
    THREAD_END=$5
    THREAD_STEP=$6
else
    echo "Enter Ndiv range (start end step) for GPU/CPU (e.g., 1 10 1): "
    read NDIV_START NDIV_END NDIV_STEP
    echo "Enter Threads range (start end step) for CPU only (e.g., 1 12 1): "
    read THREAD_START THREAD_END THREAD_STEP
fi


# === Settings ===
EXEC=./bin/meshgen
LOG_DIR=./tests
LOG_CPU_ALL=$LOG_DIR/log_cpu.txt
LOG_GPU_ALL=$LOG_DIR/log_gpu.txt
BUILD_DIR=./build

# === Create log directory ===
mkdir -p "$LOG_DIR"
rm -f "$LOG_CPU_ALL" "$LOG_GPU_ALL"

# === GPU Version ===
echo "Compiling GPU version..."
rm -rf "$BUILD_DIR"
cmake -B "$BUILD_DIR" -S . -DUSE_CUDA=ON -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build "$BUILD_DIR" --parallel > /dev/null

echo "Benchmarking GPU version..."
for (( ndiv=$NDIV_START; ndiv<=$NDIV_END; ndiv+=$NDIV_STEP )); do
    ndiv_padded=$(printf "%02d" "$ndiv")
    echo "GPU: Ndiv = $ndiv"
    log_file="$LOG_DIR/log_gpu_ndiv${ndiv_padded}.txt"
    "$EXEC" --ndiv "$ndiv" --log "$log_file"
    cat "$log_file" >> "$LOG_GPU_ALL"
done

# === CPU Version ===
echo "Compiling CPU version..."
rm -rf "$BUILD_DIR"
cmake -B "$BUILD_DIR" -S . -DUSE_CUDA=OFF -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build "$BUILD_DIR" --parallel > /dev/null

echo "Benchmarking CPU version..."
for (( threads=$THREAD_START; threads<=$THREAD_END; threads+=$THREAD_STEP )); do
  for (( ndiv=$NDIV_START; ndiv<=$NDIV_END; ndiv+=$NDIV_STEP )); do
      threads_padded=$(printf "%02d" "$threads")
      ndiv_padded=$(printf "%02d" "$ndiv")
      echo "CPU: Threads = $threads, Ndiv = $ndiv"
      log_file="$LOG_DIR/log_cpu_ndiv${ndiv_padded}_threads${threads_padded}.txt"
      "$EXEC" --ndiv "$ndiv" --threads "$threads" --log "$log_file"
      cat "$log_file" >> "$LOG_CPU_ALL"
  done
done

echo "✅ Benchmark complete. Logs saved in $LOG_DIR"
