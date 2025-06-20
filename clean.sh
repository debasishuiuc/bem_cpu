#!/bin/bash
echo "Cleaning build and output directories..."
#!/bin/bash
rm -rf build CMakeFiles CMakeCache.txt
echo "creating build and output directories..."
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug -DUSE_CUDA=OFF ..
make -j
echo "Done."
