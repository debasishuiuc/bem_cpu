#!/bin/bash
echo "Cleaning build and output directories..."
#!/bin/bash
rm -rf build CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake
echo "creating build and output directories..."
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
echo "Done."
