#!/bin/bash

echo "================================================"
echo "Compiling local CNPY library for HiSparse..."
echo "================================================"

# Create local_cnpy directory
mkdir -p local_cnpy
cd local_cnpy

# Copy source files if they don't exist
if [ ! -f "cnpy.cpp" ]; then
    echo "Copying CNPY source files..."
    if [ -f "/work/shared/common/project_build/graphblas/software/cnpy/cnpy.h" ]; then
        cp /work/shared/common/project_build/graphblas/software/cnpy/cnpy.h .
        cp /work/shared/common/project_build/graphblas/software/cnpy/cnpy.cpp .
        echo "✓ CNPY source files copied"
    else
        echo "✗ Error: CNPY source files not found at expected location"
        echo "Please ensure CNPY is available or download it from:"
        echo "https://github.com/rogersce/cnpy"
        exit 1
    fi
fi

# Compile CNPY library
echo "Compiling CNPY with current GCC version..."
g++ -fPIC -shared -std=c++14 cnpy.cpp -o libcnpy.so -lz

# Check compilation result
if [ -f "libcnpy.so" ]; then
    echo "✓ CNPY library compiled successfully"
    
    # Show library info
    echo ""
    echo "Library information:"
    echo "  File: $(pwd)/libcnpy.so"
    echo "  Size: $(ls -lh libcnpy.so | awk '{print $5}')"
    
    # Check GCC version used
    GCC_VERSION=$(strings libcnpy.so | grep "GCC:" | head -1)
    if [ -n "$GCC_VERSION" ]; then
        echo "  Compiled with: $GCC_VERSION"
    fi
    
    echo ""
    echo "✓ Local CNPY setup complete!"
    echo "Now run: source setup.sh"
else
    echo "✗ Error: Failed to compile CNPY library"
    echo "Please check for compilation errors above"
    exit 1
fi