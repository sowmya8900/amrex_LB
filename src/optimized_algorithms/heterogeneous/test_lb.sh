#!/bin/bash

# Make sure the executable is built
make -j4 USE_MPI=FALSE

# Run tests with different strategies
echo "Testing Knapsack Strategy..."
./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=knapsack

echo -e "\nTesting SFC Strategy..."
./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=sfc

echo -e "\nTesting Grouped Strategy..."
./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=grouped 