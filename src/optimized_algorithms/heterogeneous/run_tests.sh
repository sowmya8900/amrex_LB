#!/bin/bash

echo "Building tests..."
make -j4

echo -e "\nRunning heterogeneous load balancer tests..."
./main3d.gnu.x86-milan.OMP.ex test_type=test

echo -e "\nRunning example with different strategies..."
./main3d.gnu.x86-milan.OMP.ex test_type=example
# echo -e "\nKnapsack strategy:"
# ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=knapsack

# echo -e "\nSFC strategy:"
# ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=sfc

# echo -e "\nGrouped strategy:"
# ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=grouped

# echo -e "\nRij strategy:"
# ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=rij

echo -e "\nTests completed. Check output for results." 