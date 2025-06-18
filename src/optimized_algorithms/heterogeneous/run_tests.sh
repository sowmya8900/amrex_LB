#!/bin/bash

echo "Building heterogeneous load balancer tests..."
make -j4

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo -e "\n" + "="*60
echo "HETEROGENEOUS LOAD BALANCER TEST SUITE"
echo "="*60

# Create outputs directory
mkdir -p outputs

# Test 1: Basic functionality test
echo -e "\n1. Running basic functionality test..."
./main3d.gnu.x86-milan.OMP.ex test_type=test > outputs/test1_basic.txt 2>&1

# Test 2: Example demonstration  
echo -e "\n2. Running example demonstration..."
./main3d.gnu.x86-milan.OMP.ex test_type=example > outputs/test2_example.txt 2>&1

# Test 3: 4-node high heterogeneity test
echo -e "\n3. Running 4-node high heterogeneity test..."
./main3d.gnu.x86-milan.OMP.ex nnodes=4 \
    node_types="cpu0 cpu1 cpu2 cpu3" \
    perf_factors="2.0 1.0 0.5 0.25" \
    memory_caps="512.0 256.0 128.0 64.0" \
    strategy=all > outputs/test3_4node_hetero.txt 2>&1

# Test 4: 6-node balanced test
echo -e "\n4. Running 6-node balanced test..."
./main3d.gnu.x86-milan.OMP.ex nnodes=6 \
    node_types="cpu0 cpu1 cpu2 cpu3 cpu4 cpu5" \
    perf_factors="1.0 1.0 0.9 0.9 0.8 0.8" \
    memory_caps="256.0 256.0 128.0 128.0 128.0 128.0" \
    strategy=all > outputs/test4_6node_balanced.txt 2>&1

# Test 5: 8-node scalability test
echo -e "\n5. Running 8-node scalability test..."
./main3d.gnu.x86-milan.OMP.ex nnodes=8 \
    domain="256 256 256" max_grid_size="64 64 64" \
    strategy=all > outputs/test5_8node_scalability.txt 2>&1

# Test 6: Homogeneous system comparison
echo -e "\n6. Running homogeneous system test..."
./main3d.gnu.x86-milan.OMP.ex nnodes=4 \
    node_types="cpu0 cpu1 cpu2 cpu3" \
    perf_factors="1.0 1.0 1.0 1.0" \
    memory_caps="256.0 256.0 256.0 256.0" \
    strategy=all > outputs/test6_homogeneous.txt 2>&1

# Combine all outputs for analysis
echo -e "\n7. Combining outputs for analysis..."
cat outputs/test*.txt > hetero_output.txt

echo -e "\n8. Generating analysis plots and reports..."
python3 plot_hetero_output.py

echo -e "\n" + "="*60
echo "ALL TESTS COMPLETED"
echo "="*60
echo -e "\nResults saved in outputs/ directory:"
echo "• Individual test outputs: test*.txt"
echo "• Combined output: hetero_output.txt"
echo "• Analysis plots: *.png"
echo "• Performance report: performance_report.txt"
echo -e "\nCheck the outputs for detailed results and visualizations."

# #!/bin/bash

# echo "Building tests..."
# make -j4

# echo -e "\nRunning heterogeneous load balancer tests..."
# ./main3d.gnu.x86-milan.OMP.ex test_type=test

# # echo -e "\nRunning example with different strategies..."
# # ./main3d.gnu.x86-milan.OMP.ex test_type=example
# # echo -e "\nKnapsack strategy:"
# # ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=knapsack

# # echo -e "\nSFC strategy:"
# # ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=sfc

# # echo -e "\nGrouped strategy:"
# # ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=grouped

# # echo -e "\nRij strategy:"
# # ./main3d.gnu.x86-milan.OMP.ex test_type=example strategy=rij

# echo -e "\nTests completed. Check output for results." 