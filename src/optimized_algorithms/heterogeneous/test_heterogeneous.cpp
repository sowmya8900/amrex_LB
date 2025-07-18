#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include "HeterogeneousLB.H"
#include <map>
#include <iomanip>
#include <random>

using namespace amrex;

void verify_distribution(const BoxArray& ba, const DistributionMapping& dm, 
                        const Vector<ComputeNodeInfo>& nodes, const std::string& strategy,
                        const MultiFab& weights, const HeterogeneousLB& lb) {
    // Count boxes per node
    std::map<int, int> boxes_per_node;
    std::map<int, Real> load_per_node;
    
    for (int i = 0; i < ba.size(); ++i) {
        int node_id = dm[i];
        boxes_per_node[node_id]++;
        load_per_node[node_id] += ba[i].numPts();
    }
    
    // Evaluate load balance quality
    auto metrics = lb.EvaluateBalance(dm, weights);
    
    // Print distribution statistics
    amrex::Print() << "\nDistribution Statistics for " << strategy << " strategy:\n";
    amrex::Print() << std::setw(8) << "Node ID" 
                   << std::setw(15) << "CPU Type" 
                   << std::setw(12) << "Perf Factor"
                   << std::setw(12) << "Box Count"
                   << std::setw(15) << "Load"
                   << std::setw(12) << "Load %" << "\n";
    amrex::Print() << std::string(84, '-') << "\n";
    
    Real total_load = 0.0;
    for (const auto& pair : load_per_node) {
        total_load += pair.second;
    }
    
    for (const auto& node : nodes) {
        Real load_pct = (total_load > 0) ? (load_per_node[node.node_id] / total_load * 100.0) : 0.0;
        amrex::Print() << std::setw(8) << node.node_id
                       << std::setw(15) << node.node_type.substr(0, 14)
                       << std::setw(12) << std::fixed << std::setprecision(3) << node.performance_factor
                       << std::setw(12) << boxes_per_node[node.node_id]
                       << std::setw(15) << std::fixed << std::setprecision(1) << load_per_node[node.node_id]
                       << std::setw(12) << std::fixed << std::setprecision(1) << load_pct << "\n";
    }
    
    // Print quality metrics
    amrex::Print() << "\nLoad Balance Quality Metrics:\n";
    amrex::Print() << "  Makespan: " << std::fixed << std::setprecision(3) << metrics.makespan << "\n";
    amrex::Print() << "  Efficiency: " << std::fixed << std::setprecision(3) << metrics.efficiency << "\n";
    amrex::Print() << "  Load Imbalance: " << std::fixed << std::setprecision(2) 
                   << (metrics.load_imbalance * 100.0) << "%\n";
    amrex::Print() << "  Total Work: " << metrics.total_work << "\n";
}

void test_heterogeneous_lb() {
    amrex::Print() << "\n=== Testing Heterogeneous Load Balancer ===\n";
    
    // Test Case 1: Basic configuration with 3 different CPU types
    Vector<ComputeNodeInfo> nodes = {
        {0, "cpu0", 1.0, 256.0},
        {1, "cpu1", 0.9, 128.0},
        {2, "cpu2", 0.8, 128.0}
    };
    
    // Create and initialize load balancer
    HeterogeneousLB lb;
    lb.InitializeNodes(nodes);
    
    // Create test domain and boxes
    Box domain(IntVect(0), IntVect(63));  // 64^3 domain
    BoxArray ba(domain);
    ba.maxSize(16);  // Create multiple boxes
    
    // Create weights with spatial variation
    DistributionMapping dm(ba);
    MultiFab weights(ba, dm, 1, 0);
    
    // Set weights with spatial variation
    for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        IntVect center = box.smallEnd() + box.size() / 2;
        Real distance = std::sqrt(
            std::pow(center[0] - 32.0, 2) +
            std::pow(center[1] - 32.0, 2) +
            std::pow(center[2] - 32.0, 2)
        );
        weights[mfi].setVal(box.numPts() * (1.0 + std::exp(-distance/16.0)));
    }
    
    // Test all strategies
    Vector<std::string> strategies = {"rij", "knapsack", "sfc", "grouped"};
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.9, 1.1);  // 10% variation
    
    for (const auto& strategy : strategies) {
        amrex::Print() << "\nTesting " << strategy << " strategy...\n";
        try {
            auto new_dm = lb.BalanceLoad(ba, weights, strategy);
            verify_distribution(ba, new_dm, nodes, strategy, weights, lb);
            
            // Simulate some computation with random noise
            Vector<double> fake_timings;
            for (const auto& node : nodes) {
                double base_time = 1.0 / node.performance_factor;
                fake_timings.push_back(base_time * dis(gen));
            }
            lb.UpdatePerformanceMetrics(fake_timings);
            lb.PrintStats();
        } catch (const std::exception& e) {
            amrex::Print() << "Error testing " << strategy << ": " << e.what() << "\n";
        }
    }
    
    // Test Case 2: Larger system with more nodes
    amrex::Print() << "\n=== Testing with larger system configuration ===\n";
    
    nodes = {
        {0, "cpu0", 1.0, 256.0},
        {1, "cpu1", 1.0, 256.0},
        {2, "cpu2", 0.9, 128.0},
        {3, "cpu3", 0.9, 128.0},
        {4, "cpu4", 0.8, 128.0},
        {5, "cpu5", 0.8, 128.0}
    };
    
    lb.InitializeNodes(nodes);
    
    // Test with larger domain
    Box large_domain(IntVect(0), IntVect(127));  // 128^3 domain
    BoxArray large_ba(large_domain);
    large_ba.maxSize(32);  // Larger boxes
    
    DistributionMapping large_dm(large_ba);
    MultiFab large_weights(large_ba, large_dm, 1, 0);
    
    // Set weights
    for (MFIter mfi(large_weights); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        IntVect center = box.smallEnd() + box.size() / 2;
        Real distance = std::sqrt(
            std::pow(center[0] - 64.0, 2) +
            std::pow(center[1] - 64.0, 2) +
            std::pow(center[2] - 64.0, 2)
        );
        large_weights[mfi].setVal(box.numPts() * (1.0 + std::exp(-distance/32.0)));
    }
    
    for (const auto& strategy : strategies) {
        amrex::Print() << "\nTesting " << strategy << " strategy on larger system...\n";
        try {
            auto new_dm = lb.BalanceLoad(large_ba, large_weights, strategy);
            verify_distribution(large_ba, new_dm, nodes, strategy, large_weights, lb);
            
            // Simulate timings with random noise
            Vector<double> fake_timings;
            for (const auto& node : nodes) {
                double base_time = 1.0 / node.performance_factor;
                fake_timings.push_back(base_time * dis(gen));
            }
            lb.UpdatePerformanceMetrics(fake_timings);
            lb.PrintStats();
        } catch (const std::exception& e) {
            amrex::Print() << "Error testing " << strategy << " on larger system: " << e.what() << "\n";
        }
    }
    
    // Test Case 3: Performance update test
    amrex::Print() << "\n=== Testing Performance Updates ===\n";
    
    // Test individual node updates
    for (const auto& node : nodes) {
        double new_timing = dis(gen) / node.performance_factor;
        lb.UpdateSingleNodeMetrics(node.node_id, new_timing);
        amrex::Print() << "Updated node " << node.node_id << " with timing " << new_timing << "\n";
    }
    
    lb.PrintStats();
}