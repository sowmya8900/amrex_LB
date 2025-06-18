#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include "HeterogeneousLB.H"
#include <random>

using namespace amrex;

void example_main()
{
    // Set verbose output level
    amrex::SetVerbose(1);
    
    // Initialize random number generator
    InitRandom(ParallelDescriptor::MyProc() + 1);
    
    // Example heterogeneous system setup with different CPU types
    Vector<ComputeNodeInfo> nodes;
    
    // Add CPU nodes with different performance characteristics
    nodes.push_back({0, "cpu0", 1.0, 256.0});    // Base performance
    nodes.push_back({1, "cpu1", 0.9, 128.0});       // 10% slower
    nodes.push_back({2, "cpu2", 0.8, 128.0});    // 20% slower
    
    // Create and initialize the load balancer
    HeterogeneousLB lb;
    lb.InitializeNodes(nodes);
    
    // Create a sample problem with spatial variation in weights
    Box domain(IntVect(0), IntVect(31));  // 32^3 domain
    BoxArray ba(domain);
    ba.maxSize(16);  // 16^3 boxes
    
    // Create weights with spatial variation (e.g., more work in center)
    MultiFab weights(ba, DistributionMapping(ba), 1, 0);
    for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        IntVect center = box.smallEnd() + box.size() / 2;
        Real distance = std::sqrt(
            std::pow(center[0] - 16.0, 2) +
            std::pow(center[1] - 16.0, 2) +
            std::pow(center[2] - 16.0, 2)
        );
        // More weight near center
        weights[mfi].setVal(box.numPts() * (1.0 + std::exp(-distance/16.0)), box, 0, 1);
    }
    
    // Try all load balancing strategies
    amrex::Print() << "\nTesting all load balancing strategies:\n";
    
    Vector<std::string> strategies = {"rij", "knapsack", "sfc", "grouped"};
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.9, 1.1);  // 10% variation
    
    Vector<LoadBalanceMetrics> all_metrics;
    
    for (const auto& strategy : strategies) {
        amrex::Print() << "\n" << std::string(50, '=') << "\n";
        amrex::Print() << "Testing " << strategy << " balancing:\n";
        amrex::Print() << std::string(50, '=') << "\n";
        
        auto dm = lb.BalanceLoad(ba, weights, strategy);
        auto metrics = lb.EvaluateBalance(dm, weights);
        all_metrics.push_back(metrics);
        
        amrex::Print() << "Strategy: " << strategy << "\n";
        amrex::Print() << "  Makespan: " << std::fixed << std::setprecision(3) << metrics.makespan << "\n";
        amrex::Print() << "  Efficiency: " << std::fixed << std::setprecision(3) << metrics.efficiency << "\n";
        amrex::Print() << "  Load Imbalance: " << std::fixed << std::setprecision(2) 
                       << (metrics.load_imbalance * 100.0) << "%\n";
        
        lb.PrintStats();
        
        // Simulate timings with random noise
        Vector<double> fake_timings;
        for (const auto& node : nodes) {
            double base_time = 1.0 / node.performance_factor;
            fake_timings.push_back(base_time * dis(gen));
        }
        lb.UpdatePerformanceMetrics(fake_timings);
        
        amrex::Print() << "\nAfter performance update:\n";
        lb.PrintStats();
    }
    
    // Compare all strategies
    amrex::Print() << "\n" << std::string(60, '=') << "\n";
    amrex::Print() << "STRATEGY COMPARISON SUMMARY\n";
    amrex::Print() << std::string(60, '=') << "\n";
    amrex::Print() << std::setw(12) << "Strategy"
                   << std::setw(12) << "Makespan"
                   << std::setw(12) << "Efficiency"
                   << std::setw(15) << "Imbalance %" << "\n";
    amrex::Print() << std::string(51, '-') << "\n";
    
    for (size_t i = 0; i < strategies.size(); ++i) {
        amrex::Print() << std::setw(12) << strategies[i]
                       << std::setw(12) << std::fixed << std::setprecision(3) << all_metrics[i].makespan
                       << std::setw(12) << std::fixed << std::setprecision(3) << all_metrics[i].efficiency
                       << std::setw(15) << std::fixed << std::setprecision(2) 
                       << (all_metrics[i].load_imbalance * 100.0) << "\n";
    }
    
    // Find best strategy
    auto best_it = std::min_element(all_metrics.begin(), all_metrics.end(),
        [](const LoadBalanceMetrics& a, const LoadBalanceMetrics& b) {
            return a.makespan < b.makespan;
        });
    
    if (best_it != all_metrics.end()) {
        size_t best_idx = std::distance(all_metrics.begin(), best_it);
        amrex::Print() << "\nBest strategy: " << strategies[best_idx] 
                       << " (makespan: " << best_it->makespan << ")\n";
    }
}