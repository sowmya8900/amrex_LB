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
    amrex::SetVerbose(1);  // Reduce verbosity
    
    // Initialize random number generator
    InitRandom(ParallelDescriptor::MyProc() + 1);
    
    // Example heterogeneous system setup with different CPU types - using fewer nodes for local test
    Vector<ComputeNodeInfo> nodes;
    
    // Add CPU nodes with different performance characteristics
    nodes.push_back({0, "CPU_AMD_EPYC_7763", 1.0, 256.0});    // Base performance
    nodes.push_back({1, "CPU_Intel_8380", 0.9, 128.0});       // 10% slower
    nodes.push_back({2, "CPU_AMD_EPYC_7542", 0.8, 128.0});    // 20% slower
    
    // Create and initialize the load balancer
    HeterogeneousLB lb;
    lb.InitializeNodes(nodes);
    
    // Create a sample problem with spatial variation in weights - smaller domain for local test
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
    for (const auto& strategy : strategies) {
        amrex::Print() << "\nTesting " << strategy << " balancing:\n";
        auto dm = lb.BalanceLoad(ba, weights, strategy);
        lb.PrintStats();
        // Simulate timings with random noise
        Vector<double> fake_timings;
        for (const auto& node : nodes) {
            double base_time = 1.0 / node.performance_factor;
            fake_timings.push_back(base_time * dis(gen));
        }
        lb.UpdatePerformanceMetrics(fake_timings);
        amrex::Print() << "\nFinal statistics after performance update:\n";
        lb.PrintStats();
    }
} 