#ifndef AMREX_HETEROGENEOUSLB_H_
#define AMREX_HETEROGENEOUSLB_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Vector.H>
#include "../Knapsack.H"
#include "../SFC.H"
#include <map>
#include <string>
#include <vector>

namespace amrex {

struct ComputeNodeInfo {
    int node_id;
    std::string node_type;      // CPU, GPU, FPGA, etc.
    double performance_factor;   // Relative to baseline performance
    double memory_capacity;      // In GB
    double recent_performance = 0.0;  // Initialize to 0
    
    ComputeNodeInfo(int id = 0, 
                   const std::string& type = "",
                   double perf = 1.0,
                   double mem = 0.0)
        : node_id(id)
        , node_type(type)
        , performance_factor(perf)
        , memory_capacity(mem)
    {}
};

// Structure to hold timing ratio information
struct TimingRatio {
    int i;  // source node
    int j;  // target node
    double rij;  // timing ratio
    
    TimingRatio(int src = 0, int tgt = 0, double ratio = 1.0)
        : i(src), j(tgt), rij(ratio) {}
    
    bool operator<(const TimingRatio& other) const {
        return rij < other.rij;
    }
};

class HeterogeneousLB {
public:
    // Constructor
    HeterogeneousLB();
    
    // Initialize node information
    void InitializeNodes(const Vector<ComputeNodeInfo>& node_info);
    
    // Update performance metrics
    void UpdatePerformanceMetrics(const Vector<double>& timings);
    
    // Main load balancing methods
    DistributionMapping BalanceLoad(const BoxArray& ba, 
                                  const MultiFab& weights,
                                  const std::string& strategy = "rij");
    
    // Different load balancing strategies
    DistributionMapping KnapsackBalance(const BoxArray& ba, const MultiFab& weights);
    DistributionMapping SFCBalance(const BoxArray& ba, const MultiFab& weights);
    DistributionMapping GroupedRankBalance(const BoxArray& ba, const MultiFab& weights);
    
    // Utility functions
    void PrintStats() const;
    double GetImbalance() const;
    
private:
    Vector<ComputeNodeInfo> nodes_;
    std::map<int, double> performance_ratios_;  // Maps node_id to performance ratio
    Vector<Vector<double>> timing_ratios_;  // rij matrix
    Vector<TimingRatio> ordered_ratios_;    // sorted rij values
    
    // Helper functions
    double AdjustWeight(double original_weight, const ComputeNodeInfo& node) const;
    Vector<Box> SortBoxesByWeight(const BoxArray& ba, const MultiFab& weights) const;
    void UpdatePerformanceRatios();
    void UpdateTimingRatios();  // New function for rij
    
    // Convert weights to account for heterogeneous performance
    Vector<Long> AdjustWeightsForNodes(const MultiFab& weights) const;
    
    // Performance monitoring
    void TrackPerformance(const Vector<double>& timings);
    void OptimizeDistribution(Vector<int>& assignment);
    
    // New methods for rij-based load balancing
    Vector<int> GetOrderedGroupedRanks() const;
    Vector<Long> GetOrderedWeights(const MultiFab& weights) const;
    void UpdateOrderedRatios();
    
    // Load balancing strategies
    DistributionMapping RijBalance(const BoxArray& ba,
                                 const MultiFab& weights);
                                        
    // Helper methods
    double ComputeLoadImbalance(const Vector<Long>& loads) const;
    void PrintDistributionStats(const std::string& method,
                              const Vector<int>& mapping,
                              const Vector<Long>& loads) const;
};

}
#endif 