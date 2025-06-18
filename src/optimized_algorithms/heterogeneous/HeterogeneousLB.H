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
#include <cmath>

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

// Structure for load balance quality metrics
struct LoadBalanceMetrics {
    double makespan;           // Maximum completion time
    double efficiency;         // Total work / (makespan * num_nodes)
    double load_imbalance;     // (max_load - avg_load) / avg_load
    double total_work;         // Sum of all work
    
    LoadBalanceMetrics() : makespan(0.0), efficiency(0.0), load_imbalance(0.0), total_work(0.0) {}
};

// Enhanced structures for RIJ optimization
struct WeightedBox {
    int box_id;
    Long weight;
    Box box;
    IntVect center;
    double priority;
};

struct NodeLoadInfo {
    int node_idx;
    double current_load;
    int box_count;
    double performance_factor;
    Vector<int> recent_boxes;
};

class HeterogeneousLB {
public:
    // Constructor
    HeterogeneousLB();
    
    // Initialize node information
    void InitializeNodes(const Vector<ComputeNodeInfo>& node_info);
    
    // Update performance metrics
    void UpdatePerformanceMetrics(const Vector<double>& timings);
    void UpdateSingleNodeMetrics(int node_id, double timing);
    
    // Main load balancing methods
    DistributionMapping BalanceLoad(const BoxArray& ba, 
                                  const MultiFab& weights,
                                  const std::string& strategy = "rij");
    
    // Different load balancing strategies
    DistributionMapping KnapsackBalance(const BoxArray& ba, const MultiFab& weights);
    DistributionMapping SFCBalance(const BoxArray& ba, const MultiFab& weights);
    DistributionMapping GroupedRankBalance(const BoxArray& ba, const MultiFab& weights);
    DistributionMapping RijBalance(const BoxArray& ba, const MultiFab& weights);
    
    // Utility functions
    void PrintStats() const;
    double GetImbalance() const;
    LoadBalanceMetrics EvaluateBalance(const DistributionMapping& dm, 
                                     const MultiFab& weights) const;
    
private:
    Vector<ComputeNodeInfo> nodes_;
    std::map<int, int> node_id_to_index_;        // Maps node_id to array index
    std::map<int, double> performance_ratios_;   // Maps node_id to performance ratio
    Vector<Vector<double>> timing_ratios_;       // rij matrix [i][j]
    Vector<TimingRatio> ordered_ratios_;         // sorted rij values
    
    // Enhanced RIJ members
    Vector<IntVect> box_centers_;
    double spatial_locality_threshold_ = 10.0;
    
    // Helper functions
    void UpdatePerformanceRatios();
    void UpdateTimingRatios();
    void UpdateTimingRatiosForNode(int node_index);
    void UpdateOrderedRatios();
    
    // Safe node access
    int GetNodeIndex(int node_id) const;
    bool IsValidNodeIndex(int index) const;
    
    // Weight processing
    Vector<Long> AdjustWeightsForNodes(const MultiFab& weights) const;
    Vector<Long> NormalizeWeights(const Vector<Long>& weights) const;
    Vector<Long> GetOrderedWeights(const MultiFab& weights) const;
    
    // Spatial locality helpers
    uint64_t ComputeHilbertIndex(const IntVect& iv) const;
    Vector<int> CreateSpatialOrdering(const BoxArray& ba) const;
    
    // Load balancing quality
    double ComputeLoadImbalance(const Vector<Long>& loads) const;
    void PrintDistributionStats(const std::string& method,
                              const Vector<int>& mapping,
                              const Vector<Long>& loads) const;
    
    // Grouped ranking helpers
    Vector<int> GetOrderedGroupedRanks() const;
    std::map<std::string, Vector<int>> GroupNodesByType() const;
    
    // Enhanced RIJ methods
    Vector<int> EnhancedRijOptimization(const BoxArray& ba, const Vector<Long>& weights);
    Vector<int> RijAwareInitialAssignment(const BoxArray& ba, const Vector<Long>& weights);
    Vector<int> RijLocalOptimization(const Vector<int>& initial_assignment, const Vector<Long>& weights);
    Vector<int> RijGlobalRefinement(const Vector<int>& assignment, const Vector<Long>& weights);
    
    // Enhanced RIJ helper methods
    int SelectBestNodeForBox(const WeightedBox& wb, Vector<NodeLoadInfo>& node_loads, int total_boxes);
    double CalculateNodeScore(const WeightedBox& wb, const NodeLoadInfo& node_info, int node_idx, int total_boxes);
    double CalculateBoxPriority(Long weight, const IntVect& center);
    double CalculateBoxDistance(const IntVect& center1, const IntVect& center2);
    double CalculateSpatialLocalityBonus(const WeightedBox& wb, const NodeLoadInfo& node_info);
    
    int FindBetterNodeUsingRij(int box_id, int current_node_idx, Long box_weight, const Vector<int>& assignment);
    bool WouldImproveMakespan(int box_id, int from_node, int to_node, const Vector<int>& assignment, const Vector<Long>& weights);
    Vector<double> CalculateCurrentLoads(const Vector<int>& assignment, const Vector<Long>& weights);
    double CalculateLoadBalanceImprovement(int box_id, int from_node, int to_node, const Vector<double>& node_loads, Long box_weight);
    
    Vector<int> RebalanceOverloadedNodes(const Vector<int>& assignment, const Vector<Long>& weights, const Vector<double>& node_loads);
    int FindBestTargetNode(int source_node, const Vector<int>& candidate_nodes, Long box_weight);
};

}
#endif