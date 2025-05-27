#ifndef TOPOLOGY_AWARE_H
#define TOPOLOGY_AWARE_H

#include <vector>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BoxArray.H>

class TopologyAware {
public:
    // Main load balancing function
    static std::vector<int> balance(const amrex::BoxArray& ba,
        const std::vector<amrex::Long>& wgts,
        int nnodes,
        int ranks_per_node,
        amrex::Real* efficiency = nullptr);

    // Measurement functions for evaluating the load balancer
    static double calculateLoadBalanceMetric(const std::vector<amrex::Long>& wgts,
        const std::vector<int>& dmap,
        int nnodes,
        int ranks_per_node);
                                        
    static double calculateExternalCommRatio(const amrex::BoxArray& ba,
        const std::vector<int>& dmap,
        int nnodes,
        int ranks_per_node);

private:
    // Helper methods
    static void buildCommunicationGraph(const amrex::BoxArray& ba,
        std::vector<std::vector<double>>& comm_graph);

    static void buildTopologyMatrix(int nnodes,
        int ranks_per_node,
        std::vector<std::vector<double>>& distance_matrix);
        
    static int getRankNode(int rank, int ranks_per_node);
    
    static double getBoxPairCommCost(const amrex::Box& box1, 
        const amrex::Box& box2);
};

#endif

// #ifndef TOPOLOGY_AWARE_H
// #define TOPOLOGY_AWARE_H

// #include <vector>
// #include <AMReX_Utility.H>
// #include <AMReX_ParallelDescriptor.H>
// #include <AMReX_BoxArray.H>

// class TopologyAware {
// public:
//     static std::vector<int> balance(const amrex::BoxArray& ba,
//         const std::vector<amrex::Long>& wgts,
//         int nnodes,
//         int ranks_per_node,
//         amrex::Real* efficiency = nullptr);

// private:
//     static void buildCommunicationGraph(const amrex::BoxArray& ba,
//         std::vector<std::vector<double>>& comm_graph);

//     static void detectTopology(std::vector<std::vector<double>>& distance_matrix);
// };

// #endif