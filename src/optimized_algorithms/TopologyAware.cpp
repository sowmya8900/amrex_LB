#include "TopologyAware.H"
#include <algorithm>
#include <numeric>
#include <map>
#include <limits>
#include <cmath>
#include "Util.H"

using namespace amrex;

// Helper function to determine which node a rank belongs to
int TopologyAware::getRankNode(int rank, int ranks_per_node) {
    return rank / ranks_per_node;
}

// Calculate communication cost between two boxes
double TopologyAware::getBoxPairCommCost(const amrex::Box& box1, const amrex::Box& box2) {
    const Box intersect = box1 & box2;
    if (intersect.isEmpty()) return 0.0;
    
    const IntVect size = intersect.size();
    double cost = 0.0;
    
    // Communication cost proportional to touching face area
    if (size[0] == 1) { // x-face
        cost = size[1] * (AMREX_SPACEDIM > 2 ? size[2] : 1.0);
    } else if (size[1] == 1) { // y-face
        cost = size[0] * (AMREX_SPACEDIM > 2 ? size[2] : 1.0);
    } else if (AMREX_SPACEDIM > 2 && size[2] == 1) { // z-face
        cost = size[0] * size[1];
    }
    
    return cost;
}

void TopologyAware::buildCommunicationGraph(const BoxArray& ba,
    std::vector<std::vector<double>>& comm_graph) {
    BL_PROFILE("TopologyAware::buildCommunicationGraph");

    const int nboxes = ba.size();
    comm_graph.resize(nboxes, std::vector<double>(nboxes, 0.0));

    for (int i = 0; i < nboxes; ++i) {
        Box bx = ba[i]; // Create non-const copy
        bx.grow(1); // Grow the box to find neighbors
        
        // Get intersecting boxes
        auto neighbors = ba.intersections(bx);

        for (const auto& neighbor_pair : neighbors) {
            int j = neighbor_pair.first;
            const Box& neighbor = neighbor_pair.second;
            
            if (i == j) continue; // Skip self

            double cost = getBoxPairCommCost(bx, neighbor);
            
            if (cost > 0) {
                comm_graph[i][j] = cost;
                comm_graph[j][i] = cost; // Symmetric graph
            }
        }
    }
}

void TopologyAware::buildTopologyMatrix(int nnodes,
    int ranks_per_node,
    std::vector<std::vector<double>>& distance_matrix) {
    BL_PROFILE("TopologyAware::buildTopologyMatrix");

    int nranks = nnodes * ranks_per_node;
    distance_matrix.resize(nranks, std::vector<double>(nranks, 0.0));

    // Build topology distance matrix
    // 0 for same rank, 1 for same node, 10 for different nodes
    for (int i = 0; i < nranks; ++i) {
        for (int j = 0; j < nranks; ++j) {
            if (i == j) {
                distance_matrix[i][j] = 0.0; // Same rank
            } else {
                bool same_node = (getRankNode(i, ranks_per_node) == getRankNode(j, ranks_per_node));
                distance_matrix[i][j] = same_node ? 1.0 : 10.0; // Same node vs different node
            }
        }
    }
}

std::vector<int> TopologyAware::balance(const BoxArray& ba,
    const std::vector<Long>& wgts,
    int nnodes,
    int ranks_per_node,
    Real* efficiency) {
    BL_PROFILE("TopologyAware::balance");

    const int nranks = nnodes * ranks_per_node;
    const int nboxes = ba.size();

    // Guard against invalid inputs
    if (nboxes == 0 || nranks <= 0 || nnodes <= 0 || ranks_per_node <= 0) {
        if (efficiency) *efficiency = 0.0;
        return std::vector<int>(nboxes, 0); // Return default mapping
    }

    // 1. Build communication graph between boxes
    std::vector<std::vector<double>> comm_graph;
    buildCommunicationGraph(ba, comm_graph);

    // 2. Build topology distance matrix between ranks
    std::vector<std::vector<double>> distance_matrix;
    buildTopologyMatrix(nnodes, ranks_per_node, distance_matrix);

    // 3. Calculate total communication volume for each box (for sorting)
    std::vector<std::pair<double, int>> box_comm(nboxes);
    for (int i = 0; i < nboxes; ++i) {
        double total_comm = 0.0;
        for (int j = 0; j < nboxes; ++j) {
            total_comm += comm_graph[i][j];
        }
        box_comm[i] = std::make_pair(total_comm, i);
    }
    
    // Sort boxes by total communication volume (highest first)
    std::sort(box_comm.rbegin(), box_comm.rend());

    // 4. Calculate total weight
    Long total_weight = std::accumulate(wgts.begin(), wgts.end(), 0L);
    
    // Target weight per node
    double target_node_weight = static_cast<double>(total_weight) / nnodes;
    
    // Target weight per rank
    double target_rank_weight = static_cast<double>(total_weight) / nranks;
    
    // Allow imbalance up to 20%
    double imbalance_factor = 1.2;

    // Track weight assignments
    std::vector<Long> node_weights(nnodes, 0);
    std::vector<Long> rank_weights(nranks, 0);
    
    // Result distribution map
    std::vector<int> dmap(nboxes, -1);

    // 5. Greedy assignment - prioritize boxes with high communication
    for (const auto& box_pair : box_comm) {
        int box_id = box_pair.second;
        Long box_weight = wgts[box_id];
        
        // Find best placement considering both load balance and communication
        int best_rank = -1;
        double best_score = std::numeric_limits<double>::max();

        for (int rank = 0; rank < nranks; ++rank) {
            int node = getRankNode(rank, ranks_per_node);
            
            // Skip overloaded ranks and nodes
            if (rank_weights[rank] + box_weight > imbalance_factor * target_rank_weight ||
                node_weights[node] + box_weight > imbalance_factor * target_node_weight) {
                continue;
            }
            
            // Calculate communication cost if placed on this rank
            double comm_cost = 0.0;
            for (int j = 0; j < nboxes; ++j) {
                if (dmap[j] != -1 && comm_graph[box_id][j] > 0) {
                    // Use topology distance to weight communication cost
                    comm_cost += comm_graph[box_id][j] * distance_matrix[rank][dmap[j]];
                }
            }
            
            // Score combines communication cost and load balance
            // Load balance component increases as the rank gets more loaded
            double load_balance_penalty = (rank_weights[rank] / target_rank_weight) * 0.5;
            double score = comm_cost * (1.0 + load_balance_penalty);
            
            if (score < best_score) {
                best_score = score;
                best_rank = rank;
            }
        }
        
        // If no valid rank found (all are overloaded), find least loaded rank
        if (best_rank == -1) {
            Long min_weight = std::numeric_limits<Long>::max();
            for (int rank = 0; rank < nranks; ++rank) {
                if (rank_weights[rank] < min_weight) {
                    min_weight = rank_weights[rank];
                    best_rank = rank;
                }
            }
        }
        
        // Assign box to selected rank
        dmap[box_id] = best_rank;
        int node = getRankNode(best_rank, ranks_per_node);
        rank_weights[best_rank] += box_weight;
        node_weights[node] += box_weight;
    }

    // Calculate efficiency if requested
    if (efficiency) {
        *efficiency = computeEfficiency(wgts, dmap, nranks);
    }

    return dmap;
}

// Calculate load balance metric: avg load per node / max load per node (ideally 1)
double TopologyAware::calculateLoadBalanceMetric(const std::vector<Long>& wgts,
    const std::vector<int>& dmap,
    int nnodes,
    int ranks_per_node) {
    
    if (wgts.empty() || dmap.empty() || nnodes <= 0) return 0.0;
    
    const int nboxes = wgts.size();
    std::vector<Long> node_weights(nnodes, 0);
    Long total_weight = 0;
    
    // Sum weights per node
    for (int i = 0; i < nboxes; ++i) {
        int rank = dmap[i];
        if (rank >= 0) {
            int node = getRankNode(rank, ranks_per_node);
            if (node < nnodes) {
                node_weights[node] += wgts[i];
                total_weight += wgts[i];
            }
        }
    }
    
    // Find max node weight
    Long max_node_weight = *std::max_element(node_weights.begin(), node_weights.end());
    
    // Calculate average node weight
    double avg_node_weight = static_cast<double>(total_weight) / nnodes;
    
    // Return load balance metric
    return (max_node_weight > 0) ? avg_node_weight / max_node_weight : 0.0;
}

// Calculate external communication ratio: external comm / total comm (ideally 0)
double TopologyAware::calculateExternalCommRatio(const amrex::BoxArray& ba,
    const std::vector<int>& dmap,
    int nnodes,
    int ranks_per_node) {
    
    const int nboxes = ba.size();
    double total_comm = 0.0;
    double external_comm = 0.0;
    
    for (int i = 0; i < nboxes; ++i) {
        Box bx = ba[i];
        bx.grow(1);
        
        auto neighbors = ba.intersections(bx);
        
        for (const auto& neighbor_pair : neighbors) {
            int j = neighbor_pair.first;
            const Box& neighbor = neighbor_pair.second;
            
            if (i == j) continue;
            
            double comm_cost = getBoxPairCommCost(bx, neighbor);
            
            if (comm_cost > 0) {
                total_comm += comm_cost;
                
                // Check if boxes are assigned to different nodes
                int node_i = getRankNode(dmap[i], ranks_per_node);
                int node_j = getRankNode(dmap[j], ranks_per_node);
                
                if (node_i != node_j) {
                    external_comm += comm_cost;
                }
            }
        }
    }
    
    return (total_comm > 0) ? external_comm / total_comm : 0.0;
}

// #include "TopologyAware.H"
// #include <algorithm>
// #include <numeric>
// #include <map>
// #include "Util.H"

// using namespace amrex;

// void TopologyAware::buildCommunicationGraph(const BoxArray& ba,
//     std::vector<std::vector<double>>& comm_graph) {
//     BL_PROFILE("TopologyAware::buildCommunicationGraph");

//     const int nboxes = ba.size();
//     comm_graph.resize(nboxes, std::vector<double>(nboxes, 0.0));

//     for (int i = 0; i < nboxes; ++i) {
//         Box bx = ba[i]; // Create non-const copy
//         bx.grow(1); // Grow the box first
        
//         // Get intersecting boxes and their indices
//         auto neighbors = ba.intersections(bx);

//         for (const auto& neighbor_pair : neighbors) {
//             int j = neighbor_pair.first;
//             const Box& neighbor = neighbor_pair.second;
            
//             if (i == j) continue; // Skip self

//             const Box intersect = bx & neighbor;
//             const IntVect size = intersect.size();

//             // Communication cost proportional to touching face area
//             double cost = 0.0;
//             if (size[0] == 1) { // x-face
//                 cost = size[1] * size[2];
//             } else if (size[1] == 1) { // y-face
//                 cost = size[0] * size[2];
//             } else if (size[2] == 1) { // z-face
//                 cost = size[0] * size[1];
//             }

//             comm_graph[i][j] = cost;
//             comm_graph[j][i] = cost;
//         }
//     }
// }

// void TopologyAware::detectTopology(std::vector<std::vector<double>>& distance_matrix) {
//     BL_PROFILE("TopologyAware::detectTopology");

//     int nranks = ParallelDescriptor::NProcs();
//     distance_matrix.resize(nranks, std::vector<double>(nranks, 0.0));

//     // Simple topology model - penalize inter-node communication
//     for (int i = 0; i < nranks; ++i) {
//         for (int j = 0; j < nranks; ++j) {
//             if (i == j) {
//                 distance_matrix[i][j] = 0.0;
//             } else {
//                 // Assume ranks on same node are consecutive
//                 bool same_node = (i / ParallelDescriptor::NProcsPerNode()) == 
//                                 (j / ParallelDescriptor::NProcsPerNode());
//                 distance_matrix[i][j] = same_node ? 1.0 : 10.0; // 10x penalty for inter-node
//             }
//         }
//     }
// }

// std::vector<int> TopologyAware::balance(const BoxArray& ba,
//     const std::vector<Long>& wgts,
//     int nnodes,
//     int ranks_per_node,
//     Real* efficiency) {
//     BL_PROFILE("TopologyAware::balance");

//     const int nranks = nnodes * ranks_per_node;
//     const int nboxes = ba.size();

//     // 1. Build communication graph
//     std::vector<std::vector<double>> comm_graph;
//     buildCommunicationGraph(ba, comm_graph);

//     // 2. Get topology distance matrix
//     std::vector<std::vector<double>> distance_matrix;
//     detectTopology(distance_matrix);

//     // 3. Sort boxes by communication volume
//     std::vector<std::pair<double, int>> box_comm(nboxes);
//     for (int i = 0; i < nboxes; ++i) {
//         box_comm[i].first = std::accumulate(comm_graph[i].begin(), comm_graph[i].end(), 0.0);
//         box_comm[i].second = i;
//     }
//     std::sort(box_comm.rbegin(), box_comm.rend()); // Descending

//     // 4. Greedy assignment
//     std::vector<int> dmap(nboxes, -1);
//     std::vector<Long> node_weights(nnodes, 0);
//     std::vector<std::vector<Long>> rank_weights(nnodes, std::vector<Long>(ranks_per_node, 0));

//     for (const auto& [comm, box_id] : box_comm) {
//         int best_node = -1;
//         int best_rank = -1;
//         double best_score = std::numeric_limits<double>::max();
//         Long box_weight = wgts[box_id];

//         // Find best placement considering both load balance and communication
//         for (int node = 0; node < nnodes; ++node) {
//             if (node_weights[node] + box_weight > 1.2 * std::accumulate(wgts.begin(), wgts.end(), 0L) / nnodes) {
//                 continue; // Skip overloaded nodes
//             }

//             for (int rank = 0; rank < ranks_per_node; ++rank) {
//                 int global_rank = node * ranks_per_node + rank;

//                 if (rank_weights[node][rank] + box_weight > 1.2 * std::accumulate(wgts.begin(), wgts.end(), 0L) / nranks) {
//                     continue; // Skip overloaded ranks
//                 }

//                 // Calculate communication cost
//                 double comm_cost = 0.0;
//                 for (int j = 0; j < nboxes; ++j) {
//                     if (dmap[j] != -1 && comm_graph[box_id][j] > 0) {
//                         comm_cost += comm_graph[box_id][j] * distance_matrix[global_rank][dmap[j]];
//                     }
//                 }

//                 // Combine with load balance consideration
//                 double score = comm_cost + 0.1 * rank_weights[node][rank]; // Weight factor

//                 if (score < best_score) {
//                     best_score = score;
//                     best_node = node;
//                     best_rank = rank;
//                 }
//             }
//         }

//         // Assign to best found rank
//         int global_rank = best_node * ranks_per_node + best_rank;
//         dmap[box_id] = global_rank;
//         node_weights[best_node] += box_weight;
//         rank_weights[best_node][best_rank] += box_weight;
//     }

//     // Calculate efficiency if requested
//     if (efficiency) {
//         *efficiency = computeEfficiency(wgts, dmap, nranks);
//     }

//     return dmap;
// }