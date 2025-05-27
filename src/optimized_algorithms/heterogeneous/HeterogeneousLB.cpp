#include "HeterogeneousLB.H"
#include <AMReX_Print.H>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace amrex {

HeterogeneousLB::HeterogeneousLB() {}

void HeterogeneousLB::InitializeNodes(const Vector<ComputeNodeInfo>& node_info) {
    nodes_ = node_info;
    if (nodes_.empty()) {
        amrex::Abort("HeterogeneousLB::InitializeNodes: No nodes provided");
    }
    UpdatePerformanceRatios();
    UpdateTimingRatios();  // Initialize timing ratios
}

void HeterogeneousLB::UpdatePerformanceRatios() {
    BL_PROFILE("HeterogeneousLB::UpdatePerformanceRatios");
    if (nodes_.empty()) return;
    
    // Find baseline performance (minimum performance factor)
    auto min_node = std::min_element(nodes_.begin(), nodes_.end(),
        [](const ComputeNodeInfo& a, const ComputeNodeInfo& b) {
            return a.performance_factor < b.performance_factor;
        });
    
    double baseline = min_node->performance_factor;
    
    if (baseline <= 0.0) {
        amrex::Abort("HeterogeneousLB::UpdatePerformanceRatios: Invalid baseline performance factor");
    }
    
    // Calculate ratios relative to baseline
    performance_ratios_.clear();  // Clear existing ratios
    
    double max_ratio = 0.0;
    double total_ratio = 0.0;
    
    for (const auto& node : nodes_) {
        double ratio = node.performance_factor / baseline;
        
        // Validate ratio
        if (ratio <= 0.0) {
            amrex::Abort("HeterogeneousLB::UpdatePerformanceRatios: Invalid ratio for node " + 
                        std::to_string(node.node_id));
        }
        
        performance_ratios_[node.node_id] = ratio;
        max_ratio = std::max(max_ratio, ratio);
        total_ratio += ratio;
    }
    
    // Validate overall distribution
    double avg_ratio = total_ratio / nodes_.size();
    double imbalance = (max_ratio / avg_ratio) - 1.0;
    
    if (amrex::Verbose() > 1) {
        amrex::Print() << "\nPerformance Ratio Statistics:\n"
                       << "  Baseline performance: " << baseline << "\n"
                       << "  Average ratio: " << avg_ratio << "\n"
                       << "  Maximum ratio: " << max_ratio << "\n"
                       << "  Initial imbalance: " << (imbalance * 100.0) << "%\n";
        
        for (const auto& ratio : performance_ratios_) {
            amrex::Print() << "  Node " << ratio.first << " ratio: " << ratio.second << "\n";
        }
    }
    
    // Warn if imbalance is too high
    if (imbalance > 0.5) {  // More than 50% imbalance
        amrex::Print() << "\nWarning: High initial performance imbalance detected (" 
                       << (imbalance * 100.0) << "%).\n"
                       << "This may impact load balancing effectiveness.\n";
    }
}

void HeterogeneousLB::UpdateTimingRatios() {
    BL_PROFILE("HeterogeneousLB::UpdateTimingRatios");
    
    if (nodes_.empty()) {
        amrex::Print() << "Warning: No nodes available for timing ratio update\n";
        return;
    }
    
    const int n_nodes = nodes_.size();
    
    // Clear and resize the timing ratios matrix
    timing_ratios_.clear();
    timing_ratios_.resize(n_nodes, Vector<double>(n_nodes, 1.0));
    
    // Calculate timing ratios only if we have valid performance factors
    bool has_valid_factors = true;
    for (const auto& node : nodes_) {
        if (node.performance_factor <= 0.0) {
            has_valid_factors = false;
            break;
        }
    }
    
    if (!has_valid_factors) {
        amrex::Print() << "Warning: Invalid performance factors detected, using default ratios\n";
        return;
    }
    
    // Safely compute timing ratios
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            if (i == j) {
                timing_ratios_[i][j] = 1.0;  // Self-ratio is always 1
            } else {
                // Safely compute rij = ti/tj = pj/pi
                double pi = nodes_[i].performance_factor;
                double pj = nodes_[j].performance_factor;
                
                if (pi > 0.0 && pj > 0.0) {
                    timing_ratios_[i][j] = pj / pi;
                } else {
                    timing_ratios_[i][j] = 1.0;  // Default to 1 if invalid
                }
            }
        }
    }
    
    // Validate the timing ratios
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            if (timing_ratios_[i][j] <= 0.0 || !std::isfinite(timing_ratios_[i][j])) {
                amrex::Print() << "Warning: Invalid timing ratio detected at [" << i << "][" << j 
                              << "] = " << timing_ratios_[i][j] << ", setting to 1.0\n";
                timing_ratios_[i][j] = 1.0;
            }
        }
    }
    
    if (amrex::Verbose() > 1) {
        amrex::Print() << "\nTiming Ratio Matrix (rij):\n";
        for (int i = 0; i < n_nodes; ++i) {
            amrex::Print() << "Node " << nodes_[i].node_id << " (" 
                          << nodes_[i].node_type << "): ";
            for (int j = 0; j < n_nodes; ++j) {
                amrex::Print() << std::fixed << std::setprecision(3) 
                              << timing_ratios_[i][j] << " ";
            }
            amrex::Print() << "\n";
        }
    }
}

void HeterogeneousLB::UpdatePerformanceMetrics(const Vector<double>& timings) {
    if (timings.size() != nodes_.size()) return;
    
    // Update performance metrics for each node
    for (size_t i = 0; i < nodes_.size(); ++i) {
        // Update recent performance with current timing
        nodes_[i].recent_performance = timings[i];
        
        // Update performance factor (inverse of timing)
        if (timings[i] > 0.0) {
            nodes_[i].performance_factor = 1.0 / timings[i];
        }
    }
    
    // Update both performance ratios and timing ratios
    UpdatePerformanceRatios();
    UpdateTimingRatios();
}

Vector<Long> HeterogeneousLB::AdjustWeightsForNodes(const MultiFab& weights) const {
    BL_PROFILE("HeterogeneousLB::AdjustWeightsForNodes");
    const int nboxes = weights.boxArray().size();
    Vector<Long> adjusted_weights(nboxes);
    
    if (nboxes == 0) {
        amrex::Abort("HeterogeneousLB::AdjustWeightsForNodes: Empty BoxArray");
    }
    
    if (nodes_.empty() || timing_ratios_.empty()) {
        // Fallback to original weights if no timing data
        for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
            const int i = mfi.index();
            const Box& box = mfi.validbox();
            try {
                double orig_weight = weights[mfi].sum<RunOn::Host>(box, 0, 1);
                adjusted_weights[i] = static_cast<Long>(std::max(1.0, std::min(orig_weight, 
                    static_cast<double>(std::numeric_limits<Long>::max() / 2))));
            } catch (const std::exception& e) {
                amrex::Abort("Error computing box weight: " + std::string(e.what()));
            }
        }
        return adjusted_weights;
    }
    
    // Get the distribution mapping
    const DistributionMapping& dmap = weights.DistributionMap();
    
    // Track statistics
    double total_orig_weight = 0.0;
    double total_adjusted_weight = 0.0;
    
    // Safely adjust weights using timing ratios
    for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
        const int i = mfi.index();
        const int node_id = dmap[i] % nodes_.size();
        const Box& box = mfi.validbox();
        
        try {
            // Get original weight
            double orig_weight = weights[mfi].sum<RunOn::Host>(box, 0, 1);
            total_orig_weight += orig_weight;
            
            // Find maximum timing ratio for this node
            double max_ratio = 1.0;  // Default to 1.0
            if (node_id < static_cast<int>(timing_ratios_[node_id].size())) {
                for (size_t j = 0; j < timing_ratios_[node_id].size(); ++j) {
                    max_ratio = std::max(max_ratio, timing_ratios_[node_id][j]);
                }
            }
            
            // Scale weight safely
            double scaled_weight = orig_weight * max_ratio;
            scaled_weight = std::max(1.0, std::min(scaled_weight, 
                static_cast<double>(std::numeric_limits<Long>::max() / 2)));
            
            adjusted_weights[i] = static_cast<Long>(scaled_weight);
            total_adjusted_weight += adjusted_weights[i];
            
        } catch (const std::exception& e) {
            amrex::Abort("Error computing adjusted weight: " + std::string(e.what()));
        }
    }
    
    if (amrex::Verbose() > 1) {
        amrex::Print() << "Weight adjustment statistics:\n"
                       << "  Total original weight: " << total_orig_weight << "\n"
                       << "  Total adjusted weight: " << total_adjusted_weight << "\n";
    }
    
    return adjusted_weights;
}

DistributionMapping HeterogeneousLB::BalanceLoad(const BoxArray& ba, 
                                                const MultiFab& weights,
                                                const std::string& strategy) {
    if (strategy == "rij") {
        return RijBalance(ba, weights);
    } else if (strategy == "knapsack") {
        return KnapsackBalance(ba, weights);
    } else if (strategy == "sfc") {
        return SFCBalance(ba, weights);
    } else if (strategy == "grouped") {
        return GroupedRankBalance(ba, weights);
    }
    
    // Default to rij-based balancing
    return RijBalance(ba, weights);
}

DistributionMapping HeterogeneousLB::KnapsackBalance(const BoxArray& ba, const MultiFab& weights) {
    BL_PROFILE("HeterogeneousLB::KnapsackBalance");
    const int nboxes = ba.size();
    
    // Convert weights accounting for node performance
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    
    // Calculate total weight and normalize to prevent overflow
    Long total_weight = 0;
    for (const auto& w : adjusted_weights) {
        total_weight += w;
    }
    
    // Normalize weights to range [1, 10000] to prevent overflow while maintaining ratios
    Vector<Long> normalized_weights(nboxes);
    const double scale_factor = 10000.0 / static_cast<double>(total_weight);
    for (int i = 0; i < nboxes; ++i) {
        normalized_weights[i] = std::max(1L, static_cast<Long>(adjusted_weights[i] * scale_factor));
    }
    
    // Call existing knapsack implementation with normalized weights
    Real efficiency;
    auto std_proc_mapping = KnapSackDoIt(std::vector<Long>(normalized_weights.begin(), normalized_weights.end()),
                                       nodes_.size(),
                                       efficiency,
                                       true,  // do_full_knapsack
                                       std::numeric_limits<int>::max(),
                                       false, // verbose
                                       true); // sort
    
    // Convert std::vector to amrex::Vector
    Vector<int> proc_mapping(std_proc_mapping.begin(), std_proc_mapping.end());
    
    // Map processor IDs to node IDs
    for (int i = 0; i < nboxes; ++i) {
        proc_mapping[i] = nodes_[proc_mapping[i]].node_id;
    }
    
    if (amrex::Verbose() > 1) {
        // Print distribution statistics
        std::map<int, Long> weight_per_node;
        for (int i = 0; i < nboxes; ++i) {
            weight_per_node[proc_mapping[i]] += adjusted_weights[i];
        }
        
        amrex::Print() << "\nKnapsack Distribution Statistics:\n";
        for (const auto& node : nodes_) {
            amrex::Print() << "Node " << node.node_id 
                          << " (Factor: " << node.performance_factor 
                          << ") -> Weight: " << weight_per_node[node.node_id] << "\n";
        }
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

DistributionMapping HeterogeneousLB::SFCBalance(const BoxArray& ba, const MultiFab& weights) {
    BL_PROFILE("HeterogeneousLB::SFCBalance");
    amrex::Print() << "\nEntering SFCBalance...\n";
    
    const int nboxes = ba.size();
    if (nboxes == 0) {
        amrex::Abort("HeterogeneousLB::SFCBalance: Empty BoxArray");
    }
    
    // Convert weights accounting for node performance
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    
    // Create a vector to store processor assignments
    Vector<int> proc_mapping(nboxes, -1);  // Initialize with -1 to detect unassigned boxes
    
    try {
        // Calculate total available processors and their ranges
        int total_procs = 0;
        Vector<std::pair<int, int>> proc_ranges;  // (start, end) ranges for each node
        
        // First pass: calculate normalized processor counts
        double total_perf = 0.0;
        for (const auto& node : nodes_) {
            total_perf += node.performance_factor;
        }
        
        // Second pass: assign processor ranges
        int current_proc = 0;
        for (const auto& node : nodes_) {
            int node_procs = std::max(1, static_cast<int>((node.performance_factor / total_perf) * nboxes));
            proc_ranges.push_back({current_proc, current_proc + node_procs - 1});
            current_proc += node_procs;
            total_procs += node_procs;
        }
        
        // Ensure we have enough processors
        if (total_procs < nboxes) {
            total_procs = nboxes;
            // Adjust the last range to accommodate remaining boxes
            if (!proc_ranges.empty()) {
                proc_ranges.back().second = nboxes - 1;
            }
        }
        
        // Create SFC ordering of boxes
        Vector<int> sfc_order(nboxes);
        std::iota(sfc_order.begin(), sfc_order.end(), 0);
        
        // Sort boxes by their weights
        std::sort(sfc_order.begin(), sfc_order.end(),
                 [&adjusted_weights](int a, int b) {
                     return adjusted_weights[a] > adjusted_weights[b];
                 });
        
        // Assign boxes to processors based on SFC order and node ranges
        for (int i = 0; i < nboxes; ++i) {
            const int box_idx = sfc_order[i];
            const int target_proc = i % total_procs;
            
            // Find which node range contains this processor
            for (size_t n = 0; n < nodes_.size(); ++n) {
                if (target_proc >= proc_ranges[n].first && target_proc <= proc_ranges[n].second) {
                    proc_mapping[box_idx] = nodes_[n].node_id;
                    break;
                }
            }
            
            // Fallback if no valid mapping found
            if (proc_mapping[box_idx] == -1) {
                proc_mapping[box_idx] = target_proc % nodes_.size();
                amrex::Print() << "Warning: Using fallback mapping for box " << box_idx << "\n";
            }
        }
        
        // Verify all boxes were assigned
        for (int i = 0; i < nboxes; ++i) {
            if (proc_mapping[i] == -1) {
                amrex::Abort("SFCBalance: Box " + std::to_string(i) + " was not assigned to any processor");
            }
        }
        
    } catch (const std::exception& e) {
        amrex::Abort("Error in SFCBalance: " + std::string(e.what()));
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

DistributionMapping HeterogeneousLB::GroupedRankBalance(const BoxArray& ba, const MultiFab& weights) {
    // For grouped ranking, we'll use knapsack within each node type group
    std::map<std::string, Vector<ComputeNodeInfo>> node_groups;
    for (const auto& node : nodes_) {
        node_groups[node.node_type].push_back(node);
    }
    
    const int nboxes = ba.size();
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    Vector<int> proc_mapping(nboxes);
    
    // Calculate total capacity for distribution
    double total_capacity = 0.0;
    for (const auto& node : nodes_) {
        total_capacity += node.performance_factor;
    }
    
    // Distribute boxes among groups based on their total capacity
    Long current_box = 0;
    for (auto& group : node_groups) {
        double group_capacity = std::accumulate(group.second.begin(), group.second.end(), 0.0,
            [](double sum, const ComputeNodeInfo& node) {
                return sum + node.performance_factor;
            });
            
        // Calculate number of boxes for this group
        Long group_box_count = static_cast<Long>(
            (group_capacity / total_capacity) * nboxes);
            
        if (group_box_count > 0) {
            // Create subset of weights for this group
            std::vector<Long> group_weights;
            Vector<int> box_indices;
            for (Long i = current_box; i < std::min<Long>(current_box + group_box_count, adjusted_weights.size()); ++i) {
                group_weights.push_back(adjusted_weights[i]);
                box_indices.push_back(i);
            }
            
            // Balance within group using knapsack
            Real efficiency;
            auto std_group_mapping = KnapSackDoIt(group_weights,
                                               group.second.size(),
                                               efficiency,
                                               true,
                                               std::numeric_limits<int>::max(),
                                               false,
                                               true);
            
            // Convert std::vector to amrex::Vector
            Vector<int> group_mapping(std_group_mapping.begin(), std_group_mapping.end());
            
            // Map assignments back to global indices
            for (size_t i = 0; i < group_mapping.size(); ++i) {
                proc_mapping[box_indices[i]] = group.second[group_mapping[i]].node_id;
            }
            
            current_box += group_box_count;
        }
    }
    
    // Assign any remaining boxes
    if (current_box < nboxes) {
        Real efficiency;
        std::vector<Long> remaining_weights(adjusted_weights.begin() + current_box, adjusted_weights.end());
        auto std_remaining_mapping = KnapSackDoIt(remaining_weights,
                                               nodes_.size(),
                                               efficiency,
                                               true,
                                               std::numeric_limits<int>::max(),
                                               false,
                                               true);
                                               
        // Convert std::vector to amrex::Vector
        Vector<int> remaining_mapping(std_remaining_mapping.begin(), std_remaining_mapping.end());
        
        for (size_t i = 0; i < remaining_mapping.size(); ++i) {
            proc_mapping[current_box + i] = nodes_[remaining_mapping[i]].node_id;
        }
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

double HeterogeneousLB::GetImbalance() const {
    if (nodes_.empty()) return 0.0;
    
    Vector<double> adjusted_loads(nodes_.size(), 0.0);
    double total_load = 0.0;
    
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i].recent_performance > 0.0) {
            auto ratio_it = performance_ratios_.find(nodes_[i].node_id);
            if (ratio_it != performance_ratios_.end()) {
                adjusted_loads[i] = nodes_[i].recent_performance * ratio_it->second;
                total_load += adjusted_loads[i];
            }
        }
    }
    
    if (total_load <= 0.0) return 0.0;
    
    double avg_load = total_load / nodes_.size();
    double max_load = *std::max_element(adjusted_loads.begin(), adjusted_loads.end());
    
    return (max_load / avg_load) - 1.0;
}

void HeterogeneousLB::PrintStats() const {
    amrex::Print() << "\nHeterogeneous Load Balancer Statistics:\n";
    amrex::Print() << "Number of compute nodes: " << nodes_.size() << "\n";
    
    // Print current load imbalance
    if (std::any_of(nodes_.begin(), nodes_.end(),
                    [](const ComputeNodeInfo& node) { return node.recent_performance > 0.0; })) {
        double max_time = 0.0;
        double total_time = 0.0;
        for (const auto& node : nodes_) {
            if (node.recent_performance > 0.0) {
                max_time = std::max(max_time, node.recent_performance);
                total_time += node.recent_performance;
            }
        }
        double avg_time = total_time / nodes_.size();
        double imbalance = ((max_time / avg_time) - 1.0) * 100.0;
        amrex::Print() << "Current load imbalance: " << imbalance << "%\n";
    } else {
        amrex::Print() << "Current load imbalance: No performance data available\n";
    }
    
    // Print node info
    amrex::Print() << "\nNode Performance Breakdown:\n";
    for (const auto& node : nodes_) {
        amrex::Print() << "Node " << node.node_id << " (" << node.node_type << "):\n"
                       << "  Performance Factor: " << node.performance_factor << "\n";
        if (node.recent_performance > 0.0) {
            amrex::Print() << "  Recent Performance: " << node.recent_performance << " s\n";
        }
    }
}

void HeterogeneousLB::UpdateOrderedRatios() {
    ordered_ratios_.clear();
    const int n_nodes = nodes_.size();
    
    // Safety check
    if (timing_ratios_.empty() || timing_ratios_.size() != n_nodes) {
        UpdateTimingRatios();  // Re-initialize if not properly initialized
    }
    
    // Create ordered list of timing ratios
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            if (i != j && i < timing_ratios_.size() && j < timing_ratios_[i].size()) {
                ordered_ratios_.push_back(TimingRatio(i, j, timing_ratios_[i][j]));
            }
        }
    }
    
    // Sort ratios in ascending order
    std::sort(ordered_ratios_.begin(), ordered_ratios_.end());
    
    if (amrex::Verbose() > 1) {
        amrex::Print() << "\nOrdered Timing Ratios (rij):\n";
        for (const auto& ratio : ordered_ratios_) {
            amrex::Print() << "Node " << ratio.i << " -> Node " << ratio.j 
                          << ": " << ratio.rij << "\n";
        }
    }
}

Vector<int> HeterogeneousLB::GetOrderedGroupedRanks() const {
    const int n_nodes = nodes_.size();
    Vector<int> ranks(n_nodes);
    
    // First, group nodes by type
    std::map<std::string, Vector<int>> type_groups;
    for (int i = 0; i < n_nodes; ++i) {
        type_groups[nodes_[i].node_type].push_back(i);
    }
    
    // Sort groups by average performance
    Vector<std::pair<std::string, double>> sorted_groups;
    for (const auto& group : type_groups) {
        double avg_perf = 0.0;
        for (int node_id : group.second) {
            avg_perf += nodes_[node_id].performance_factor;
        }
        avg_perf /= group.second.size();
        sorted_groups.push_back({group.first, avg_perf});
    }
    
    std::sort(sorted_groups.begin(), sorted_groups.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Create ordered rank list
    int current_rank = 0;
    for (const auto& group : sorted_groups) {
        for (int node_id : type_groups[group.first]) {
            ranks[node_id] = current_rank++;
        }
    }
    
    return ranks;
}

Vector<Long> HeterogeneousLB::GetOrderedWeights(const MultiFab& weights) const {
    const int nboxes = weights.boxArray().size();
    Vector<Long> ordered_weights(nboxes);
    
    // Get base weights
    for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
        const int i = mfi.index();
        const Box& box = mfi.validbox();
        ordered_weights[i] = static_cast<Long>(weights[mfi].sum<RunOn::Host>(box, 0, 1));
    }
    
    // Sort weights based on timing ratios
    std::sort(ordered_weights.begin(), ordered_weights.end(), std::greater<Long>());
    
    return ordered_weights;
}

DistributionMapping HeterogeneousLB::RijBalance(const BoxArray& ba, const MultiFab& weights) {
    BL_PROFILE("HeterogeneousLB::RijBalance");
    
    const int nboxes = ba.size();
    if (nboxes == 0) return DistributionMapping();
    
    // Update timing ratios and get ordered lists
    UpdateOrderedRatios();
    Vector<int> ordered_ranks = GetOrderedGroupedRanks();
    Vector<Long> ordered_weights = GetOrderedWeights(weights);
    
    // Create two different mappings
    Vector<int> rank_based_mapping(nboxes);
    Vector<int> weight_based_mapping(nboxes);
    
    // Method 1: Assign based on ordered ranks
    for (int i = 0; i < nboxes; ++i) {
        int rank = ordered_ranks[i % nodes_.size()];
        rank_based_mapping[i] = nodes_[rank].node_id;  // Map rank to actual node ID
    }
    
    // Method 2: Assign based on ordered weights
    Vector<Long> node_loads(nodes_.size(), 0);
    for (int i = 0; i < nboxes; ++i) {
        // Find least loaded node considering timing ratios
        int target_node = 0;
        double min_relative_load = std::numeric_limits<double>::max();
        
        for (int j = 0; j < static_cast<int>(nodes_.size()); ++j) {
            // Calculate relative load by considering both source and target ratios
            double relative_load = static_cast<double>(node_loads[j]);
            for (const auto& ratio : ordered_ratios_) {
                if (ratio.i == j) {
                    // This node is source, multiply by rij
                    relative_load *= ratio.rij;
                } else if (ratio.j == j) {
                    // This node is target, divide by rij
                    relative_load /= ratio.rij;
                }
            }
            
            if (relative_load < min_relative_load) {
                min_relative_load = relative_load;
                target_node = j;
            }
        }
        
        weight_based_mapping[i] = nodes_[target_node].node_id;  // Map to actual node ID
        node_loads[target_node] += ordered_weights[i];
    }
    
    // Compare the two methods
    Vector<Long> rank_loads(nodes_.size(), 0);
    Vector<Long> weight_loads(nodes_.size(), 0);
    
    for (int i = 0; i < nboxes; ++i) {
        rank_loads[rank_based_mapping[i]] += ordered_weights[i];
        weight_loads[weight_based_mapping[i]] += ordered_weights[i];
    }
    
    double rank_imbalance = ComputeLoadImbalance(rank_loads);
    double weight_imbalance = ComputeLoadImbalance(weight_loads);
    
    if (amrex::Verbose() > 0) {
        PrintDistributionStats("Rank-based", rank_based_mapping, rank_loads);
        PrintDistributionStats("Weight-based", weight_based_mapping, weight_loads);
    }
    
    // Choose the better method
    if (rank_imbalance <= weight_imbalance) {
        return DistributionMapping(std::move(rank_based_mapping));
    } else {
        return DistributionMapping(std::move(weight_based_mapping));
    }
}

double HeterogeneousLB::ComputeLoadImbalance(const Vector<Long>& loads) const {
    if (loads.empty()) return 0.0;
    
    Long total_load = 0;
    Long max_load = 0;
    for (Long load : loads) {
        total_load += load;
        max_load = std::max(max_load, load);
    }
    
    double avg_load = static_cast<double>(total_load) / loads.size();
    return (max_load / avg_load) - 1.0;
}

void HeterogeneousLB::PrintDistributionStats(const std::string& method,
                                           const Vector<int>& mapping,
                                           const Vector<Long>& loads) const {
    amrex::Print() << "\n" << method << " Distribution Statistics:\n";
    amrex::Print() << std::setw(15) << "Node ID" 
                   << std::setw(15) << "Load"
                   << std::setw(15) << "Box Count\n";
    amrex::Print() << std::string(45, '-') << "\n";
    
    Vector<int> box_counts(nodes_.size(), 0);
    for (int i = 0; i < mapping.size(); ++i) {
        box_counts[mapping[i]]++;
    }
    
    for (size_t i = 0; i < nodes_.size(); ++i) {
        amrex::Print() << std::setw(15) << nodes_[i].node_id
                       << std::setw(15) << loads[i]
                       << std::setw(15) << box_counts[i] << "\n";
    }
    
    amrex::Print() << "Load imbalance: " 
                   << (ComputeLoadImbalance(loads) * 100.0) << "%\n";
}

} // namespace amrex 