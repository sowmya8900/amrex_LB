#include "HeterogeneousLB.H"
#include <AMReX_Print.H>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>

namespace amrex {

HeterogeneousLB::HeterogeneousLB() {
    spatial_locality_threshold_ = 10.0;
}

void HeterogeneousLB::InitializeNodes(const Vector<ComputeNodeInfo>& node_info) {
    nodes_ = node_info;
    if (nodes_.empty()) {
        amrex::Abort("HeterogeneousLB::InitializeNodes: No nodes provided");
    }
    
    // Build node_id to index mapping
    node_id_to_index_.clear();
    for (size_t i = 0; i < nodes_.size(); ++i) {
        node_id_to_index_[nodes_[i].node_id] = static_cast<int>(i);
    }
    
    UpdatePerformanceRatios();
    UpdateTimingRatios();
}

int HeterogeneousLB::GetNodeIndex(int node_id) const {
    auto it = node_id_to_index_.find(node_id);
    return (it != node_id_to_index_.end()) ? it->second : -1;
}

bool HeterogeneousLB::IsValidNodeIndex(int index) const {
    return index >= 0 && index < static_cast<int>(nodes_.size());
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
    performance_ratios_.clear();
    
    double max_ratio = 0.0;
    double total_ratio = 0.0;
    
    for (const auto& node : nodes_) {
        double ratio = node.performance_factor / baseline;
        
        if (ratio <= 0.0) {
            amrex::Abort("HeterogeneousLB::UpdatePerformanceRatios: Invalid ratio for node " + 
                        std::to_string(node.node_id));
        }
        
        performance_ratios_[node.node_id] = ratio;
        max_ratio = std::max(max_ratio, ratio);
        total_ratio += ratio;
    }
    
    if (amrex::Verbose() > 1) {
        double avg_ratio = total_ratio / nodes_.size();
        double imbalance = (max_ratio / avg_ratio) - 1.0;
        
        amrex::Print() << "\nPerformance Ratio Statistics:\n"
                       << "  Baseline performance: " << baseline << "\n"
                       << "  Average ratio: " << avg_ratio << "\n"
                       << "  Maximum ratio: " << max_ratio << "\n"
                       << "  Initial imbalance: " << (imbalance * 100.0) << "%\n";
        
        for (const auto& ratio : performance_ratios_) {
            amrex::Print() << "  Node " << ratio.first << " ratio: " << ratio.second << "\n";
        }
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
    
    // Calculate timing ratios: rij = tj/ti = pi/pj
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            if (i == j) {
                timing_ratios_[i][j] = 1.0;  // Self-ratio is always 1
            } else {
                double pi = nodes_[i].performance_factor;
                double pj = nodes_[j].performance_factor;
                
                if (pi > 0.0 && pj > 0.0) {
                    // rij = time_on_j / time_on_i = perf_i / perf_j
                    timing_ratios_[i][j] = pi / pj;
                } else {
                    timing_ratios_[i][j] = 1.0;
                }
            }
        }
    }
    
    // Validate timing ratios
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            if (timing_ratios_[i][j] <= 0.0 || !std::isfinite(timing_ratios_[i][j])) {
                amrex::Print() << "Warning: Invalid timing ratio at [" << i << "][" << j 
                              << "] = " << timing_ratios_[i][j] << ", setting to 1.0\n";
                timing_ratios_[i][j] = 1.0;
            }
        }
    }
    
    if (amrex::Verbose() > 1) {
        amrex::Print() << "\nTiming Ratio Matrix (rij = ti/tj):\n";
        amrex::Print() << std::setw(8) << " ";
        for (int j = 0; j < n_nodes; ++j) {
            amrex::Print() << std::setw(8) << ("N" + std::to_string(nodes_[j].node_id));
        }
        amrex::Print() << "\n";
        
        for (int i = 0; i < n_nodes; ++i) {
            amrex::Print() << std::setw(8) << ("N" + std::to_string(nodes_[i].node_id));
            for (int j = 0; j < n_nodes; ++j) {
                amrex::Print() << std::setw(8) << std::fixed << std::setprecision(2) 
                              << timing_ratios_[i][j];
            }
            amrex::Print() << "\n";
        }
    }
}

void HeterogeneousLB::UpdateTimingRatiosForNode(int node_index) {
    if (!IsValidNodeIndex(node_index) || timing_ratios_.empty()) {
        UpdateTimingRatios();  // Full update if invalid
        return;
    }
    
    const int n_nodes = nodes_.size();
    
    // Update row i (how work from node_index performs on other nodes)
    for (int j = 0; j < n_nodes; ++j) {
        if (node_index == j) {
            timing_ratios_[node_index][j] = 1.0;
        } else {
            double pi = nodes_[node_index].performance_factor;
            double pj = nodes_[j].performance_factor;
            if (pi > 0.0 && pj > 0.0) {
                timing_ratios_[node_index][j] = pi / pj;
            }
        }
    }
    
    // Update column j (how work from other nodes performs on node_index)
    for (int i = 0; i < n_nodes; ++i) {
        if (i == node_index) continue;
        
        double pi = nodes_[i].performance_factor;
        double pj = nodes_[node_index].performance_factor;
        if (pi > 0.0 && pj > 0.0) {
            timing_ratios_[i][node_index] = pi / pj;
        }
    }
}

void HeterogeneousLB::UpdatePerformanceMetrics(const Vector<double>& timings) {
    if (timings.size() != nodes_.size()) {
        amrex::Print() << "Warning: Timing vector size mismatch\n";
        return;
    }
    
    // Update performance metrics for each node
    for (size_t i = 0; i < nodes_.size(); ++i) {
        nodes_[i].recent_performance = timings[i];
        
        // Update performance factor (inverse of timing)
        if (timings[i] > 0.0) {
            nodes_[i].performance_factor = 1.0 / timings[i];
        }
    }
    
    UpdatePerformanceRatios();
    UpdateTimingRatios();
}

void HeterogeneousLB::UpdateSingleNodeMetrics(int node_id, double timing) {
    int index = GetNodeIndex(node_id);
    if (!IsValidNodeIndex(index)) {
        amrex::Print() << "Warning: Invalid node_id " << node_id << "\n";
        return;
    }
    
    nodes_[index].recent_performance = timing;
    if (timing > 0.0) {
        nodes_[index].performance_factor = 1.0 / timing;
    }
    
    // Update only affected ratios
    UpdateTimingRatiosForNode(index);
    UpdatePerformanceRatios();
}

Vector<Long> HeterogeneousLB::NormalizeWeights(const Vector<Long>& weights) const {
    if (weights.empty()) return weights;
    
    auto [min_it, max_it] = std::minmax_element(weights.begin(), weights.end());
    Long min_weight = *min_it;
    Long max_weight = *max_it;
    
    if (max_weight == min_weight) {
        return Vector<Long>(weights.size(), 1000);  // All equal weights
    }
    
    // Scale to range [1, 10000] preserving ratios
    Vector<Long> normalized(weights.size());
    double range = static_cast<double>(max_weight - min_weight);
    
    for (size_t i = 0; i < weights.size(); ++i) {
        double ratio = static_cast<double>(weights[i] - min_weight) / range;
        normalized[i] = 1 + static_cast<Long>(ratio * 9999);
    }
    
    return normalized;
}

Vector<Long> HeterogeneousLB::AdjustWeightsForNodes(const MultiFab& weights) const {
    BL_PROFILE("HeterogeneousLB::AdjustWeightsForNodes");
    const int nboxes = weights.boxArray().size();
    Vector<Long> adjusted_weights(nboxes);
    
    if (nboxes == 0) {
        amrex::Abort("HeterogeneousLB::AdjustWeightsForNodes: Empty BoxArray");
    }
    
    if (nodes_.empty() || timing_ratios_.empty()) {
        // Fallback to original weights with reasonable scaling
        for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
            const int i = mfi.index();
            const Box& box = mfi.validbox();
            double orig_weight = weights[mfi].sum<RunOn::Host>(box, 0, 1);
            // Scale by box size to get reasonable values
            adjusted_weights[i] = static_cast<Long>(orig_weight * box.numPts());
        }
        return adjusted_weights;
    }
    
    // Get the distribution mapping
    const DistributionMapping& dmap = weights.DistributionMap();
    
    // Better weight adjustment strategy
    for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
        const int i = mfi.index();
        const int proc_id = dmap[i];
        const int node_index = proc_id % nodes_.size();
        const Box& box = mfi.validbox();
        
        // Get original weight
        double orig_weight = weights[mfi].sum<RunOn::Host>(box, 0, 1);
        
        // Scale by computational complexity (box size * weight factor)
        double computational_weight = orig_weight * box.numPts();
        
        // Apply performance adjustment
        if (IsValidNodeIndex(node_index) && 
            node_index < static_cast<int>(timing_ratios_.size())) {
            
            // Find average timing ratio for this node
            double avg_ratio = 1.0;
            int valid_ratios = 0;
            for (size_t j = 0; j < timing_ratios_[node_index].size(); ++j) {
                if (timing_ratios_[node_index][j] > 0.0 && 
                    std::isfinite(timing_ratios_[node_index][j])) {
                    avg_ratio += timing_ratios_[node_index][j];
                    valid_ratios++;
                }
            }
            if (valid_ratios > 0) {
                avg_ratio /= (valid_ratios + 1);  // +1 for the initial 1.0
            }
            
            computational_weight *= avg_ratio;
        }
        
        // Ensure reasonable range [1, 1000000]
        adjusted_weights[i] = static_cast<Long>(
            std::max(1.0, std::min(computational_weight, 1000000.0))
        );
    }
    
    return adjusted_weights;
}

LoadBalanceMetrics HeterogeneousLB::EvaluateBalance(const DistributionMapping& dm, 
                                                   const MultiFab& weights) const {
    LoadBalanceMetrics metrics;
    
    if (nodes_.empty()) return metrics;
    
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    Vector<double> node_loads(nodes_.size(), 0.0);
    Vector<Long> node_work(nodes_.size(), 0);
    
    // Calculate load per node
    for (int i = 0; i < dm.size(); ++i) {
        int node_idx = GetNodeIndex(dm[i]);
        if (IsValidNodeIndex(node_idx)) {
            Long work = adjusted_weights[i];
            node_work[node_idx] += work;
            // Proper execution time calculation
            double execution_time = static_cast<double>(work) / 
                                  (nodes_[node_idx].performance_factor * 1000000.0); // Normalize
            node_loads[node_idx] += execution_time;
            metrics.total_work += work;
        }
    }
    
    // Calculate metrics
    if (!node_loads.empty()) {
        metrics.makespan = *std::max_element(node_loads.begin(), node_loads.end());
        
        if (metrics.makespan > 0.0) {
            double total_execution_time = std::accumulate(node_loads.begin(), node_loads.end(), 0.0);
            metrics.efficiency = total_execution_time / (metrics.makespan * nodes_.size());
        }
        
        double avg_load = std::accumulate(node_loads.begin(), node_loads.end(), 0.0) / node_loads.size();
        if (avg_load > 0.0) {
            metrics.load_imbalance = (metrics.makespan - avg_load) / avg_load;
        }
    }
    
    return metrics;
}

uint64_t HeterogeneousLB::ComputeHilbertIndex(const IntVect& iv) const {
    // Simple Morton code implementation (can be replaced with true Hilbert)
    uint64_t result = 0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        uint32_t coord = static_cast<uint32_t>(iv[d]);
        for (int bit = 0; bit < 21; ++bit) {  // 21 bits per dimension
            if (coord & (1u << bit)) {
                result |= (1ull << (bit * AMREX_SPACEDIM + d));
            }
        }
    }
    return result;
}

Vector<int> HeterogeneousLB::CreateSpatialOrdering(const BoxArray& ba) const {
    Vector<std::pair<uint64_t, int>> spatial_pairs;
    spatial_pairs.reserve(ba.size());
    
    for (int i = 0; i < ba.size(); ++i) {
        Box box = ba[i];
        uint64_t spatial_index = ComputeHilbertIndex(box.smallEnd());
        spatial_pairs.push_back({spatial_index, i});
    }
    
    // Sort by spatial index
    std::sort(spatial_pairs.begin(), spatial_pairs.end());
    
    Vector<int> ordering;
    ordering.reserve(ba.size());
    for (const auto& pair : spatial_pairs) {
        ordering.push_back(pair.second);
    }
    
    return ordering;
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
    
    if (nboxes == 0 || nodes_.empty()) {
        return DistributionMapping(ba, nodes_.size());
    }
    
    // Better weight adjustment strategy
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    
    // Improved normalization to prevent overflow and maintain precision
    Vector<Long> normalized_weights(nboxes);
    if (!adjusted_weights.empty()) {
        auto [min_it, max_it] = std::minmax_element(adjusted_weights.begin(), adjusted_weights.end());
        Long min_weight = *min_it;
        Long max_weight = *max_it;
        
        if (max_weight > min_weight) {
            // Scale to range [100, 10000] to maintain precision
            double scale = 9900.0 / (max_weight - min_weight);
            for (int i = 0; i < nboxes; ++i) {
                normalized_weights[i] = 100 + static_cast<Long>((adjusted_weights[i] - min_weight) * scale);
            }
        } else {
            // All weights equal
            std::fill(normalized_weights.begin(), normalized_weights.end(), 1000L);
        }
    }
    
    // Multiple knapsack attempts with different parameters
    std::vector<std::vector<int>> candidate_assignments;
    std::vector<double> candidate_efficiencies;
    
    // Try different knapsack configurations
    std::vector<bool> knapsack_modes = {true, false}; // full vs partial knapsack
    
    for (bool full_knapsack : knapsack_modes) {
        try {
            Real efficiency;
            std::vector<Long> std_weights(normalized_weights.begin(), normalized_weights.end());
            auto assignment = KnapSackDoIt(std_weights,
                                         nodes_.size(),
                                         efficiency,
                                         full_knapsack,
                                         std::numeric_limits<int>::max(),
                                         false,
                                         true);
            
            if (!assignment.empty() && assignment.size() == nboxes) {
                candidate_assignments.push_back(assignment);
                candidate_efficiencies.push_back(efficiency);
            }
        } catch (...) {
            // Skip failed attempts
            continue;
        }
    }
    
    // Choose best assignment based on load balance
    int best_assignment_idx = -1;
    double best_imbalance = std::numeric_limits<double>::max();
    
    for (size_t a = 0; a < candidate_assignments.size(); ++a) {
        // Calculate load imbalance for this assignment
        Vector<double> node_loads(nodes_.size(), 0.0);
        for (int i = 0; i < nboxes; ++i) {
            int node_idx = candidate_assignments[a][i];
            if (node_idx >= 0 && node_idx < static_cast<int>(nodes_.size())) {
                node_loads[node_idx] += static_cast<double>(adjusted_weights[i]) / nodes_[node_idx].performance_factor;
            }
        }
        
        double max_load = *std::max_element(node_loads.begin(), node_loads.end());
        double avg_load = std::accumulate(node_loads.begin(), node_loads.end(), 0.0) / nodes_.size();
        double imbalance = (avg_load > 0) ? (max_load - avg_load) / avg_load : 0.0;
        
        if (imbalance < best_imbalance) {
            best_imbalance = imbalance;
            best_assignment_idx = a;
        }
    }
    
    // Use best assignment or fallback to greedy
    Vector<int> proc_mapping(nboxes);
    if (best_assignment_idx >= 0) {
        auto& best_assignment = candidate_assignments[best_assignment_idx];
        for (int i = 0; i < nboxes; ++i) {
            int knapsack_node = best_assignment[i];
            if (knapsack_node >= 0 && knapsack_node < static_cast<int>(nodes_.size())) {
                proc_mapping[i] = nodes_[knapsack_node].node_id;
            } else {
                // Fallback assignment
                proc_mapping[i] = nodes_[i % nodes_.size()].node_id;
            }
        }
        
        if (amrex::Verbose() > 1) {
            amrex::Print() << "\nKnapsack Distribution (efficiency=" 
                          << candidate_efficiencies[best_assignment_idx] 
                          << ", imbalance=" << (best_imbalance * 100.0) << "%):\n";
        }
    } else {
        // Fallback to simple round-robin if knapsack fails
        amrex::Print() << "Warning: Knapsack failed, using round-robin fallback\n";
        for (int i = 0; i < nboxes; ++i) {
            proc_mapping[i] = nodes_[i % nodes_.size()].node_id;
        }
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

DistributionMapping HeterogeneousLB::SFCBalance(const BoxArray& ba, const MultiFab& weights) {
    BL_PROFILE("HeterogeneousLB::SFCBalance");
    const int nboxes = ba.size();
    
    if (nboxes == 0 || nodes_.empty()) {
        return DistributionMapping(ba, nodes_.size());
    }
    
    // Create spatial ordering using Hilbert-like curve
    Vector<int> spatial_order = CreateSpatialOrdering(ba);
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    
    // Calculate cumulative performance for load balancing
    Vector<double> cumulative_perf(nodes_.size() + 1, 0.0);
    for (size_t i = 0; i < nodes_.size(); ++i) {
        cumulative_perf[i + 1] = cumulative_perf[i] + nodes_[i].performance_factor;
    }
    double total_perf = cumulative_perf.back();
    
    // Calculate total work
    Long total_work = 0;
    for (Long w : adjusted_weights) {
        total_work += w;
    }
    
    // Distribute boxes maintaining spatial locality
    Vector<int> proc_mapping(nboxes);
    Long current_work = 0;
    int current_node = 0;
    
    for (int i = 0; i < nboxes; ++i) {
        int box_idx = spatial_order[i];
        current_work += adjusted_weights[box_idx];
        
        // Check if we should move to next node
        double work_fraction = static_cast<double>(current_work) / total_work;
        double target_perf_fraction = static_cast<double>(cumulative_perf[current_node + 1]) / total_perf;
        
        if (work_fraction > target_perf_fraction && current_node < static_cast<int>(nodes_.size()) - 1) {
            current_node++;
        }
        
        proc_mapping[box_idx] = nodes_[current_node].node_id;
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

std::map<std::string, Vector<int>> HeterogeneousLB::GroupNodesByType() const {
    std::map<std::string, Vector<int>> groups;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        groups[nodes_[i].node_type].push_back(i);
    }
    return groups;
}

DistributionMapping HeterogeneousLB::GroupedRankBalance(const BoxArray& ba, const MultiFab& weights) {
    BL_PROFILE("HeterogeneousLB::GroupedRankBalance");
    const int nboxes = ba.size();
    
    if (nboxes == 0 || nodes_.empty()) {
        return DistributionMapping(ba, nodes_.size());
    }
    
    auto node_groups = GroupNodesByType();
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    Vector<int> proc_mapping(nboxes);
    
    // Calculate total capacity
    double total_capacity = 0.0;
    for (const auto& node : nodes_) {
        total_capacity += node.performance_factor;
    }
    
    // Distribute work among groups
    int assigned_boxes = 0;
    for (auto& group : node_groups) {
        // Calculate group capacity
        double group_capacity = 0.0;
        for (int node_idx : group.second) {
            group_capacity += nodes_[node_idx].performance_factor;
        }
        
        // Calculate boxes for this group
        int group_box_count = static_cast<int>(
            std::round((group_capacity / total_capacity) * nboxes));
        
        // Ensure we don't exceed total boxes
        group_box_count = std::min(group_box_count, nboxes - assigned_boxes);
        
        if (group_box_count > 0) {
            // Create weights for this group
            std::vector<Long> group_weights;
            for (int i = assigned_boxes; i < assigned_boxes + group_box_count; ++i) {
                group_weights.push_back(adjusted_weights[i]);
            }
            
            // Balance within group using knapsack
            Real efficiency;
            auto group_mapping = KnapSackDoIt(group_weights,
                                            group.second.size(),
                                            efficiency,
                                            true,
                                            std::numeric_limits<int>::max(),
                                            false,
                                            true);
            
            // Map back to global indices
            for (int i = 0; i < group_box_count; ++i) {
                int global_box_idx = assigned_boxes + i;
                int group_node_idx = group_mapping[i];
                if (group_node_idx >= 0 && group_node_idx < static_cast<int>(group.second.size())) {
                    int actual_node_idx = group.second[group_node_idx];
                    proc_mapping[global_box_idx] = nodes_[actual_node_idx].node_id;
                } else {
                    // Fallback
                    proc_mapping[global_box_idx] = nodes_[group.second[0]].node_id;
                }
            }
            
            assigned_boxes += group_box_count;
        }
    }
    
    // Assign any remaining boxes to first available nodes
    for (int i = assigned_boxes; i < nboxes; ++i) {
        proc_mapping[i] = nodes_[i % nodes_.size()].node_id;
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

void HeterogeneousLB::UpdateOrderedRatios() {
    ordered_ratios_.clear();
    const int n_nodes = nodes_.size();
    
    if (timing_ratios_.empty() || timing_ratios_.size() != n_nodes) {
        UpdateTimingRatios();
    }
    
    // Create ordered list of timing ratios
    for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
            if (i != j && i < static_cast<int>(timing_ratios_.size()) && 
                j < static_cast<int>(timing_ratios_[i].size())) {
                ordered_ratios_.push_back(TimingRatio(i, j, timing_ratios_[i][j]));
            }
        }
    }
    
    // Sort ratios in ascending order
    std::sort(ordered_ratios_.begin(), ordered_ratios_.end());
}

Vector<int> HeterogeneousLB::GetOrderedGroupedRanks() const {
    const int n_nodes = nodes_.size();
    Vector<int> ranks(n_nodes);
    
    // Group nodes by type and sort by performance
    auto type_groups = GroupNodesByType();
    
    Vector<std::pair<std::string, double>> sorted_groups;
    for (const auto& group : type_groups) {
        double avg_perf = 0.0;
        for (int node_idx : group.second) {
            avg_perf += nodes_[node_idx].performance_factor;
        }
        avg_perf /= group.second.size();
        sorted_groups.push_back({group.first, avg_perf});
    }
    
    std::sort(sorted_groups.begin(), sorted_groups.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Create ordered rank list
    int current_rank = 0;
    for (const auto& group : sorted_groups) {
        for (int node_idx : type_groups.at(group.first)) {
            ranks[node_idx] = current_rank++;
        }
    }
    
    return ranks;
}

Vector<Long> HeterogeneousLB::GetOrderedWeights(const MultiFab& weights) const {
    Vector<Long> ordered_weights = AdjustWeightsForNodes(weights);
    std::sort(ordered_weights.begin(), ordered_weights.end(), std::greater<Long>());
    return ordered_weights;
}

// ENHANCED RIJ IMPLEMENTATION

DistributionMapping HeterogeneousLB::RijBalance(const BoxArray& ba, const MultiFab& weights) {
    BL_PROFILE("HeterogeneousLB::RijBalance");
    
    const int nboxes = ba.size();
    if (nboxes == 0 || nodes_.empty()) {
        return DistributionMapping(ba, nodes_.size());
    }
    
    // ENHANCED: Update timing ratios with better precision
    UpdateOrderedRatios();
    Vector<Long> adjusted_weights = AdjustWeightsForNodes(weights);
    
    // Initialize box centers for spatial locality
    box_centers_.clear();
    box_centers_.reserve(nboxes);
    for (int i = 0; i < nboxes; ++i) {
        Box box = ba[i];
        box_centers_.push_back(box.smallEnd() + box.size() / 2);
    }
    
    // ENHANCED: Multi-phase RIJ optimization
    Vector<int> proc_mapping = EnhancedRijOptimization(ba, adjusted_weights);
    
    if (amrex::Verbose() > 1) {
        Vector<Long> node_weights(nodes_.size(), 0);
        for (int i = 0; i < nboxes; ++i) {
            int node_idx = GetNodeIndex(proc_mapping[i]);
            if (IsValidNodeIndex(node_idx)) {
                node_weights[node_idx] += adjusted_weights[i];
            }
        }
        PrintDistributionStats("Enhanced Rij-based", proc_mapping, node_weights);
    }
    
    return DistributionMapping(std::move(proc_mapping));
}

Vector<int> HeterogeneousLB::EnhancedRijOptimization(const BoxArray& ba, const Vector<Long>& weights) {
    BL_PROFILE("HeterogeneousLB::EnhancedRijOptimization");
    
    // Phase 1: Intelligent initial assignment using RIJ ratios
    Vector<int> initial_assignment = RijAwareInitialAssignment(ba, weights);
    
    // Phase 2: Local optimization using timing ratios
    Vector<int> optimized_assignment = RijLocalOptimization(initial_assignment, weights);
    
    // Phase 3: Global refinement with load balancing
    Vector<int> final_assignment = RijGlobalRefinement(optimized_assignment, weights);
    
    return final_assignment;
}

Vector<int> HeterogeneousLB::RijAwareInitialAssignment(const BoxArray& ba, const Vector<Long>& weights) {
    BL_PROFILE("HeterogeneousLB::RijAwareInitialAssignment");
    
    const int nboxes = ba.size();
    Vector<int> assignment(nboxes);
    
    // Create weighted boxes with spatial information
    Vector<WeightedBox> weighted_boxes;
    weighted_boxes.reserve(nboxes);
    
    for (int i = 0; i < nboxes; ++i) {
        WeightedBox wb;
        wb.box_id = i;
        wb.weight = weights[i];
        wb.box = ba[i];
        wb.center = wb.box.smallEnd() + wb.box.size() / 2;
        wb.priority = CalculateBoxPriority(wb.weight, wb.center);
        weighted_boxes.push_back(wb);
    }
    
    // Sort by priority (heaviest and most central first)
    std::sort(weighted_boxes.begin(), weighted_boxes.end(),
              [](const WeightedBox& a, const WeightedBox& b) {
                  return a.priority > b.priority;
              });
    
    // Enhanced node load tracking
    Vector<NodeLoadInfo> node_loads(nodes_.size());
    for (size_t i = 0; i < nodes_.size(); ++i) {
        node_loads[i].node_idx = i;
        node_loads[i].current_load = 0.0;
        node_loads[i].box_count = 0;
        node_loads[i].performance_factor = nodes_[i].performance_factor;
        node_loads[i].recent_boxes.clear();
    }
    
    // Assign boxes using enhanced RIJ logic
    for (const auto& wb : weighted_boxes) {
        int best_node = SelectBestNodeForBox(wb, node_loads, nboxes);
        assignment[wb.box_id] = nodes_[best_node].node_id;
        
        // Update node load
        double execution_time = static_cast<double>(wb.weight) / nodes_[best_node].performance_factor;
        node_loads[best_node].current_load += execution_time;
        node_loads[best_node].box_count++;
        node_loads[best_node].recent_boxes.push_back(wb.box_id);
        
        // Keep recent_boxes list manageable
        if (node_loads[best_node].recent_boxes.size() > 10) {
            node_loads[best_node].recent_boxes.erase(node_loads[best_node].recent_boxes.begin());
        }
    }
    
    return assignment;
}

int HeterogeneousLB::SelectBestNodeForBox(const WeightedBox& wb, Vector<NodeLoadInfo>& node_loads, int total_boxes) {
    int best_node = 0;
    double best_score = std::numeric_limits<double>::max();
    
    for (size_t node_idx = 0; node_idx < nodes_.size(); ++node_idx) {
        double score = CalculateNodeScore(wb, node_loads[node_idx], node_idx, total_boxes);
        
        if (score < best_score) {
            best_score = score;
            best_node = node_idx;
        }
    }
    
    return best_node;
}

double HeterogeneousLB::CalculateNodeScore(const WeightedBox& wb, const NodeLoadInfo& node_info, int node_idx, int total_boxes) {
    // Multi-factor scoring function
    double score = 0.0;
    
    // Factor 1: Projected completion time (primary)
    double execution_time = static_cast<double>(wb.weight) / node_info.performance_factor;
    double completion_time = node_info.current_load + execution_time;
    score += completion_time * 1.0;  // Weight: 1.0
    
    // Factor 2: Load imbalance penalty
    double total_current_load = 0.0;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        total_current_load += nodes_[i].recent_performance;
    }
    double avg_load = total_current_load / nodes_.size();
    
    if (avg_load > 0) {
        double imbalance_penalty = std::abs(completion_time - avg_load) / avg_load;
        score += imbalance_penalty * 0.3;  // Weight: 0.3
    }
    
    // Factor 3: RIJ timing ratio advantage
    if (!timing_ratios_.empty() && node_idx < static_cast<int>(timing_ratios_.size())) {
        double min_ratio = 1.0;
        for (size_t j = 0; j < timing_ratios_[node_idx].size(); ++j) {
            if (j != node_idx && timing_ratios_[node_idx][j] > 0) {
                min_ratio = std::min(min_ratio, timing_ratios_[node_idx][j]);
            }
        }
        
        double ratio_bonus = (1.0 - min_ratio) * 0.2;  // Weight: 0.2
        score -= ratio_bonus;  // Subtract because lower score is better
    }
    
    // Factor 4: Spatial locality bonus
    double locality_bonus = CalculateSpatialLocalityBonus(wb, node_info);
    score -= locality_bonus * 0.1;  // Weight: 0.1
    
    // Factor 5: Box count balancing
    double expected_boxes_per_node = static_cast<double>(total_boxes) / nodes_.size();
    double box_count_penalty = std::abs(static_cast<double>(node_info.box_count) - expected_boxes_per_node) / expected_boxes_per_node;
    score += box_count_penalty * 0.2;  // Weight: 0.2
    
    return score;
}

double HeterogeneousLB::CalculateSpatialLocalityBonus(const WeightedBox& wb, const NodeLoadInfo& node_info) {
    if (node_info.recent_boxes.empty()) return 0.0;
    
    double bonus = 0.0;
    int nearby_boxes = 0;
    
    // Check spatial proximity to recently assigned boxes on this node
    for (int recent_box_id : node_info.recent_boxes) {
        if (recent_box_id >= 0 && recent_box_id < static_cast<int>(box_centers_.size())) {
            IntVect recent_center = box_centers_[recent_box_id];
            double distance = CalculateBoxDistance(wb.center, recent_center);
            
            if (distance < spatial_locality_threshold_) {
                bonus += 1.0 / (1.0 + distance);  // Closer boxes get higher bonus
                nearby_boxes++;
            }
        }
    }
    
    return (nearby_boxes > 0) ? bonus / nearby_boxes : 0.0;
}

Vector<int> HeterogeneousLB::RijLocalOptimization(const Vector<int>& initial_assignment, const Vector<Long>& weights) {
    BL_PROFILE("HeterogeneousLB::RijLocalOptimization");
    
    Vector<int> assignment = initial_assignment;
    const int nboxes = assignment.size();
    bool improved = true;
    int iteration = 0;
    const int max_iterations = 5;
    
    while (improved && iteration < max_iterations) {
        improved = false;
        iteration++;
        
        // Try to improve assignment using RIJ ratios
        for (int box_id = 0; box_id < nboxes; ++box_id) {
            int current_node_idx = GetNodeIndex(assignment[box_id]);
            if (!IsValidNodeIndex(current_node_idx)) continue;
            
            // Find better node using timing ratios
            int better_node_idx = FindBetterNodeUsingRij(box_id, current_node_idx, weights[box_id], assignment);
            
            if (better_node_idx != current_node_idx && IsValidNodeIndex(better_node_idx)) {
                // Check if this move improves overall balance
                if (WouldImproveMakespan(box_id, current_node_idx, better_node_idx, assignment, weights)) {
                    assignment[box_id] = nodes_[better_node_idx].node_id;
                    improved = true;
                }
            }
        }
    }
    
    return assignment;
}

int HeterogeneousLB::FindBetterNodeUsingRij(int box_id, int current_node_idx, Long box_weight, 
                                           const Vector<int>& assignment) {
    if (timing_ratios_.empty() || current_node_idx >= static_cast<int>(timing_ratios_.size())) {
        return current_node_idx;
    }
    
    int best_node_idx = current_node_idx;
    double best_improvement = 0.0;
    
    // Check all other nodes using RIJ ratios
    for (size_t target_node_idx = 0; target_node_idx < nodes_.size(); ++target_node_idx) {
        if (target_node_idx == current_node_idx) continue;
        
        // Calculate timing ratio advantage
        if (target_node_idx >= timing_ratios_[current_node_idx].size()) continue;
        
        double rij = timing_ratios_[current_node_idx][target_node_idx];
        if (rij <= 0.0 || !std::isfinite(rij)) continue;
        
        // Calculate potential improvement
        double current_time = static_cast<double>(box_weight) / nodes_[current_node_idx].performance_factor;
        double target_time = static_cast<double>(box_weight) / nodes_[target_node_idx].performance_factor;
        
        double time_improvement = current_time - target_time;
        
        // Consider load balance impact
        Vector<double> node_loads = CalculateCurrentLoads(assignment, Vector<Long>(assignment.size(), box_weight));
        double load_balance_improvement = CalculateLoadBalanceImprovement(
            box_id, current_node_idx, target_node_idx, node_loads, box_weight);
        
        double total_improvement = time_improvement + load_balance_improvement * 0.5;
        
        if (total_improvement > best_improvement) {
            best_improvement = total_improvement;
            best_node_idx = target_node_idx;
        }
    }
    
    return best_node_idx;
}

Vector<int> HeterogeneousLB::RijGlobalRefinement(const Vector<int>& assignment, const Vector<Long>& weights) {
    BL_PROFILE("HeterogeneousLB::RijGlobalRefinement");
    
    Vector<int> refined_assignment = assignment;
    
    // Global load balancing pass
    Vector<double> node_loads = CalculateCurrentLoads(assignment, weights);
    double max_load = *std::max_element(node_loads.begin(), node_loads.end());
    double avg_load = std::accumulate(node_loads.begin(), node_loads.end(), 0.0) / nodes_.size();
    
    if (max_load > avg_load * 1.1) {  // If imbalance > 10%
        refined_assignment = RebalanceOverloadedNodes(refined_assignment, weights, node_loads);
    }
    
    return refined_assignment;
}

Vector<int> HeterogeneousLB::RebalanceOverloadedNodes(const Vector<int>& assignment, 
                                                     const Vector<Long>& weights,
                                                     const Vector<double>& node_loads) {
    Vector<int> rebalanced = assignment;
    const int nboxes = assignment.size();
    
    // Find overloaded and underloaded nodes
    double avg_load = std::accumulate(node_loads.begin(), node_loads.end(), 0.0) / nodes_.size();
    
    Vector<int> overloaded_nodes, underloaded_nodes;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (node_loads[i] > avg_load * 1.2) {
            overloaded_nodes.push_back(i);
        } else if (node_loads[i] < avg_load * 0.8) {
            underloaded_nodes.push_back(i);
        }
    }
    
    // Move boxes from overloaded to underloaded nodes
    Vector<double> current_loads = node_loads;  // Make a copy we can modify
    
    for (int overloaded_idx : overloaded_nodes) {
        Vector<std::pair<int, Long>> boxes_on_node;
        
        // Find boxes on overloaded node
        for (int box_id = 0; box_id < nboxes; ++box_id) {
            int node_idx = GetNodeIndex(rebalanced[box_id]);
            if (node_idx == overloaded_idx) {
                boxes_on_node.push_back({box_id, weights[box_id]});
            }
        }
        
        // Sort by weight (move lighter boxes first for stability)
        std::sort(boxes_on_node.begin(), boxes_on_node.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
        
        // Move boxes to underloaded nodes
        for (const auto& [box_id, box_weight] : boxes_on_node) {
            if (underloaded_nodes.empty()) break;
            
            // Find best underloaded node using RIJ ratios
            int best_target = FindBestTargetNode(overloaded_idx, underloaded_nodes, box_weight);
            
            if (best_target >= 0) {
                rebalanced[box_id] = nodes_[best_target].node_id;
                
                // Update loads and check if target is no longer underloaded
                double moved_time = static_cast<double>(box_weight) / nodes_[best_target].performance_factor;
                current_loads[best_target] += moved_time;
                
                if (current_loads[best_target] > avg_load * 0.9) {
                    underloaded_nodes.erase(
                        std::remove(underloaded_nodes.begin(), underloaded_nodes.end(), best_target),
                        underloaded_nodes.end());
                }
                
                // Check if source is no longer overloaded
                double removed_time = static_cast<double>(box_weight) / nodes_[overloaded_idx].performance_factor;
                current_loads[overloaded_idx] -= removed_time;
                
                if (current_loads[overloaded_idx] < avg_load * 1.1) {
                    break;  // This node is balanced enough
                }
            }
        }
    }
    
    return rebalanced;
}

// HELPER METHOD IMPLEMENTATIONS

double HeterogeneousLB::CalculateBoxPriority(Long weight, const IntVect& center) {
    // Higher priority for heavier boxes and boxes closer to domain center
    double weight_factor = static_cast<double>(weight);
    
    // Calculate distance from domain center (assuming domain is roughly centered at origin)
    double distance = std::sqrt(static_cast<double>(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]));
    double locality_factor = 1.0 / (1.0 + distance * 0.01);  // Closer to center = higher priority
    
    return weight_factor * (1.0 + locality_factor);
}

double HeterogeneousLB::CalculateBoxDistance(const IntVect& center1, const IntVect& center2) {
    double dist_sq = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        double diff = static_cast<double>(center1[d] - center2[d]);
        dist_sq += diff * diff;
    }
    return std::sqrt(dist_sq);
}

bool HeterogeneousLB::WouldImproveMakespan(int box_id, int from_node, int to_node, 
                                          const Vector<int>& assignment, const Vector<Long>& weights) {
    Vector<double> current_loads = CalculateCurrentLoads(assignment, weights);
    double current_makespan = *std::max_element(current_loads.begin(), current_loads.end());
    
    // Simulate the move
    Long box_weight = weights[box_id];
    double from_time = static_cast<double>(box_weight) / nodes_[from_node].performance_factor;
    double to_time = static_cast<double>(box_weight) / nodes_[to_node].performance_factor;
    
    current_loads[from_node] -= from_time;
    current_loads[to_node] += to_time;
    
    double new_makespan = *std::max_element(current_loads.begin(), current_loads.end());
    
    return new_makespan < current_makespan;
}

Vector<double> HeterogeneousLB::CalculateCurrentLoads(const Vector<int>& assignment, const Vector<Long>& weights) {
    Vector<double> loads(nodes_.size(), 0.0);
    
    for (size_t i = 0; i < assignment.size() && i < weights.size(); ++i) {
        int node_idx = GetNodeIndex(assignment[i]);
        if (IsValidNodeIndex(node_idx)) {
            double execution_time = static_cast<double>(weights[i]) / nodes_[node_idx].performance_factor;
            loads[node_idx] += execution_time;
        }
    }
    
    return loads;
}

double HeterogeneousLB::CalculateLoadBalanceImprovement(int box_id, int from_node, int to_node, 
                                                       const Vector<double>& node_loads, Long box_weight) {
    double from_time = static_cast<double>(box_weight) / nodes_[from_node].performance_factor;
    double to_time = static_cast<double>(box_weight) / nodes_[to_node].performance_factor;
    
    double current_imbalance = *std::max_element(node_loads.begin(), node_loads.end()) - 
                              *std::min_element(node_loads.begin(), node_loads.end());
    
    Vector<double> new_loads = node_loads;
    new_loads[from_node] -= from_time;
    new_loads[to_node] += to_time;
    
    double new_imbalance = *std::max_element(new_loads.begin(), new_loads.end()) - 
                          *std::min_element(new_loads.begin(), new_loads.end());
    
    return current_imbalance - new_imbalance;  // Positive if improvement
}

int HeterogeneousLB::FindBestTargetNode(int source_node, const Vector<int>& candidate_nodes, Long box_weight) {
    if (candidate_nodes.empty()) return -1;
    
    int best_node = candidate_nodes[0];
    double best_score = std::numeric_limits<double>::max();
    
    for (int target_node : candidate_nodes) {
        double score = 0.0;
        
        // Factor 1: Performance advantage
        double target_time = static_cast<double>(box_weight) / nodes_[target_node].performance_factor;
        score += target_time;
        
        // Factor 2: RIJ timing ratio
        if (!timing_ratios_.empty() && 
            source_node < static_cast<int>(timing_ratios_.size()) && 
            target_node < static_cast<int>(timing_ratios_[source_node].size())) {
            double rij = timing_ratios_[source_node][target_node];
            if (rij > 0.0 && std::isfinite(rij)) {
                score *= rij;  // Lower ratio = better target
            }
        }
        
        if (score < best_score) {
            best_score = score;
            best_node = target_node;
        }
    }
    
    return best_node;
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
        int active_nodes = 0;
        
        for (const auto& node : nodes_) {
            if (node.recent_performance > 0.0) {
                max_time = std::max(max_time, node.recent_performance);
                total_time += node.recent_performance;
                active_nodes++;
            }
        }
        
        if (active_nodes > 0) {
            double avg_time = total_time / active_nodes;
            double imbalance = ((max_time / avg_time) - 1.0) * 100.0;
            amrex::Print() << "Current load imbalance: " << std::fixed << std::setprecision(2) 
                          << imbalance << "%\n";
        }
    } else {
        amrex::Print() << "Current load imbalance: No performance data available\n";
    }
    
    // Print node info
    amrex::Print() << "\nNode Performance Breakdown:\n";
    amrex::Print() << std::setw(8) << "Node ID" 
                   << std::setw(15) << "Type"
                   << std::setw(12) << "Perf Factor"
                   << std::setw(15) << "Recent Time"
                   << std::setw(12) << "Memory GB" << "\n";
    amrex::Print() << std::string(62, '-') << "\n";
    
    for (const auto& node : nodes_) {
        amrex::Print() << std::setw(8) << node.node_id
                       << std::setw(15) << node.node_type.substr(0, 14)
                       << std::setw(12) << std::fixed << std::setprecision(3) << node.performance_factor;
        
        if (node.recent_performance > 0.0) {
            amrex::Print() << std::setw(15) << std::fixed << std::setprecision(3) << node.recent_performance;
        } else {
            amrex::Print() << std::setw(15) << "N/A";
        }
        
        amrex::Print() << std::setw(12) << std::fixed << std::setprecision(1) << node.memory_capacity << "\n";
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
    
    if (total_load <= 0) return 0.0;
    
    double avg_load = static_cast<double>(total_load) / loads.size();
    return (static_cast<double>(max_load) / avg_load) - 1.0;
}

void HeterogeneousLB::PrintDistributionStats(const std::string& method,
                                           const Vector<int>& mapping,
                                           const Vector<Long>& loads) const {
    amrex::Print() << "\n" << method << " Distribution Statistics:\n";
    amrex::Print() << std::setw(8) << "Node ID" 
                   << std::setw(15) << "Load"
                   << std::setw(12) << "Box Count"
                   << std::setw(15) << "Efficiency" << "\n";
    amrex::Print() << std::string(50, '-') << "\n";
    
    Vector<int> box_counts(nodes_.size(), 0);
    for (size_t i = 0; i < mapping.size(); ++i) {
        int node_idx = GetNodeIndex(mapping[i]);
        if (IsValidNodeIndex(node_idx)) {
            box_counts[node_idx]++;
        }
    }
    
    Long total_load = std::accumulate(loads.begin(), loads.end(), 0L);
    
    for (size_t i = 0; i < nodes_.size(); ++i) {
        double efficiency = (total_load > 0) ? 
            (static_cast<double>(loads[i]) / total_load) / (1.0 / nodes_.size()) : 0.0;
        
        amrex::Print() << std::setw(8) << nodes_[i].node_id
                       << std::setw(15) << loads[i]
                       << std::setw(12) << box_counts[i]
                       << std::setw(15) << std::fixed << std::setprecision(3) << efficiency << "\n";
    }
    
    amrex::Print() << "Load imbalance: " 
                   << std::fixed << std::setprecision(2) 
                   << (ComputeLoadImbalance(loads) * 100.0) << "%\n";
}

} // namespace amrex