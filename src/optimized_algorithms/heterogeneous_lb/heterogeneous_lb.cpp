#include "heterogeneous_lb.h"
#include "../Knapsack.H"
#include <iostream>
#include <vector>
#include <random>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <AMReX_INT.H>
#include <numeric>
#include <fstream>

struct HeterogeneityResult {
    double heterogeneity_factor;
    double homo_efficiency;
    double without_rij_efficiency;
    double with_rij_efficiency;
    double homo_makespan;
    double without_rij_makespan;
    double with_rij_makespan;
};

// Compute rij matrix: rij[i][j] = perf[j] / perf[i]
std::vector<std::vector<double>> compute_rij(const std::vector<Node>& nodes) {
    size_t n = nodes.size();
    std::vector<std::vector<double>> rij(n, std::vector<double>(n, 1.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            rij[i][j] = nodes[j].performance_factor / nodes[i].performance_factor;
        }
    }
    return rij;
}

// Print rij matrix
void print_rij_matrix(const std::vector<std::vector<double>>& rij, const std::vector<Node>& nodes) {
    std::cout << "rij matrix (rij[i][j] = perf[j]/perf[i]):\n";
    std::cout << "     ";
    for (size_t j = 0; j < nodes.size(); ++j) {
        std::cout << std::setw(6) << "N" + std::to_string(j);
    }
    std::cout << "\n";
    
    for (size_t i = 0; i < rij.size(); ++i) {
        std::cout << "N" << i << "   ";
        for (size_t j = 0; j < rij[i].size(); ++j) {
            std::cout << std::setw(6) << std::setprecision(2) << std::fixed << rij[i][j];
        }
        std::cout << " (perf=" << nodes[i].performance_factor << ")\n";
    }
    std::cout << "\n";
}

// Homogeneous knapsack - all nodes have performance_factor = 1.0
std::vector<int> homogeneous_knapsack_assignment(const std::vector<Task>& tasks, int num_nodes) {
    std::vector<amrex::Long> weights;
    weights.reserve(tasks.size());
    for (const auto& task : tasks) {
        weights.push_back(static_cast<amrex::Long>(task.base_time * 1000));
    }
    amrex::Real efficiency = 0.0;
    return KnapSackDoIt(weights, num_nodes, efficiency, true,
                        std::numeric_limits<int>::max(), false, false,
                        std::vector<amrex::Long>());
}

// Heterogeneous knapsack WITHOUT rij - basic performance-aware approach
std::vector<int> knapsack_without_rij_assignment(const std::vector<Task>& tasks, const std::vector<Node>& nodes) {
    std::vector<int> assignments(tasks.size());
    std::vector<double> node_loads(nodes.size(), 0.0);
    
    // Sort tasks by size (largest first) for better balancing
    std::vector<std::pair<double, int>> task_indices;
    for (size_t i = 0; i < tasks.size(); ++i) {
        task_indices.push_back({tasks[i].base_time, i});
    }
    std::sort(task_indices.rbegin(), task_indices.rend());
    
    // Assign each task using a naive performance-aware approach
    for (const auto& task_pair : task_indices) {
        int task_idx = task_pair.second;
        const Task& task = tasks[task_idx];
        
        // assign to the node with the lowest current load
        // scaled by performance factor (to account for heterogeneity)
        int best_node = 0;
        double min_scaled_load = std::numeric_limits<double>::max();
        
        for (size_t candidate = 0; candidate < nodes.size(); ++candidate) {
            // Scale current load by inverse of performance factor
            // This gives preference to faster nodes but doesn't account for task size
            double scaled_load = node_loads[candidate] + (task.base_time / nodes[candidate].performance_factor);
            
            if (scaled_load < min_scaled_load) {
                min_scaled_load = scaled_load;
                best_node = candidate;
            }
        }
        
        // Assign task and update load
        assignments[task_idx] = best_node;
        node_loads[best_node] += task.base_time / nodes[best_node].performance_factor;
    }
    
    return assignments;
}

std::vector<int> knapsack_with_rij_assignment(const std::vector<Task>& tasks, const std::vector<Node>& nodes, const std::vector<std::vector<double>>& rij) {
    std::vector<int> assignments(tasks.size());
    std::vector<double> node_loads(nodes.size(), 0.0);

    // Calculate total work
    double total_work = 0.0;
    for (const auto& task : tasks) {
        total_work += task.base_time;
    }

    // Sort tasks by size (largest first)
    std::vector<std::pair<double, int>> task_indices;
    for (size_t i = 0; i < tasks.size(); ++i) {
        task_indices.push_back({tasks[i].base_time, i});
    }
    std::sort(task_indices.rbegin(), task_indices.rend());
    
    for (const auto& task_pair : task_indices) {
        int task_idx = task_pair.second;
        const Task& task = tasks[task_idx];
        
        int best_node = 0;
        double best_score = std::numeric_limits<double>::max();
        
        for (size_t candidate = 0; candidate < nodes.size(); ++candidate) {
            double exec_time = task.base_time / nodes[candidate].performance_factor;
            double completion_time = node_loads[candidate] + exec_time;
            
            // Use rij matrix for load balancing awareness
            double load_balance_score = completion_time;
            
            // Factor 1: Consider how this assignment affects overall balance
            // Look at load distribution after this assignment
            double projected_load = node_loads[candidate] + exec_time;
            double balance_penalty = 0.0;

            // Calculate load variance among nodes
            double mean_load = std::accumulate(node_loads.begin(), node_loads.end(), 0.0) / nodes.size();
            double load_variance = 0.0;
            for (double load : node_loads) {
                load_variance += (load - mean_load) * (load - mean_load);
            }
            load_variance /= nodes.size();

            for (size_t other = 0; other < nodes.size(); ++other) {
                if (other != candidate) {
                    // Use rij to understand relative capacity
                    double relative_capacity = rij[other][candidate]; // How much faster 'candidate' is vs 'other'
                    double other_relative_load = node_loads[other] * relative_capacity;

                    // Dynamic imbalance threshold based on current loads and variance
                    double imbalance_threshold = 1.0 + (node_loads[candidate] / (total_work / nodes.size())) + (load_variance * 0.1);
                    if (projected_load > other_relative_load * imbalance_threshold) {
                        balance_penalty += (projected_load - other_relative_load) * 0.3;
                    }
                }
            }

            // Factor 2: Prefer nodes that can handle load redistribution well
            // Nodes with good rij relationships to other nodes are preferred
            double redistribution_bonus = 0.0;
            for (size_t other = 0; other < nodes.size(); ++other) {
                if (other != candidate) {
                    double rij_val = rij[candidate][other];
                    // Bonus for nodes that have good load transfer potential
                    if (rij_val >= 0.8 && rij_val <= 1.25) { // Similar performance range
                        redistribution_bonus += 0.15 * exec_time;
                    }
                }
            }

            // Factor 3: Consider node capacity
            double node_capacity = nodes[candidate].performance_factor;
            double capacity_bonus = node_capacity * 0.1 * exec_time;

            load_balance_score = completion_time + balance_penalty - redistribution_bonus + capacity_bonus;

            if (load_balance_score < best_score) {
                best_score = load_balance_score;
                best_node = candidate;
            }
        }

        assignments[task_idx] = best_node;
        node_loads[best_node] += task.base_time / nodes[best_node].performance_factor;
    }

    return assignments;
}

std::vector<HeterogeneityResult> test_heterogeneity_levels(const std::vector<Task>& tasks) {
    std::vector<HeterogeneityResult> results;
    
    // Different heterogeneity levels
    std::vector<std::pair<double, std::vector<double>>> heterogeneity_configs = {
        {1.0, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}},  // Homogeneous
        {1.5, {1.2, 1.2, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8}},  // Low heterogeneity
        {2.14, {1.5, 1.5, 1.2, 1.2, 1.0, 1.0, 0.7, 0.7}},  // Medium-low
        {4.0, {2.0, 2.0, 1.5, 1.5, 1.0, 1.0, 0.5, 0.5}},  // Medium (my current)
        {6.25, {2.5, 2.5, 2.0, 2.0, 1.0, 1.0, 0.4, 0.4}},  // Medium-high
        {9.09, {3.0, 3.0, 2.0, 2.0, 1.0, 1.0, 0.33, 0.33}}, // High
        {16.0, {4.0, 4.0, 2.0, 2.0, 1.0, 1.0, 0.25, 0.25}}  // Very high
    };
    
    for (const auto& config : heterogeneity_configs) {
        double het_factor = config.first; // max_perf / min_perf
        const auto& perf_factors = config.second;
        
        // Create nodes for this configuration
        std::vector<Node> nodes;
        for (size_t i = 0; i < perf_factors.size(); ++i) {
            nodes.push_back(Node(i, "cpu" + std::to_string(i), perf_factors[i], 0.0));
        }
        
        // Homogeneous nodes (all perf = 1.0)
        std::vector<Node> homo_nodes;
        for (size_t i = 0; i < perf_factors.size(); ++i) {
            homo_nodes.push_back(Node(i, "cpu" + std::to_string(i), 1.0, 0.0));
        }
        
        // Compute rij matrix
        auto rij = compute_rij(nodes);
        
        // Run all three algorithms
        auto homo_assign = homogeneous_knapsack_assignment(tasks, homo_nodes.size());
        auto without_assign = knapsack_without_rij_assignment(tasks, nodes);
        auto with_assign = knapsack_with_rij_assignment(tasks, nodes, rij);
        
        // Calculate loads and metrics
        std::vector<double> homo_loads(nodes.size(), 0.0);
        std::vector<double> without_loads(nodes.size(), 0.0);
        std::vector<double> with_loads(nodes.size(), 0.0);
        
        double total_work = 0.0;
        for (const auto& task : tasks) {
            total_work += task.base_time;
        }
        
        for (size_t i = 0; i < tasks.size(); ++i) {
            homo_loads[homo_assign[i]] += tasks[i].base_time / 1.0;
            without_loads[without_assign[i]] += tasks[i].base_time / nodes[without_assign[i]].performance_factor;
            with_loads[with_assign[i]] += tasks[i].base_time / nodes[with_assign[i]].performance_factor;
        }
        
        double homo_makespan = *std::max_element(homo_loads.begin(), homo_loads.end());
        double without_makespan = *std::max_element(without_loads.begin(), without_loads.end());
        double with_makespan = *std::max_element(with_loads.begin(), with_loads.end());
        
        double total_capacity = std::accumulate(perf_factors.begin(), perf_factors.end(), 0.0);
        double ideal_makespan = total_work / total_capacity;
        
        HeterogeneityResult result;
        result.heterogeneity_factor = *std::max_element(perf_factors.begin(), perf_factors.end()) / *std::min_element(perf_factors.begin(), perf_factors.end());
        result.homo_efficiency = ideal_makespan / homo_makespan;
        result.without_rij_efficiency = ideal_makespan / without_makespan;
        result.with_rij_efficiency = ideal_makespan / with_makespan;
        result.homo_makespan = homo_makespan;
        result.without_rij_makespan = without_makespan;
        result.with_rij_makespan = with_makespan;
        
        results.push_back(result);
    }
    
    return results;
}

int main() {
    std::cout << "=== RIJ MATRIX ADVANTAGE DEMONSTRATION ===\n\n";
    
    // Create a more heterogeneous system - wider performance gap
    std::vector<double> perf_factors = {2.0, 2.0, 1.5, 1.5, 1.0, 1.0, 0.5, 0.5};
    std::vector<Node> nodes;
    
    // std::cout << "Machine groups:\n";
    // std::cout << "  Group 1 (Very Fast): Nodes 0,1 with performance 2.0\n";
    // std::cout << "  Group 2 (Fast):      Nodes 2,3 with performance 1.5\n";
    // std::cout << "  Group 3 (Medium):    Nodes 4,5 with performance 1.0\n";
    // std::cout << "  Group 4 (Slow):      Nodes 6,7 with performance 0.5\n\n";
    
    for (size_t i = 0; i < perf_factors.size(); ++i) {
        nodes.push_back(Node(i, "cpu" + std::to_string(i), perf_factors[i], 0.0));
    }
    
    // Compute and print rij matrix
    auto rij = compute_rij(nodes);
    print_rij_matrix(rij, nodes);
    
    // Generate tasks with higher variance
    std::vector<Task> tasks;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> d(5.0, 50.0); // Wider range of task sizes
    double total_work = 0.0;
    
    for (int i = 0; i < 100; ++i) {
        double task_time = d(gen);
        tasks.push_back(Task(i, task_time));
        total_work += task_time;
    }
    
    std::cout << "Generated " << tasks.size() << " tasks, total work: " << total_work << "\n\n";
    
    // Homogeneous nodes: all perf = 1.0
    std::vector<Node> homo_nodes;
    for (size_t i = 0; i < nodes.size(); ++i) {
        homo_nodes.push_back(Node(i, "cpu" + std::to_string(i), 1.0, 0.0));
    }
    auto without_rij_nodes = nodes;
    auto with_rij_nodes = nodes;
    
    // Assignment
    auto homo_assign = homogeneous_knapsack_assignment(tasks, homo_nodes.size());
    auto without_assign = knapsack_without_rij_assignment(tasks, without_rij_nodes);
    auto with_assign = knapsack_with_rij_assignment(tasks, with_rij_nodes, rij);
    
    // Calculate loads for each approach
    std::vector<double> homo_loads(homo_nodes.size(), 0.0);
    std::vector<double> without_loads(nodes.size(), 0.0);
    std::vector<double> with_loads(nodes.size(), 0.0);
    
    // For homogeneous, use perf=1.0 for all nodes; for others, use actual perf factors
    for (size_t i = 0; i < tasks.size(); ++i) {
        // Homogeneous: calculate load as if all nodes have performance 1.0
        homo_loads[homo_assign[i]] += tasks[i].base_time / 1.0;
        // Without Rij and With Rij: use actual performance factors
        without_loads[without_assign[i]] += tasks[i].base_time / nodes[without_assign[i]].performance_factor;
        with_loads[with_assign[i]] += tasks[i].base_time / nodes[with_assign[i]].performance_factor;
    }
    
    // Calculate metrics
    double homo_makespan = *std::max_element(homo_loads.begin(), homo_loads.end());
    double without_makespan = *std::max_element(without_loads.begin(), without_loads.end());
    double with_makespan = *std::max_element(with_loads.begin(), with_loads.end());
    
    double total_capacity = std::accumulate(perf_factors.begin(), perf_factors.end(), 0.0);
    double ideal_makespan = total_work / total_capacity;
    
    double homo_efficiency = ideal_makespan / homo_makespan;
    double without_efficiency = ideal_makespan / without_makespan;
    double with_efficiency = ideal_makespan / with_makespan;

    // // Debug: Print first 10 assignments
    // std::cout << "First 10 task assignments:\n";
    // std::cout << "Task | Homo | Without | With\n";
    // for (int i = 0; i < std::min(10, (int)tasks.size()); ++i) {
    //     std::cout << std::setw(4) << i << " | " 
    //             << std::setw(4) << homo_assign[i] << " | "
    //             << std::setw(7) << without_assign[i] << " | "
    //             << std::setw(4) << with_assign[i] << "\n";
    // }
    
    // Print comparison
    std::cout << "=== RESULTS COMPARISON ===\n";
    std::cout << "Ideal makespan (perfect balance): " << std::setprecision(2) << std::fixed << ideal_makespan << "\n\n";
    
    std::cout << "Approach                | Makespan | Efficiency | Improvement\n";
    std::cout << "------------------------|----------|------------|------------\n";
    std::cout << "1. Homogeneous Knapsack | " << std::setw(8) << homo_makespan 
              << " | " << std::setw(10) << std::setprecision(3) << homo_efficiency << " | baseline\n";
    std::cout << "2. Without Rij          | " << std::setw(8) << without_makespan 
              << " | " << std::setw(10) << without_efficiency 
              << " | " << std::setprecision(1) << ((without_efficiency - homo_efficiency) / homo_efficiency * 100) << "%\n";
    std::cout << "3. With Rij             | " << std::setw(8) << with_makespan 
              << " | " << std::setw(10) << with_efficiency 
              << " | " << std::setprecision(1) << ((with_efficiency - homo_efficiency) / homo_efficiency * 100) << "%\n\n";
    
    std::cout << "Rij advantage over basic heterogeneous: " 
              << std::setprecision(1) << ((with_efficiency - without_efficiency) / without_efficiency * 100) << "%\n\n";
    
    // Show load distribution
    std::cout << "Load Distribution:\n";
    std::cout << "Node | Perf | Homogeneous | Without Rij |  With Rij   \n";
    std::cout << "-----|------|-------------|-------------|-------------\n";
    
    for (size_t i = 0; i < nodes.size(); ++i) {
        std::cout << std::setw(4) << i << " | " 
                  << std::setw(4) << std::setprecision(1) << perf_factors[i] << " | "
                  << std::setw(11) << std::setprecision(2) << homo_loads[i] << " | "
                  << std::setw(11) << without_loads[i] << " | "
                  << std::setw(11) << with_loads[i] << "\n";
    }
    
    // Generate CSV data for plotting
    std::system("mkdir -p plots");
    
    // CSV for load distribution comparison
    std::ofstream csv_file("plots/load_distribution.csv");
    csv_file << "Node,Performance,Homogeneous,Without_Rij,With_Rij\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        csv_file << i << "," << perf_factors[i] << "," 
                 << homo_loads[i] << "," << without_loads[i] << "," << with_loads[i] << "\n";
    }
    csv_file.close();
    
    // CSV for efficiency comparison
    std::ofstream eff_file("plots/efficiency_comparison.csv");
    eff_file << "Approach,Makespan,Efficiency\n";
    eff_file << "Homogeneous," << homo_makespan << "," << homo_efficiency << "\n";
    eff_file << "Without_Rij," << without_makespan << "," << without_efficiency << "\n";
    eff_file << "With_Rij," << with_makespan << "," << with_efficiency << "\n";
    eff_file.close();
    
    std::cout << "\nCSV files generated in plots/ directory\n";

    std::cout << "\n=== HETEROGENEITY LEVEL ANALYSIS ===\n";
    auto het_results = test_heterogeneity_levels(tasks);

    // Save heterogeneity results to CSV
    std::ofstream het_file("plots/heterogeneity_analysis.csv");
    het_file << "Heterogeneity_Factor,Homo_Efficiency,Without_Rij_Efficiency,With_Rij_Efficiency,";
    het_file << "Homo_Makespan,Without_Rij_Makespan,With_Rij_Makespan\n";

    for (const auto& result : het_results) {
        het_file << result.heterogeneity_factor << ","
                << result.homo_efficiency << ","
                << result.without_rij_efficiency << ","
                << result.with_rij_efficiency << ","
                << result.homo_makespan << ","
                << result.without_rij_makespan << ","
                << result.with_rij_makespan << "\n";
    }
    het_file.close();

    std::cout << "Heterogeneity analysis saved to plots/heterogeneity_analysis.csv\n";

        return 0;
}

